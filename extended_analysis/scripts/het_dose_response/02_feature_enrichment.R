#!/usr/bin/env Rscript
# =============================================================================
# het_dose_response/02_feature_enrichment.R — peak features that distinguish dose-response classes
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Given that peaks are haploinsufficient vs buffered vs linear, what peak-level
# features predict which class they belong to? Features examined:
#   1. Baseline WT accessibility (VST mean, DESeq2 baseMean)
#   2. Peak width
#   3. Genomic annotation (ChIPseeker: promoter / intron / intergenic / ...)
#   4. Distance to nearest TSS
#   5. D5 WT CUT&RUN overlap: ARID1A, H3K27ac, Tbet, BATF, ETS1 (also H3K27me3)
#   6. TF co-occupancy count (0..4 of {ARID1A, Tbet, BATF, ETS1})
#
# Motif enrichment is not run here: the script writes per-class BEDs and a
# HOMER driver, run_homer_dose_response.sh (haploinsufficient-lost vs
# buffered-lost background, per subset; requires HOMER findMotifsGenome.pl and
# the mm39 genome FASTA). 03_known_motif_summary.R and 04_denovo_motif_summary.R
# summarize its output.
#
# Inputs:  results/extended_analysis/het_dose_response/dose_classes_combined.csv (01_dose_classes.R)
#          results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/
#            consensus_peaks.mRp.clN.featureCounts.txt
#          results/atac/differential/consensus_peaks_annotated.csv (core/02_atacseq_analysis.R)
#          results/cutrun/03_peak_calling/04_called_peaks/macs2/*_D5_WT_R*.macs2_peaks.narrowPeak
# Outputs: results/extended_analysis/het_dose_response/
#            feature_summary_by_class.csv, feature_enrichment_vs_insensitive.csv,
#            beds/, run_homer_dose_response.sh, het_feature_enrichment.{RData,log}
#          figures/extended_analysis/het_dose_response/
# Usage:   Rscript extended_analysis/scripts/het_dose_response/02_feature_enrichment.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(GenomicFeatures)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})
select <- dplyr::select
filter <- dplyr::filter
count  <- dplyr::count

outdir <- file.path(paths$ext_results, "het_dose_response")
figdir <- file.path(paths$ext_figures, "het_dose_response")
stopifnot("Run het_dose_response/01_dose_classes.R first" = dir.exists(outdir))

logfile <- file.path(outdir, "het_feature_enrichment.log")
sink(logfile, split = TRUE)
on.exit(sink(), add = TRUE)

message("=== het_dose_response/02_feature_enrichment.R started: ", Sys.time(), " ===")

# =============================================================================
# 1. Load dose-classes table + peak coordinates
# =============================================================================

message("\n=== Section 1: Load inputs ===")

combined <- read.csv(file.path(outdir, "dose_classes_combined.csv"),
                     stringsAsFactors = FALSE)
combined$class <- factor(combined$class, levels = c(
  "insensitive", "buffered", "linear", "haploinsufficient",
  "nonmonotonic", "other_responsive"))

# featureCounts peak coordinates (same file used in 01_dose_classes.R)
fc_file <- file.path(paths$atac,
  "bowtie2/merged_replicate/macs2/narrow_peak/consensus",
  "consensus_peaks.mRp.clN.featureCounts.txt")
raw <- read.delim(fc_file, comment.char = "#", check.names = FALSE)
peak_info <- raw[, 1:6]
colnames(peak_info) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

peak_gr <- GRanges(seqnames = peak_info$Chr,
                   ranges   = IRanges(peak_info$Start, peak_info$End),
                   peak_id  = peak_info$Geneid,
                   width_bp = peak_info$Length)
names(peak_gr) <- peak_info$Geneid
message(sprintf("  Loaded %d consensus peaks", length(peak_gr)))

# =============================================================================
# 2. ChIPseeker genomic annotation (once, across all peaks)
# =============================================================================

message("\n=== Section 2: Genomic annotation ===")

anno_csv <- file.path(paths$results,
  "atac/differential/consensus_peaks_annotated.csv")
if (file.exists(anno_csv)) {
  anno_df <- read.csv(anno_csv, stringsAsFactors = FALSE)
  # Harmonize column name: existing file uses 'peak_id'
  stopifnot("peak_id" %in% colnames(anno_df))
  message(sprintf("  Re-using %s (%d rows)", anno_csv, nrow(anno_df)))
} else {
  # This file is a manuscript input (Fig 1/2, 4, 5 read it); only core/02_atacseq_analysis.R writes it
  stop("consensus_peaks_annotated.csv not found — run scripts/core/02_atacseq_analysis.R first: ",
       anno_csv)
}

# Simplify annotation into a few categories
simplify_anno <- function(a) {
  case_when(
    grepl("^Promoter", a)                       ~ "Promoter",
    grepl("^5' UTR|5UTR", a)                    ~ "5'UTR",
    grepl("^3' UTR|3UTR", a)                    ~ "3'UTR",
    grepl("Intron", a)                          ~ "Intron",
    grepl("Exon",   a)                          ~ "Exon",
    grepl("Intergenic|Downstream|Distal", a)    ~ "Intergenic",
    TRUE                                        ~ "Other"
  )
}
anno_df$feature <- simplify_anno(anno_df$annotation)

# =============================================================================
# 3. D5 WT CUT&RUN overlap features
# =============================================================================

message("\n=== Section 3: CUT&RUN D5 WT overlap features ===")

cutrun_peak_dir <- file.path(paths$cutrun,
  "03_peak_calling/04_called_peaks/macs2")

cutrun_tracks <- list(
  ARID1A    = c("ARID1A_D5_WT_R1", "ARID1A_D5_WT_R2"),
  H3K27ac   = "H3K27ac_D5_WT_R1",
  H3K27me3  = NULL,   # not called at D5 WT
  Tbet      = "Tbet_D5_WT_R1",
  BATF      = "BATF_D5_WT_R1",
  ETS1      = "ETS1_D5_WT_R1"
)

load_cutrun <- function(names_vec) {
  if (is.null(names_vec)) return(GRanges())
  files <- file.path(cutrun_peak_dir,
                     paste0(names_vec, ".macs2_peaks.narrowPeak"))
  files <- files[file.exists(files)]
  if (length(files) == 0) return(GRanges())
  grs <- lapply(files, import, format = "narrowPeak")
  reduce(do.call(c, grs))
}

cutrun_gr <- lapply(cutrun_tracks, load_cutrun)
for (tf in names(cutrun_gr)) {
  message(sprintf("  %-9s : %s peaks",
                  tf,
                  if (length(cutrun_gr[[tf]])) length(cutrun_gr[[tf]]) else "none"))
}

# Harmonize seqlevels (featureCounts uses chr1..chrY style already)
peak_gr_chr <- peak_gr
overlap_flags <- sapply(cutrun_gr, function(gr) {
  if (length(gr) == 0) return(rep(NA, length(peak_gr_chr)))
  seqlevelsStyle(gr) <- seqlevelsStyle(peak_gr_chr)[1]
  overlapsAny(peak_gr_chr, gr)
})
overlap_df <- data.frame(peak_id = names(peak_gr_chr), overlap_flags,
                         stringsAsFactors = FALSE,
                         check.names = FALSE)
colnames(overlap_df)[-1] <- paste0("has_", colnames(overlap_df)[-1])

# TF co-occupancy count (of present tracks among ARID1A/Tbet/BATF/ETS1)
tf_cols <- intersect(c("has_ARID1A", "has_Tbet", "has_BATF", "has_ETS1"),
                     colnames(overlap_df))
overlap_df$tf_cooccupancy <- rowSums(overlap_df[, tf_cols, drop = FALSE],
                                     na.rm = TRUE)

# =============================================================================
# 4. Assemble feature table and join with class
# =============================================================================

message("\n=== Section 4: Assemble feature table ===")

feat_df <- peak_info %>%
  select(peak_id = Geneid, chr = Chr, start = Start, end = End,
         peak_width = Length) %>%
  left_join(anno_df %>% select(peak_id, feature,
                               distance_to_TSS = distanceToTSS),
            by = "peak_id") %>%
  left_join(overlap_df, by = "peak_id")

feat_class <- combined %>%
  select(peak_id, subset, class, direction,
         vst_WT, vst_Het, vst_KO,
         lfc_HetWT, lfc_KOWT, baseMean) %>%
  left_join(feat_df, by = "peak_id")

# Restrict downstream comparisons to D8 "lost" peaks (the biology we care about)
lost_classes <- c("buffered", "linear", "haploinsufficient",
                  "nonmonotonic", "other_responsive", "insensitive")
feat_lost <- feat_class %>%
  filter(direction == "lost" | class == "insensitive") %>%
  mutate(class = factor(class, levels = lost_classes))

# =============================================================================
# 5. Summary statistics per class (per subset)
# =============================================================================

message("\n=== Section 5: Feature summaries per class ===")

summarize_feat <- function(df) {
  df %>%
    group_by(subset, class) %>%
    summarize(
      n_peaks           = n(),
      med_peak_width    = median(peak_width,      na.rm = TRUE),
      med_baseWT_VST    = median(vst_WT,          na.rm = TRUE),
      med_baseMean      = median(baseMean,        na.rm = TRUE),
      med_dist_TSS      = median(abs(distance_to_TSS), na.rm = TRUE),
      pct_promoter      = 100 * mean(feature == "Promoter",    na.rm = TRUE),
      pct_intergenic    = 100 * mean(feature == "Intergenic",  na.rm = TRUE),
      pct_intron        = 100 * mean(feature == "Intron",      na.rm = TRUE),
      pct_ARID1A        = 100 * mean(has_ARID1A,  na.rm = TRUE),
      pct_H3K27ac       = 100 * mean(has_H3K27ac, na.rm = TRUE),
      pct_Tbet          = 100 * mean(has_Tbet,    na.rm = TRUE),
      pct_BATF          = 100 * mean(has_BATF,    na.rm = TRUE),
      pct_ETS1          = 100 * mean(has_ETS1,    na.rm = TRUE),
      mean_cooccupancy  = mean(tf_cooccupancy,   na.rm = TRUE),
      .groups = "drop"
    )
}

feat_summary <- summarize_feat(feat_lost)
write.csv(feat_summary,
          file.path(outdir, "feature_summary_by_class.csv"),
          row.names = FALSE)
message("  Wrote feature_summary_by_class.csv")
print(as.data.frame(feat_summary))

# Enrichment vs insensitive: fold-change + Fisher p-value, TF overlap per class
fisher_feature <- function(df, feature_col, subset_, class_) {
  bg <- df$class == "insensitive" & df$subset == subset_
  fg <- df$class == class_        & df$subset == subset_
  if (sum(fg) < 10 || sum(bg) < 10) return(NULL)
  a <- sum(df[[feature_col]][fg], na.rm = TRUE)
  b <- sum(fg) - a
  c <- sum(df[[feature_col]][bg], na.rm = TRUE)
  d <- sum(bg) - c
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = FALSE)
  ft  <- fisher.test(mat)
  data.frame(subset = subset_, class = class_, feature = feature_col,
             n_class = sum(fg), frac_class = a / sum(fg),
             frac_bg = c / sum(bg),
             odds_ratio = unname(ft$estimate),
             pvalue = ft$p.value,
             stringsAsFactors = FALSE)
}

feat_lost$is_promoter <- feat_lost$feature == "Promoter"
feat_lost$is_distal   <- feat_lost$feature == "Intergenic" |
                        feat_lost$feature == "Intron"

classes_to_test <- c("buffered", "linear", "haploinsufficient", "nonmonotonic")
features_to_test <- c("has_ARID1A", "has_H3K27ac", "has_Tbet",
                      "has_BATF", "has_ETS1", "is_promoter", "is_distal")

enrich_rows <- list()
for (sub in levels(factor(feat_lost$subset))) {
  for (cls in classes_to_test) {
    for (ft in features_to_test) {
      enrich_rows[[length(enrich_rows) + 1]] <-
        fisher_feature(feat_lost, ft, sub, cls)
    }
  }
}
enrich_df <- bind_rows(enrich_rows)
enrich_df$padj <- p.adjust(enrich_df$pvalue, method = "BH")
write.csv(enrich_df, file.path(outdir, "feature_enrichment_vs_insensitive.csv"),
          row.names = FALSE)
message("  Wrote feature_enrichment_vs_insensitive.csv")

# =============================================================================
# 6. Figures: class-level feature profiles
# =============================================================================

message("\n=== Section 6: Figures ===")

class_colors <- c(
  "insensitive"        = "#BBBBBB",
  "buffered"           = "#4DAF4A",
  "linear"             = "#377EB8",
  "haploinsufficient"  = "#E41A1C",
  "nonmonotonic"       = "#984EA3",
  "other_responsive"   = "#999999"
)

focus_classes <- c("insensitive", "buffered", "linear",
                   "haploinsufficient", "nonmonotonic")

# Baseline WT accessibility by class, per subset
plot_df <- feat_lost %>%
  filter(class %in% focus_classes) %>%
  mutate(class = factor(class, levels = focus_classes))

p_baseWT <- ggplot(plot_df, aes(x = class, y = vst_WT, fill = class)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ subset) +
  scale_fill_manual(values = class_colors) +
  labs(title = "Baseline WT accessibility by dose-response class",
       y = "VST (WT)", x = NULL) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
save_figure(p_baseWT, "baseline_WT_by_class", width = 10, height = 4, dir = figdir)

# Peak width
p_width <- ggplot(plot_df, aes(x = class, y = peak_width, fill = class)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ subset) +
  scale_fill_manual(values = class_colors) +
  scale_y_log10() +
  labs(title = "Peak width by class", y = "peak width (bp, log10)", x = NULL) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
save_figure(p_width, "peak_width_by_class", width = 10, height = 4, dir = figdir)

# Genomic annotation breakdown
anno_plot <- plot_df %>%
  count(subset, class, feature) %>%
  group_by(subset, class) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

p_anno <- ggplot(anno_plot, aes(x = class, y = frac, fill = feature)) +
  geom_col() +
  facet_wrap(~ subset) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Genomic annotation by class",
       y = "Fraction", x = NULL, fill = "Feature") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_figure(p_anno, "annotation_by_class", width = 12, height = 4, dir = figdir)

# CUT&RUN overlap % (long form)
cutrun_cols <- c("has_ARID1A", "has_H3K27ac",
                 "has_Tbet", "has_BATF", "has_ETS1")
cutrun_cols <- intersect(cutrun_cols, colnames(plot_df))
cr_long <- plot_df %>%
  select(subset, class, all_of(cutrun_cols)) %>%
  pivot_longer(all_of(cutrun_cols), names_to = "track", values_to = "hit") %>%
  mutate(track = sub("^has_", "", track)) %>%
  group_by(subset, class, track) %>%
  summarize(pct = 100 * mean(hit, na.rm = TRUE), .groups = "drop")

p_cr <- ggplot(cr_long, aes(x = class, y = pct, fill = class)) +
  geom_col() +
  facet_grid(track ~ subset) +
  scale_fill_manual(values = class_colors) +
  labs(title = "D5 WT CUT&RUN overlap by class",
       y = "% peaks with overlap", x = NULL) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
save_figure(p_cr, "cutrun_overlap_by_class", width = 10, height = 10, dir = figdir)

# TF cooccupancy distribution
if ("tf_cooccupancy" %in% colnames(plot_df)) {
  cooc_plot <- plot_df %>%
    count(subset, class, tf_cooccupancy) %>%
    group_by(subset, class) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()

  p_cooc <- ggplot(cooc_plot, aes(x = class, y = frac,
                                  fill = factor(tf_cooccupancy))) +
    geom_col() +
    facet_wrap(~ subset) +
    scale_fill_brewer(palette = "YlOrRd",
                      name = "# TFs bound (ARID1A/Tbet/BATF/ETS1)") +
    labs(title = "TF co-occupancy depth by class",
         y = "Fraction", x = NULL) +
    theme_paper +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_figure(p_cooc, "tf_cooccupancy_by_class",
              width = 12, height = 4, dir = figdir)
}

# =============================================================================
# 7. Export BEDs + write HOMER motif script
# =============================================================================

message("\n=== Section 7: Export BEDs for HOMER ===")

bed_dir <- file.path(outdir, "beds")
dir.create(bed_dir, showWarnings = FALSE)

write_class_bed <- function(df, sub, cls, direction_ = "lost") {
  rows <- df %>%
    filter(subset == sub, class == cls, (direction == direction_ |
                                         cls == "insensitive"))
  if (nrow(rows) < 20) return(NA_character_)
  bed <- rows %>% select(chr, start, end, peak_id) %>%
    mutate(score = 0, strand = "+")
  path <- file.path(bed_dir, sprintf("%s_%s_%s.bed", sub, cls, direction_))
  write.table(bed, path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE)
  path
}

bed_paths <- list()
for (sub in levels(factor(feat_lost$subset))) {
  for (cls in c("buffered", "linear", "haploinsufficient",
                "insensitive")) {
    bed_paths[[paste(sub, cls, sep = "_")]] <-
      write_class_bed(feat_lost, sub, cls, "lost")
  }
}
message(sprintf("  Wrote %d class BEDs to %s",
                sum(!is.na(unlist(bed_paths))), bed_dir))

# HOMER: haploinsufficient-lost vs buffered-lost (the key comparison),
# and vs insensitive (as a looser comparator), per subset.
homer_outdir_base <- file.path(outdir, "homer_motifs")
dir.create(homer_outdir_base, showWarnings = FALSE)

homer_cmds <- c("#!/bin/bash",
  "# HOMER motif enrichment: what separates ARID1A-fragile peaks?",
  "set -euo pipefail",
  sprintf("export PATH=%s:$PATH", "/usr/bin/homer/bin"),
  "")

add_cmd <- function(fg, bg, tag) {
  if (is.na(fg) || (!is.null(bg) && is.na(bg))) return(invisible())
  out <- file.path(homer_outdir_base, tag)
  bg_arg <- if (is.null(bg)) "" else sprintf("-bg %s", bg)
  cmd <- sprintf(
    "findMotifsGenome.pl %s %s %s -size 200 -mask -p 8 %s > %s/homer.log 2>&1",
    fg, genome$fasta, out, bg_arg, out)
  homer_cmds <<- c(homer_cmds, sprintf("mkdir -p %s", out), cmd, "")
}

for (sub in levels(factor(feat_lost$subset))) {
  haplo <- bed_paths[[paste(sub, "haploinsufficient", sep = "_")]]
  buff  <- bed_paths[[paste(sub, "buffered",          sep = "_")]]
  ins   <- bed_paths[[paste(sub, "insensitive",       sep = "_")]]
  lin   <- bed_paths[[paste(sub, "linear",            sep = "_")]]

  # Key contrast: haploinsufficient vs buffered (what makes a site fragile?)
  add_cmd(haplo, buff, sprintf("%s_haplo_vs_buffered", sub))
  # Looser: haploinsufficient vs insensitive
  add_cmd(haplo, ins,  sprintf("%s_haplo_vs_insensitive", sub))
  # Buffered vs insensitive (what's dose-robust about these)
  add_cmd(buff,  ins,  sprintf("%s_buffered_vs_insensitive", sub))
  # Linear vs insensitive
  add_cmd(lin,   ins,  sprintf("%s_linear_vs_insensitive", sub))
}

homer_script <- file.path(outdir, "run_homer_dose_response.sh")
writeLines(homer_cmds, homer_script)
Sys.chmod(homer_script, "0755")
message("  Wrote ", homer_script)
message("  To run: bash ", homer_script)

# =============================================================================
# 8. Save RData checkpoint
# =============================================================================

message("\n=== Section 8: Save checkpoint ===")

save(feat_class, feat_lost, feat_summary, enrich_df,
     cutrun_gr, peak_gr,
     file = file.path(outdir, "het_feature_enrichment.RData"))

message("\n=== Done: ", Sys.time(), " ===")
