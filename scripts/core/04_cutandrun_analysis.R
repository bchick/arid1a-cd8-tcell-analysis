#!/usr/bin/env Rscript
# =============================================================================
# core/04_cutandrun_analysis.R — CUT&RUN peak-level analysis
# McDonald, Chick et al. 2023 Immunity 56:1303 — core analysis
#
# Antibodies: ARID1A, H3K27ac, H3K27me3, T-bet, BATF, ETS1 (GSE228380).
#   1. ARID1A binding across activation (Naive, 48h, D5)
#   2. TF / histone mark peak overlap, WT vs KO at D5 (n = 1 per group)
#   3. ARID1A-dependent vs -independent TF binding at D5
#   4. Genomic annotation of peaks (ChIPseeker, GENCODE vM35)
#   5. Manifest of D8 subset BAMs (WT TE/EEC/MP) for signal-level analyses
# Signal heatmaps are drawn by core/06_deeptools_heatmaps.sh.
#
# Inputs:  results/cutrun/02_alignment/bowtie2/target/markdup/*.bam
#          results/cutrun/03_peak_calling/04_called_peaks/macs2/*.narrowPeak
#          data/reference/gencode.vM35.primary_assembly.annotation.gtf
# Outputs: results/cutrun/differential/ (sample manifest, D5 WT/KO overlap BEDs
#            and tables, ARID1A temporal BEDs, TF_ARID1A_dependency.csv,
#            annotation tables, cutandrun_analysis.RData)
#          figures/cutrun/
# Usage:   Rscript scripts/core/04_cutandrun_analysis.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(GenomicFeatures)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggvenn)
})
select <- dplyr::select

outdir <- file.path(paths$results, "cutrun/differential")
figdir <- file.path(paths$figures, "cutrun")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Build sample manifest from actual files on disk
# =============================================================================

message("=== Building CUT&RUN sample manifest ===")

bam_dir  <- file.path(paths$cutrun, "02_alignment/bowtie2/target/markdup")
peak_dir <- file.path(paths$cutrun, "03_peak_calling/04_called_peaks/macs2")

# Map known files to metadata
samples <- tribble(
  ~sample_id,                  ~antibody,   ~timepoint, ~genotype, ~cell_subset, ~treatment,
  "ARID1A_Naive_WT_R1",        "ARID1A",    "Naive",    "WT",      NA,           NA,
  "ARID1A_48h_WT_R1",          "ARID1A",    "48h",      "WT",      NA,           NA,
  "ARID1A_48h_WT_R2",          "ARID1A",    "48h",      "WT",      NA,           NA,
  "ARID1A_D5_WT_R1",           "ARID1A",    "D5",       "WT",      NA,           NA,
  "ARID1A_D5_WT_R2",           "ARID1A",    "D5",       "WT",      NA,           NA,
  "ARID1A_D5_KO_R1",           "ARID1A",    "D5",       "KO",      NA,           NA,
  "ARID1A_D5_TbetKO_R1",       "ARID1A",    "D5",       "TbetKO",  NA,           NA,
  "ARID1A_D8_WT_TE_R1",        "ARID1A",    "D8",       "WT",      "TE",         NA,
  "ARID1A_D8_WT_EEC_R1",       "ARID1A",    "D8",       "WT",      "EEC",        NA,
  "ARID1A_D8_WT_MP_R1",        "ARID1A",    "D8",       "WT",      "MP",         NA,
  "BATF_D5_WT_R1",             "BATF",      "D5",       "WT",      NA,           NA,
  "BATF_D5_KO_R1",             "BATF",      "D5",       "KO",      NA,           NA,
  "ETS1_D5_WT_R1",             "ETS1",      "D5",       "WT",      NA,           NA,
  "ETS1_D5_KO_R1",             "ETS1",      "D5",       "KO",      NA,           NA,
  "Tbet_D5_WT_R1",             "Tbet",      "D5",       "WT",      NA,           NA,
  "Tbet_D5_KO_R1",             "Tbet",      "D5",       "KO",      NA,           NA,
  "Tbet_D5_WT_TbetOE_R1",      "Tbet",      "D5",       "WT",      NA,           "TbetOE",
  "Tbet_D5_KO_TbetOE_R1",      "Tbet",      "D5",       "KO",      NA,           "TbetOE",
  "H3K27ac_D5_WT_R1",          "H3K27ac",   "D5",       "WT",      NA,           NA,
  "H3K27ac_D5_KO_R1",          "H3K27ac",   "D5",       "KO",      NA,           NA,
  "H3K27ac_D8_WT_TE_R1",       "H3K27ac",   "D8",       "WT",      "TE",         NA,
  "H3K27ac_D8_WT_EEC_R1",      "H3K27ac",   "D8",       "WT",      "EEC",        NA,
  "H3K27ac_D8_WT_MP_R1",       "H3K27ac",   "D8",       "WT",      "MP",         NA,
  "H3K27me3_D8_WT_TE_R1",      "H3K27me3",  "D8",       "WT",      "TE",         NA,
  "H3K27me3_D8_WT_EEC_R1",     "H3K27me3",  "D8",       "WT",      "EEC",        NA,
  "H3K27me3_D8_WT_MP_R1",      "H3K27me3",  "D8",       "WT",      "MP",         NA,
  "IgG_D5_WT_R1",              "IgG",       "D5",       "WT",      NA,           NA,
  "IgG_D5_KO_R1",              "IgG",       "D5",       "KO",      NA,           NA,
) %>%
  mutate(
    bam  = file.path(bam_dir,  paste0(sample_id, ".target.markdup.sorted.bam")),
    peak = file.path(peak_dir, paste0(sample_id, ".macs2_peaks.narrowPeak"))
  )

# Check file existence
samples <- samples %>%
  mutate(
    bam_exists  = file.exists(bam),
    peak_exists = file.exists(peak)
  )

message(sprintf("  BAMs:  %d/%d found", sum(samples$bam_exists),  nrow(samples)))
message(sprintf("  Peaks: %d/%d found", sum(samples$peak_exists), nrow(samples)))

# Samples with peaks (D5 antibodies only; D8 and Naive BAMs exist but weren't peak-called individually)
peak_samples <- samples %>% filter(peak_exists)
message("  Peak samples: ", paste(peak_samples$sample_id, collapse = ", "))

write_csv(samples, file.path(outdir, "cutrun_sample_manifest.csv"))

# =============================================================================
# 2. Load all peaks into GRanges
# =============================================================================

message("=== Loading peaks as GRanges ===")

load_narrowpeak <- function(path) {
  if (!file.exists(path)) return(GRanges())
  tryCatch(
    import(path, format = "narrowPeak"),
    error = function(e) {
      message("  Warning: could not load ", basename(path), ": ", e$message)
      GRanges()
    }
  )
}

peaks_gr <- setNames(
  lapply(peak_samples$peak, load_narrowpeak),
  peak_samples$sample_id
)

peak_counts <- sapply(peaks_gr, length)
message("Peak counts:")
for (nm in names(peak_counts)) message(sprintf("  %s: %d peaks", nm, peak_counts[nm]))

# Safe accessor — returns GRanges() for samples with no peaks (e.g. Naive/48h/D8 BAM-only)
get_gr <- function(id) { x <- peaks_gr[[id]]; if (is.null(x)) GRanges() else x }

# =============================================================================
# 3. Build TxDb for annotation
# =============================================================================

message("=== Building TxDb ===")
txdb <- get_txdb()

annotate_peaks <- function(gr, txdb, label) {
  if (length(gr) == 0) return(NULL)
  anno <- annotatePeak(gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
  df <- as.data.frame(anno) %>% as_tibble() %>% mutate(sample = label)
  # Map gene symbols via AnnotationDbi (ChIPseeker misses them with custom TxDb)
  if ("geneId" %in% colnames(df) && any(!is.na(df$geneId))) {
    gene_ids <- gsub("\\.\\d+$", "", df$geneId)
    sym_map <- tryCatch(
      AnnotationDbi::select(org.Mm.eg.db, keys = gene_ids,
                            keytype = "ENSEMBL", columns = "SYMBOL"),
      error = function(e) data.frame(ENSEMBL = character(), SYMBOL = character())
    )
    df <- df %>%
      mutate(ensembl_clean = gsub("\\.\\d+$", "", geneId)) %>%
      left_join(sym_map, by = c("ensembl_clean" = "ENSEMBL")) %>%
      mutate(SYMBOL = coalesce(SYMBOL, geneId))
  }
  df
}

# =============================================================================
# 4. WT vs KO peak overlap (D5 antibodies — n=1 per group, use bedtools logic)
# =============================================================================

message("=== D5 WT vs KO peak overlap analysis ===")

d5_antibodies <- c("ARID1A", "BATF", "ETS1", "Tbet", "H3K27ac")

overlap_results <- list()

for (ab in d5_antibodies) {
  wt_ids <- peak_samples %>% filter(antibody == ab, genotype == "WT", timepoint == "D5",
                                    is.na(treatment)) %>% pull(sample_id)
  ko_ids <- peak_samples %>% filter(antibody == ab, genotype == "KO", timepoint == "D5",
                                    is.na(treatment)) %>% pull(sample_id)

  if (length(wt_ids) == 0 || length(ko_ids) == 0) {
    message(sprintf("  %s: missing WT or KO peaks, skipping", ab))
    next
  }

  # Merge multiple WT reps with reduce()
  wt_gr <- Reduce(c, peaks_gr[wt_ids]) %>% reduce()
  ko_gr <- peaks_gr[[ko_ids[1]]]

  # Overlap: shared = reciprocal ≥1bp overlap
  ov_wt_in_ko <- countOverlaps(wt_gr, ko_gr) > 0
  ov_ko_in_wt <- countOverlaps(ko_gr, wt_gr) > 0

  wt_only     <- wt_gr[!ov_wt_in_ko]
  ko_only     <- ko_gr[!ov_ko_in_wt]
  shared_wt   <- wt_gr[ov_wt_in_ko]

  overlap_results[[ab]] <- tibble(
    antibody   = ab,
    WT_only    = length(wt_only),
    KO_only    = length(ko_only),
    Shared     = length(shared_wt),
    WT_total   = length(wt_gr),
    KO_total   = length(ko_gr),
    pct_WT_lost = length(wt_only) / length(wt_gr) * 100
  )

  message(sprintf("  %s: WT=%d, KO=%d, Shared=%d, WT-only=%d (%.1f%% lost in KO)",
                  ab, length(wt_gr), length(ko_gr),
                  length(shared_wt), length(wt_only),
                  length(wt_only) / length(wt_gr) * 100))

  # Export BED files
  export(wt_only,   file.path(outdir, paste0(ab, "_D5_WT_only.bed")))
  export(ko_only,   file.path(outdir, paste0(ab, "_D5_KO_only.bed")))
  export(shared_wt, file.path(outdir, paste0(ab, "_D5_shared.bed")))
}

overlap_df <- bind_rows(overlap_results)
write_csv(overlap_df, file.path(outdir, "D5_wt_ko_peak_overlap.csv"))
message(overlap_df)

# Bar plot: peak categories per antibody
overlap_long <- overlap_df %>%
  select(antibody, WT_only, Shared, KO_only) %>%
  pivot_longer(-antibody, names_to = "category", values_to = "n_peaks") %>%
  mutate(category = factor(category, levels = c("WT_only", "Shared", "KO_only")))

p_overlap <- ggplot(overlap_long, aes(x = antibody, y = n_peaks, fill = category)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(WT_only = "#1F77B4", Shared = "#AEC7E8", KO_only = "#2CA02C"),
                    labels = c("WT-specific", "Shared", "KO-specific")) +
  labs(title = "CUT&RUN peak overlap: WT vs KO at D5",
       y = "Number of peaks", x = "Antibody", fill = "") +
  theme_paper

save_figure(p_overlap, "D5_peak_overlap_stacked", width = 7, height = 5, dir = figdir)

# % WT peaks retained bar chart
p_pct <- ggplot(overlap_df, aes(x = antibody, y = 100 - pct_WT_lost, fill = antibody)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.0f%%", 100 - pct_WT_lost)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(ARID1A = "#B0B0B0", BATF = "#CC6677",
                                ETS1 = "#DDCC77", Tbet = "#AA4499", H3K27ac = "#2CA02C")) +
  scale_y_continuous(limits = c(0, 110)) +
  labs(title = "% WT peaks retained in ARID1A KO (D5)",
       y = "% WT peaks overlapping KO", x = "Antibody") +
  theme_paper + theme(legend.position = "none")

save_figure(p_pct, "D5_pct_peaks_retained_KO", width = 6, height = 5, dir = figdir)

# =============================================================================
# 5. ARID1A temporal binding: Naive → 48h → D5
# =============================================================================

message("=== ARID1A temporal binding across timepoints ===")

arid1a_tp <- list(
  Naive = reduce(get_gr("ARID1A_Naive_WT_R1")),
  h48   = reduce(c(get_gr("ARID1A_48h_WT_R1"), get_gr("ARID1A_48h_WT_R2"))),
  D5    = reduce(c(get_gr("ARID1A_D5_WT_R1"),  get_gr("ARID1A_D5_WT_R2")))
)

# 3-way overlap summary
n_naive <- length(arid1a_tp$Naive)
n_48h   <- length(arid1a_tp$h48)
n_d5    <- length(arid1a_tp$D5)
n_naive_48h  <- sum(countOverlaps(arid1a_tp$Naive, arid1a_tp$h48) > 0)
n_naive_d5   <- sum(countOverlaps(arid1a_tp$Naive, arid1a_tp$D5)  > 0)
n_48h_d5     <- sum(countOverlaps(arid1a_tp$h48,   arid1a_tp$D5)  > 0)
n_all        <- sum(countOverlaps(arid1a_tp$Naive, arid1a_tp$h48) > 0 &
                    countOverlaps(arid1a_tp$Naive, arid1a_tp$D5)  > 0)

arid1a_temporal_summary <- tibble(
  timepoint = c("Naive", "48h", "D5"),
  n_peaks   = c(n_naive, n_48h, n_d5)
)
write_csv(arid1a_temporal_summary, file.path(outdir, "ARID1A_temporal_peak_counts.csv"))
message("ARID1A peaks: Naive=", n_naive, " 48h=", n_48h, " D5=", n_d5)

# Bar plot: ARID1A peak count over time
p_arid1a_tp <- ggplot(arid1a_temporal_summary,
                       aes(x = factor(timepoint, c("Naive","48h","D5")), y = n_peaks)) +
  geom_col(fill = "#B0B0B0", color = "black", width = 0.6) +
  geom_text(aes(label = scales::comma(n_peaks)), vjust = -0.3, size = 3.5) +
  labs(title = "ARID1A CUT&RUN peaks across activation (WT)",
       x = "Timepoint", y = "Number of peaks") +
  theme_paper

save_figure(p_arid1a_tp, "ARID1A_temporal_peak_counts", width = 5, height = 5, dir = figdir)

# ARID1A-specific categories
arid1a_naive_only <- arid1a_tp$Naive[countOverlaps(arid1a_tp$Naive, arid1a_tp$D5) == 0]
arid1a_d5_only    <- arid1a_tp$D5[countOverlaps(arid1a_tp$D5, arid1a_tp$Naive) == 0]
arid1a_conserved  <- arid1a_tp$D5[countOverlaps(arid1a_tp$D5, arid1a_tp$Naive) > 0]

export(arid1a_naive_only, file.path(outdir, "ARID1A_Naive_specific.bed"))
export(arid1a_d5_only,    file.path(outdir, "ARID1A_D5_gained.bed"))
export(arid1a_conserved,  file.path(outdir, "ARID1A_conserved.bed"))

message(sprintf("  ARID1A: Naive-specific=%d, D5-gained=%d, conserved=%d",
                length(arid1a_naive_only), length(arid1a_d5_only), length(arid1a_conserved)))

# =============================================================================
# 6. ARID1A-dependent vs independent TF/mark classification at D5
# =============================================================================

message("=== Classifying TF peaks by ARID1A dependency ===")

# ARID1A-dependent = WT-only peaks (lost in KO)
arid1a_dep   <- import(file.path(outdir, "ARID1A_D5_WT_only.bed"))
arid1a_indep <- import(file.path(outdir, "ARID1A_D5_shared.bed"))

tf_classification <- list()

for (ab in c("BATF", "ETS1", "Tbet", "H3K27ac")) {
  wt_id <- peak_samples %>%
    filter(antibody == ab, genotype == "WT", timepoint == "D5", is.na(treatment)) %>%
    pull(sample_id)
  if (length(wt_id) == 0) next

  tf_gr <- peaks_gr[[wt_id[1]]]
  if (length(tf_gr) == 0) next

  pct_dep   <- mean(countOverlaps(tf_gr, arid1a_dep)   > 0) * 100
  pct_indep <- mean(countOverlaps(tf_gr, arid1a_indep) > 0) * 100
  pct_none  <- 100 - pct_dep - pct_indep

  tf_classification[[ab]] <- tibble(
    antibody       = ab,
    pct_dep        = pct_dep,
    pct_indep      = pct_indep,
    pct_no_arid1a  = max(0, pct_none)
  )
}

tf_class_df <- bind_rows(tf_classification)
write_csv(tf_class_df, file.path(outdir, "TF_ARID1A_dependency.csv"))
message(tf_class_df)

# Stacked bar: % TF peaks at ARID1A-dep vs indep sites
tf_class_long <- tf_class_df %>%
  pivot_longer(-antibody, names_to = "category", values_to = "pct") %>%
  mutate(category = recode(category,
    pct_dep        = "ARID1A-dependent",
    pct_indep      = "ARID1A-independent",
    pct_no_arid1a  = "No ARID1A"
  ),
  category = factor(category, c("ARID1A-dependent", "ARID1A-independent", "No ARID1A")))

p_tf_dep <- ggplot(tf_class_long, aes(x = antibody, y = pct, fill = category)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(
    "ARID1A-dependent"   = "#D62728",
    "ARID1A-independent" = "#1F77B4",
    "No ARID1A"          = "#C0C0C0"
  )) +
  labs(title = "TF peaks at ARID1A-dependent vs independent sites (D5)",
       y = "% of TF peaks", x = "Factor", fill = "") +
  theme_paper

save_figure(p_tf_dep, "TF_ARID1A_dependency", width = 7, height = 5, dir = figdir)

# =============================================================================
# 7. Peak annotation: genomic context
# =============================================================================

message("=== Annotating peaks (genomic context) ===")

anno_sets <- list(
  "ARID1A_D5_WT"    = reduce(c(peaks_gr[["ARID1A_D5_WT_R1"]], peaks_gr[["ARID1A_D5_WT_R2"]])),
  "ARID1A_D5_KO"    = peaks_gr[["ARID1A_D5_KO_R1"]],
  "ARID1A_D5_WT_only"   = arid1a_dep,
  "ARID1A_D5_shared"    = arid1a_indep,
  "BATF_D5_WT"      = peaks_gr[["BATF_D5_WT_R1"]],
  "ETS1_D5_WT"      = peaks_gr[["ETS1_D5_WT_R1"]],
  "Tbet_D5_WT"      = peaks_gr[["Tbet_D5_WT_R1"]],
  "H3K27ac_D5_WT"   = peaks_gr[["H3K27ac_D5_WT_R1"]]
)

anno_list <- lapply(names(anno_sets), function(nm) {
  gr <- anno_sets[[nm]]
  if (length(gr) == 0) return(NULL)
  annotatePeak(gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
})
names(anno_list) <- names(anno_sets)
anno_list <- Filter(Negate(is.null), anno_list)

# Combined annotation bar plot
p_anno <- plotAnnoBar(anno_list) + theme_paper +
  ggtitle("Genomic annotation of CUT&RUN peaks") +
  theme(axis.text.y = element_text(size = 8))

save_figure(p_anno, "peak_annotation_bar", width = 9, height = 6, dir = figdir)

# Export annotated ARID1A D5 WT table
if (!is.null(anno_list[["ARID1A_D5_WT"]])) {
  arid1a_anno_df <- as.data.frame(anno_list[["ARID1A_D5_WT"]]) %>%
    as_tibble() %>%
    mutate(ensembl_clean = gsub("\\.\\d+$", "", geneId)) %>%
    left_join(
      tryCatch(
        AnnotationDbi::select(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", .$geneId),
                              keytype = "ENSEMBL", columns = "SYMBOL"),
        error = function(e) data.frame(ENSEMBL = character(), SYMBOL = character())
      ),
      by = c("ensembl_clean" = "ENSEMBL")
    )
  write_csv(arid1a_anno_df, file.path(outdir, "ARID1A_D5_WT_annotated.csv"))
}

# =============================================================================
# 8. H3K27ac at D8 subsets (WT only — chromatin state landscape)
# =============================================================================

message("=== H3K27ac at D8 subsets ===")

# Note: D8 samples don't have individual peak calls from the pipeline
# Use ARID1A D5 WT peak universe as anchor for comparison
# (deepTools heatmaps will be the primary viz — see core/06_deeptools_heatmaps.sh)

# Summarize what data exists
d8_bam_samples <- samples %>%
  filter(timepoint == "D8", bam_exists) %>%
  select(sample_id, antibody, genotype, cell_subset)

message("  D8 BAM files available:")
message(capture.output(print(d8_bam_samples)))

write_csv(d8_bam_samples, file.path(outdir, "D8_bam_manifest.csv"))

# =============================================================================
# Save session
# =============================================================================

message("=== Saving R session ===")
save.image(file.path(outdir, "cutandrun_analysis.RData"))
message("CUT&RUN analysis complete.")
message("  Results: ", outdir)
message("  Figures: ", figdir)
