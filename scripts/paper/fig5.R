#!/usr/bin/env Rscript
# =============================================================================
# paper/fig5.R — Figure 5 + Figure S5A: cBAF is required for targeting of T-bet
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# Defines every OCR set used in Figure 5, draws the count / annotation panels,
# and writes BED files for the signal panels (scripts/paper/fig4_fig5_signal.sh)
# and the HOMER motif runs (scripts/paper/fig5a_homer.sh -> fig5a_motifs.R).
#
#   5A  BED inputs: OCRs lost / gained in Arid1a KO at d5 and in d8 TE/EEC/MP
#   5B  NOT reproduced: needs public in vitro TF ChIP-seq (GSE192390)
#   5C  ARID1A-dependent vs -independent OCRs at d5 (annotation pies + BEDs)
#   5D-F ACBI1-dependent vs -independent T-bet-bound OCRs in vitro (pies + BEDs)
#   5G  Overlap of OCRs lost in Arid1a KO vs Tbx21 KO (TE, MP)
#   5H  BEDs of the 5G groups for the ARID1A CUT&RUN WT vs Tbx21 KO heatmap
#   5I  wet-lab (flow cytometry)
#   5J  uses the published Fig 1A clusters directly (see fig4_fig5_signal.sh)
#   S5A OCRs lost / gained in Tbx21 KO vs WT d8 subsets (2-fold, padj < 0.01)
#
# Set definitions (see PANEL_MAP.md):
#   ARID1A-bound OCR     consensus OCR overlapping an ARID1A CUT&RUN peak in
#                        d5 WT (MACS2, either replicate)
#   ARID1A-dependent     ARID1A-bound OCR lost in d5 Arid1a KO (2-fold, padj<0.05)
#   ARID1A-independent   ARID1A-bound OCR tested at d5 and not lost
#   T-bet-bound OCR      consensus OCR overlapping a T-bet ChIP peak (DMSO + IL-12)
#   ACBI1-dependent      T-bet-bound OCR lost in ACBI1 + IL-12 vs DMSO + IL-12
#                        (DESeq2 on the IL-12 / ACBI1 / BRM014 libraries,
#                        2-fold, padj < 0.05)
#   ACBI1-independent    T-bet-bound OCR tested and not ACBI1-dependent
#
# 5D-F normalization. The inhibitor libraries have lower FRiP (reads in
# consensus OCRs / all counted reads) than DMSO + IL-12, so DESeq2's
# median-of-ratios size factors absorb a global loss of accessibility. The paper
# used HOMER, which normalizes to total tags, so the PRIMARY split here uses
# total-read size factors (all reads counted by featureCounts); median-of-ratios
# is written alongside as a SENSITIVITY analysis (fig5d_normalization_summary.csv,
# fig5_norm_sensitivity.*). Caveats:
#   - the paper's 7,487 ACBI1-dependent OCRs are not reached under either scheme;
#   - the global FRiP drop is consistent with BAF disruption but is confounded
#     with library quality (FRiP ~0.1 overall; BRM014 REP1 0.088);
#   - two replicates per condition.
# The paper's two sets sum to 18,502, which matches the ~23.7k T-bet-bound tested
# OCRs here rather than all ~72.5k tested OCRs, so the T-bet-bound denominator
# is used.
#
# Inputs:  results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/ (featureCounts)
#          results/atac/differential/{da_*.csv, consensus_peaks_annotated.csv}
#          results/cutrun/03_peak_calling/04_called_peaks/macs2/ (ARID1A d5 WT)
#          results/chipseq/bowtie2/merged_library/macs3/narrow_peak/ (T-bet)
# Outputs: results/paper/fig5_beds/*.bed (for fig4_fig5_signal.sh, fig5a_homer.sh)
#          figures/paper/fig5*_*, figS5a_*.{pdf,png}; results/paper/fig5*_*.csv
# Usage:   Rscript scripts/paper/fig5.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggforce)
})
select <- dplyr::select
filter <- dplyr::filter

consensus_dir <- file.path(paths$atac, "bowtie2/merged_replicate/macs2/narrow_peak/consensus")
cutrun_peaks  <- file.path(paths$cutrun, "03_peak_calling/04_called_peaks/macs2")
chip_peaks    <- file.path(paths$chipseq, "bowtie2/merged_library/macs3/narrow_peak")
bed_dir       <- file.path(paths$paper_tab, "fig5_beds")
dir.create(bed_dir, showWarnings = FALSE)

saf <- read_tsv(file.path(consensus_dir, "consensus_peaks.mRp.clN.saf"), show_col_types = FALSE)
names(saf)[1:4] <- c("peak_id", "chr", "start", "end")
ocr_gr <- GRanges(saf$chr, IRanges(saf$start, saf$end), peak_id = saf$peak_id)

read_narrowpeak <- function(f) {
  np <- read_tsv(f, col_names = FALSE, col_types = "cii-------", progress = FALSE)
  GRanges(np$X1, IRanges(np$X2 + 1L, np$X3))
}

#' OCR ids overlapping any of the given peak sets
bound_ocrs <- function(peak_files) {
  pk <- reduce(do.call(c, unname(lapply(peak_files, read_narrowpeak))))
  ocr_gr$peak_id[overlapsAny(ocr_gr, pk)]
}

#' Write OCR ids as a BED (0-based) for deepTools / HOMER
write_ocr_bed <- function(ids, name) {
  f <- file.path(bed_dir, paste0(name, ".bed"))
  saf |> filter(peak_id %in% ids) |>
    transmute(chr, start = start - 1L, end, peak_id, score = 0, strand = "+") |>
    arrange(chr, start) |>
    write_tsv(f, col_names = FALSE)
  message(sprintf("  %-40s %6d OCRs", basename(f), length(unique(ids))))
  invisible(f)
}

# Genomic annotation, collapsed to the paper's five categories
anno <- read_csv(file.path(paths$atac, "differential/consensus_peaks_annotated.csv"),
                 show_col_types = FALSE) |>
  transmute(peak_id,
            category = case_when(grepl("^Promoter", annotation) ~ "Promoter",
                                 grepl("Intergenic", annotation) ~ "Intergenic",
                                 grepl("^Intron", annotation)    ~ "Intron",
                                 grepl("^Exon", annotation)      ~ "Exon",
                                 TRUE                            ~ "Other"),
            category = factor(category, levels = names(pal_annotation)))

annotation_pies <- function(sets, id) {
  tab <- purrr::imap(sets, ~ tibble(set = .y, peak_id = .x)) |> bind_rows() |>
    inner_join(anno, by = "peak_id") |>
    count(set, category, .drop = FALSE) |>
    group_by(set) |> arrange(category, .by_group = TRUE) |>
    # ggplot stacks the first factor level on top: label at each slice's midpoint
    mutate(pct = 100 * n / sum(n), total = sum(n), ypos = 100 - (cumsum(pct) - pct / 2)) |>
    ungroup() |>
    mutate(set = factor(set, levels = names(sets)))
  write_panel_table(tab, id)
  p <- ggplot(tab, aes(x = 1, y = pct, fill = category)) +
    geom_col(width = 1, color = "white", linewidth = 0.3) +
    geom_text(data = ~ filter(.x, pct >= 5),
              aes(y = ypos, label = sprintf("%.1f%%", pct),
                  color = ifelse(category %in% c("Promoter", "Intergenic"), "white", "black")),
              size = 2.8) +
    scale_color_identity() +
    coord_polar(theta = "y") +
    facet_wrap(~ set, ncol = 1,
               labeller = as_labeller(setNames(sprintf("%s\n(n=%s)", levels(tab$set),
                 format(tapply(tab$n, tab$set, sum), big.mark = ",", trim = TRUE)), levels(tab$set)))) +
    scale_fill_manual(values = pal_annotation, drop = FALSE) +
    labs(fill = NULL) +
    theme_void() + theme(strip.text = element_text(size = 9))
  save_panel(p, id, width = 2.8, height = 4.6)
}

# =============================================================================
# 5A. OCR sets for HOMER known-motif enrichment
# =============================================================================

message("=== 5A: HOMER input sets ===")

contrasts_5a <- c(d5 = "KO_vs_WT_D5", TE = "KO_vs_WT_D8_TE",
                  EEC = "KO_vs_WT_D8_EEC", MP = "KO_vs_WT_D8_MP")
da5 <- purrr::map(contrasts_5a, read_da)
for (s in names(da5)) for (dir in c("lost", "gained"))
  write_ocr_bed(da5[[s]]$peak_id[da5[[s]]$direction == dir], sprintf("fig5a_%s_%s", s, dir))

# =============================================================================
# 5C. ARID1A-dependent vs -independent OCRs (day 5)
# =============================================================================

message("=== 5C: ARID1A-dependent OCRs ===")

arid1a_bound <- bound_ocrs(file.path(cutrun_peaks,
  c("ARID1A_D5_WT_R1.macs2_peaks.narrowPeak", "ARID1A_D5_WT_R2.macs2_peaks.narrowPeak")))
d5 <- da5$d5
sets_5c <- list(
  "ARID1A-dependent OCRs"   = intersect(arid1a_bound, d5$peak_id[d5$direction == "lost"]),
  "ARID1A-independent OCRs" = intersect(arid1a_bound, d5$peak_id[d5$direction != "lost"])
)
write_ocr_bed(sets_5c[[1]], "fig5c_arid1a_dependent")
write_ocr_bed(sets_5c[[2]], "fig5c_arid1a_independent")
annotation_pies(sets_5c, "fig5c_annotation_pies")

# =============================================================================
# 5D-F. ACBI1-dependent vs -independent T-bet-bound OCRs (in vitro)
# =============================================================================

message("=== 5D-F: ACBI1-dependent OCRs ===")

fc <- read.delim(file.path(consensus_dir, "consensus_peaks.mRp.clN.featureCounts.txt"),
                 comment.char = "#", check.names = FALSE)
cnt <- as.matrix(fc[, -(1:6)])
rownames(cnt) <- fc$Geneid
colnames(cnt) <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(cnt))
fc_sum <- read.delim(file.path(consensus_dir, "consensus_peaks.mRp.clN.featureCounts.txt.summary"),
                     check.names = FALSE)
total_reads <- colSums(fc_sum[, -1])
names(total_reads) <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", names(total_reads))

inh_libs <- grep("^ATAC_IL-12(_ACBI1|_BRM014)?_REP", colnames(cnt), value = TRUE)
inh_cd <- data.frame(
  treatment = factor(ifelse(grepl("ACBI1", inh_libs), "ACBI1",
                            ifelse(grepl("BRM014", inh_libs), "BRM014", "IL12")),
                     levels = c("IL12", "ACBI1", "BRM014")),
  row.names = inh_libs)
inh_x <- cnt[, inh_libs]
inh_x <- inh_x[rowSums(inh_x >= 10) >= 2, ]

libs_tab <- tibble(library = inh_libs, treatment = as.character(inh_cd$treatment),
                   total_reads = total_reads[inh_libs],
                   reads_in_peaks = colSums(cnt[, inh_libs]),
                   frip = reads_in_peaks / total_reads)

inhibitor_da <- function(norm) {
  dds <- DESeq2::DESeqDataSetFromMatrix(inh_x, inh_cd, design = ~ treatment)
  if (norm == "totalreads") {
    tr <- total_reads[inh_libs]
    DESeq2::sizeFactors(dds) <- tr / exp(mean(log(tr)))
  } else {
    dds <- DESeq2::estimateSizeFactors(dds)
  }
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  purrr::map_dfr(c("ACBI1", "BRM014"), function(drug)
    DESeq2::results(dds, contrast = c("treatment", drug, "IL12")) |>
      as.data.frame() |> tibble::rownames_to_column("peak_id") |>
      mutate(drug = drug, normalization = norm,
             size_factor_set = paste(round(DESeq2::sizeFactors(dds), 3), collapse = ";")))
}

tbet_bound <- bound_ocrs(file.path(chip_peaks, "ChIP_IL-12_T-bet_REP1_peaks.narrowPeak"))
tbet_tested <- intersect(tbet_bound, rownames(inh_x))

acbi1_sets <- list(); norm_summary <- list()
for (norm in c("totalreads", "medianratio")) {
  res <- inhibitor_da(norm) |>
    mutate(direction = case_when(
      !is.na(padj) & padj < PAPER$atac_padj & log2FoldChange <= -PAPER$atac_lfc ~ "lost",
      !is.na(padj) & padj < PAPER$atac_padj & log2FoldChange >=  PAPER$atac_lfc ~ "gained",
      TRUE ~ "ns"))
  norm_summary[[norm]] <- res |>
    group_by(normalization, drug) |>
    summarise(lost = sum(direction == "lost"), gained = sum(direction == "gained"),
              lost_tbet_bound = sum(direction == "lost" & peak_id %in% tbet_tested),
              median_log2fc = median(log2FoldChange, na.rm = TRUE), .groups = "drop")
  lost_acbi1 <- res$peak_id[res$drug == "ACBI1" & res$direction == "lost"]
  acbi1_sets[[norm]] <- tibble(peak_id = tbet_tested,
                               set = ifelse(tbet_tested %in% lost_acbi1,
                                            "ACBI1-dependent", "ACBI1-independent"),
                               normalization = norm)
  write_panel_table(acbi1_sets[[norm]], paste0("fig5d_acbi1_sets_", norm))
}
norm_summary <- bind_rows(norm_summary)
write_panel_table(bind_rows(
  norm_summary |> mutate(record = "differential"),
  libs_tab |> mutate(record = "library")), "fig5d_normalization_summary")
print(norm_summary); print(libs_tab)

# Panels 5D-F use the primary (total-read) split
prim <- acbi1_sets$totalreads
sets_5f <- list(
  "ACBI1-dependent"   = prim$peak_id[prim$set == "ACBI1-dependent"],
  "ACBI1-independent" = prim$peak_id[prim$set == "ACBI1-independent"]
)
write_ocr_bed(sets_5f[[1]], "fig5d_acbi1_dependent")
write_ocr_bed(sets_5f[[2]], "fig5d_acbi1_independent")
annotation_pies(sets_5f, "fig5f_annotation_pies")

# Supplementary: sensitivity of the lost-OCR count to normalization, with FRiP
p_ns1 <- norm_summary |>
  mutate(normalization = recode(normalization, totalreads = "Total reads (primary)",
                                medianratio = "Median-of-ratios")) |>
  ggplot(aes(drug, lost, fill = normalization)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = lost), position = position_dodge(0.8), vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = c("Total reads (primary)" = "grey20", "Median-of-ratios" = "grey65")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "OCRs lost vs DMSO + IL-12", fill = NULL) +
  theme(legend.position = "top")
p_ns2 <- libs_tab |>
  mutate(library = sub("^ATAC_", "", library),
         treatment = factor(treatment, levels = c("IL12", "ACBI1", "BRM014"))) |>
  ggplot(aes(library, frip, fill = treatment)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", frip)), vjust = -0.3, size = 2.4) +
  scale_fill_manual(values = c(IL12 = unname(pal_treatment["DMSO_IL12"]),
                               ACBI1 = unname(pal_treatment["ACBI1_IL12"]),
                               BRM014 = unname(pal_treatment["BRM014_IL12"])), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "FRiP (consensus OCRs)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(patchwork::wrap_plots(p_ns1, p_ns2, nrow = 1, widths = c(1, 1.4)),
           "fig5_norm_sensitivity", width = 6.5, height = 3.2)

# =============================================================================
# 5G / 5H. OCRs lost in Arid1a KO vs Tbx21 KO (TE, MP)
# =============================================================================

message("=== 5G/5H: Arid1a KO vs Tbx21 KO ===")

#' Two-circle area-proportional Euler diagram (circle areas = set sizes)
euler2 <- function(a_only, both, b_only) {
  ra <- sqrt((a_only + both) / pi); rb <- sqrt((b_only + both) / pi)
  lens <- function(d) {                       # intersection area of two circles
    if (d >= ra + rb) return(0)
    if (d <= abs(ra - rb)) return(pi * min(ra, rb)^2)
    ra^2 * acos((d^2 + ra^2 - rb^2) / (2 * d * ra)) +
      rb^2 * acos((d^2 + rb^2 - ra^2) / (2 * d * rb)) -
      0.5 * sqrt((-d + ra + rb) * (d + ra - rb) * (d - ra + rb) * (d + ra + rb))
  }
  d <- if (both >= pi * min(ra, rb)^2) abs(ra - rb) else
    uniroot(function(d) lens(d) - both, c(abs(ra - rb) + 1e-9, ra + rb))$root
  tibble(set = c("Arid1a", "Tbx21"), x0 = c(0, d), y0 = 0, r = c(ra, rb))
}

g_counts <- list(); g_plots <- list()
for (s in c("TE", "MP")) {
  ko   <- read_da(paste0("KO_vs_WT_D8_", s))
  tbko <- read_da(paste0("TbetKO_vs_WT_D8_", s))
  lost_ko <- ko$peak_id[ko$direction == "lost"]
  lost_tb <- tbko$peak_id[tbko$direction == "lost"]
  grp <- list(ko_only = setdiff(lost_ko, lost_tb), both = intersect(lost_ko, lost_tb),
              tbx21_only = setdiff(lost_tb, lost_ko))
  write_ocr_bed(grp$both,       sprintf("fig5h_%s_lost_both", s))
  write_ocr_bed(grp$tbx21_only, sprintf("fig5h_%s_lost_tbx21ko_only", s))
  write_ocr_bed(grp$ko_only,    sprintf("fig5h_%s_lost_arid1ako_only", s))
  g_counts[[s]] <- tibble(subset = s, lost_arid1a_ko = length(lost_ko),
                          lost_tbx21_ko = length(lost_tb), lost_both = length(grp$both),
                          arid1a_only = length(grp$ko_only), tbx21_only = length(grp$tbx21_only))
  circ <- euler2(length(grp$ko_only), length(grp$both), length(grp$tbx21_only))
  g_plots[[s]] <- ggplot(circ) +
    geom_circle(aes(x0 = x0, y0 = y0, r = r, fill = set), color = "black", alpha = 0.75) +
    annotate("text", x = circ$x0[1] - 0.45 * circ$r[1], y = 0,
             label = sprintf("Lost in\nArid1a cKO\n(n=%d)", length(lost_ko)), size = 2.8) +
    annotate("text", x = circ$x0[2], y = 0,
             label = sprintf("Both\n(n=%d)", length(grp$both)), size = 2.6) +
    annotate("text", x = circ$x0[2] + circ$r[2], y = -circ$r[1] * 0.9, hjust = 0,
             label = sprintf("Lost in Tbx21 KO\n(n=%d)", length(lost_tb)), size = 2.8) +
    scale_fill_manual(values = c(Arid1a = unname(pal_genotype["KO"]), Tbx21 = "#A6BDDB"), guide = "none") +
    coord_fixed(clip = "off") + labs(title = s) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", color = pal_subset[[s]]),
          plot.margin = margin(5, 60, 5, 5))
}
write_panel_table(bind_rows(g_counts), "fig5g_arid1a_vs_tbx21_lost")
print(bind_rows(g_counts))
save_panel(patchwork::wrap_plots(g_plots, nrow = 1), "fig5g_arid1a_vs_tbx21_venn", width = 6.5, height = 3)

# =============================================================================
# S5A. Tbx21 KO vs WT differential accessibility per d8 subset
# =============================================================================

message("=== S5A: Tbx21 KO DA counts ===")

s5a <- purrr::map(c("TE", "EEC", "MP"), function(s) {
  d <- read_da(paste0("TbetKO_vs_WT_D8_", s), lfc = PAPER$tbet_lfc, padj = PAPER$tbet_padj)
  tibble(subset = s, lost = sum(d$direction == "lost"), gained = sum(d$direction == "gained"))
}) |> bind_rows()
write_panel_table(s5a, "figS5a_tbx21ko_da_counts")
print(s5a)

ps5a <- s5a |>
  pivot_longer(c(lost, gained), names_to = "direction", values_to = "n") |>
  mutate(subset = factor(subset, levels = c("TE", "EEC", "MP")),
         direction = factor(direction, levels = c("lost", "gained"),
                            labels = c("Lost in Tbx21 KO", "Gained in Tbx21 KO"))) |>
  ggplot(aes(subset, n, fill = direction)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("grey30", "grey70")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "# differentially accessible OCRs", fill = NULL,
       title = "Tbx21 KO vs WT, day 8") +
  theme(legend.position = "top")
save_panel(ps5a, "figS5a_tbx21ko_da_counts", width = 3.6, height = 3.2)

message("=== fig5.R done ===")
