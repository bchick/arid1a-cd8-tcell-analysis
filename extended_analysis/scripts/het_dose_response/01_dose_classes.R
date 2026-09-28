#!/usr/bin/env Rscript
# =============================================================================
# het_dose_response/01_dose_classes.R — ARID1A dose-response classes of ATAC peaks
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Design: D8 effector subsets (TE, EEC, MP), WT/Het/KO, Exp2-only
#   n = 3 WT, 3 Het, 2 KO per subset (balanced, single batch)
#
# Classifies each OCR by its response to ARID1A dosage, per subset (DESeq2):
#   - insensitive       : no change in Het or KO
#   - linear            : monotonic with dose, |Het lfc| < |KO lfc|
#   - haploinsufficient : Het ~= KO, both differ from WT  (one copy not enough)
#   - buffered          : Het ~= WT, only KO differs       (one copy sufficient)
#   - nonmonotonic      : Het and KO change in opposite directions
#
# Inputs:  results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/
#            consensus_peaks.mRp.clN.featureCounts.txt (nf-core/atacseq)
# Outputs: results/extended_analysis/het_dose_response/
#            dose_classes_<subset>.csv  (per-peak class + stats)
#            dose_classes_combined.csv  (all subsets, long form)
#            class_summary.csv          (counts per class x subset)
#            het_dose_response.RData    (DESeq2 objects, VST, results)
#            het_dose_response.log
#          figures/extended_analysis/het_dose_response/ (PCA, class counts, dose curves)
# Usage:   Rscript extended_analysis/scripts/het_dose_response/01_dose_classes.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
  library(ComplexHeatmap)
  library(circlize)
})
select <- dplyr::select
filter <- dplyr::filter
count  <- dplyr::count

outdir <- file.path(paths$ext_results, "het_dose_response")
figdir <- file.path(paths$ext_figures, "het_dose_response")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(outdir, "het_dose_response.log")
sink(logfile, split = TRUE)
on.exit(sink(), add = TRUE)

message("=== het_dose_response/01_dose_classes.R started: ", Sys.time(), " ===")

# =============================================================================
# 1. Load featureCounts and subset to Exp2 D8 main samples
# =============================================================================

message("\n=== Section 1: Load featureCounts and build metadata ===")

fc_file <- file.path(paths$atac,
  "bowtie2/merged_replicate/macs2/narrow_peak/consensus",
  "consensus_peaks.mRp.clN.featureCounts.txt")
stopifnot("featureCounts file not found" = file.exists(fc_file))

raw <- read.delim(fc_file, comment.char = "#", check.names = FALSE)
peak_info <- raw[, 1:6]
colnames(peak_info) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

count_mat <- as.matrix(raw[, -(1:6)])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- peak_info$Geneid
colnames(count_mat) <- gsub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(count_mat))
colnames(count_mat) <- gsub("REP(\\d)", "Rep\\1", colnames(count_mat))

message(sprintf("  Full count matrix: %d peaks x %d samples",
                nrow(count_mat), ncol(count_mat)))

# Restrict to D8 main effector subsets, Exp2 only, genotypes WT/Het/KO
keep_pat <- "^D8_(WT|Het|KO)_(TE|EEC|MP)_Exp2_Rep[0-9]+$"
keep_samples <- grep(keep_pat, colnames(count_mat), value = TRUE)
count_mat <- count_mat[, keep_samples]

# Parse metadata from sample names
parse_meta <- function(nm) {
  m <- regmatches(nm, regexec(
    "^D8_(WT|Het|KO)_(TE|EEC|MP)_(Exp[12])_Rep([0-9]+)$", nm))[[1]]
  data.frame(
    sample_name = nm,
    genotype    = m[2],
    cell_subset = m[3],
    experiment  = m[4],
    replicate   = paste0("Rep", m[5]),
    stringsAsFactors = FALSE
  )
}
meta <- do.call(rbind, lapply(keep_samples, parse_meta))
rownames(meta) <- meta$sample_name

# Ordered factor: WT(2) > Het(1) > KO(0) copies of ARID1A
meta$genotype    <- factor(meta$genotype, levels = c("WT", "Het", "KO"))
meta$cell_subset <- factor(meta$cell_subset, levels = c("TE", "EEC", "MP"))
meta$dose        <- c("WT" = 2L, "Het" = 1L, "KO" = 0L)[as.character(meta$genotype)]

message("  Sample design (Exp2 D8 only):")
print(table(meta$genotype, meta$cell_subset))

stopifnot(all(colnames(count_mat) == rownames(meta)))

# =============================================================================
# 2. Per-subset DESeq2: LRT + Wald contrasts
# =============================================================================

message("\n=== Section 2: Per-subset DESeq2 ===")

subsets <- levels(meta$cell_subset)
dds_list  <- list()
vsd_list  <- list()
res_tabs  <- list()   # holds long-form results per subset

for (sub in subsets) {
  message(sprintf("\n--- Subset: %s ---", sub))
  s_meta  <- meta[meta$cell_subset == sub, ]
  s_meta$genotype <- droplevels(s_meta$genotype)
  s_counts <- count_mat[, rownames(s_meta)]

  dds <- DESeqDataSetFromMatrix(countData = s_counts,
                                colData   = s_meta,
                                design    = ~ genotype)

  # Pre-filter: >= 10 counts in >= min(replicates) samples
  min_rep <- min(table(s_meta$genotype))
  keep <- rowSums(counts(dds) >= 10) >= min_rep
  dds <- dds[keep, ]
  message(sprintf("  Kept %d / %d peaks (>=10 counts in >=%d samples)",
                  sum(keep), length(keep), min_rep))

  # LRT: any genotype effect
  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
  lrt_res <- results(dds_lrt, independentFiltering = FALSE)

  # Wald contrasts for effect sizes and per-contrast significance
  dds_wald <- DESeq(dds, quiet = TRUE)
  het_res <- results(dds_wald, contrast = c("genotype", "Het", "WT"),
                     independentFiltering = FALSE)
  ko_res  <- results(dds_wald, contrast = c("genotype", "KO",  "WT"),
                     independentFiltering = FALSE)
  koVhet  <- results(dds_wald, contrast = c("genotype", "KO",  "Het"),
                     independentFiltering = FALSE)

  vsd <- vst(dds_wald, blind = FALSE)

  # Per-genotype mean VST (used for dose-curve ordering)
  vmat <- assay(vsd)
  mean_by <- function(g) rowMeans(vmat[, s_meta$genotype == g, drop = FALSE])
  wt_mean  <- mean_by("WT")
  het_mean <- mean_by("Het")
  ko_mean  <- mean_by("KO")

  tab <- data.frame(
    peak_id    = rownames(dds_wald),
    subset     = sub,
    baseMean   = lrt_res$baseMean,
    lrt_padj   = lrt_res$padj,
    lrt_stat   = lrt_res$stat,
    lfc_HetWT  = het_res$log2FoldChange,
    padj_HetWT = het_res$padj,
    lfc_KOWT   = ko_res$log2FoldChange,
    padj_KOWT  = ko_res$padj,
    lfc_KOHet  = koVhet$log2FoldChange,
    padj_KOHet = koVhet$padj,
    vst_WT     = wt_mean,
    vst_Het    = het_mean,
    vst_KO     = ko_mean,
    stringsAsFactors = FALSE
  )

  dds_list[[sub]] <- dds_wald
  vsd_list[[sub]] <- vsd
  res_tabs[[sub]] <- tab
}

# =============================================================================
# 3. Dose-response classification
# =============================================================================

message("\n=== Section 3: Classify dose response ===")

# Thresholds
ALPHA   <- 0.05    # BH-adj significance
LFC_MIN <- 0.5     # minimum |log2FC| to call a change
HET_KO_RATIO <- 0.5  # for linear: |Het lfc| must be between HET_KO_RATIO * |KO lfc| and (1 - ...)

classify <- function(df) {
  # Use LRT padj to gate "any response"; refine with pairwise Wald LFCs.
  responds <- !is.na(df$lrt_padj) & df$lrt_padj < ALPHA

  # Helpers
  sig_HetWT <- !is.na(df$padj_HetWT) & df$padj_HetWT < ALPHA
  sig_KOWT  <- !is.na(df$padj_KOWT)  & df$padj_KOWT  < ALPHA
  sig_KOHet <- !is.na(df$padj_KOHet) & df$padj_KOHet < ALPHA

  big_HetWT <- abs(df$lfc_HetWT) >= LFC_MIN
  big_KOWT  <- abs(df$lfc_KOWT)  >= LFC_MIN

  same_sign <- sign(df$lfc_HetWT) == sign(df$lfc_KOWT)

  cls <- rep("insensitive", nrow(df))

  # Only classify peaks that pass LRT
  # 1. nonmonotonic: Het and KO move in opposite directions, both with meaningful lfc
  is_nonmono <- responds & !same_sign & big_HetWT & big_KOWT &
                (sig_HetWT | sig_KOWT)
  cls[is_nonmono] <- "nonmonotonic"

  # 2. linear (dose-dependent): Het sig, KO sig, same sign, |Het lfc| < |KO lfc|,
  #    and KO vs Het itself significant (continues to drop with second allele lost)
  is_linear <- responds & same_sign &
               sig_HetWT & sig_KOWT & sig_KOHet &
               abs(df$lfc_HetWT) >= HET_KO_RATIO * abs(df$lfc_KOWT) &
               abs(df$lfc_HetWT) <  abs(df$lfc_KOWT)
  cls[is_linear] <- "linear"

  # 3. haploinsufficient: Het sig, KO sig, same sign, |Het lfc| ~ |KO lfc|,
  #    KO vs Het NOT significant (one-copy loss already produces full effect)
  is_haplo <- responds & same_sign &
              sig_HetWT & sig_KOWT & !sig_KOHet &
              big_HetWT & big_KOWT
  cls[is_haplo] <- "haploinsufficient"

  # 4. buffered / recessive: Het NOT sig (and small lfc), KO sig
  is_buffered <- responds & sig_KOWT & big_KOWT &
                 (!sig_HetWT) & abs(df$lfc_HetWT) < LFC_MIN
  cls[is_buffered] <- "buffered"

  # 5. other responders (sig LRT but don't fit clean classes)
  is_other <- responds & cls == "insensitive"
  cls[is_other] <- "other_responsive"

  # Direction (gained vs lost in KO)
  direction <- ifelse(df$lfc_KOWT > 0, "gained", "lost")
  direction[!(sig_KOWT | sig_HetWT)] <- NA_character_

  df$class     <- factor(cls, levels = c("insensitive", "buffered",
                                         "linear", "haploinsufficient",
                                         "nonmonotonic", "other_responsive"))
  df$direction <- direction
  df
}

res_tabs <- lapply(res_tabs, classify)

# Write per-subset tables
for (sub in names(res_tabs)) {
  out_csv <- file.path(outdir, sprintf("dose_classes_%s.csv", sub))
  write.csv(res_tabs[[sub]], out_csv, row.names = FALSE)
  message(sprintf("  Wrote %s (%d peaks)", out_csv, nrow(res_tabs[[sub]])))
}

combined <- bind_rows(res_tabs)
write.csv(combined, file.path(outdir, "dose_classes_combined.csv"), row.names = FALSE)

# Class summary table
class_summary <- combined %>%
  count(subset, class, direction) %>%
  arrange(subset, class, direction)
write.csv(class_summary, file.path(outdir, "class_summary.csv"), row.names = FALSE)

message("\n  Class counts (subset x class x direction):")
print(as.data.frame(class_summary))

# =============================================================================
# 4. QC: per-subset PCA
# =============================================================================

message("\n=== Section 4: QC plots ===")

pca_frames <- list()
for (sub in subsets) {
  vsd <- vsd_list[[sub]]
  vmat <- assay(vsd)
  rv <- rowVars(vmat)
  top <- order(rv, decreasing = TRUE)[seq_len(min(5000, length(rv)))]
  pr <- prcomp(t(vmat[top, ]), center = TRUE, scale. = FALSE)
  pct <- round(100 * summary(pr)$importance[2, 1:2], 1)
  pca_frames[[sub]] <- data.frame(
    PC1 = pr$x[, 1], PC2 = pr$x[, 2],
    sample = colnames(vsd),
    genotype = colData(vsd)$genotype,
    subset = sub,
    replicate = colData(vsd)$replicate,
    pc1_var = pct[1], pc2_var = pct[2]
  )
}
pca_df <- bind_rows(pca_frames)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = genotype, label = replicate)) +
  geom_point(size = 3) +
  facet_wrap(~ subset, scales = "free") +
  scale_color_manual(values = pal_genotype) +
  labs(title = "ATAC Het dose-response: PCA per subset (Exp2, D8)") +
  theme_paper
save_figure(p_pca, "pca_per_subset", width = 10, height = 4, dir = figdir)

# =============================================================================
# 5. Class summary plots
# =============================================================================

message("\n=== Section 5: Class summary plots ===")

# Stacked bar: proportion of responsive peaks in each class, per subset
class_bar_df <- combined %>%
  filter(class != "insensitive") %>%
  count(subset, class, direction) %>%
  mutate(class_dir = paste0(class, " (", direction, ")"))

p_bar <- ggplot(class_bar_df,
                aes(x = subset, y = n, fill = class)) +
  geom_col(position = "stack") +
  facet_wrap(~ direction) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Dose-response classes per subset",
       y = "# peaks", x = NULL, fill = "Class") +
  theme_paper
save_figure(p_bar, "class_counts_bar", width = 8, height = 4, dir = figdir)

# Dose curves: z-scored VST means per class, subset TE as example
plot_dose_curves <- function(sub) {
  df <- res_tabs[[sub]]
  curves <- df %>%
    filter(class %in% c("linear", "haploinsufficient",
                        "buffered", "nonmonotonic")) %>%
    mutate(
      z_WT  = (vst_WT  - (vst_WT + vst_Het + vst_KO)/3),
      z_Het = (vst_Het - (vst_WT + vst_Het + vst_KO)/3),
      z_KO  = (vst_KO  - (vst_WT + vst_Het + vst_KO)/3)
    ) %>%
    select(peak_id, class, direction, z_WT, z_Het, z_KO) %>%
    pivot_longer(z_WT:z_KO, names_to = "genotype", values_to = "z") %>%
    mutate(genotype = factor(sub("^z_", "", genotype),
                             levels = c("WT", "Het", "KO")))

  ggplot(curves, aes(x = genotype, y = z, group = peak_id)) +
    geom_line(alpha = 0.03, color = "grey40") +
    stat_summary(aes(group = direction, color = direction),
                 fun = median, geom = "line", linewidth = 1.2) +
    facet_grid(direction ~ class) +
    labs(title = sprintf("Dose-response curves — %s", sub),
         y = "VST (centered)") +
    theme_paper +
    theme(legend.position = "none")
}

for (sub in subsets) {
  p <- plot_dose_curves(sub)
  save_figure(p, sprintf("dose_curves_%s", sub),
              width = 9, height = 5, dir = figdir)
}

# =============================================================================
# 6. Save RData checkpoint
# =============================================================================

message("\n=== Section 6: Save checkpoint ===")

save(dds_list, vsd_list, res_tabs, combined, class_summary, meta,
     file = file.path(outdir, "het_dose_response.RData"))

message("\n=== Done: ", Sys.time(), " ===")
