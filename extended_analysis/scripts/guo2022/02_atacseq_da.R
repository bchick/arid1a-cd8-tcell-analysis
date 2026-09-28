#!/usr/bin/env Rscript
# =============================================================================
# guo2022/02_atacseq_da.R — Guo et al. 2022 ATAC-seq differential accessibility
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Cross-study comparison with Guo et al. 2022 (Nature): DESeq2 (apeglm
# shrinkage), ChIPseeker, clusterProfiler, ComplexHeatmap on the nf-core/atacseq
# reprocessing of Guo's ATAC-seq (mm39 / GENCODE vM35): merged-replicate
# consensus featureCounts, 182,052 peaks x 30 samples.
#
# Mirrors core/02_atacseq_analysis.R (§1-5,7-9), adapted for Guo:
#   - No master sample sheet: metadata parsed from Guo_<gse>_<cond>_REP<N> BAM
#     column names. run_da() reused ~verbatim (it already has a group_col param).
#   - Contrasts are WITHIN-GSE only (GSE confounded with condition); the single
#     shared 182k-peak consensus across all Guo groups is what makes the
#     cross-arm LFC comparison valid downstream.
#   - GSE183618 bundles TWO experiments under one "WT" GEO label. Its six WT
#     samples split into the ARID1A-arm controls and the MYC-arm controls, and
#     the paper (Guo ED Fig. 6a) compares 3 vs 3 within an arm — never pooled.
#     See the relabelling block in Section 1 for the evidence and mapping.
#   - No batch_col: after that relabelling every contrast is again clean
#     replicates of two conditions from one arm (unlike McDonald D8).
#   - No OCR-dynamics timecourse (core/02_atacseq_analysis.R §6) — Guo has no
#     WT Naive->D8 series.
#   - ChIPseeker GENCODE version-strip fix reused verbatim from
#     core/02_atacseq_analysis.R.
#
# Key test: does ARID1A loss in Guo replicate the ETS dose-buffered enhancer
# program? The KO_vs_WT / vehKO_vs_vehWT / inhib_vs_vehWT DA tables carry it,
# and feed the cross-study fold-change correlation.
#
# Inputs:  results/meta_analysis/guo2022/atacseq/bowtie2/merged_replicate/macs2/
#            narrow_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt
# Outputs: results/extended_analysis/guo2022/atacseq/differential/
#            da_<contrast>.csv, da_summary.csv, consensus_peaks_annotated.csv,
#            go_da_<dir>_<contrast>.csv, checkpoint_sections1to3.RData,
#            guo_atacseq_analysis.RData
#          figures/extended_analysis/guo2022/atacseq/
# Usage:   Rscript extended_analysis/scripts/guo2022/02_atacseq_da.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(ChIPseeker)
  library(rtracklayer)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(viridisLite)
  library(AnnotationDbi)
  library(matrixStats)
  library(txdbmaker)
})

# Resolve namespace conflicts: Bioconductor packages mask dplyr verbs
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate
slice  <- dplyr::slice

# Output directories
outdir <- file.path(paths$ext_results, "guo2022/atacseq/differential")
figdir <- file.path(paths$ext_figures, "guo2022/atacseq")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# Guo-specific palettes (conditions span three arms; no genotype/timepoint schema)
pal_guo_cond <- c(
  "WT"    = "#000000",  # g183618 ARID1A-arm control
  "KO"    = "#2CA02C",  # g183618 Arid1a knockout
  "MycWT" = "#3B3B3B",  # g183618 MYC-arm control (separate experiment/batch)
  "MycKO" = "#1F77B4",  # g183618 c-Myc knockout
  "vehWT" = "#666666",  # g198894 vehicle WT
  "vehKO" = "#8C564B",  # g198894 vehicle KO (independent genetic replication)
  "inhib" = "#D62728",  # g198894 cBAF inhibitor
  "MycLo" = "#9EDAE5",  # g183616 Myc-low
  "MycHi" = "#FF7F0E",  # g183616 Myc-high
  "Naive" = "#7F7F7F"   # g183616 naive
)
pal_guo_gse <- c("g183618" = "#4A90D9", "g198894" = "#DDCC77", "g183616" = "#AA4499")

# Within-GSE, unconfounded contrasts (num vs denom). KO_vs_WT is the primary
# ARID1A/ETS-replication test on the ATAC side.
contrast_specs <- list(
  list(name = "KO_vs_WT",       gse = "g183618", num = "KO",    denom = "WT"),     # 3 vs 3 — PRIMARY
  list(name = "MycKO_vs_MycWT", gse = "g183618", num = "MycKO", denom = "MycWT"),  # 3 vs 3 — c-Myc
  list(name = "vehKO_vs_vehWT", gse = "g198894", num = "vehKO", denom = "vehWT"),  # 3 vs 3 — independent genetic
  list(name = "inhib_vs_vehWT", gse = "g198894", num = "inhib", denom = "vehWT"),  # 3 vs 3 — pharmacologic
  list(name = "MycHi_vs_MycLo", gse = "g183616", num = "MycHi", denom = "MycLo"),  # 3 vs 3 — Myc dose
  list(name = "MycHi_vs_Naive", gse = "g183616", num = "MycHi", denom = "Naive")   # 3 vs 3 — activation baseline
)

# Checkpoint from a previous run (skip Sections 1-3 if available)
checkpoint_file <- file.path(outdir, "checkpoint_sections1to3.RData")
if (file.exists(checkpoint_file)) {
  message("=== Loading checkpoint (Sections 1-3 already complete) ===")
  load(checkpoint_file)
  message(sprintf("  Loaded: %d DA results, %d peaks, %d samples",
                  length(da_results), nrow(peak_info_filt), ncol(vsd)))
} else {

# =============================================================================
# 1. Import featureCounts consensus peak matrix + parse metadata
# =============================================================================

message("=== Section 1: Importing featureCounts consensus peak matrix ===")

fc_file <- file.path(paths$meta,
  "atacseq/bowtie2/merged_replicate/macs2/narrow_peak/consensus",
  "consensus_peaks.mRp.clN.featureCounts.txt")
stopifnot("featureCounts file not found" = file.exists(fc_file))

raw <- read.delim(fc_file, comment.char = "#", check.names = FALSE)
message(sprintf("  Loaded %d peaks x %d columns", nrow(raw), ncol(raw)))

# Peak coordinates
peak_info <- raw[, 1:6]
colnames(peak_info) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

# Count matrix (columns 7+)
count_mat <- as.matrix(raw[, -(1:6)])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- peak_info$Geneid

# Clean column names: strip the BAM suffix
clean_names <- gsub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(count_mat))
colnames(count_mat) <- clean_names
message(sprintf("  Count matrix: %d peaks x %d samples", nrow(count_mat), ncol(count_mat)))

# --- Parse metadata from Guo_<gse>_<cond>_REP<N> column names ---
parsed <- stringr::str_match(clean_names, "^Guo_(g\\d+)_(.+)_REP(\\d+)$")
stopifnot("Unparseable Guo ATAC sample name(s)" = !any(is.na(parsed[, 1])))

meta <- data.frame(
  sample    = clean_names,
  gse       = parsed[, 2],
  condition = parsed[, 3],
  replicate = parsed[, 4],
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample

# --- Split the g183618 "WT" label into its two source experiments -------------
# GSE183618 bundles two experiments that share a single "WT" GEO condition
# label: an ARID1A arm (Arid1a KO vs WT) and a MYC arm (Myc KO vs WT). Guo
# compares 3 vs 3 *within* an arm (ED Fig. 6a, "n = 3 biological replicates");
# pooling all six WT mixes two batches and inflates dispersion.
#
# nf-core replicate numbers are assigned by samplesheet row order and do NOT
# follow GEO's rep numbering, so the arms must be assigned explicitly. Mapping
# verified two independent ways:
#   (1) samplesheet SRR -> GEO GSM (GSE183618 filelist: GSM5563555/56/57 =
#       WT.rep1/2/3 sequenced with the Arid1a KOs; GSM5563561/62/63 =
#       WT.rep4/5/6 sequenced with the Myc KOs), and
#   (2) a reads-per-GB batch signature that separates the two arms with a clean
#       2x gap (ARID1A arm 5.4-7.5 M/GB; MYC arm 15.4-16.8 M/GB) and groups each
#       WT with its own KOs.
#     WT_REP1 = GSM5563557   WT_REP5 = GSM5563555   WT_REP6 = GSM5563556  -> ARID1A arm
#     WT_REP2 = GSM5563561   WT_REP3 = GSM5563562   WT_REP4 = GSM5563563  -> MYC arm
myc_arm_wt <- paste0("Guo_g183618_WT_REP", c(2, 3, 4))
stopifnot(
  "Expected g183618 WT samples not found — check nf-core replicate numbering" =
    all(myc_arm_wt %in% meta$sample),
  "Expected 6 WT samples in g183618" =
    sum(meta$gse == "g183618" & meta$condition == "WT") == 6
)
meta$condition[meta$sample %in% myc_arm_wt] <- "MycWT"
message(sprintf("  Relabelled %d g183618 WT samples as MycWT (MYC-arm controls): %s",
                length(myc_arm_wt), paste(myc_arm_wt, collapse = ", ")))

meta$gse       <- factor(meta$gse)
meta$condition <- factor(meta$condition)
meta$replicate <- factor(meta$replicate)

message(sprintf("  Metadata: %d samples", nrow(meta)))
message("  GSE x condition breakdown:")
print(table(meta$gse, meta$condition))

# DESeqDataSet with minimal design (subset per contrast)
dds_full <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ 1
)

# Pre-filter: peaks with >= 10 counts in >= 3 samples
keep <- rowSums(counts(dds_full) >= 10) >= 3
dds_full <- dds_full[keep, ]
message(sprintf("  Kept %d / %d peaks after pre-filtering", sum(keep), length(keep)))

# Update peak_info to match filtered peaks
peak_info_filt <- peak_info[peak_info$Geneid %in% rownames(dds_full), ]

# VST for visualization
vsd <- vst(dds_full, blind = FALSE)
message("  VST transformation complete")

# =============================================================================
# 2. QC Visualization (PCA + sample-correlation heatmap)
# =============================================================================

message("=== Section 2: QC Visualization ===")

# --- PCA all samples ---
pca_mat <- assay(vsd)
rv <- rowVars(pca_mat)
top_var <- order(rv, decreasing = TRUE)[seq_len(min(5000, length(rv)))]
pca_res <- prcomp(t(pca_mat[top_var, ]), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  condition = colData(vsd)$condition,
  gse       = colData(vsd)$gse,
  sample    = colnames(vsd)
)
pct_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

p_pca_all <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = gse)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = pal_guo_cond) +
  labs(
    title = "Guo ATAC-seq PCA: All Samples",
    x = sprintf("PC1 (%.1f%%)", pct_var[1]),
    y = sprintf("PC2 (%.1f%%)", pct_var[2]),
    color = "Condition", shape = "GSE"
  ) +
  theme_paper
save_figure(p_pca_all, "guo_atacseq_pca_all", width = 8, height = 6, dir = figdir)

# --- PCA g183618 (genetic Arid1a arm) only ---
gen_vsd <- vsd[, colData(vsd)$gse == "g183618"]
if (ncol(gen_vsd) > 5) {
  pca_mat_g <- assay(gen_vsd)
  rv_g <- rowVars(pca_mat_g)
  top_g <- order(rv_g, decreasing = TRUE)[seq_len(min(5000, length(rv_g)))]
  pca_g <- prcomp(t(pca_mat_g[top_g, ]), center = TRUE, scale. = FALSE)
  pca_g_df <- data.frame(
    PC1       = pca_g$x[, 1],
    PC2       = pca_g$x[, 2],
    condition = colData(gen_vsd)$condition,
    sample    = colnames(gen_vsd)
  )
  pct_g <- round(100 * summary(pca_g)$importance[2, 1:2], 1)

  p_pca_g <- ggplot(pca_g_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 4, alpha = 0.85) +
    scale_color_manual(values = pal_guo_cond) +
    labs(
      title = "Guo ATAC-seq PCA: g183618 (Arid1a genetic)",
      x = sprintf("PC1 (%.1f%%)", pct_g[1]),
      y = sprintf("PC2 (%.1f%%)", pct_g[2]),
      color = "Condition"
    ) +
    theme_paper
  save_figure(p_pca_g, "guo_atacseq_pca_g183618", width = 8, height = 6, dir = figdir)
}

# --- Sample-sample correlation heatmap ---
cor_mat <- cor(assay(vsd), method = "pearson")

ha_col <- HeatmapAnnotation(
  Condition = colData(vsd)$condition,
  GSE       = colData(vsd)$gse,
  col = list(
    Condition = pal_guo_cond[levels(droplevels(colData(vsd)$condition))],
    GSE       = pal_guo_gse[levels(droplevels(colData(vsd)$gse))]
  ),
  annotation_name_side = "left"
)

col_cor <- colorRamp2(c(0.7, 0.85, 1), c("#2166AC", "white", "#B2182B"))
ht_cor <- Heatmap(cor_mat,
  name = "Pearson r",
  col = col_cor,
  top_annotation = ha_col,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = "Sample-sample Correlation (Guo ATAC-seq)"
)
save_heatmap(ht_cor, "guo_atacseq_sample_correlation", width = 10, height = 9, dir = figdir)

# =============================================================================
# 3. Differential Accessibility (DESeq2) — within-GSE, two-level contrasts
# =============================================================================

message("=== Section 3: Differential Accessibility ===")

# Helper: run DESeq2 on a subset for a two-level contrast (reused ~verbatim from
# core/02_atacseq_analysis.R, which already carries a group_col param). batch_col unused for Guo.
run_da <- function(count_matrix, col_data, contrast_name,
                   num, denom, group_col = "condition",
                   batch_col = NULL) {

  cd <- col_data[col_data[[group_col]] %in% c(num, denom), , drop = FALSE]
  cd[[group_col]] <- droplevels(factor(cd[[group_col]], levels = c(denom, num)))

  n_num   <- sum(cd[[group_col]] == num)
  n_denom <- sum(cd[[group_col]] == denom)
  if (n_num < 2 || n_denom < 2) {
    message(sprintf("    Skipping %s: only %d vs %d samples", contrast_name, n_num, n_denom))
    return(NULL)
  }

  cm <- count_matrix[, rownames(cd), drop = FALSE]

  if (!is.null(batch_col) && length(unique(cd[[batch_col]])) > 1) {
    design_formula <- as.formula(paste("~", batch_col, "+", group_col))
  } else {
    design_formula <- as.formula(paste("~", group_col))
  }
  message(sprintf("    Design: %s (%d vs %d)", deparse(design_formula), n_num, n_denom))

  dds_sub <- DESeqDataSetFromMatrix(countData = cm, colData = cd, design = design_formula)
  keep_sub <- rowSums(counts(dds_sub) >= 10) >= 2
  dds_sub <- dds_sub[keep_sub, ]
  dds_sub <- DESeq(dds_sub)

  coef_name <- paste0(group_col, "_", num, "_vs_", denom)
  if (!(coef_name %in% resultsNames(dds_sub))) {
    coefs <- grep(group_col, resultsNames(dds_sub), value = TRUE)
    coef_name <- coefs[length(coefs)]
  }

  res <- tryCatch(
    lfcShrink(dds_sub, coef = coef_name, type = "apeglm"),
    error = function(e) {
      message(sprintf("    apeglm failed for %s, using normal shrinkage", contrast_name))
      lfcShrink(dds_sub, coef = coef_name, type = "normal")
    }
  )

  res_unshrunk <- results(dds_sub, name = coef_name)

  res_df <- as.data.frame(res) %>%
    rownames_to_column("peak_id") %>%
    as_tibble() %>%
    mutate(
      stat = res_unshrunk$stat[match(peak_id, rownames(res_unshrunk))],
      # Unshrunken LFC retained alongside the apeglm estimate: Guo used DiffBind
      # /DESeq2 defaults (no shrinkage), so published-count comparisons need it.
      # apeglm stays the primary estimate for consistency with
      # core/02_atacseq_analysis.R (and hence the cross-study LFC correlation).
      log2FoldChange_unshrunk =
        res_unshrunk$log2FoldChange[match(peak_id, rownames(res_unshrunk))],
      comparison = contrast_name,
      sig = case_when(
        padj < 0.05 & log2FoldChange > 1  ~ "Gained",
        padj < 0.05 & log2FoldChange < -1 ~ "Lost",
        TRUE ~ "NS"
      )
    ) %>%
    arrange(padj)

  n_gained <- sum(res_df$sig == "Gained", na.rm = TRUE)
  n_lost   <- sum(res_df$sig == "Lost", na.rm = TRUE)
  message(sprintf("    %s: %d gained, %d lost (|LFC|>1, padj<0.05)",
                  contrast_name, n_gained, n_lost))

  return(res_df)
}

all_counts <- counts(dds_full)

da_results <- list()
for (spec in contrast_specs) {
  message(sprintf("  %s (%s: %s vs %s)", spec$name, spec$gse, spec$num, spec$denom))
  gse_meta <- meta[meta$gse == spec$gse, ]
  da_results[[spec$name]] <- run_da(
    all_counts, gse_meta, spec$name, spec$num, spec$denom, group_col = "condition"
  )
}
da_results <- compact(da_results)  # drop any NULLs (skipped contrasts)

# Save checkpoint after Section 3 (avoids re-running DA on restart)
save(dds_full, vsd, meta, peak_info, peak_info_filt, all_counts,
     da_results, contrast_specs,
     file = checkpoint_file)
message(sprintf("  Saved checkpoint: %s", checkpoint_file))

}  # end else (Sections 1-3)

# =============================================================================
# 4. Peak Annotation (ChIPseeker)
# =============================================================================

message("=== Section 4: Peak Annotation ===")

message("  Building TxDb from GENCODE GTF (this may take a few minutes)...")
txdb <- get_txdb()

consensus_gr <- GRanges(
  seqnames = peak_info_filt$Chr,
  ranges   = IRanges(start = peak_info_filt$Start, end = peak_info_filt$End),
  strand   = peak_info_filt$Strand
)
names(consensus_gr) <- peak_info_filt$Geneid

message("  Annotating consensus peaks...")
peak_anno <- annotatePeak(consensus_gr,
  TxDb = txdb,
  annoDb = "org.Mm.eg.db",
  tssRegion = c(-3000, 3000),
  level = "gene"
)

anno_df <- as.data.frame(peak_anno)
anno_df$peak_id <- names(consensus_gr)[match(
  paste(anno_df$seqnames, anno_df$start, anno_df$end),
  paste(peak_info_filt$Chr, peak_info_filt$Start, peak_info_filt$End)
)]
message(sprintf("  annotatePeak returned %d / %d peaks (some scaffolds may be absent from TxDb)",
                nrow(anno_df), length(consensus_gr)))

write_csv(anno_df, file.path(outdir, "consensus_peaks_annotated.csv"))

# ChIPseeker SYMBOL/ENTREZID are NA because GENCODE uses versioned Ensembl IDs
# (e.g. ENSMUSG00000079800.3) that org.Mm.eg.db can't match. Strip + map manually.
anno_df$ensembl_clean <- gsub("\\.\\d+$", "", anno_df$geneId)
gene_map <- AnnotationDbi::select(org.Mm.eg.db,
  keys    = unique(na.omit(anno_df$ensembl_clean)),
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)
anno_df$SYMBOL   <- gene_map$SYMBOL[match(anno_df$ensembl_clean, gene_map$ENSEMBL)]
anno_df$ENTREZID <- gene_map$ENTREZID[match(anno_df$ensembl_clean, gene_map$ENSEMBL)]
message(sprintf("  Gene mapping: %d / %d peaks got SYMBOL, %d got ENTREZID",
                sum(!is.na(anno_df$SYMBOL)), nrow(anno_df),
                sum(!is.na(anno_df$ENTREZID))))

peak_to_gene <- anno_df %>%
  select(peak_id, SYMBOL, ENTREZID, annotation, distanceToTSS) %>%
  rename(gene_name = SYMBOL, entrez_id = ENTREZID)

# Merge gene names into DA results and re-save
for (name in names(da_results)) {
  da_results[[name]] <- da_results[[name]] %>%
    left_join(peak_to_gene, by = "peak_id")
  write_csv(da_results[[name]], file.path(outdir, paste0("da_", name, ".csv")))
}
message("  Merged gene annotations into all DA results")

# --- Annotation pie chart ---
pdf(file.path(figdir, "peak_annotation_pie.pdf"), width = 7, height = 5)
plotAnnoPie(peak_anno)
dev.off()
png(file.path(figdir, "peak_annotation_pie.png"), width = 7, height = 5,
    units = "in", res = 300)
plotAnnoPie(peak_anno)
dev.off()

# --- Annotation comparison for the primary genetic contrast (KO vs WT) ---
if ("KO_vs_WT" %in% names(da_results)) {
  res_pr <- da_results[["KO_vs_WT"]]
  lost_ids   <- res_pr$peak_id[res_pr$sig == "Lost"]
  gained_ids <- res_pr$peak_id[res_pr$sig == "Gained"]

  peak_lists <- list(
    "All consensus" = consensus_gr,
    "Lost in KO"    = consensus_gr[names(consensus_gr) %in% lost_ids],
    "Gained in KO"  = consensus_gr[names(consensus_gr) %in% gained_ids]
  )
  peak_lists <- peak_lists[sapply(peak_lists, length) > 10]

  if (length(peak_lists) > 1) {
    anno_list <- lapply(peak_lists, function(x) {
      annotatePeak(x, TxDb = txdb, tssRegion = c(-3000, 3000), level = "gene")
    })
    p_anno <- plotAnnoBar(anno_list) + theme_paper +
      ggtitle("Genomic annotation: DA peaks in Guo KO vs WT (g183618)")
    save_figure(p_anno, "guo_annotation_comparison_KO_vs_WT",
                width = 10, height = 5, dir = figdir)
  }
}

# =============================================================================
# 5. Volcano Plots
# =============================================================================

message("=== Section 5: Volcano Plots ===")

make_volcano <- function(res, title, filename) {
  top_labs <- res %>%
    filter(sig != "NS", !is.na(gene_name), gene_name != "") %>%
    slice_min(padj, n = 20) %>%
    pull(gene_name)

  n_gained <- sum(res$sig == "Gained", na.rm = TRUE)
  n_lost   <- sum(res$sig == "Lost", na.rm = TRUE)

  lab_col <- ifelse(is.na(res$gene_name) | res$gene_name == "",
                    res$peak_id, res$gene_name)

  p <- EnhancedVolcano(res,
    lab = lab_col,
    x = "log2FoldChange",
    y = "padj",
    title = title,
    subtitle = sprintf("%d gained, %d lost (|LFC|>1, padj<0.05)", n_gained, n_lost),
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 1.5,
    labSize = 3,
    col = c("grey80", "#2CA02C", "#4A90D9", "#D62728"),
    colAlpha = 0.7,
    drawConnectors = TRUE,
    widthConnectors = 0.3,
    maxoverlapsConnectors = 20,
    selectLab = top_labs
  ) + theme_paper

  save_figure(p, filename, width = 8, height = 7, dir = figdir)
  return(invisible(p))
}

walk2(da_results, names(da_results), function(res, name) {
  make_volcano(res, paste("Guo:", gsub("_", " ", name)), paste0("da_volcano_", name))
})

# =============================================================================
# 7. DA Heatmaps — top 50 DA peaks per contrast
# =============================================================================

message("=== Section 7: DA Heatmaps ===")

make_da_heatmap <- function(res, vsd_obj, comparison_name, n_peaks = 50) {
  top_peaks <- res %>%
    filter(sig != "NS") %>%
    slice_min(padj, n = n_peaks) %>%
    pull(peak_id)

  available <- intersect(top_peaks, rownames(assay(vsd_obj)))
  if (length(available) < 5) {
    message(sprintf("  %s: only %d DA peaks available, skipping heatmap",
                    comparison_name, length(available)))
    return(invisible(NULL))
  }

  mat <- assay(vsd_obj)[available, , drop = FALSE]
  mat_z <- t(scale(t(mat)))

  row_labels <- peak_to_gene$gene_name[match(available, peak_to_gene$peak_id)]
  row_labels[is.na(row_labels) | row_labels == ""] <- available[is.na(row_labels) | row_labels == ""]

  cd <- as.data.frame(colData(vsd_obj))
  cd$condition <- droplevels(factor(cd$condition))

  ha <- HeatmapAnnotation(
    Condition = cd$condition,
    col = list(Condition = pal_guo_cond[levels(cd$condition)]),
    annotation_name_side = "left"
  )

  ht <- Heatmap(mat_z,
    name = "z-score",
    col = col_expression,
    top_annotation = ha,
    column_split = cd$condition,
    row_labels = row_labels,
    show_row_names = length(available) <= 50,
    row_names_gp = gpar(fontsize = 7),
    show_column_names = FALSE,
    cluster_columns = FALSE,
    column_title = comparison_name,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    use_raster = TRUE,
    raster_quality = 3
  )

  save_heatmap(ht, paste0("guo_da_heatmap_", comparison_name),
               width = max(8, ncol(mat) * 0.3),
               height = min(14, max(6, length(available) * 0.2)),
               dir = figdir)
}

for (spec in contrast_specs) {
  if (!(spec$name %in% names(da_results))) next
  hm_idx <- which(colData(vsd)$gse == spec$gse &
                  colData(vsd)$condition %in% c(spec$num, spec$denom))
  if (length(hm_idx) > 2) {
    make_da_heatmap(da_results[[spec$name]], vsd[, hm_idx], spec$name)
  }
}

# =============================================================================
# 8. GO Enrichment
# =============================================================================

message("=== Section 8: GO Enrichment ===")

run_go_enrichment <- function(res, contrast_name, direction = "both") {
  results_list <- list()

  for (dir in c("Gained", "Lost")) {
    if (direction != "both" && direction != dir) next

    gene_ids <- res %>%
      filter(sig == dir, !is.na(entrez_id), entrez_id != "") %>%
      distinct(entrez_id) %>%
      pull(entrez_id)

    if (length(gene_ids) < 10) {
      message(sprintf("    %s %s: only %d genes, skipping GO", contrast_name, dir, length(gene_ids)))
      next
    }

    go_res <- enrichGO(
      gene     = gene_ids,
      OrgDb    = org.Mm.eg.db,
      keyType  = "ENTREZID",
      ont      = "BP",
      pvalueCutoff = 0.05,
      readable = TRUE
    )

    if (!is.null(go_res) && nrow(go_res@result) > 0) {
      go_df <- as.data.frame(go_res)
      write_csv(go_df, file.path(outdir,
        paste0("go_da_", tolower(dir), "_", contrast_name, ".csv")))
      message(sprintf("    %s %s: %d GO terms (padj<0.05)",
                      contrast_name, dir, sum(go_df$p.adjust < 0.05)))
      results_list[[dir]] <- go_df
    }
  }

  return(results_list)
}

go_results <- list()
for (name in names(da_results)) {
  message(sprintf("  GO enrichment: %s", name))
  go_results[[name]] <- run_go_enrichment(da_results[[name]], name)
}

# =============================================================================
# 9. DA Summary & Session Save
# =============================================================================

message("=== Section 9: DA Summary ===")

da_summary <- bind_rows(lapply(names(da_results), function(name) {
  res <- da_results[[name]]
  tibble(
    contrast = name,
    n_tested = nrow(res),
    n_gained = sum(res$sig == "Gained", na.rm = TRUE),
    n_lost   = sum(res$sig == "Lost", na.rm = TRUE),
    n_ns     = sum(res$sig == "NS", na.rm = TRUE)
  )
}))
write_csv(da_summary, file.path(outdir, "da_summary.csv"))
message("DA summary:")
print(as.data.frame(da_summary), row.names = FALSE)

da_plot_df <- da_summary %>%
  select(contrast, n_gained, n_lost) %>%
  pivot_longer(cols = c(n_gained, n_lost), names_to = "direction", values_to = "count") %>%
  mutate(
    direction = ifelse(direction == "n_gained", "Gained", "Lost"),
    count_signed = ifelse(direction == "Lost", -count, count)
  )

p_summary <- ggplot(da_plot_df, aes(x = reorder(contrast, abs(count_signed)),
                                     y = count_signed, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Gained" = "#D62728", "Lost" = "#2166AC")) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(title = "Guo: Differentially Accessible Peaks per Contrast",
       subtitle = "|LFC| > 1, padj < 0.05",
       x = "", y = "Number of DA peaks", fill = "Direction") +
  theme_paper
save_figure(p_summary, "guo_barplot_da_summary", width = 10, height = 5, dir = figdir)

message("=== Saving R session ===")
save.image(file.path(outdir, "guo_atacseq_analysis.RData"))

message("\n============================================")
message("Guo ATAC-seq analysis complete!")
message("  Tables: ", outdir)
message("  Figures: ", figdir)
message("============================================")
