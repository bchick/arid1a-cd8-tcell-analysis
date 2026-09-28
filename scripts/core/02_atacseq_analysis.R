#!/usr/bin/env Rscript
# =============================================================================
# core/02_atacseq_analysis.R — ATAC-seq differential accessibility analysis
# McDonald, Chick et al. 2023 Immunity 56:1303 — core analysis
#
# DESeq2 on the nf-core/atacseq featureCounts consensus peak matrix
# (129,314 peaks x 62 samples; faster than DiffBind BAM counting). Contrasts:
# KO vs WT at d3/d5; KO, Het and TbetKO vs WT per d8 subset (experiment batch
# as a covariate where it varies) and in d8 pseudobulk (TE + EEC + MP); and
# the in vitro BAF inhibitor series (IL-12 vs untreated; ACBI1 and BRM014 vs
# IL-12). apeglm shrinkage. Peaks are annotated with ChIPseeker against GENCODE vM35; WT
# time-course OCR dynamics are clustered with k-means; GO enrichment with
# clusterProfiler. Sections 1-3 are checkpointed so a restart skips DA.
#
# Inputs:  results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/
#            consensus_peaks.mRp.clN.featureCounts.txt
#          data/reference/gencode.vM35.primary_assembly.annotation.gtf
# Outputs: results/atac/differential/ (da_*.csv, consensus_peaks_annotated.csv,
#            ocr_clusters_kmeans.csv, ocr_cluster_means.csv, ocr_cluster_*.bed,
#            GO tables, da_summary.csv, atacseq_analysis.RData)
#          figures/atacseq/
# Usage:   Rscript scripts/core/02_atacseq_analysis.R   (from the repository root)
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
outdir <- file.path(paths$results, "atac/differential")
figdir <- file.path(paths$figures, "atacseq")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# Check for checkpoint from previous run (skip Sections 1-3 if available)
checkpoint_file <- file.path(outdir, "checkpoint_sections1to3.RData")
if (file.exists(checkpoint_file)) {
  message("=== Loading checkpoint (Sections 1-3 already complete) ===")
  load(checkpoint_file)
  message(sprintf("  Loaded: %d DA results, %d peaks, %d samples",
                  length(da_results), nrow(peak_info_filt), ncol(vsd)))
} else {

# =============================================================================
# 1. Import featureCounts consensus peak matrix
# =============================================================================

message("=== Section 1: Importing featureCounts consensus peak matrix ===")

fc_file <- file.path(paths$atac,
  "bowtie2/merged_replicate/macs2/narrow_peak/consensus",
  "consensus_peaks.mRp.clN.featureCounts.txt")
stopifnot("featureCounts file not found" = file.exists(fc_file))

# Read featureCounts (skip comment line starting with #)
raw <- read.delim(fc_file, comment.char = "#", check.names = FALSE)
message(sprintf("  Loaded %d peaks x %d columns", nrow(raw), ncol(raw)))

# Extract peak coordinates
peak_info <- raw[, 1:6]
colnames(peak_info) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

# Extract count matrix (columns 7+)
count_mat <- as.matrix(raw[, -(1:6)])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- peak_info$Geneid

# Clean column names: strip BAM suffix, fix REP → Rep
clean_names <- colnames(count_mat)
clean_names <- gsub("\\.mLb\\.clN\\.sorted\\.bam$", "", clean_names)
clean_names <- gsub("REP(\\d)", "Rep\\1", clean_names)
colnames(count_mat) <- clean_names

message(sprintf("  Count matrix: %d peaks x %d samples", nrow(count_mat), ncol(count_mat)))

# Build sample metadata from master sample sheet
master <- load_sample_sheet() %>%
  filter(assay %in% c("atacseq", "atacseq_inhibitors"))

# Verify all featureCounts samples have metadata
missing <- setdiff(colnames(count_mat), master$sample_name)
if (length(missing) > 0) {
  warning(sprintf("  %d samples not in master sheet: %s",
                  length(missing), paste(missing, collapse = ", ")))
}
stopifnot("Sample name mismatch" = all(colnames(count_mat) %in% master$sample_name))

# Order metadata to match count matrix columns and convert to data.frame
meta <- as.data.frame(master[match(colnames(count_mat), master$sample_name), ])
rownames(meta) <- meta$sample_name

# Parse factors
meta$genotype  <- factor(meta$genotype, levels = c("WT", "Het", "KO", "TbetKO"))
meta$timepoint <- factor(meta$timepoint, levels = c("Naive", "48h", "D3", "D5", "D8"))
meta$cell_subset <- ifelse(is.na(meta$cell_subset) | meta$cell_subset == "",
                           "total", meta$cell_subset)
meta$cell_subset <- factor(meta$cell_subset, levels = c("total", "TE", "EEC", "MP"))
meta$experiment  <- ifelse(is.na(meta$experiment) | meta$experiment == "",
                           "Exp0", meta$experiment)
meta$experiment  <- factor(meta$experiment)
meta$replicate   <- factor(meta$replicate)
meta$is_inhibitor <- meta$assay == "atacseq_inhibitors"

# Clean treatment names for inhibitor samples (R-safe factor levels)
meta$treatment <- ifelse(is.na(meta$treatment) | meta$treatment == "",
                         "none", meta$treatment)
meta$treatment_clean <- gsub("[^A-Za-z0-9]", "_", meta$treatment)
meta$treatment_clean <- gsub("_+", "_", meta$treatment_clean)
meta$treatment_clean <- gsub("^_|_$", "", meta$treatment_clean)

# Create group labels
meta$group <- paste(meta$timepoint, meta$genotype, meta$cell_subset, sep = "_")

message(sprintf("  Metadata: %d samples (%d main, %d inhibitor)",
                nrow(meta), sum(!meta$is_inhibitor), sum(meta$is_inhibitor)))
message("  Timepoint x Genotype breakdown (main):")
print(table(meta$timepoint[!meta$is_inhibitor], meta$genotype[!meta$is_inhibitor]))

# Separate main vs inhibitor sample names
main_samples  <- meta$sample_name[!meta$is_inhibitor]
inhib_samples <- meta$sample_name[meta$is_inhibitor]

# Create DESeqDataSet with minimal design (subset per contrast)
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
# 2. QC Visualization
# =============================================================================

message("=== Section 2: QC Visualization ===")

# --- PCA all samples (main experiment only) ---
main_vsd <- vsd[, !colData(vsd)$is_inhibitor]
pca_mat <- assay(main_vsd)
rv <- rowVars(pca_mat)
top_var <- order(rv, decreasing = TRUE)[seq_len(min(5000, length(rv)))]
pca_res <- prcomp(t(pca_mat[top_var, ]), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  genotype    = colData(main_vsd)$genotype,
  timepoint   = colData(main_vsd)$timepoint,
  cell_subset = colData(main_vsd)$cell_subset,
  sample = colnames(main_vsd)
)
pct_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

p_pca_all <- ggplot(pca_df, aes(x = PC1, y = PC2, color = genotype, shape = timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = pal_genotype) +
  labs(
    title = "ATAC-seq PCA: All Samples",
    x = sprintf("PC1 (%.1f%%)", pct_var[1]),
    y = sprintf("PC2 (%.1f%%)", pct_var[2]),
    color = "Genotype", shape = "Timepoint"
  ) +
  theme_paper
save_figure(p_pca_all, "atacseq_pca_all", width = 8, height = 6, dir = figdir)

# --- PCA D8 only ---
d8_idx <- which(colData(main_vsd)$timepoint == "D8")
if (length(d8_idx) > 5) {
  d8_vsd <- main_vsd[, d8_idx]
  pca_mat_d8 <- assay(d8_vsd)
  rv_d8 <- rowVars(pca_mat_d8)
  top_d8 <- order(rv_d8, decreasing = TRUE)[seq_len(min(5000, length(rv_d8)))]
  pca_d8 <- prcomp(t(pca_mat_d8[top_d8, ]), center = TRUE, scale. = FALSE)
  pca_d8_df <- data.frame(
    PC1 = pca_d8$x[, 1],
    PC2 = pca_d8$x[, 2],
    genotype    = colData(d8_vsd)$genotype,
    cell_subset = colData(d8_vsd)$cell_subset,
    experiment  = colData(d8_vsd)$experiment,
    sample = colnames(d8_vsd)
  )
  pct_d8 <- round(100 * summary(pca_d8)$importance[2, 1:2], 1)

  p_pca_d8 <- ggplot(pca_d8_df, aes(x = PC1, y = PC2, color = genotype, shape = cell_subset)) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_color_manual(values = pal_genotype) +
    labs(
      title = "ATAC-seq PCA: D8 Subsets",
      x = sprintf("PC1 (%.1f%%)", pct_d8[1]),
      y = sprintf("PC2 (%.1f%%)", pct_d8[2]),
      color = "Genotype", shape = "Subset"
    ) +
    theme_paper
  save_figure(p_pca_d8, "atacseq_pca_d8", width = 8, height = 6, dir = figdir)
}

# --- PCA inhibitor samples ---
inhib_vsd <- vsd[, colData(vsd)$is_inhibitor]
if (ncol(inhib_vsd) > 3) {
  pca_mat_inh <- assay(inhib_vsd)
  rv_inh <- rowVars(pca_mat_inh)
  top_inh <- order(rv_inh, decreasing = TRUE)[seq_len(min(5000, length(rv_inh)))]
  pca_inh <- prcomp(t(pca_mat_inh[top_inh, ]), center = TRUE, scale. = FALSE)
  pca_inh_df <- data.frame(
    PC1 = pca_inh$x[, 1],
    PC2 = pca_inh$x[, 2],
    treatment = colData(inhib_vsd)$treatment,
    sample = colnames(inhib_vsd)
  )
  pct_inh <- round(100 * summary(pca_inh)$importance[2, 1:2], 1)

  pal_inhib <- c("Untreated" = "#000000", "IL-12" = "#D62728",
                 "IL-12 + ACBI1" = "#2CA02C", "IL-12 + BRM014" = "#9467BD")

  p_pca_inh <- ggplot(pca_inh_df, aes(x = PC1, y = PC2, color = treatment)) +
    geom_point(size = 4, alpha = 0.85) +
    scale_color_manual(values = pal_inhib) +
    labs(
      title = "ATAC-seq PCA: BAF Inhibitor Samples",
      x = sprintf("PC1 (%.1f%%)", pct_inh[1]),
      y = sprintf("PC2 (%.1f%%)", pct_inh[2]),
      color = "Treatment"
    ) +
    theme_paper
  save_figure(p_pca_inh, "inhibitor_pca", width = 8, height = 6, dir = figdir)
}

# --- Sample-sample correlation heatmap (main samples only) ---
cor_mat <- cor(assay(main_vsd), method = "pearson")

ha_col <- HeatmapAnnotation(
  Genotype  = colData(main_vsd)$genotype,
  Timepoint = colData(main_vsd)$timepoint,
  Subset    = colData(main_vsd)$cell_subset,
  col = list(
    Genotype  = pal_genotype[levels(droplevels(colData(main_vsd)$genotype))],
    Timepoint = pal_timepoint[levels(droplevels(colData(main_vsd)$timepoint))]
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
  column_title = "Sample-sample Correlation (ATAC-seq)"
)
save_heatmap(ht_cor, "atacseq_sample_correlation", width = 10, height = 9, dir = figdir)

# =============================================================================
# 3. Differential Accessibility (DESeq2)
# =============================================================================

message("=== Section 3: Differential Accessibility ===")

# Helper: run DESeq2 on a subset of samples with optional batch covariate
run_da <- function(count_matrix, col_data, contrast_name,
                   num, denom, group_col = "genotype",
                   batch_col = NULL) {

  # Subset to relevant groups
  cd <- col_data[col_data[[group_col]] %in% c(num, denom), , drop = FALSE]
  cd[[group_col]] <- droplevels(factor(cd[[group_col]], levels = c(denom, num)))

  # Check minimum sample sizes
  n_num   <- sum(cd[[group_col]] == num)
  n_denom <- sum(cd[[group_col]] == denom)
  if (n_num < 2 || n_denom < 2) {
    message(sprintf("    Skipping %s: only %d vs %d samples", contrast_name, n_num, n_denom))
    return(NULL)
  }

  # Build count matrix subset
  cm <- count_matrix[, rownames(cd), drop = FALSE]

  # Build design formula (include batch if variable and not confounded)
  if (!is.null(batch_col) && length(unique(cd[[batch_col]])) > 1) {
    design_formula <- as.formula(paste("~", batch_col, "+", group_col))
  } else {
    design_formula <- as.formula(paste("~", group_col))
  }
  message(sprintf("    Design: %s (%d vs %d)", deparse(design_formula), n_num, n_denom))

  dds_sub <- DESeqDataSetFromMatrix(countData = cm, colData = cd, design = design_formula)

  # Pre-filter for this subset
  keep_sub <- rowSums(counts(dds_sub) >= 10) >= 2
  dds_sub <- dds_sub[keep_sub, ]

  dds_sub <- DESeq(dds_sub)

  # Find coefficient name for the contrast
  coef_name <- paste0(group_col, "_", num, "_vs_", denom)
  if (!(coef_name %in% resultsNames(dds_sub))) {
    coefs <- grep(group_col, resultsNames(dds_sub), value = TRUE)
    coef_name <- coefs[length(coefs)]
  }

  # apeglm shrinkage with fallback to normal if it fails
  res <- tryCatch(
    lfcShrink(dds_sub, coef = coef_name, type = "apeglm"),
    error = function(e) {
      message(sprintf("    apeglm failed for %s, using normal shrinkage", contrast_name))
      lfcShrink(dds_sub, coef = coef_name, type = "normal")
    }
  )

  # Get Wald stat from unshrunken results (for ranking)
  res_unshrunk <- results(dds_sub, name = coef_name)

  # Convert to tibble
  res_df <- as.data.frame(res) %>%
    rownames_to_column("peak_id") %>%
    as_tibble() %>%
    mutate(
      stat = res_unshrunk$stat[match(peak_id, rownames(res_unshrunk))],
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

# Count matrices for main and inhibitor (pre-filtered peaks)
main_counts  <- counts(dds_full)[, main_samples]
main_meta    <- meta[main_samples, ]
inhib_counts <- counts(dds_full)[, inhib_samples]
inhib_meta   <- meta[inhib_samples, ]

da_results <- list()

# --- D3 KO vs WT ---
message("  D3 KO vs WT")
d3_meta <- main_meta[main_meta$timepoint == "D3", ]
da_results[["KO_vs_WT_D3"]] <- run_da(
  main_counts, d3_meta, "KO_vs_WT_D3", "KO", "WT"
)

# --- D5 KO vs WT ---
message("  D5 KO vs WT")
d5_meta <- main_meta[main_meta$timepoint == "D5", ]
da_results[["KO_vs_WT_D5"]] <- run_da(
  main_counts, d5_meta, "KO_vs_WT_D5", "KO", "WT"
)

# --- D8 per-subset contrasts ---
for (sub in c("TE", "EEC", "MP")) {
  sub_meta <- main_meta[main_meta$timepoint == "D8" & main_meta$cell_subset == sub, ]

  # KO vs WT
  message(sprintf("  D8 %s KO vs WT", sub))
  da_results[[paste0("KO_vs_WT_D8_", sub)]] <- run_da(
    main_counts, sub_meta, paste0("KO_vs_WT_D8_", sub), "KO", "WT",
    batch_col = "experiment"
  )

  # Het vs WT
  message(sprintf("  D8 %s Het vs WT", sub))
  da_results[[paste0("Het_vs_WT_D8_", sub)]] <- run_da(
    main_counts, sub_meta, paste0("Het_vs_WT_D8_", sub), "Het", "WT",
    batch_col = "experiment"
  )

  # TbetKO vs WT
  message(sprintf("  D8 %s TbetKO vs WT", sub))
  da_results[[paste0("TbetKO_vs_WT_D8_", sub)]] <- run_da(
    main_counts, sub_meta, paste0("TbetKO_vs_WT_D8_", sub), "TbetKO", "WT",
    batch_col = "experiment"
  )
}

# --- D8 pseudobulk (sum TE + EEC + MP per genotype + experiment + replicate) ---
message("  Creating D8 pseudobulk samples...")
d8_meta <- main_meta[main_meta$timepoint == "D8", ]

pb_candidates <- d8_meta %>%
  group_by(genotype, experiment, replicate) %>%
  summarize(
    subsets = list(sort(as.character(cell_subset))),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n == 3, sapply(subsets, function(x) all(c("EEC", "MP", "TE") %in% x)))

pb_counts_list <- list()
pb_meta_list <- list()

for (i in seq_len(nrow(pb_candidates))) {
  geno <- as.character(pb_candidates$genotype[i])
  exp  <- as.character(pb_candidates$experiment[i])
  rep  <- as.character(pb_candidates$replicate[i])

  idx <- which(d8_meta$genotype == geno & d8_meta$experiment == exp &
               d8_meta$replicate == rep)
  sample_ids <- d8_meta$sample_name[idx]

  pb_name <- paste0("D8_", geno, "_pb_", exp, "_", rep)
  pb_counts_list[[pb_name]] <- rowSums(main_counts[, sample_ids, drop = FALSE])

  pb_meta_list[[pb_name]] <- data.frame(
    sample_name = pb_name,
    genotype    = geno,
    experiment  = exp,
    replicate   = rep,
    stringsAsFactors = FALSE
  )
}

pb_count_mat <- do.call(cbind, pb_counts_list)
storage.mode(pb_count_mat) <- "integer"
pb_meta_df <- bind_rows(pb_meta_list)
rownames(pb_meta_df) <- pb_meta_df$sample_name
pb_meta_df$genotype   <- factor(pb_meta_df$genotype, levels = c("WT", "Het", "KO", "TbetKO"))
pb_meta_df$experiment <- factor(pb_meta_df$experiment)
pb_meta_df$replicate  <- factor(pb_meta_df$replicate)

message(sprintf("  Created %d pseudobulk samples", ncol(pb_count_mat)))
print(table(pb_meta_df$genotype, pb_meta_df$experiment))

# Pseudobulk contrasts
message("  D8 pseudobulk KO vs WT")
da_results[["KO_vs_WT_D8_pseudobulk"]] <- run_da(
  pb_count_mat, pb_meta_df, "KO_vs_WT_D8_pseudobulk", "KO", "WT",
  batch_col = "experiment"
)

message("  D8 pseudobulk Het vs WT")
da_results[["Het_vs_WT_D8_pseudobulk"]] <- run_da(
  pb_count_mat, pb_meta_df, "Het_vs_WT_D8_pseudobulk", "Het", "WT",
  batch_col = "experiment"
)

message("  D8 pseudobulk TbetKO vs WT")
da_results[["TbetKO_vs_WT_D8_pseudobulk"]] <- run_da(
  pb_count_mat, pb_meta_df, "TbetKO_vs_WT_D8_pseudobulk", "TbetKO", "WT",
  batch_col = "experiment"
)

# --- Inhibitor contrasts ---
message("  Inhibitor contrasts")

# IL-12 vs Untreated
inh_meta_il12_ut <- inhib_meta[inhib_meta$treatment %in% c("Untreated", "IL-12"), ]
da_results[["IL12_vs_Untreated"]] <- run_da(
  inhib_counts, inh_meta_il12_ut, "IL12_vs_Untreated",
  "IL_12", "Untreated", group_col = "treatment_clean"
)

# ACBI1+IL-12 vs IL-12
inh_meta_acbi1 <- inhib_meta[inhib_meta$treatment %in% c("IL-12", "IL-12 + ACBI1"), ]
da_results[["ACBI1_vs_IL12"]] <- run_da(
  inhib_counts, inh_meta_acbi1, "ACBI1_vs_IL12",
  "IL_12_ACBI1", "IL_12", group_col = "treatment_clean"
)

# BRM014+IL-12 vs IL-12
inh_meta_brm <- inhib_meta[inhib_meta$treatment %in% c("IL-12", "IL-12 + BRM014"), ]
da_results[["BRM014_vs_IL12"]] <- run_da(
  inhib_counts, inh_meta_brm, "BRM014_vs_IL12",
  "IL_12_BRM014", "IL_12", group_col = "treatment_clean"
)

# Remove NULL results (skipped contrasts)
da_results <- da_results[!sapply(da_results, is.null)]

# Save all DA results
walk2(da_results, names(da_results), function(res, name) {
  write_csv(res, file.path(outdir, paste0("da_", name, ".csv")))
})
message(sprintf("  Saved %d DA result tables", length(da_results)))

# Save checkpoint after Section 3 (avoids re-running DA on restart)
checkpoint_file <- file.path(outdir, "checkpoint_sections1to3.RData")
save(dds_full, vsd, main_vsd, meta, main_meta, inhib_meta,
     main_samples, inhib_samples, main_counts, inhib_counts,
     peak_info, peak_info_filt, da_results, pb_count_mat, pb_meta_df,
     file = checkpoint_file)
message(sprintf("  Saved checkpoint: %s", checkpoint_file))

}  # end else (Sections 1-3)

# =============================================================================
# 4. Peak Annotation (ChIPseeker)
# =============================================================================

message("=== Section 4: Peak Annotation ===")

# Build TxDb from GENCODE vM35 GTF
message("  Building TxDb from GENCODE GTF (this may take a few minutes)...")
txdb <- get_txdb()

# Build GRanges of filtered consensus peaks
consensus_gr <- GRanges(
  seqnames = peak_info_filt$Chr,
  ranges   = IRanges(start = peak_info_filt$Start, end = peak_info_filt$End),
  strand   = peak_info_filt$Strand
)
names(consensus_gr) <- peak_info_filt$Geneid

# Annotate all consensus peaks
message("  Annotating consensus peaks...")
peak_anno <- annotatePeak(consensus_gr,
  TxDb = txdb,
  annoDb = "org.Mm.eg.db",
  tssRegion = c(-3000, 3000),
  level = "gene"
)

# Extract annotation as data frame
anno_df <- as.data.frame(peak_anno)
anno_df$peak_id <- names(consensus_gr)[match(
  paste(anno_df$seqnames, anno_df$start, anno_df$end),
  paste(peak_info_filt$Chr, peak_info_filt$Start, peak_info_filt$End)
)]
message(sprintf("  annotatePeak returned %d / %d peaks (some scaffolds may be absent from TxDb)",
                nrow(anno_df), length(consensus_gr)))

# Build peak_id → gene symbol lookup
# ChIPseeker SYMBOL/ENTREZID are NA because GENCODE uses versioned Ensembl IDs
# (e.g., ENSMUSG00000079800.3) that org.Mm.eg.db can't match.
# Fix: strip version suffix and map manually.
anno_df$ensembl_clean <- gsub("\\.\\d+$", "", anno_df$geneId)
gene_map <- AnnotationDbi::select(org.Mm.eg.db,
  keys    = unique(na.omit(anno_df$ensembl_clean)),
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)  # one mapping per Ensembl ID
anno_df$SYMBOL   <- gene_map$SYMBOL[match(anno_df$ensembl_clean, gene_map$ENSEMBL)]
anno_df$ENTREZID <- gene_map$ENTREZID[match(anno_df$ensembl_clean, gene_map$ENSEMBL)]
message(sprintf("  Gene mapping: %d / %d peaks got SYMBOL, %d got ENTREZID",
                sum(!is.na(anno_df$SYMBOL)), nrow(anno_df),
                sum(!is.na(anno_df$ENTREZID))))

# written after the Ensembl -> SYMBOL/ENTREZID mapping so the CSV carries ensembl_clean
write_csv(anno_df, file.path(outdir, "consensus_peaks_annotated.csv"))
message(sprintf("  Annotated %d peaks, saved consensus_peaks_annotated.csv", nrow(anno_df)))

peak_to_gene <- anno_df %>%
  select(peak_id, SYMBOL, ENTREZID, annotation, distanceToTSS) %>%
  rename(gene_name = SYMBOL, entrez_id = ENTREZID)

# Merge gene names into DA results
for (name in names(da_results)) {
  da_results[[name]] <- da_results[[name]] %>%
    left_join(peak_to_gene, by = "peak_id")
  # Re-save with gene annotations
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
message("  Saved peak_annotation_pie")

# --- Annotation comparison for D8 KO vs WT (pseudobulk, main contrast) ---
if ("KO_vs_WT_D8_pseudobulk" %in% names(da_results)) {
  res_pb <- da_results[["KO_vs_WT_D8_pseudobulk"]]
  lost_ids   <- res_pb$peak_id[res_pb$sig == "Lost"]
  gained_ids <- res_pb$peak_id[res_pb$sig == "Gained"]

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
      ggtitle("Genomic annotation: DA peaks in D8 KO vs WT (pseudobulk)")
    save_figure(p_anno, "annotation_comparison_D8_KO_vs_WT",
                width = 10, height = 5, dir = figdir)
  }
}

# =============================================================================
# 5. Volcano Plots
# =============================================================================

message("=== Section 5: Volcano Plots ===")

make_volcano <- function(res, title, filename) {
  # Top peaks to label (by padj among significant, labeled by nearest gene)
  top_labs <- res %>%
    filter(sig != "NS", !is.na(gene_name), gene_name != "") %>%
    slice_min(padj, n = 20) %>%
    pull(gene_name)

  n_gained <- sum(res$sig == "Gained", na.rm = TRUE)
  n_lost   <- sum(res$sig == "Lost", na.rm = TRUE)

  # Use gene_name for labeling, fall back to peak_id
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
  make_volcano(res, gsub("_", " ", name), paste0("da_volcano_", name))
})

# =============================================================================
# 6. OCR Dynamics (WT time course, k-means clustering)
# =============================================================================

message("=== Section 6: OCR Dynamics ===")

# Use WT samples only across all timepoints
wt_idx <- which(colData(main_vsd)$genotype == "WT")
wt_vsd <- main_vsd[, wt_idx]
wt_meta_info <- as.data.frame(colData(wt_vsd))

# Define timepoint groups (D8 subsets averaged together)
wt_meta_info$tp_group <- as.character(wt_meta_info$timepoint)

# Average VST signal per timepoint group
tp_levels <- c("Naive", "48h", "D3", "D5", "D8")
wt_avg <- sapply(tp_levels, function(tp) {
  samps <- rownames(wt_meta_info)[wt_meta_info$tp_group == tp]
  if (length(samps) == 1) return(assay(wt_vsd)[, samps])
  rowMeans(assay(wt_vsd)[, samps])
})
colnames(wt_avg) <- tp_levels

# Z-score normalize each peak across timepoints
wt_avg_z <- t(scale(t(wt_avg)))
wt_avg_z[is.nan(wt_avg_z)] <- 0

# k-means clustering (k=5 to match paper)
set.seed(42)
km <- kmeans(wt_avg_z, centers = 5, nstart = 25, iter.max = 100)

message(sprintf("  k-means: %s peaks per cluster",
                paste(table(km$cluster), collapse = ", ")))

# Assign cluster names based on peak timepoint of cluster centers
cluster_centers <- km$centers
peak_tp_idx <- apply(cluster_centers, 1, which.max)
# Sort clusters by peak timepoint index for natural ordering
cluster_order <- order(peak_tp_idx, apply(cluster_centers, 1, max))

# Name mapping: assign biological names based on temporal profile
cluster_names_map <- c(
  "Conserved", "Naive", "Early Activation", "Activation", "Late Activation"
)
# Create ordered cluster labels
ordered_labels <- character(5)
for (i in seq_along(cluster_order)) {
  ordered_labels[cluster_order[i]] <- cluster_names_map[i]
}

cluster_df <- tibble(
  peak_id = rownames(wt_avg_z),
  cluster_num = km$cluster,
  cluster_name = ordered_labels[km$cluster]
)
write_csv(cluster_df, file.path(outdir, "ocr_clusters_kmeans.csv"))

# Export cluster center profiles
centers_df <- as.data.frame(cluster_centers)
centers_df$cluster_num <- 1:5
centers_df$cluster_name <- ordered_labels
write_csv(centers_df, file.path(outdir, "ocr_cluster_means.csv"))

# --- OCR dynamics heatmap ---
row_order <- order(km$cluster)

# mako palette (viridisLite)
col_mako <- colorRamp2(
  seq(-2, 2, length.out = 256),
  rev(mako(256))
)

# Cluster annotation sidebar with distinct colors
cluster_cols <- setNames(
  c("#1A1A1A", "#9ECAE1", "#FF7F0E", "#2CA02C", "#9467BD"),
  cluster_names_map
)

cluster_anno <- rowAnnotation(
  Cluster = factor(ordered_labels[km$cluster[row_order]], levels = cluster_names_map),
  col = list(Cluster = cluster_cols),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = 9)
)

ht_dynamics <- Heatmap(
  wt_avg_z[row_order, ],
  name = "z-score",
  col = col_mako,
  left_annotation = cluster_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 10),
  column_title = "OCR dynamics across timepoints (WT)",
  use_raster = TRUE,
  raster_quality = 3
)
save_heatmap(ht_dynamics, "ocr_dynamics_heatmap", width = 8, height = 12, dir = figdir)

# Export cluster BED files for downstream motif analysis
for (cname in cluster_names_map) {
  cpeaks <- cluster_df %>% filter(cluster_name == cname) %>% pull(peak_id)
  cbed <- peak_info_filt %>%
    filter(Geneid %in% cpeaks) %>%
    select(Chr, Start, End, Geneid)
  write_tsv(cbed, file.path(outdir, paste0("ocr_cluster_",
            gsub(" ", "_", cname), ".bed")), col_names = FALSE)
}
message("  Exported cluster BED files")

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

  # Map peak_ids to gene names for row labels
  row_labels <- peak_to_gene$gene_name[match(available, peak_to_gene$peak_id)]
  row_labels[is.na(row_labels) | row_labels == ""] <- available[is.na(row_labels) | row_labels == ""]

  cd <- as.data.frame(colData(vsd_obj))

  ha <- HeatmapAnnotation(
    Genotype = cd$genotype,
    col = list(
      Genotype = pal_genotype[levels(droplevels(cd$genotype))]
    ),
    annotation_name_side = "left"
  )

  col_split <- factor(
    paste0(cd$timepoint, "\n", cd$cell_subset),
    levels = unique(paste0(cd$timepoint, "\n", cd$cell_subset))
  )

  ht <- Heatmap(mat_z,
    name = "z-score",
    col = col_expression,
    top_annotation = ha,
    column_split = col_split,
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

  save_heatmap(ht, paste0("da_heatmap_", comparison_name),
               width = max(8, ncol(mat) * 0.3),
               height = min(14, max(6, length(available) * 0.2)),
               dir = figdir)
}

# Heatmaps for key contrasts using the main experiment VST
key_contrasts_main <- grep("^(KO_vs_WT|Het_vs_WT)", names(da_results), value = TRUE)
for (name in key_contrasts_main) {
  # Determine which samples to show in heatmap
  if (grepl("pseudobulk", name)) next  # skip pseudobulk for heatmaps
  if (grepl("D3", name)) {
    hm_idx <- which(colData(main_vsd)$timepoint == "D3")
  } else if (grepl("D5", name)) {
    hm_idx <- which(colData(main_vsd)$timepoint == "D5")
  } else if (grepl("D8_(TE|EEC|MP)", name)) {
    sub <- sub(".*D8_(TE|EEC|MP)$", "\\1", name)
    hm_idx <- which(colData(main_vsd)$timepoint == "D8" &
                    colData(main_vsd)$cell_subset == sub)
  } else {
    next
  }
  if (length(hm_idx) > 2) {
    make_da_heatmap(da_results[[name]], main_vsd[, hm_idx], name)
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

# Run GO for key contrasts
go_results <- list()
for (name in names(da_results)) {
  message(sprintf("  GO enrichment: %s", name))
  go_results[[name]] <- run_go_enrichment(da_results[[name]], name)
}

# =============================================================================
# 9. DA Summary & Session Save
# =============================================================================

message("=== Section 9: DA Summary ===")

# Summary table
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

# Summary barplot
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
  labs(title = "Differentially Accessible Peaks per Contrast",
       subtitle = "|LFC| > 1, padj < 0.05",
       x = "", y = "Number of DA peaks", fill = "Direction") +
  theme_paper
save_figure(p_summary, "da_summary_barplot", width = 10, height = 7, dir = figdir)

# Save R session
message("=== Saving R session ===")
save.image(file.path(outdir, "atacseq_analysis.RData"))
message("ATAC-seq analysis complete.")
message("  Results: ", outdir)
message("  Figures: ", figdir)
