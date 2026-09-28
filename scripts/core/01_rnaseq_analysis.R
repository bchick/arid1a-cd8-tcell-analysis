#!/usr/bin/env Rscript
# =============================================================================
# core/01_rnaseq_analysis.R — RNA-seq differential expression analysis
# McDonald, Chick et al. 2023 Immunity 56:1303 — core analysis
#
# Differential expression between WT, Arid1a cHet and cKO CD8+ T cells
# (GSE227634; d3, d5 and d8 TE/EEC/MP) from the nf-core/rnaseq STAR+Salmon
# gene counts. D8 subsets are also pseudobulked (TE + EEC + MP summed per
# genotype and replicate). Each contrast is fit separately with DESeq2
# (~ genotype) and log2 fold changes are shrunk with apeglm; the Wald statistic
# from the unshrunken fit is kept for GSEA ranking. Also produces QC PCA and
# sample correlation, volcano plots, DE heatmaps, T cell signature gene plots,
# fgsea (Hallmark + ImmuneSigDB), GO over-representation (clusterProfiler) and
# a DE summary.
#
# Inputs:  results/rnaseq/star_salmon/salmon.merged.gene_counts.tsv
# Outputs: results/rnaseq/differential/ (normalized_counts.csv, de_*.csv,
#            gsea_*.csv, go_{up,down}_*.csv, de_summary.csv, rnaseq_analysis.RData)
#          figures/rnaseq/
# Usage:   Rscript scripts/core/01_rnaseq_analysis.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(fgsea)
  library(msigdbr)
  library(org.Mm.eg.db)
  library(viridisLite)
})

# Resolve namespace conflicts: Bioconductor packages mask dplyr verbs
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate
slice  <- dplyr::slice

# Output directories
outdir <- file.path(paths$results, "rnaseq/differential")
figdir <- file.path(paths$figures, "rnaseq")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Import merged Salmon gene counts
# =============================================================================

message("=== Section 1: Importing merged Salmon gene counts ===")

# Load the merged count matrix produced by nf-core/rnaseq
counts_file <- file.path(paths$rnaseq, "star_salmon/salmon.merged.gene_counts.tsv")
stopifnot("Merged gene counts file not found" = file.exists(counts_file))

raw <- read.delim(counts_file, check.names = FALSE)
message(sprintf("  Loaded %d genes x %d sample columns", nrow(raw), ncol(raw) - 2))

# Build gene_id -> gene_name mapping
gene_map <- data.frame(
  gene_id   = raw$gene_id,
  gene_name = raw$gene_name,
  stringsAsFactors = FALSE
)
# Handle duplicate gene names by appending gene_id suffix
dupes <- duplicated(gene_map$gene_name) | duplicated(gene_map$gene_name, fromLast = TRUE)
gene_map$gene_name_unique <- ifelse(
  dupes & gene_map$gene_name != "",
  paste0(gene_map$gene_name, "_", gsub("\\..*", "", gene_map$gene_id)),
  gene_map$gene_name
)
# For empty gene names, use the Ensembl ID
gene_map$gene_name_unique[gene_map$gene_name_unique == ""] <- gene_map$gene_id[gene_map$gene_name_unique == ""]

# Extract count matrix (round to integers for DESeq2)
count_mat <- as.matrix(raw[, -(1:2)])
rownames(count_mat) <- raw$gene_id
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

# Build sample metadata from master sample sheet
master <- load_sample_sheet() %>% filter(assay == "rnaseq")

# Ensure all count matrix columns have metadata
sample_names <- colnames(count_mat)
stopifnot("Sample name mismatch between counts and metadata" =
            all(sample_names %in% master$sample_name))

# Order metadata to match count matrix columns and convert to data.frame
meta <- as.data.frame(master[match(sample_names, master$sample_name), ])
rownames(meta) <- meta$sample_name

# Parse factors
meta$genotype  <- factor(meta$genotype, levels = c("WT", "Het", "KO"))
meta$timepoint <- factor(meta$timepoint, levels = c("D3", "D5", "D8"))
meta$cell_subset <- factor(meta$cell_subset)
meta$replicate <- factor(meta$replicate)

# Create a combined group label for convenience
meta$group <- paste(meta$timepoint, meta$genotype, meta$cell_subset, sep = "_")

message(sprintf("  Sample metadata: %d samples", nrow(meta)))
message("  Timepoint x Genotype breakdown:")
print(table(meta$timepoint, meta$genotype))

# Build DESeqDataSet with minimal design (we subset per contrast)
dds_full <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ 1
)

# Pre-filter: genes with >= 10 counts in >= 3 samples
keep <- rowSums(counts(dds_full) >= 10) >= 3
dds_full <- dds_full[keep, ]
message(sprintf("  Kept %d / %d genes after pre-filtering", sum(keep), length(keep)))

# VST for visualization (using all samples)
vsd <- vst(dds_full, blind = FALSE)
message("  VST transformation complete")

# Save normalized counts
dds_norm <- DESeq(dds_full)
norm_counts <- counts(dds_norm, normalized = TRUE)
# Add gene names
norm_out <- data.frame(
  gene_id   = rownames(norm_counts),
  gene_name = gene_map$gene_name[match(rownames(norm_counts), gene_map$gene_id)],
  norm_counts,
  check.names = FALSE
)
write.csv(norm_out, file.path(outdir, "normalized_counts.csv"), row.names = FALSE)
message("  Saved normalized_counts.csv")

# =============================================================================
# 2. Pseudo-bulk D8 samples (sum TE + EEC + MP per genotype + replicate)
# =============================================================================

message("=== Section 2: Creating pseudo-bulk D8 samples ===")

d8_meta <- meta %>% filter(timepoint == "D8")

# Find replicates that have all 3 subsets (TE, EEC, MP)
pb_candidates <- d8_meta %>%
  group_by(genotype, replicate) %>%
  summarize(subsets = list(sort(as.character(cell_subset))), n = n(), .groups = "drop") %>%
  filter(n == 3, sapply(subsets, function(x) all(c("EEC", "MP", "TE") %in% x)))

message(sprintf("  Pseudo-bulk replicates found: %d", nrow(pb_candidates)))
print(pb_candidates %>% select(genotype, replicate, n))

# Sum counts across subsets for each genotype+replicate
pb_counts_list <- list()
pb_meta_list <- list()

for (i in seq_len(nrow(pb_candidates))) {
  geno <- as.character(pb_candidates$genotype[i])
  rep  <- as.character(pb_candidates$replicate[i])

  # Get the 3 samples (TE, EEC, MP) for this genotype+replicate
  idx <- which(d8_meta$genotype == geno & d8_meta$replicate == rep)
  sample_ids <- d8_meta$sample_name[idx]

  # Sum raw counts
  pb_name <- paste0("D8_", geno, "_pseudobulk_", rep)
  pb_counts_list[[pb_name]] <- rowSums(count_mat[rownames(dds_full), sample_ids, drop = FALSE])

  pb_meta_list[[pb_name]] <- data.frame(
    sample_name = pb_name,
    genotype    = geno,
    timepoint   = "D8",
    cell_subset = "pseudobulk",
    replicate   = rep,
    stringsAsFactors = FALSE
  )
}

pb_count_mat <- do.call(cbind, pb_counts_list)
storage.mode(pb_count_mat) <- "integer"
pb_meta <- bind_rows(pb_meta_list)
rownames(pb_meta) <- pb_meta$sample_name
pb_meta$genotype  <- factor(pb_meta$genotype, levels = c("WT", "Het", "KO"))
pb_meta$timepoint <- factor(pb_meta$timepoint)
pb_meta$replicate <- factor(pb_meta$replicate)

message(sprintf("  Created %d pseudo-bulk samples", ncol(pb_count_mat)))
print(table(pb_meta$genotype))

# =============================================================================
# 3. QC Visualization
# =============================================================================

message("=== Section 3: QC Visualization ===")

# --- PCA all samples ---
pca_mat <- assay(vsd)
rv <- rowVars(pca_mat)
top_var <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca_res <- prcomp(t(pca_mat[top_var, ]), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  genotype  = colData(vsd)$genotype,
  timepoint = colData(vsd)$timepoint,
  cell_subset = colData(vsd)$cell_subset,
  sample = colnames(vsd)
)
pct_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

p_pca_all <- ggplot(pca_df, aes(x = PC1, y = PC2, color = genotype, shape = timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = pal_genotype) +
  labs(
    title = "RNA-seq PCA: All Samples",
    x = sprintf("PC1 (%.1f%%)", pct_var[1]),
    y = sprintf("PC2 (%.1f%%)", pct_var[2]),
    color = "Genotype", shape = "Timepoint"
  ) +
  theme_paper
save_figure(p_pca_all, "rnaseq_pca_all", width = 8, height = 6, dir = figdir)

# --- PCA D8 subsets only ---
d8_idx <- which(pca_df$timepoint == "D8")
if (length(d8_idx) > 5) {
  d8_vsd <- vsd[, colData(vsd)$timepoint == "D8"]
  pca_mat_d8 <- assay(d8_vsd)
  rv_d8 <- rowVars(pca_mat_d8)
  top_d8 <- order(rv_d8, decreasing = TRUE)[seq_len(min(500, length(rv_d8)))]
  pca_d8 <- prcomp(t(pca_mat_d8[top_d8, ]), center = TRUE, scale. = FALSE)
  pca_d8_df <- data.frame(
    PC1 = pca_d8$x[, 1],
    PC2 = pca_d8$x[, 2],
    genotype    = colData(d8_vsd)$genotype,
    cell_subset = colData(d8_vsd)$cell_subset,
    sample = colnames(d8_vsd)
  )
  pct_d8 <- round(100 * summary(pca_d8)$importance[2, 1:2], 1)

  p_pca_d8 <- ggplot(pca_d8_df, aes(x = PC1, y = PC2, color = genotype, shape = cell_subset)) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_color_manual(values = pal_genotype) +
    labs(
      title = "RNA-seq PCA: D8 Subsets",
      x = sprintf("PC1 (%.1f%%)", pct_d8[1]),
      y = sprintf("PC2 (%.1f%%)", pct_d8[2]),
      color = "Genotype", shape = "Subset"
    ) +
    theme_paper
  save_figure(p_pca_d8, "rnaseq_pca_d8_subsets", width = 8, height = 6, dir = figdir)
}

# --- Sample-sample correlation heatmap ---
cor_mat <- cor(assay(vsd), method = "pearson")

ha_col <- HeatmapAnnotation(
  Genotype  = colData(vsd)$genotype,
  Timepoint = colData(vsd)$timepoint,
  Subset    = colData(vsd)$cell_subset,
  col = list(
    Genotype  = pal_genotype[levels(colData(vsd)$genotype)],
    Timepoint = pal_timepoint[as.character(unique(colData(vsd)$timepoint))]
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
  column_title = "Sample-sample Correlation (RNA-seq)"
)
save_heatmap(ht_cor, "rnaseq_sample_correlation", width = 10, height = 9, dir = figdir)

# =============================================================================
# 4. Differential Expression (DESeq2) — per-contrast subsetting
# =============================================================================

message("=== Section 4: Differential Expression ===")

# Helper: run DESeq2 on a subset of samples, return annotated results
run_de <- function(count_matrix, col_data, contrast_name, num = "KO", denom = "WT") {
  # Subset to relevant genotypes
  cd <- col_data[col_data$genotype %in% c(num, denom), , drop = FALSE]
  cd$genotype <- droplevels(factor(cd$genotype, levels = c(denom, num)))

  # Select matching columns from the count matrix by sample name
  gene_rows <- intersect(rownames(dds_full), rownames(count_matrix))
  cm <- count_matrix[gene_rows, rownames(cd), drop = FALSE]

  dds_sub <- DESeqDataSetFromMatrix(countData = cm, colData = cd, design = ~ genotype)
  keep_sub <- rowSums(counts(dds_sub) >= 10) >= 3
  dds_sub <- dds_sub[keep_sub, ]
  dds_sub <- DESeq(dds_sub)

  # apeglm shrinkage
  coef_name <- paste0("genotype_", num, "_vs_", denom)
  if (!(coef_name %in% resultsNames(dds_sub))) {
    # Fallback: find the genotype coefficient
    coef_name <- grep("genotype", resultsNames(dds_sub), value = TRUE)[1]
  }
  res <- lfcShrink(dds_sub, coef = coef_name, type = "apeglm")

  # Get Wald stat from unshrunken results (needed for GSEA ranking)
  res_unshrunk <- results(dds_sub, name = coef_name)

  # Convert to tibble with gene names
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    mutate(stat = res_unshrunk$stat[match(gene_id, rownames(res_unshrunk))]) %>%
    left_join(gene_map %>% select(gene_id, gene_name), by = "gene_id") %>%
    mutate(
      comparison = contrast_name,
      sig = case_when(
        padj < 0.05 & log2FoldChange > 1  ~ "Up",
        padj < 0.05 & log2FoldChange < -1 ~ "Down",
        TRUE ~ "NS"
      )
    ) %>%
    arrange(padj)

  return(res_df)
}

de_results <- list()

# --- D3 bulk: KO vs WT ---
message("  Running D3 KO vs WT...")
d3_idx <- meta$timepoint == "D3"
de_results[["KO_vs_WT_D3"]] <- run_de(
  count_mat, meta[d3_idx, ], "KO_vs_WT_D3", "KO", "WT"
)

# --- D8 per-subset contrasts ---
for (sub in c("TE", "EEC", "MP")) {
  sub_idx <- meta$timepoint == "D8" & meta$cell_subset == sub

  # KO vs WT
  message(sprintf("  Running D8 %s KO vs WT...", sub))
  de_results[[paste0("KO_vs_WT_D8_", sub)]] <- run_de(
    count_mat, meta[sub_idx, ], paste0("KO_vs_WT_D8_", sub), "KO", "WT"
  )

  # Het vs WT
  message(sprintf("  Running D8 %s Het vs WT...", sub))
  de_results[[paste0("Het_vs_WT_D8_", sub)]] <- run_de(
    count_mat, meta[sub_idx, ], paste0("Het_vs_WT_D8_", sub), "Het", "WT"
  )
}

# --- D8 pseudo-bulk contrasts ---
message("  Running D8 pseudo-bulk KO vs WT...")
de_results[["KO_vs_WT_D8_pseudobulk"]] <- run_de(
  pb_count_mat, pb_meta, "KO_vs_WT_D8_pseudobulk", "KO", "WT"
)

message("  Running D8 pseudo-bulk Het vs WT...")
de_results[["Het_vs_WT_D8_pseudobulk"]] <- run_de(
  pb_count_mat, pb_meta, "Het_vs_WT_D8_pseudobulk", "Het", "WT"
)

# Save all DE results
walk2(de_results, names(de_results), function(res, name) {
  write_csv(res, file.path(outdir, paste0("de_", name, ".csv")))
  n_up   <- sum(res$sig == "Up", na.rm = TRUE)
  n_down <- sum(res$sig == "Down", na.rm = TRUE)
  message(sprintf("  %s: %d up, %d down (|LFC|>1, padj<0.05)", name, n_up, n_down))
})

# =============================================================================
# 5. Volcano Plots
# =============================================================================

message("=== Section 5: Volcano Plots ===")

make_volcano <- function(res, title, filename) {
  # Top genes to label (by padj among significant)
  top_labs <- res %>%
    filter(sig != "NS", !is.na(gene_name), gene_name != "") %>%
    slice_min(padj, n = 20) %>%
    pull(gene_name)

  n_up   <- sum(res$sig == "Up", na.rm = TRUE)
  n_down <- sum(res$sig == "Down", na.rm = TRUE)

  p <- EnhancedVolcano(res,
    lab = ifelse(is.na(res$gene_name) | res$gene_name == "", res$gene_id, res$gene_name),
    x = "log2FoldChange",
    y = "padj",
    title = title,
    subtitle = sprintf("%d up, %d down (|LFC|>1, padj<0.05)", n_up, n_down),
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

walk2(de_results, names(de_results), function(res, name) {
  make_volcano(res, gsub("_", " ", name), paste0("volcano_", name))
})

# =============================================================================
# 6. DE Heatmaps — top 50 DE genes per contrast
# =============================================================================

message("=== Section 6: DE Heatmaps ===")

# Helper: deterministic sample ordering for heatmaps
# Groups by timepoint -> cell_subset -> genotype (WT, Het, KO) -> replicate
order_samples <- function(col_data) {
  cd <- as.data.frame(col_data)
  cd$.idx <- seq_len(nrow(cd))

  tp_order   <- c("D3", "D5", "D8")
  sub_order  <- c("total", "TE", "EEC", "MP")
  geno_order <- c("WT", "Het", "KO")

  cd$tp_rank   <- match(as.character(cd$timepoint), tp_order)
  cd$sub_rank  <- match(as.character(cd$cell_subset), sub_order)
  cd$geno_rank <- match(as.character(cd$genotype), geno_order)
  cd$rep_rank  <- as.integer(cd$replicate)

  # Replace NAs with high values so they sort last
  cd$tp_rank[is.na(cd$tp_rank)]   <- 99
  cd$sub_rank[is.na(cd$sub_rank)] <- 99
  cd$geno_rank[is.na(cd$geno_rank)] <- 99
  cd$rep_rank[is.na(cd$rep_rank)]   <- 99

  ord <- order(cd$tp_rank, cd$sub_rank, cd$geno_rank, cd$rep_rank)
  return(ord)
}

# Build consistent top annotation + column split for ordered samples
build_heatmap_anno <- function(col_data_ordered) {
  cd <- as.data.frame(col_data_ordered)

  # Genotype bar only — timepoint/subset encoded in split labels
  ha <- HeatmapAnnotation(
    Genotype  = cd$genotype,
    col = list(
      Genotype = pal_genotype[intersect(levels(cd$genotype), names(pal_genotype))]
    ),
    show_legend = FALSE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
  )

  # Two-level column split: group + genotype
  # Labels like "D3 total\nWT", "D8 TE\nKO"
  group_label <- paste0(cd$timepoint, " ", cd$cell_subset, "\n", cd$genotype)
  col_split <- factor(group_label, levels = unique(group_label))

  return(list(annotation = ha, split = col_split))
}

make_de_heatmap <- function(res, vsd_obj, comparison_name, n_genes = 50) {
  top_genes <- res %>%
    filter(sig != "NS") %>%
    slice_min(padj, n = n_genes) %>%
    pull(gene_id)

  # Intersect with genes in vsd
  available <- intersect(top_genes, rownames(assay(vsd_obj)))

  if (length(available) < 5) {
    message(sprintf("  %s: only %d DE genes available, skipping heatmap", comparison_name, length(available)))
    return(invisible(NULL))
  }

  # Order columns by timepoint -> subset -> genotype (no hierarchical clustering)
  col_order <- order_samples(colData(vsd_obj))

  mat <- assay(vsd_obj)[available, col_order, drop = FALSE]
  mat_z <- t(scale(t(mat)))

  # Map rownames to gene symbols
  row_labels <- gene_map$gene_name[match(available, gene_map$gene_id)]
  row_labels[is.na(row_labels) | row_labels == ""] <- available[is.na(row_labels) | row_labels == ""]

  cd <- colData(vsd_obj)[col_order, ]
  col_split <- factor(
    paste0(cd$timepoint, "\n", cd$cell_subset),
    levels = unique(paste0(cd$timepoint, "\n", cd$cell_subset))
  )

  ha <- HeatmapAnnotation(
    Genotype  = cd$genotype,
    col = list(
      Genotype = pal_genotype[levels(colData(vsd_obj)$genotype)]
    ),
    annotation_name_side = "left"
  )

  ht <- Heatmap(mat_z,
    name = "z-score",
    col = col_expression,
    top_annotation = ha,
    cluster_columns = FALSE,
    column_split = col_split,
    cluster_column_slices = FALSE,
    cluster_rows = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_labels = row_labels,
    row_names_gp = gpar(fontsize = 7),
    column_title = gsub("_", " ", comparison_name),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_gap = unit(2, "mm"),
    row_title = sprintf("Top %d DE genes", length(available))
  )

  save_heatmap(ht, paste0("heatmap_top_DE_", comparison_name),
               width = 10, height = max(6, length(available) * 0.2), dir = figdir)
}

walk2(de_results, names(de_results), function(res, name) {
  make_de_heatmap(res, vsd, name)
})

# =============================================================================
# 7. T Cell Signature Gene Analysis
# =============================================================================

message("=== Section 7: T Cell Signature Gene Analysis ===")

# Curated gene lists from the paper
te_signature <- c("Tbx21", "Zeb2", "Cx3cr1", "Klrg1", "S1pr5", "Gzmb",
                  "Gzma", "Prf1", "Id2", "Bhlhe40", "Runx3")
mp_signature <- c("Tcf7", "Id3", "Bach2", "Bcl2", "Il7r", "Sell",
                  "Ccr7", "Lef1", "Myb", "Foxo1")
trm_signature <- c("Zfp683", "Itgae", "Itga1", "Cd69", "Bhlhe40",
                    "Runx3", "Prdm1", "Nr4a1")
effector_tfs <- c("Tbx21", "Eomes", "Batf", "Irf4", "Runx3",
                   "Bhlhe40", "Zeb2", "Id2", "Prdm1")

gene_signatures <- list(
  "TE_signature"  = te_signature,
  "MP_signature"  = mp_signature,
  "Trm_signature" = trm_signature,
  "Effector_TFs"  = effector_tfs
)

# Map gene symbols -> Ensembl IDs using our gene_map
symbol_to_ensembl <- function(symbols) {
  matched <- gene_map %>%
    filter(gene_name %in% symbols) %>%
    distinct(gene_name, .keep_all = TRUE)
  # Only keep genes present in our filtered dataset
  matched <- matched %>% filter(gene_id %in% rownames(assay(vsd)))
  return(matched)
}

make_signature_heatmap <- function(gene_list, sig_name) {
  matched <- symbol_to_ensembl(gene_list)

  if (nrow(matched) < 3) {
    message(sprintf("  %s: only %d/%d genes found, skipping", sig_name,
                    nrow(matched), length(gene_list)))
    return(invisible(NULL))
  }

  # Order columns by timepoint -> subset -> genotype
  col_order <- order_samples(colData(vsd))

  mat <- assay(vsd)[matched$gene_id, col_order, drop = FALSE]
  mat_z <- t(scale(t(mat)))

  cd <- colData(vsd)[col_order, ]
  col_split <- factor(
    paste0(cd$timepoint, "\n", cd$cell_subset),
    levels = unique(paste0(cd$timepoint, "\n", cd$cell_subset))
  )

  ha <- HeatmapAnnotation(
    Genotype  = cd$genotype,
    col = list(
      Genotype = pal_genotype[levels(colData(vsd)$genotype)]
    )
  )

  ht <- Heatmap(mat_z,
    name = "z-score",
    col = col_expression,
    top_annotation = ha,
    cluster_columns = FALSE,
    column_split = col_split,
    cluster_column_slices = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    row_labels = matched$gene_name,
    row_names_gp = gpar(fontsize = 10, fontface = "italic"),
    column_title = gsub("_", " ", sig_name),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_gap = unit(2, "mm")
  )

  save_heatmap(ht, paste0("heatmap_", sig_name),
               width = 10, height = max(4, nrow(matched) * 0.4), dir = figdir)

  message(sprintf("  %s: plotted %d/%d genes", sig_name, nrow(matched), length(gene_list)))
}

walk2(gene_signatures, names(gene_signatures), make_signature_heatmap)

# --- Figure 2H: Curated gene expression heatmap (D3 KO vs WT) ---
message("  Generating Figure 2H heatmap...")

# Curated gene list from Figure 2H, in exact column order with group assignments
fig2h_genes <- data.frame(
  gene_name = c(
    # Transcription factors
    "Bhlhe40", "Eomes", "Ar", "Rxra", "Atf3", "Zeb2", "Id2", "Tcf3",
    "Tbx21", "Batf", "Batf3", "Runx3", "Irf4", "Maf", "Myb", "Tcf7",
    "Klf3", "Foxn3", "Zfp395", "Bhlhe41", "Tox2", "Irf6", "Prdm1",
    "Runx1", "Chd3", "Tnfaip3",
    # Cell cycle
    "Top2a", "Cdkn2d", "Cdkn1b", "Cdkn2a", "Mki67",
    # Cytokine / TNF receptors
    "Tnfrsf14", "Il12rb2", "Il2rb", "Il7r", "Il18r1",
    # Chemokine receptors / migration
    "Cx3cr1", "Cxcr3", "S1pr1", "Ccr7", "Ccr5", "Ccr2", "Ccr9",
    # Integrins / adhesion
    "Itga1", "Itga2", "Itgax", "Itgam", "Itga4", "Itga6", "Cd69",
    # Effector molecules
    "Gzma", "Gzmb", "Gzmk", "Gzmm", "Ccl3", "Ccl5",
    # Metabolism / proliferation
    "Myc", "Srm", "Dusp2", "Tfrc", "Slc7a5"
  ),
  group = c(
    rep("TFs", 26),
    rep("Cell cycle", 5),
    rep("Receptors", 5),
    rep("Chemokine R", 7),
    rep("Integrins", 7),
    rep("Effector", 6),
    rep("Metabolism", 5)
  ),
  stringsAsFactors = FALSE
)
fig2h_genes$order <- seq_len(nrow(fig2h_genes))

# Map to Ensembl IDs
fig2h_matched <- fig2h_genes %>%
  left_join(gene_map %>% distinct(gene_name, .keep_all = TRUE), by = "gene_name") %>%
  filter(!is.na(gene_id), gene_id %in% rownames(assay(vsd)))

message(sprintf("  Figure 2H: matched %d/%d genes", nrow(fig2h_matched), nrow(fig2h_genes)))

if (nrow(fig2h_matched) >= 30) {
  # Subset to D3 samples only
  d3_vsd <- vsd[, colData(vsd)$timepoint == "D3"]
  # Order: WT first, then KO
  d3_cd <- as.data.frame(colData(d3_vsd))
  d3_order <- order(match(d3_cd$genotype, c("WT", "KO")), d3_cd$replicate)
  d3_vsd <- d3_vsd[, d3_order]
  d3_cd <- as.data.frame(colData(d3_vsd))

  # Extract expression matrix: samples (rows) x genes (columns) — transposed
  mat <- t(assay(d3_vsd)[fig2h_matched$gene_id, , drop = FALSE])
  colnames(mat) <- fig2h_matched$gene_name

  # Z-score per gene (column)
  mat_z <- scale(mat)
  mat_z[mat_z > 2] <- 2
  mat_z[mat_z < -2] <- -2

  # DEG annotation
  d3_de <- de_results[["KO_vs_WT_D3"]]
  is_deg <- fig2h_matched$gene_name %in% (d3_de %>% filter(sig != "NS") %>% pull(gene_name))

  # Column split by gene group
  col_split <- factor(fig2h_matched$group,
    levels = c("TFs", "Cell cycle", "Receptors", "Chemokine R",
               "Integrins", "Effector", "Metabolism"))

  # Row split
  row_split <- factor(d3_cd$genotype, levels = c("WT", "KO"))

  # Define color palettes to generate
  fig2h_palettes <- list(
    "blured" = colorRamp2(
      c(-2, -1, 0, 1, 2),
      c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
    ),
    "mako" = colorRamp2(seq(-2, 2, length.out = 9), viridisLite::mako(9)),
    "rocket" = colorRamp2(seq(-2, 2, length.out = 9), viridisLite::rocket(9)),
    "viridis" = colorRamp2(seq(-2, 2, length.out = 9), viridisLite::viridis(9)),
    "magma" = colorRamp2(seq(-2, 2, length.out = 9), viridisLite::magma(9)),
    "turbo" = colorRamp2(seq(-2, 2, length.out = 9), viridisLite::turbo(9))
  )

  # Function to build and save the heatmap for a given palette
  make_fig2h <- function(col_fun, palette_name) {
    # Top annotation: DEG bar with border
    deg_status <- ifelse(is_deg, "Yes", "No")
    ha_top <- HeatmapAnnotation(
      DEG = deg_status,
      col = list(DEG = c("Yes" = "#7B2D8E", "No" = "grey80")),
      border = TRUE,
      simple_anno_size = unit(3, "mm"),
      show_annotation_name = TRUE,
      annotation_name_side = "left",
      annotation_name_gp = gpar(fontsize = 8)
    )

    # Left annotation: genotype as colored blocks with text inside
    geno_colors <- c("WT" = "black", "KO" = "#2CA02C")
    ha_row <- rowAnnotation(
      Genotype = anno_block(
        gp = gpar(fill = geno_colors[levels(row_split)], col = "black"),
        labels = levels(row_split),
        labels_gp = gpar(col = "white", fontsize = 9, fontface = "bold"),
        width = unit(8, "mm")
      )
    )

    ht <- Heatmap(mat_z,
      name = "z-score",
      col = col_fun,
      top_annotation = ha_top,
      left_annotation = ha_row,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_split = col_split,
      cluster_column_slices = FALSE,
      show_row_names = FALSE,
      show_column_names = TRUE,
      column_names_gp = gpar(fontsize = 7, fontface = "italic"),
      column_names_rot = 90,
      column_gap = unit(1.5, "mm"),
      column_title_gp = gpar(fontsize = 8),
      border = TRUE,
      row_split = row_split,
      row_gap = unit(1, "mm"),
      row_title = NULL,
      show_row_dend = FALSE,
      heatmap_legend_param = list(
        title = "z-score",
        at = c(-2, -1, 0, 1, 2),
        legend_height = unit(3, "cm")
      ),
      show_heatmap_legend = TRUE
    )

    suffix <- if (palette_name == "blured") "" else paste0("_", palette_name)
    fname <- paste0("fig2H_D3_gene_heatmap", suffix)
    base <- tools::file_path_sans_ext(fname)
    w <- 14; h <- 3.5
    pdf(file.path(figdir, paste0(base, ".pdf")), width = w, height = h)
    draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
    dev.off()
    png(file.path(figdir, paste0(base, ".png")), width = w, height = h,
        units = "in", res = 300)
    draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
    dev.off()
    message(sprintf("  Saved Figure 2H (%s)", palette_name))
  }

  # Generate all palette versions
  for (pal_name in names(fig2h_palettes)) {
    make_fig2h(fig2h_palettes[[pal_name]], pal_name)
  }
}

# --- Boxplots of key individual genes ---
message("  Generating key gene boxplots...")

key_genes <- c("Arid1a", "Tbx21", "Tcf7", "Prdm1", "Klrg1",
               "Il7r", "Gzmb", "Id2", "Id3", "Bach2")

key_matched <- symbol_to_ensembl(key_genes)

if (nrow(key_matched) >= 3) {
  # Get normalized counts for these genes
  key_expr <- as.data.frame(assay(vsd)[key_matched$gene_id, , drop = FALSE])
  key_expr$gene_id <- rownames(key_expr)
  key_expr <- key_expr %>%
    left_join(key_matched %>% select(gene_id, gene_name), by = "gene_id") %>%
    pivot_longer(cols = -c(gene_id, gene_name), names_to = "sample", values_to = "vst_expr") %>%
    left_join(
      as.data.frame(colData(vsd)) %>%
        rownames_to_column("sample") %>%
        select(sample, genotype, timepoint, cell_subset),
      by = "sample"
    ) %>%
    mutate(gene_name = factor(gene_name, levels = key_genes[key_genes %in% key_matched$gene_name]))

  p_box <- ggplot(key_expr, aes(x = genotype, y = vst_expr, fill = genotype)) +
    geom_boxplot(outlier.size = 0.5, width = 0.7) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
    facet_grid(gene_name ~ timepoint + cell_subset, scales = "free_y") +
    scale_fill_manual(values = pal_genotype_fill) +
    labs(title = "Key Gene Expression", y = "VST expression", x = NULL) +
    theme_paper +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text.y = element_text(face = "italic", size = 8))
  save_figure(p_box, "boxplot_key_genes", width = 16, height = 20, dir = figdir)
}

# =============================================================================
# 8. GSEA (fgsea) — Hallmark + ImmuneSigDB
# =============================================================================

message("=== Section 8: GSEA ===")

# Load MSigDB gene sets for mouse
# msigdbr v10+ uses collection/subcollection instead of category/subcategory
msig_hallmark <- msigdbr(species = "Mus musculus", collection = "H") %>%
  select(gs_name, gene_symbol)
msig_immunesig <- msigdbr(species = "Mus musculus", collection = "C7", subcollection = "IMMUNESIGDB") %>%
  select(gs_name, gene_symbol)

pathways_hallmark <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)
pathways_immune  <- split(msig_immunesig$gene_symbol, msig_immunesig$gs_name)

run_gsea <- function(res, name) {
  # Build ranked gene list using DESeq2 stat, mapped to gene symbols
  ranked_df <- res %>%
    filter(!is.na(stat), !is.na(gene_name), gene_name != "") %>%
    group_by(gene_name) %>%
    slice_max(abs(stat), n = 1) %>%
    ungroup() %>%
    arrange(desc(stat))

  ranks <- setNames(ranked_df$stat, ranked_df$gene_name)

  if (length(ranks) < 100) {
    message(sprintf("  %s: only %d ranked genes, skipping GSEA", name, length(ranks)))
    return(NULL)
  }

  # fgsea on Hallmark
  gsea_h <- fgsea(pathways = pathways_hallmark, stats = ranks,
                   nPermSimple = 10000, eps = 0) %>%
    as_tibble() %>%
    mutate(collection = "Hallmark", comparison = name)

  # fgsea on ImmuneSigDB
  gsea_i <- fgsea(pathways = pathways_immune, stats = ranks,
                   nPermSimple = 10000, eps = 0) %>%
    as_tibble() %>%
    mutate(collection = "ImmuneSigDB", comparison = name)

  gsea_all <- bind_rows(gsea_h, gsea_i) %>% arrange(padj)
  write_csv(gsea_all, file.path(outdir, paste0("gsea_", name, ".csv")))

  n_sig <- sum(gsea_all$padj < 0.05, na.rm = TRUE)
  message(sprintf("  %s: %d significant gene sets (padj<0.05)", name, n_sig))

  return(gsea_all)
}

gsea_results <- map2(de_results, names(de_results), run_gsea)
gsea_results <- compact(gsea_results)  # remove NULLs

# =============================================================================
# 9. GSEA Summary Heatmap — top Hallmark pathways across contrasts
# =============================================================================

message("=== Section 9: GSEA Hallmark Heatmap ===")

if (length(gsea_results) > 0) {
  gsea_all <- bind_rows(gsea_results)

  # Top Hallmark pathways: significant in >= 1 contrast
  top_hallmark <- gsea_all %>%
    filter(collection == "Hallmark", padj < 0.05) %>%
    distinct(pathway) %>%
    pull(pathway)

  if (length(top_hallmark) > 0) {
    # Limit to top 30 by mean |NES|
    if (length(top_hallmark) > 30) {
      top_hallmark <- gsea_all %>%
        filter(collection == "Hallmark", pathway %in% top_hallmark) %>%
        group_by(pathway) %>%
        summarize(mean_nes = mean(abs(NES), na.rm = TRUE), .groups = "drop") %>%
        slice_max(mean_nes, n = 30) %>%
        pull(pathway)
    }

    nes_mat <- gsea_all %>%
      filter(collection == "Hallmark", pathway %in% top_hallmark) %>%
      select(pathway, comparison, NES) %>%
      pivot_wider(names_from = comparison, values_from = NES, values_fill = 0) %>%
      column_to_rownames("pathway") %>%
      as.matrix()

    pval_mat <- gsea_all %>%
      filter(collection == "Hallmark", pathway %in% top_hallmark) %>%
      select(pathway, comparison, padj) %>%
      pivot_wider(names_from = comparison, values_from = padj, values_fill = 1) %>%
      column_to_rownames("pathway") %>%
      as.matrix()

    # Ensure same row/col order
    nes_mat <- nes_mat[rownames(pval_mat), colnames(pval_mat), drop = FALSE]

    sig_mat <- ifelse(pval_mat < 0.001, "***",
               ifelse(pval_mat < 0.01, "**",
               ifelse(pval_mat < 0.05, "*", "")))

    # Clean pathway names for display
    row_labels_gsea <- gsub("^HALLMARK_", "", rownames(nes_mat))
    row_labels_gsea <- gsub("_", " ", row_labels_gsea)

    ht_gsea <- Heatmap(nes_mat,
      name = "NES",
      col = col_nes,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 8))
      },
      cluster_columns = FALSE,
      row_labels = row_labels_gsea,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 9),
      column_title = "GSEA: Hallmark Pathways",
      row_names_max_width = unit(15, "cm")
    )

    save_heatmap(ht_gsea, "gsea_hallmark_heatmap", width = 12, height = max(8, length(top_hallmark) * 0.35), dir = figdir)
  }
}

# =============================================================================
# 10. GO Over-representation Analysis
# =============================================================================

message("=== Section 10: GO Enrichment ===")

run_go <- function(res, name) {
  # Strip Ensembl version for org.Mm.eg.db lookup
  up_genes   <- res %>% filter(sig == "Up") %>% pull(gene_id) %>% gsub("\\..*", "", .)
  down_genes <- res %>% filter(sig == "Down") %>% pull(gene_id) %>% gsub("\\..*", "", .)
  bg_genes   <- res %>% pull(gene_id) %>% gsub("\\..*", "", .)

  if (length(up_genes) >= 5) {
    go_up <- tryCatch(
      enrichGO(gene = up_genes, universe = bg_genes,
               OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
               ont = "BP", pvalueCutoff = 0.05, readable = TRUE),
      error = function(e) { message("  GO up error: ", e$message); NULL }
    )
    if (!is.null(go_up) && nrow(go_up@result) > 0) {
      write_csv(as.data.frame(go_up), file.path(outdir, paste0("go_up_", name, ".csv")))
      message(sprintf("  %s GO up: %d terms", name, sum(go_up@result$p.adjust < 0.05)))
    }
  }

  if (length(down_genes) >= 5) {
    go_down <- tryCatch(
      enrichGO(gene = down_genes, universe = bg_genes,
               OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
               ont = "BP", pvalueCutoff = 0.05, readable = TRUE),
      error = function(e) { message("  GO down error: ", e$message); NULL }
    )
    if (!is.null(go_down) && nrow(go_down@result) > 0) {
      write_csv(as.data.frame(go_down), file.path(outdir, paste0("go_down_", name, ".csv")))
      message(sprintf("  %s GO down: %d terms", name, sum(go_down@result$p.adjust < 0.05)))
    }
  }
}

walk2(de_results, names(de_results), run_go)

# =============================================================================
# 11. DE Summary Table & Barplot
# =============================================================================

message("=== Section 11: DE Summary ===")

de_summary <- map_dfr(names(de_results), function(name) {
  res <- de_results[[name]]
  tibble(
    contrast = name,
    n_up     = sum(res$sig == "Up", na.rm = TRUE),
    n_down   = sum(res$sig == "Down", na.rm = TRUE),
    n_total  = sum(res$sig != "NS", na.rm = TRUE),
    n_tested = sum(!is.na(res$padj))
  )
})

write_csv(de_summary, file.path(outdir, "de_summary.csv"))
message("  DE summary:")
print(de_summary)

# Barplot of DEG counts
de_bar <- de_summary %>%
  select(contrast, Up = n_up, Down = n_down) %>%
  pivot_longer(cols = c(Up, Down), names_to = "direction", values_to = "count") %>%
  mutate(
    count_signed = ifelse(direction == "Down", -count, count),
    contrast = factor(contrast, levels = rev(de_summary$contrast))
  )

p_bar <- ggplot(de_bar, aes(x = contrast, y = count_signed, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "#D62728", "Down" = "#2166AC")) +
  labs(title = "Differentially Expressed Genes per Contrast",
       subtitle = "|LFC| > 1, padj < 0.05",
       x = NULL, y = "Number of DEGs", fill = "Direction") +
  theme_paper
save_figure(p_bar, "barplot_de_summary", width = 10, height = 6, dir = figdir)

# =============================================================================
# Save session
# =============================================================================

message("=== Saving R session ===")
save.image(file.path(outdir, "rnaseq_analysis.RData"))

message("\n============================================")
message("RNA-seq analysis complete!")
message("  Tables: ", outdir)
message("  Figures: ", figdir)
message("============================================")
