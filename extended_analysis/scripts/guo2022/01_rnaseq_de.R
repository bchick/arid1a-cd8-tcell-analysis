#!/usr/bin/env Rscript
# =============================================================================
# guo2022/01_rnaseq_de.R — Guo et al. 2022 RNA-seq differential expression
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Cross-study comparison with Guo et al. 2022 (Nature): DESeq2 (apeglm
# shrinkage), fgsea (MSigDB Hallmark + ImmuneSigDB), EnhancedVolcano, ComplexHeatmap on the
# nf-core/rnaseq reprocessing of Guo's RNA-seq (mm39 / GENCODE vM35).
#
# Mirrors core/01_rnaseq_analysis.R (§1,3,4,5,8,11), adapted for Guo:
#   - No master sample sheet: metadata parsed from the Guo_<gse>_<cond>_r<N>
#     column names produced by build_meta_nfcore_samplesheets.py.
#   - Contrasts are WITHIN-GSE only (GSE is confounded with condition; the two
#     WT arms — g183615 genetic control vs g199184 vehicle — are NEVER pooled).
#   - run_de() generalized with a group_col param (the core version hardcodes
#     `genotype`).
#   - No pseudo-bulk section (Guo RNA is already bulk).
#
# Key test: does ARID1A loss in Guo replicate the ETS dose-buffered enhancer
# program found in McDonald? These KO-vs-WT DE tables feed the cross-study
# fold-change correlation.
#
# Inputs:  results/meta_analysis/guo2022/rnaseq/star_salmon/salmon.merged.gene_counts.tsv
# Outputs: results/extended_analysis/guo2022/rnaseq/differential/
#            de_<contrast>.csv, gsea_<contrast>.csv, de_summary.csv,
#            guo_normalized_counts.csv, guo_rnaseq_analysis.RData
#          figures/extended_analysis/guo2022/rnaseq/
# Usage:   Rscript extended_analysis/scripts/guo2022/01_rnaseq_de.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
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
outdir <- file.path(paths$ext_results, "guo2022/rnaseq/differential")
figdir <- file.path(paths$ext_figures, "guo2022/rnaseq")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# Guo-specific palettes (local; conditions differ from the main project schema)
pal_guo_cond <- c(
  "WT"    = "#000000",  # genetic control (g183615) / vehicle control (g199184)
  "KO"    = "#2CA02C",  # Arid1a knockout
  "MycWT" = "#7F7F7F",  # c-Myc arm control
  "MycKO" = "#1F77B4",  # c-Myc knockout
  "inhib" = "#D62728"   # cBAF inhibitor
)
pal_guo_gse <- c("g183615" = "#4A90D9", "g199184" = "#DDCC77")

# =============================================================================
# 1. Import merged Salmon gene counts + parse metadata from column names
# =============================================================================

message("=== Section 1: Importing Guo merged Salmon gene counts ===")

counts_file <- file.path(paths$meta, "rnaseq/star_salmon/salmon.merged.gene_counts.tsv")
stopifnot("Merged gene counts file not found" = file.exists(counts_file))

raw <- read.delim(counts_file, check.names = FALSE)
message(sprintf("  Loaded %d genes x %d sample columns", nrow(raw), ncol(raw) - 2))

# Build gene_id -> gene_name mapping (same handling as core/01_rnaseq_analysis.R)
gene_map <- data.frame(
  gene_id   = raw$gene_id,
  gene_name = raw$gene_name,
  stringsAsFactors = FALSE
)
dupes <- duplicated(gene_map$gene_name) | duplicated(gene_map$gene_name, fromLast = TRUE)
gene_map$gene_name_unique <- ifelse(
  dupes & gene_map$gene_name != "",
  paste0(gene_map$gene_name, "_", gsub("\\..*", "", gene_map$gene_id)),
  gene_map$gene_name
)
gene_map$gene_name_unique[gene_map$gene_name_unique == ""] <-
  gene_map$gene_id[gene_map$gene_name_unique == ""]

# Extract count matrix (round to integers for DESeq2)
count_mat <- as.matrix(raw[, -(1:2)])
rownames(count_mat) <- raw$gene_id
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

# --- Parse sample metadata from Guo_<gse>_<cond>_r<N> column names ---
sample_names <- colnames(count_mat)
parsed <- stringr::str_match(sample_names, "^Guo_(g\\d+)_(.+)_r(\\d+)$")
stopifnot("Unparseable Guo sample name(s)" = !any(is.na(parsed[, 1])))

meta <- data.frame(
  sample    = sample_names,
  gse       = parsed[, 2],
  condition = parsed[, 3],
  replicate = parsed[, 4],
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample
meta$gse       <- factor(meta$gse)
meta$condition <- factor(meta$condition)
meta$replicate <- factor(meta$replicate)
# Distinct arm label (keeps the two WT arms separable in QC)
meta$arm <- factor(paste(meta$gse, meta$condition, sep = "_"))

message(sprintf("  Sample metadata: %d samples", nrow(meta)))
message("  GSE x condition breakdown:")
print(table(meta$gse, meta$condition))

# Build DESeqDataSet with minimal design (we subset per contrast)
dds_full <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ 1
)

# Pre-filter: genes with >= 10 counts in >= 3 samples (as in core/01_rnaseq_analysis.R)
keep <- rowSums(counts(dds_full) >= 10) >= 3
dds_full <- dds_full[keep, ]
message(sprintf("  Kept %d / %d genes after pre-filtering", sum(keep), length(keep)))

# VST for visualization
vsd <- vst(dds_full, blind = FALSE)
message("  VST transformation complete")

# Save normalized counts
dds_norm <- DESeq(dds_full)
norm_counts <- counts(dds_norm, normalized = TRUE)
norm_out <- data.frame(
  gene_id   = rownames(norm_counts),
  gene_name = gene_map$gene_name[match(rownames(norm_counts), gene_map$gene_id)],
  norm_counts,
  check.names = FALSE
)
write_csv(norm_out, file.path(outdir, "guo_normalized_counts.csv"))

# =============================================================================
# 2. QC Visualization (PCA + sample-correlation heatmap)
# =============================================================================

message("=== Section 2: QC Visualization ===")

# --- PCA all samples ---
pca_mat <- assay(vsd)
rv <- rowVars(pca_mat)
top_var <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
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
    title = "Guo RNA-seq PCA: All Samples",
    x = sprintf("PC1 (%.1f%%)", pct_var[1]),
    y = sprintf("PC2 (%.1f%%)", pct_var[2]),
    color = "Condition", shape = "GSE"
  ) +
  theme_paper
save_figure(p_pca_all, "guo_rnaseq_pca_all", width = 8, height = 6, dir = figdir)

# --- PCA g183615 (genetic Arid1a + c-Myc arm) only ---
gen_vsd <- vsd[, colData(vsd)$gse == "g183615"]
if (ncol(gen_vsd) > 5) {
  pca_mat_g <- assay(gen_vsd)
  rv_g <- rowVars(pca_mat_g)
  top_g <- order(rv_g, decreasing = TRUE)[seq_len(min(500, length(rv_g)))]
  pca_g <- prcomp(t(pca_mat_g[top_g, ]), center = TRUE, scale. = FALSE)
  pca_g_df <- data.frame(
    PC1       = pca_g$x[, 1],
    PC2       = pca_g$x[, 2],
    condition = colData(gen_vsd)$condition,
    sample    = colnames(gen_vsd)
  )
  pct_g <- round(100 * summary(pca_g)$importance[2, 1:2], 1)

  p_pca_g <- ggplot(pca_g_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_color_manual(values = pal_guo_cond) +
    labs(
      title = "Guo RNA-seq PCA: g183615 (Arid1a + c-Myc)",
      x = sprintf("PC1 (%.1f%%)", pct_g[1]),
      y = sprintf("PC2 (%.1f%%)", pct_g[2]),
      color = "Condition"
    ) +
    theme_paper
  save_figure(p_pca_g, "guo_rnaseq_pca_g183615", width = 8, height = 6, dir = figdir)
}

# --- Sample-sample correlation heatmap ---
cor_mat <- cor(assay(vsd), method = "pearson")

ha_col <- HeatmapAnnotation(
  Condition = colData(vsd)$condition,
  GSE       = colData(vsd)$gse,
  col = list(
    Condition = pal_guo_cond[levels(colData(vsd)$condition)],
    GSE       = pal_guo_gse[levels(colData(vsd)$gse)]
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
  column_title = "Sample-sample Correlation (Guo RNA-seq)"
)
save_heatmap(ht_cor, "guo_rnaseq_sample_correlation", width = 10, height = 9, dir = figdir)

# =============================================================================
# 3. Differential Expression (DESeq2) — within-GSE, two-level contrasts
# =============================================================================

message("=== Section 3: Differential Expression ===")

# Helper: run DESeq2 on a subset of samples for a two-level contrast on an
# arbitrary grouping column. Generalizes run_de() in core/01_rnaseq_analysis.R (which hardcodes `genotype`).
run_de <- function(count_matrix, col_data, contrast_name, num, denom,
                   group_col = "condition") {
  # Subset to the two relevant levels of the grouping column
  cd <- col_data[col_data[[group_col]] %in% c(num, denom), , drop = FALSE]
  cd[[group_col]] <- droplevels(factor(cd[[group_col]], levels = c(denom, num)))

  # Restrict to globally pre-filtered genes, matching columns by sample name
  gene_rows <- intersect(rownames(dds_full), rownames(count_matrix))
  cm <- count_matrix[gene_rows, rownames(cd), drop = FALSE]

  design_formula <- as.formula(paste("~", group_col))
  dds_sub <- DESeqDataSetFromMatrix(countData = cm, colData = cd, design = design_formula)
  keep_sub <- rowSums(counts(dds_sub) >= 10) >= 3
  dds_sub <- dds_sub[keep_sub, ]
  dds_sub <- DESeq(dds_sub)

  # apeglm shrinkage on the num-vs-denom coefficient
  coef_name <- paste0(group_col, "_", num, "_vs_", denom)
  if (!(coef_name %in% resultsNames(dds_sub))) {
    coef_name <- grep(group_col, resultsNames(dds_sub), value = TRUE)[1]
  }
  res <- lfcShrink(dds_sub, coef = coef_name, type = "apeglm")

  # Wald stat from unshrunken results (apeglm drops `stat`; needed for GSEA)
  res_unshrunk <- results(dds_sub, name = coef_name)

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

# --- Arid1a KO vs WT (g183615, 5v5) — the ETS-replication test, primary ---
message("  Running g183615 KO vs WT (Arid1a)...")
g183615 <- meta[meta$gse == "g183615", ]
de_results[["KO_vs_WT"]] <- run_de(
  count_mat, g183615, "KO_vs_WT", num = "KO", denom = "WT"
)

# --- c-Myc KO vs WT (g183615, 10v10) — for the Myc-vs-Arid1a contrast ---
message("  Running g183615 MycKO vs MycWT (c-Myc)...")
de_results[["MycKO_vs_MycWT"]] <- run_de(
  count_mat, g183615, "MycKO_vs_MycWT", num = "MycKO", denom = "MycWT"
)

# --- Inhibitor vs vehicle (g199184, 4v4) — pharmacologic replication ---
message("  Running g199184 inhib vs WT (cBAF inhibitor)...")
g199184 <- meta[meta$gse == "g199184", ]
de_results[["inhib_vs_WT"]] <- run_de(
  count_mat, g199184, "inhib_vs_WT", num = "inhib", denom = "WT"
)

# Save all DE results (same schema as the main project -> drop-in for the cross-study comparison)
walk2(de_results, names(de_results), function(res, name) {
  write_csv(res, file.path(outdir, paste0("de_", name, ".csv")))
  n_up   <- sum(res$sig == "Up", na.rm = TRUE)
  n_down <- sum(res$sig == "Down", na.rm = TRUE)
  message(sprintf("  %s: %d up, %d down (|LFC|>1, padj<0.05)", name, n_up, n_down))
})

# =============================================================================
# 4. Volcano Plots
# =============================================================================

message("=== Section 4: Volcano Plots ===")

make_volcano <- function(res, title, filename) {
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
  make_volcano(res, paste("Guo:", gsub("_", " ", name)), paste0("volcano_", name))
})

# =============================================================================
# 5. GSEA (fgsea) — Hallmark + ImmuneSigDB
# =============================================================================

message("=== Section 5: GSEA ===")

# msigdbr v10+ uses collection/subcollection instead of category/subcategory
msig_hallmark <- msigdbr(species = "Mus musculus", collection = "H") %>%
  select(gs_name, gene_symbol)
msig_immunesig <- msigdbr(species = "Mus musculus", collection = "C7",
                          subcollection = "IMMUNESIGDB") %>%
  select(gs_name, gene_symbol)

pathways_hallmark <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)
pathways_immune   <- split(msig_immunesig$gene_symbol, msig_immunesig$gs_name)

run_gsea <- function(res, name) {
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

  gsea_h <- fgsea(pathways = pathways_hallmark, stats = ranks,
                  nPermSimple = 10000, eps = 0) %>%
    as_tibble() %>%
    mutate(collection = "Hallmark", comparison = name)

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
gsea_results <- compact(gsea_results)

# --- GSEA Hallmark summary heatmap across contrasts ---
if (length(gsea_results) > 0) {
  gsea_all <- bind_rows(gsea_results)

  top_hallmark <- gsea_all %>%
    filter(collection == "Hallmark", padj < 0.05) %>%
    distinct(pathway) %>%
    pull(pathway)

  if (length(top_hallmark) > 0) {
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

    nes_mat <- nes_mat[rownames(pval_mat), colnames(pval_mat), drop = FALSE]

    sig_mat <- ifelse(pval_mat < 0.001, "***",
               ifelse(pval_mat < 0.01, "**",
               ifelse(pval_mat < 0.05, "*", "")))

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
      column_title = "Guo GSEA: Hallmark Pathways",
      row_names_max_width = unit(15, "cm")
    )

    save_heatmap(ht_gsea, "guo_gsea_hallmark_heatmap",
                 width = 10, height = max(8, length(top_hallmark) * 0.35), dir = figdir)
  }
}

# =============================================================================
# 6. DE Summary Table & Barplot
# =============================================================================

message("=== Section 6: DE Summary ===")

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
  labs(title = "Guo: Differentially Expressed Genes per Contrast",
       subtitle = "|LFC| > 1, padj < 0.05",
       x = NULL, y = "Number of DEGs", fill = "Direction") +
  theme_paper
save_figure(p_bar, "guo_barplot_de_summary", width = 10, height = 5, dir = figdir)

# =============================================================================
# Save session
# =============================================================================

message("=== Saving R session ===")
save.image(file.path(outdir, "guo_rnaseq_analysis.RData"))

message("\n============================================")
message("Guo RNA-seq analysis complete!")
message("  Tables: ", outdir)
message("  Figures: ", figdir)
message("============================================")
