#!/usr/bin/env Rscript
# =============================================================================
# paper/rna_helpers.R — Shared RNA-seq inputs for the paper panels (Fig 2G-I, Fig 3)
# McDonald, Chick et al. 2023 Immunity 56:1303 — shared helpers
#
# Reads only compact tables, so the panels do not depend on the large
# rnaseq_analysis.RData session. Source after utils.R and utils_paper.R.
#
# Inputs:  results/rnaseq/star_salmon/salmon.merged.gene_counts.tsv  (nf-core/rnaseq)
#          results/rnaseq/differential/de_<contrast>.csv             (core/01_rnaseq_analysis.R)
# Outputs: none (returns data frames / matrices to the caller)
# Sourced by: scripts/paper/fig2g_2i.R, scripts/paper/fig3.R
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename

rna_dir <- file.path(paths$rnaseq, "differential")

#' VST matrix of all RNA-seq samples, built exactly as in core/01_rnaseq_analysis.R
#' (Salmon gene counts, >= 10 counts in >= 3 samples, vst(blind = FALSE) on ~1).
#' Returns list(vst = matrix, meta = data.frame, gene_map = data.frame).
load_rna_vst <- function() {
  f <- file.path(paths$rnaseq, "star_salmon/salmon.merged.gene_counts.tsv")
  stopifnot("Salmon gene counts not found" = file.exists(f))
  raw <- read.delim(f, check.names = FALSE)
  cm <- round(as.matrix(raw[, -(1:2)]))
  rownames(cm) <- raw$gene_id
  storage.mode(cm) <- "integer"

  meta <- load_sample_sheet() |>
    filter(assay == "rnaseq") |>
    as.data.frame()
  meta <- meta[match(colnames(cm), meta$sample_name), ]
  rownames(meta) <- meta$sample_name
  meta$genotype    <- factor(meta$genotype, levels = c("WT", "Het", "KO"))
  meta$timepoint   <- factor(meta$timepoint, levels = c("D3", "D5", "D8"))
  meta$cell_subset <- factor(meta$cell_subset, levels = c("total", "TE", "EEC", "MP"))

  dds <- DESeqDataSetFromMatrix(cm, meta, design = ~ 1)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  vsd <- vst(dds, blind = FALSE)

  list(vst = assay(vsd), meta = meta,
       gene_map = data.frame(gene_id = raw$gene_id, gene_name = raw$gene_name))
}

#' DE table from core/01_rnaseq_analysis.R (apeglm-shrunken LFC; Wald stat unshrunken)
read_de <- function(contrast) {
  readr::read_csv(file.path(rna_dir, paste0("de_", contrast, ".csv")),
                  show_col_types = FALSE)
}

#' Count DEGs at the legend threshold (2-fold) and the STAR Methods threshold
#' (log2FC 0.585), padj < 0.05.
count_degs <- function(res, contrast) {
  purrr::map_dfr(c(legend_2fold = 1, methods_1.5fold = 0.585), function(l) {
    tibble::tibble(
      contrast = contrast, lfc_threshold = l,
      up   = sum(res$padj < PAPER$rna_padj & res$log2FoldChange >=  l, na.rm = TRUE),
      down = sum(res$padj < PAPER$rna_padj & res$log2FoldChange <= -l, na.rm = TRUE))
  }, .id = "threshold")
}

#' Ranked gene vector (Wald stat, one entry per gene symbol) for fgsea
rank_stat <- function(res) {
  r <- res |>
    filter(!is.na(stat), !is.na(gene_name), gene_name != "") |>
    group_by(gene_name) |>
    slice_max(abs(stat), n = 1, with_ties = FALSE) |>
    ungroup()
  sort(setNames(r$stat, r$gene_name), decreasing = TRUE)
}

#' MSigDB gene sets as a named list (msigdbr >= 10 argument names)
msig_sets <- function(collection, subcollection = NULL, pattern = NULL) {
  m <- msigdbr(species = "Mus musculus", collection = collection,
               subcollection = subcollection)
  if (!is.null(pattern)) m <- filter(m, grepl(pattern, gs_name))
  split(m$gene_symbol, m$gs_name)
}

#' fgsea with the paper's 10,000 permutations (seeded for reproducibility)
run_fgsea <- function(res, pathways, seed = 42) {
  set.seed(seed)
  fgsea(pathways = pathways, stats = rank_stat(res), nPermSimple = 10000, eps = 0) |>
    tibble::as_tibble()
}

#' Display name for a Hallmark set: HALLMARK_MYC_TARGETS_V1 -> MYC_TARGETS_V1
hallmark_label <- function(x) sub("^HALLMARK_", "", x)
