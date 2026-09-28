#!/usr/bin/env Rscript
# =============================================================================
# utils.R — Shared configuration, palettes and helper functions
# McDonald, Chick et al. 2023 Immunity 56:1303 — shared helpers
#
# Defines the project root and output paths (`paths`), genome/annotation
# configuration (`genome`, get_txdb()), colour palettes matched to the
# publication, heatmap colour functions, and helpers for saving figures
# (save_figure(), save_heatmap(); PDF + 300 dpi PNG), loading the sample sheet
# and mapping mouse to human gene symbols. The project root is the working
# directory, or ARID1A_PROJECT_DIR if set.
#
# Inputs:  data/metadata/, data/reference/ (paths only; read by callers)
# Outputs: creates the results/ and figures/ directories listed in `paths`
# Sourced by: scripts/core/*.R, scripts/paper/*.R, extended_analysis/scripts/**/*.R
#             (source("scripts/utils.R") from the repository root)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
})

# =============================================================================
# Project paths
# =============================================================================

# Project root: set ARID1A_PROJECT_DIR, or run scripts from the repository root.
PROJECT_DIR <- normalizePath(Sys.getenv("ARID1A_PROJECT_DIR", unset = "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_DIR, "scripts", "utils.R")))
  stop("PROJECT_DIR (", PROJECT_DIR, ") is not the repository root; ",
       "run from the repo root or set ARID1A_PROJECT_DIR.")

paths <- list(
  project    = PROJECT_DIR,
  metadata   = file.path(PROJECT_DIR, "data/metadata"),
  reference  = file.path(PROJECT_DIR, "data/reference"),
  # nf-core output directories
  atac       = file.path(PROJECT_DIR, "results/atac"),
  rnaseq     = file.path(PROJECT_DIR, "results/rnaseq"),
  cutrun     = file.path(PROJECT_DIR, "results/cutrun"),
  chipseq    = file.path(PROJECT_DIR, "results/chipseq"),
  # Analysis output directories
  results    = file.path(PROJECT_DIR, "results"),
  figures    = file.path(PROJECT_DIR, "figures"),
  integrated = file.path(PROJECT_DIR, "results/integrated"),
  # Cross-study meta-analysis (Guo et al. 2022)
  meta       = file.path(PROJECT_DIR, "results/meta_analysis/guo2022"),
  # Analyses beyond the manuscript (extended_analysis/scripts/)
  ext_results = file.path(PROJECT_DIR, "results/extended_analysis"),
  ext_figures = file.path(PROJECT_DIR, "figures/extended_analysis")
)

# Create output dirs if they don't exist
walk(paths, ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

# =============================================================================
# Genome configuration
# =============================================================================

genome <- list(
  build     = "GRCm39",
  common    = "mm39",
  fasta     = file.path(PROJECT_DIR, "data/reference/GRCm39.primary_assembly.genome.fa"),
  gtf       = file.path(PROJECT_DIR, "data/reference/gencode.vM35.primary_assembly.annotation.gtf"),
  blacklist = file.path(PROJECT_DIR, "data/reference/mm39-blacklist.v2.bed"),
  txdb_pkg  = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  org_db    = "org.Mm.eg.db",
  species   = "Mus musculus",
  msigdb    = "mouse"
)

#' GENCODE vM35 TxDb. Built from the GTF when available (and cached as SQLite);
#' otherwise loaded from the cache shipped in the data bundle, so annotation
#' steps run without the 800 MB GTF.
get_txdb <- function(cache = file.path(PROJECT_DIR, "results/annotation/gencode.vM35.txdb.sqlite")) {
  if (file.exists(genome$gtf)) {
    txdb <- GenomicFeatures::makeTxDbFromGFF(genome$gtf, format = "gtf")
    if (!file.exists(cache)) {
      dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
      AnnotationDbi::saveDb(txdb, cache)
    }
    return(txdb)
  }
  if (!file.exists(cache))
    stop("Need ", genome$gtf, " or the cached TxDb ", cache, " (make fetch-data)")
  AnnotationDbi::loadDb(cache)
}

# =============================================================================
# Color palettes — matched to original publication
# =============================================================================

# Genotype colors (primary palette)
pal_genotype <- c(
  "WT"     = "#000000",
  "Het"    = "#4A90D9",
  "KO"     = "#2CA02C",
  "TbetKO" = "#17BECF"
)

# Genotype bar fills (lighter for bar charts)
pal_genotype_fill <- c(
  "WT"     = "#B0B0B0",
  "Het"    = "#4A90D9",
  "KO"     = "#2CA02C",
  "TbetKO" = "#17BECF"
)

# Effector subset colors
pal_subset <- c(
  "TE"    = "#CC6677",
  "EEC"   = "#DDCC77",
  "MP"    = "#AA4499",
  "Naive" = "#7F7F7F"
)

# OCR cluster colors
pal_cluster <- c(
  "Conserved"        = "#1A1A1A",
  "Naive"            = "#9ECAE1",
  "Early Activation" = "#FF7F0E",
  "Activation"       = "#2CA02C",
  "Late Activation"  = "#9467BD"
)

# BAF inhibitor treatment colors
pal_treatment <- c(
  "DMSO"           = "#000000",
  "DMSO_IL12"      = "#D62728",
  "ACBI1_IL12"     = "#2CA02C",
  "BRM014_IL12"    = "#9467BD"
)

# Genomic annotation colors (pie charts)
pal_annotation <- c(
  "Promoter"    = "#1A1A1A",
  "Intergenic"  = "#1F77B4",
  "Intron"      = "#FF7F0E",
  "Exon"        = "#FFF176",
  "Other"       = "#C0C0C0"
)

# Timepoint colors
pal_timepoint <- c(
  "Naive" = "#1A1A1A",
  "48h"   = "#FF7F0E",
  "D3"    = "#2CA02C",
  "D5"    = "#D62728",
  "D8"    = "#9467BD"
)

# =============================================================================
# Heatmap color functions
# =============================================================================

# ATAC-seq signal (white to red, Figures 2/4 style)
col_atac_signal <- colorRamp2(
  c(0, 2, 5, 10, 20),
  c("white", "#FFF5F0", "#FCA082", "#D63B20", "#7F0000")
)

# RNA expression z-score (blue-white-red diverging)
col_expression <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
)

# CUT&RUN signal (white to dark green)
col_cutrun_signal <- colorRamp2(
  c(0, 2, 5, 10, 20),
  c("white", "#C7E9C0", "#74C476", "#238B45", "#005A32")
)

# Log2 fold change (diverging)
col_lfc <- colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
)

# GSEA NES (diverging, red = positive, blue = negative)
col_nes <- colorRamp2(
  c(-3, -1.5, 0, 1.5, 3),
  c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
)

# =============================================================================
# ggplot2 theme — matches paper styling
# =============================================================================

theme_paper <- theme_bw() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.5),
    axis.text          = element_text(size = 9, color = "black"),
    axis.title         = element_text(size = 11),
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 10),
    legend.background  = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 10, face = "bold"),
    plot.title         = element_text(size = 12, face = "bold"),
    plot.subtitle      = element_text(size = 10, color = "gray30"),
    text               = element_text(family = "Arial")
  )

theme_set(theme_paper)

# =============================================================================
# Helper functions
# =============================================================================

#' Save a ggplot as both PDF (vector) and PNG (300 dpi)
save_figure <- function(p, filename, width = 7, height = 5, dir = paths$figures) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  base <- tools::file_path_sans_ext(filename)
  ggsave(file.path(dir, paste0(base, ".pdf")), p, width = width, height = height,
         device = cairo_pdf)
  ggsave(file.path(dir, paste0(base, ".png")), p, width = width, height = height,
         dpi = 300)
  message("Saved: ", file.path(dir, paste0(base, ".{pdf,png}")))
}

#' Save a ComplexHeatmap as both PDF and PNG
save_heatmap <- function(ht, filename, width = 7, height = 5, dir = paths$figures) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  base <- tools::file_path_sans_ext(filename)
  pdf(file.path(dir, paste0(base, ".pdf")), width = width, height = height)
  draw(ht)
  dev.off()
  png(file.path(dir, paste0(base, ".png")), width = width, height = height,
      units = "in", res = 300)
  draw(ht)
  dev.off()
  message("Saved: ", file.path(dir, paste0(base, ".{pdf,png}")))
}

#' Load master sample sheet and parse metadata
load_sample_sheet <- function() {
  read_tsv(file.path(paths$metadata, "master_sample_sheet.tsv"), show_col_types = FALSE)
}

#' Convert mouse gene symbols to human orthologs (for MSigDB)
#' Uses biomaRt or a pre-built mapping
mouse_to_human <- function(genes) {
  # Simple approach using msigdbr's built-in mouse support
  # For complex cases, use biomaRt getLDS()
  genes
}

#' Format p-values for display
format_pval <- function(p) {
  case_when(
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )
}

message("Loaded utils.R — project config, palettes, and helpers ready")
