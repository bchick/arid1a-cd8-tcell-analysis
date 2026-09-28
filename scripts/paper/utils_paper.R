#!/usr/bin/env Rscript
# =============================================================================
# paper/utils_paper.R — Helpers for the panel-by-panel paper reproduction
# McDonald, Chick et al. 2023 Immunity 56:1303 — shared helpers
#
# Paper thresholds, output paths and helpers shared by scripts/paper/fig*.R
# (save_panel(), write_panel_table(), require_inputs(), read_da(),
# read_ocr_clusters()). Every panel script writes figures to figures/paper/
# and the numbers behind each panel to results/paper/, so panels can be
# checked without re-plotting.
#
# Inputs:  results/atac/differential/da_*.csv, results/paper/fig1a_ocr_clusters.csv
#          (read on demand by the helpers)
# Outputs: creates figures/paper/ and results/paper/
# Sourced by: scripts/paper/fig*.R, after utils.R, from the repository root:
#               source("scripts/utils.R"); source("scripts/paper/utils_paper.R")
# =============================================================================

stopifnot(exists("paths"), exists("save_figure"), exists("save_heatmap"))

paths$paper_fig <- file.path(PROJECT_DIR, "figures/paper")
paths$paper_tab <- file.path(PROJECT_DIR, "results/paper")
invisible(lapply(c(paths$paper_fig, paths$paper_tab), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

# -----------------------------------------------------------------------------
# Thresholds as stated in the paper's legends / STAR Methods
# -----------------------------------------------------------------------------
PAPER <- list(
  atac_lfc  = 1,       # ATAC DA: 2-fold change
  atac_padj = 0.05,    # FDR < 0.05 (Benjamini-Hochberg)
  rna_lfc   = 1,       # Fig 2G/3F legends: >2-fold change
  rna_padj  = 0.05,
  tbet_lfc  = 1,       # Fig S5A: fold change >= 2 ...
  tbet_padj = 0.01,    # ... adjusted p < 0.01
  motif_window = 200,  # HOMER findMotifsGenome.pl -size 200 (100 bp of centre)
  profile_window = 1000  # deepTools +/- 1 kb around peak centres
)

# Fig 1 OCR cluster order used by every downstream panel
OCR_CLUSTERS <- c("Conserved", "Naive", "Early Activation", "Activation", "Late Activation")

# -----------------------------------------------------------------------------
# I/O helpers
# -----------------------------------------------------------------------------

#' Save a ggplot panel as figures/paper/<id>.{pdf,png}
save_panel <- function(p, id, width = 4, height = 4) {
  save_figure(p, id, width = width, height = height, dir = paths$paper_fig)
}

#' Save a ComplexHeatmap panel as figures/paper/<id>.{pdf,png}
save_panel_heatmap <- function(ht, id, width = 4, height = 4) {
  save_heatmap(ht, id, width = width, height = height, dir = paths$paper_fig)
}

#' Write the table behind a panel to results/paper/<id>.csv
write_panel_table <- function(df, id) {
  f <- file.path(paths$paper_tab, paste0(id, ".csv"))
  readr::write_csv(df, f)
  message("Wrote: ", f)
  invisible(f)
}

#' Stop with a clear message when an input that is not in the data bundle
#' (BAM, bigWig, genome FASTA) is missing. Returns TRUE when all exist.
require_inputs <- function(files, what = "input") {
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    message("Skipping: ", what, " not found (needs the upstream/raw-data mode):\n  ",
            paste(missing, collapse = "\n  "))
    return(FALSE)
  }
  TRUE
}

#' Read a DA table (results/atac/differential/da_<contrast>.csv) and flag
#' lost/gained OCRs with the paper's thresholds.
read_da <- function(contrast, lfc = PAPER$atac_lfc, padj = PAPER$atac_padj) {
  f <- file.path(paths$atac, "differential", paste0("da_", contrast, ".csv"))
  readr::read_csv(f, show_col_types = FALSE) |>
    dplyr::mutate(direction = dplyr::case_when(
      !is.na(padj) & padj < !!padj & log2FoldChange <= -lfc ~ "lost",
      !is.na(padj) & padj < !!padj & log2FoldChange >=  lfc ~ "gained",
      TRUE ~ "ns"))
}

#' The Fig 1A OCR cluster assignment (written by scripts/paper/fig1_ocr_clusters.R)
read_ocr_clusters <- function() {
  f <- file.path(paths$paper_tab, "fig1a_ocr_clusters.csv")
  if (!file.exists(f)) stop("Run scripts/paper/fig1_ocr_clusters.R first (", f, ")")
  readr::read_csv(f, show_col_types = FALSE) |>
    dplyr::mutate(cluster = factor(cluster, levels = OCR_CLUSTERS))
}
