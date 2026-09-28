#!/usr/bin/env Rscript
# =============================================================================
# utils_figures.R — shared figure helpers for the extended analyses
# McDonald, Chick et al. 2023 reanalysis
#
# Layered ON TOP of utils.R. Always source in this order, and re-source BOTH
# after any load(*.RData) (RData sessions overwrite utils.R palettes/helpers):
#
#   source("scripts/utils.R")
#   source("extended_analysis/scripts/utils_figures.R")
#
# Provides: theme_ext, save_ext_figure / save_ext_heatmap (PDF + 300 dpi PNG to
# the extended-analysis summary figure directory), ordered-factor helpers, the
# dose-response class palette, and small composition helpers.
# =============================================================================

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(ggrepel)
})

stopifnot(exists("theme_paper"), exists("save_figure"), exists("paths"))  # utils.R must be loaded first

# =============================================================================
# Font.
# Arial is not installed on this server. Use "Helvetica": it is resolvable by
# BOTH the base PostScript pdf() device (used by ComplexHeatmap via
# save_heatmap) and cairo/fontconfig (used by ggplot's cairo_pdf), where it
# maps to Nimbus Sans (URW Helvetica clone). A single family name keeps ggplot
# and heatmap text visually identical and avoids "font not found" errors that
# occur when a fontconfig-only family (e.g. "Nimbus Sans") hits the PS device.
# =============================================================================

FIG_FONT <- "Helvetica"
message("utils_figures.R: figure font family = '", FIG_FONT, "'")

#' gpar with the project figure font (for ComplexHeatmap text)
gp_font <- function(size = 10, ...) grid::gpar(fontfamily = FIG_FONT, fontsize = size, ...)

# =============================================================================
# Output directories
# =============================================================================

paths$summary_figures <- file.path(PROJECT_DIR, "figures/extended_analysis")
dir.create(paths$summary_figures, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# Themes
# =============================================================================

# The publication theme, with the font fixed to FIG_FONT
# (theme_paper hardcodes "Arial"; child elements inherit family from `text`).
theme_ext <- theme_paper + theme(text = element_text(family = FIG_FONT))


theme_set(theme_ext)

# =============================================================================
# Savers — PDF + 300 dpi PNG
# =============================================================================

#' Save a ggplot to the summary figure directory
save_ext_figure <- function(p, filename, width = 7, height = 5) {
  save_figure(p, filename, width = width, height = height, dir = paths$summary_figures)
}

#' Save a ComplexHeatmap to the summary figure directory
save_ext_heatmap <- function(ht, filename, width = 7, height = 5) {
  save_heatmap(ht, filename, width = width, height = height, dir = paths$summary_figures)
}


# =============================================================================
# Ordered-factor helpers (consistent axis/legend ordering everywhere)
# =============================================================================

gfac <- function(x) factor(as.character(x), levels = c("WT", "Het", "KO", "TbetKO"))
sfac <- function(x) factor(as.character(x), levels = c("Naive", "TE", "EEC", "MP"))
tfac <- function(x) factor(as.character(x), levels = c("Naive", "48h", "D3", "D5", "D8"))

# =============================================================================
# Dose-response class palette & ordering (shared by all dose-response figures)
# Okabe-Ito based, colorblind-safe. Story emphasis:
#   buffered = dose-tolerant (blue), haploinsufficient = dose-sensitive (vermilion),
#   linear = rare direct targets (highlight purple).
# =============================================================================

pal_doseclass <- c(
  "insensitive"       = "#BBBBBB",
  "buffered"          = "#0072B2",
  "haploinsufficient" = "#D55E00",
  "linear"            = "#CC79A7",
  "nonmonotonic"      = "#F0E442",
  "other_responsive"  = "#56B4E9"
)
dose_class_levels <- names(pal_doseclass)
dfac <- function(x) factor(as.character(x), levels = dose_class_levels)

# Direction of accessibility change (matches col_lfc orientation)
pal_direction <- c("lost" = "#2166AC", "gained" = "#B2182B", "ns" = "#BBBBBB")

# TF-family palette for motif figures (stable colors across the dose-response and TF-binding figures)
pal_tf_family <- c(
  "ETS"       = "#0072B2", "RUNX"     = "#009E73", "T-box"   = "#56B4E9",
  "AP-1/bZIP" = "#D55E00", "NFkB"     = "#CC79A7", "KLF/SP"  = "#E69F00",
  "TCF/LEF"   = "#117733", "bHLH"     = "#882255", "GATA"    = "#44AA99",
  "IRF/STAT"  = "#999933", "EGR"      = "#AA4499", "NR"      = "#DDCC77",
  "Zf"        = "#888888", "Homeo"    = "#661100", "other"   = "#DDDDDD"
)

# =============================================================================
# Composition helpers
# =============================================================================

#' Add bold A/B/C... panel tags to a patchwork composition
tag_panels <- function(pw, levels = "A") {
  pw + plot_annotation(tag_levels = levels) &
    theme(plot.tag = element_text(face = "bold", size = 14, family = FIG_FONT))
}


message("Loaded utils_figures.R — themes, savers, palettes ready")
