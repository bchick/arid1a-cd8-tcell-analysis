#!/usr/bin/env Rscript
# =============================================================================
# paper/fig4_fig5_profiles.R — average signal profiles for Figures 4C, 5C, 5E and 5J
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# Reads the deepTools matrices written by scripts/paper/fig4_fig5_signal.sh
# and draws mean-coverage profiles in the paper's layout; the per-bin means are
# written to results/paper/ for checking.
#
#   4C  ATAC coverage, WT / Het / KO, one panel per subset (all lost OCRs)
#   5C  CUT&RUN coverage, d5 WT vs Arid1a KO, at ARID1A-dependent OCRs
#   5E  ATAC and T-bet ChIP coverage per treatment at ACBI1-dependent / -independent OCRs
#   5J  T-bet CUT&RUN, WT / Arid1a KO x EV / T-bet-OE, published Activation / Late Activation OCRs
#
# Inputs:  results/paper/signal/*.gz
# Outputs: figures/paper/fig4c_atac_profiles, fig5c_cutrun_profiles,
#          fig5e_inhibitor_profiles, fig5j_tbet_cutrun_profiles (.pdf/.png)
#          results/paper/<panel>.csv
# Usage:   Rscript scripts/paper/fig4_fig5_profiles.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")

suppressPackageStartupMessages(library(jsonlite))
select <- dplyr::select
filter <- dplyr::filter

mat_dir <- file.path(paths$paper_tab, "signal")

#' Mean signal per bin for each (sample, region group) of a deepTools matrix
read_profile <- function(name) {
  f <- file.path(mat_dir, paste0(name, ".gz"))
  if (!require_inputs(f, paste("deepTools matrix", name))) return(NULL)
  hdr <- fromJSON(sub("^@", "", readLines(f, n = 1)))
  m <- as.matrix(read.delim(f, header = FALSE, skip = 1)[, -(1:6)])
  bins <- hdr$sample_boundaries
  groups <- hdr$group_boundaries
  bin_mid <- seq(-hdr$upstream[1], hdr$downstream[1] - hdr$`bin size`[1],
                 by = hdr$`bin size`[1]) + hdr$`bin size`[1] / 2
  purrr::map_dfr(seq_along(hdr$sample_labels), function(i) {
    cols <- (bins[i] + 1):bins[i + 1]
    purrr::map_dfr(seq_along(hdr$group_labels), function(j) {
      rows <- (groups[j] + 1):groups[j + 1]
      tibble(sample = hdr$sample_labels[i], group = sub("\\.bed$", "", hdr$group_labels[j]),
             n_regions = length(rows), position = bin_mid,
             mean = colMeans(m[rows, cols, drop = FALSE], na.rm = TRUE))
    })
  })
}

profile_theme <- theme(legend.position = "right", strip.text = element_text(size = 9),
                       legend.key.width = unit(0.6, "cm"))
kb_axis <- scale_x_continuous(breaks = c(-1000, 0, 1000), labels = c("-1", "0", "+1"))

# =============================================================================
# 4C profiles: WT / Het / KO per subset over all lost OCR groups combined
# =============================================================================

p4c <- read_profile("fig4c_atac_heatmap")
if (!is.null(p4c)) {
  prof <- p4c |>
    separate(sample, c("subset", "genotype"), sep = "_", remove = FALSE) |>
    group_by(subset, genotype, position) |>
    summarise(mean = weighted.mean(mean, n_regions), .groups = "drop") |>
    mutate(subset = factor(subset, levels = c("TE", "EEC", "MP")),
           genotype = factor(genotype, levels = c("WT", "Het", "KO")))
  write_panel_table(prof, "fig4c_atac_profiles")
  p <- ggplot(prof, aes(position, mean, color = genotype)) +
    geom_line(linewidth = 0.7) + facet_wrap(~ subset, nrow = 1) + kb_axis +
    scale_color_manual(values = pal_genotype[c("WT", "Het", "KO")],
                       labels = c("WT", "Arid1a cHet", "Arid1a cKO")) +
    labs(x = "Distance from ATAC peak (kb)", y = "ATAC coverage", color = NULL) + profile_theme
  save_panel(p, "fig4c_atac_profiles", width = 7, height = 2.4)
}

# =============================================================================
# 5C profiles: WT vs KO per antibody at ARID1A-dependent OCRs
# =============================================================================

p5c <- read_profile("fig5c_cutrun_heatmap")
if (!is.null(p5c)) {
  prof <- p5c |>
    separate(sample, c("antibody", "genotype"), sep = "_", remove = FALSE) |>
    mutate(antibody = factor(recode(antibody, Tbet = "T-bet"),
                             levels = c("ARID1A", "BATF", "ETS1", "T-bet")),
           genotype = factor(genotype, levels = c("WT", "KO")))
  write_panel_table(prof, "fig5c_cutrun_profiles")
  p <- ggplot(filter(prof, group == "fig5c_arid1a_dependent"), aes(position, mean, color = genotype)) +
    geom_line(linewidth = 0.7) + facet_wrap(~ antibody, nrow = 1, scales = "free_y") + kb_axis +
    scale_color_manual(values = pal_genotype[c("WT", "KO")], labels = c("d5 WT", "d5 Arid1a KO")) +
    labs(x = "Distance from OCR centre (kb) at ARID1A-dependent sites",
         y = "CUT&RUN coverage", color = NULL) + profile_theme
  save_panel(p, "fig5c_cutrun_profiles", width = 7.5, height = 2.4)
}

# =============================================================================
# 5E histograms: ATAC and T-bet ChIP per treatment
# =============================================================================

trt_levels <- c("DMSO", "DMSO+IL-12", "ACBI1+IL-12", "BRM014+IL-12")
pal_trt <- setNames(unname(pal_treatment[c("DMSO", "DMSO_IL12", "ACBI1_IL12", "BRM014_IL12")]),
                    trt_levels)
p5e <- purrr::compact(list(ATAC = read_profile("fig5d_atac_heatmap"),
                           `T-bet ChIP` = read_profile("fig5d_tbet_chip_heatmap")))
if (length(p5e)) {
  # The nf-core bigWigs are depth-scaled (1e6 / mapped fragments), i.e. the same
  # total-read normalization as the primary 5D-F split in fig5.R, so they are
  # plotted as is (a median-of-ratios rescaling would hide the global FRiP drop).
  prof <- bind_rows(p5e, .id = "assay") |>
    mutate(sample = factor(sample, levels = trt_levels),
           group = factor(recode(group, fig5d_acbi1_dependent = "ACBI1-dependent",
                                 fig5d_acbi1_independent = "ACBI1-independent"),
                          levels = c("ACBI1-dependent", "ACBI1-independent")))
  write_panel_table(prof, "fig5e_inhibitor_profiles")
  p <- ggplot(prof, aes(position, mean, color = sample)) +
    geom_line(linewidth = 0.7) +
    facet_grid(group ~ assay, scales = "free_y") + kb_axis +
    scale_color_manual(values = pal_trt) +
    labs(x = "Distance from OCR centre (kb)", y = "Coverage", color = NULL) + profile_theme
  save_panel(p, "fig5e_inhibitor_profiles", width = 5.5, height = 4)
}

# =============================================================================
# 5J: T-bet CUT&RUN at Activation / Late Activation OCRs
# =============================================================================

p5j <- read_profile("fig5j_tbet_cutrun")
if (!is.null(p5j)) {
  prof <- p5j |>
    mutate(genotype = ifelse(grepl("_KO_", sample), "Arid1a cKO", "WT"),
           vector   = ifelse(grepl("TbetOE", sample), "Tbet-OE", "EV"),
           condition = factor(paste(genotype, "+", vector),
                              levels = c("WT + EV", "Arid1a cKO + EV", "WT + Tbet-OE", "Arid1a cKO + Tbet-OE")),
           group = factor(recode(group, activation.specific.sig = "Activation ATAC cluster",
                                 late.activation.specific.sig = "Late Activation ATAC cluster")))
  write_panel_table(prof, "fig5j_tbet_cutrun_profiles")
  p <- ggplot(prof, aes(position, mean, color = condition, linetype = condition)) +
    geom_line(linewidth = 0.7) + facet_wrap(~ group, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = c(-1000, 0, 1000)) +
    scale_color_manual(values = c("black", pal_genotype[["KO"]], "black", pal_genotype[["KO"]])) +
    scale_linetype_manual(values = c("solid", "solid", "dashed", "dashed")) +
    labs(x = "Distance from ATAC peak centre (bp)", y = "T-bet C&R coverage",
         color = NULL, linetype = NULL) + profile_theme
  save_panel(p, "fig5j_tbet_cutrun_profiles", width = 6.5, height = 2.6)
}

message("=== fig4_fig5_profiles.R done ===")
