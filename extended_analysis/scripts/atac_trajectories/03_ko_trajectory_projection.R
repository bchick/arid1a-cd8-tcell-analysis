#!/usr/bin/env Rscript
# =============================================================================
# atac_trajectories/03_ko_trajectory_projection.R — WT trajectory classes, and the same peaks in KO
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Projects the ARID1A-KO trajectory of each classified peak onto the WT class
# centroids, for two sets of WT classes from atac_trajectories/01_timecourse_patterns.sh:
#   all_dynamic    pass 1, every peak dynamic in the WT timecourse
#   lost_in_ko     pass 2, the peaks dynamic in WT but not in KO
#                  (atac_trajectories/02_dynamic_site_sets.R)
#
#   * WT and KO libraries are variance-stabilized together (blind VST), so
#     both trajectories are on one scale. Naive WT is the shared D0 (there are
#     no Naive KO libraries); D8 is TE only (Exp1 + Exp2 in both genotypes).
#   * Per-timepoint means are z-scored within genotype and correlated with the
#     WT class centroids; the best match is the KO class.
#   * KO range (max - min of timepoint means) below the workflow's dynamic gate
#     (0.5, VST) is "Static": the dynamics are gone rather than reshaped. Best
#     correlation < 0.6 (the workflow's assignment floor) is "Unassigned".
#   * Control: WT profiles re-assigned the same way should recover the
#     workflow's labels; the agreement rate is reported per set.
#
# Inputs:  results/atac/.../consensus_peaks.mRp.clN.featureCounts.txt (nf-core consensus counts)
#          results/extended_analysis/atac_trajectories/{pass1_wt_ko,pass2_lost_in_ko}/results/clusters/WT_clusters.tsv
# Outputs: results/extended_analysis/atac_trajectories/ko_projection/<set>/
#          figures/extended_analysis/atac_trajectories/<set>/
# Usage:   Rscript extended_analysis/scripts/atac_trajectories/03_ko_trajectory_projection.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("extended_analysis/scripts/utils_trajectories.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggalluvial)
})

MIN_RANGE <- 0.5   # timecourse-patterns differential.min_range (VST scale)
MIN_COR   <- 0.6   # timecourse-patterns cluster.assign_min_cor

sets <- list(
  all_dynamic = list(clusters = file.path(traj_dir, "pass1_wt_ko/results/clusters/WT_clusters.tsv"),
                     what = "WT-dynamic peaks"),
  lost_in_ko  = list(clusters = file.path(traj_dir, "pass2_lost_in_ko/results/clusters/WT_clusters.tsv"),
                     what = "peaks dynamic in WT but not in KO")
)

# =============================================================================
# 1. Joint VST of the WT and KO trajectory libraries (shared by both sets)
# =============================================================================

count_mat <- read_atac_consensus_counts()
libs <- atac_trajectory_libraries(colnames(count_mat))
print(count(libs, group, time))

std_chr <- paste0("chr", c(1:19, "X", "Y"))
cnt <- count_mat[attr(count_mat, "chr") %in% std_chr, libs$sample]
cnt <- cnt[rowMeans(cnt) >= 10, ]
dds <- DESeqDataSetFromMatrix(cnt, data.frame(row.names = libs$sample, time = factor(libs$time)),
                              design = ~ 1)
vst_mat <- assay(vst(dds, blind = TRUE))
message(sprintf("  Joint VST: %d peaks x %d libraries", nrow(vst_mat), ncol(vst_mat)))

tp_mean <- function(ids, geno) {
  m <- sapply(c(0, 3, 5, 8), function(t) {
    s <- libs$sample[libs$time == t & libs$group == (if (t == 0) "shared" else geno)]
    rowMeans(vst_mat[ids, s, drop = FALSE])
  })
  dimnames(m) <- list(ids, paste0("D", c(0, 3, 5, 8))); m
}

fate_pal <- c(Retained = "#000000", Redistributed = "#4A90D9",
              `Static in KO` = "#2CA02C", Unassigned = "grey85")

# =============================================================================
# 2. Projection + figures, per set
# =============================================================================

project_set <- function(set_name, spec) {
  message(sprintf("\n=== %s: %s ===", set_name, spec$what))
  outdir <- file.path(traj_dir, "ko_projection", set_name)
  figdir <- file.path(paths$ext_figures, "atac_trajectories", set_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  wt_cl <- read_tsv(spec$clusters, show_col_types = FALSE) |>
    filter(supercluster_label != "Unassigned")
  missing <- setdiff(wt_cl$feature_id, rownames(vst_mat))
  if (length(missing))
    message(sprintf("  %d classified peaks fail the joint count filter; dropped", length(missing)))
  wt_cl <- filter(wt_cl, feature_id %in% rownames(vst_mat))
  class_levels <- names(sort(table(wt_cl$supercluster_label), decreasing = TRUE))
  print(table(factor(wt_cl$supercluster_label, class_levels)))

  traj <- list(WT = tp_mean(wt_cl$feature_id, "WT"), KO = tp_mean(wt_cl$feature_id, "KO"))
  zs   <- lapply(traj, function(m) t(scale(t(m))))
  centroids <- sapply(class_levels, function(k)
    colMeans(zs$WT[wt_cl$supercluster_label == k, , drop = FALSE]))  # 4 x classes

  assign_class <- function(z, m) {
    r <- cor(t(z), centroids)
    best   <- class_levels[max.col(r, ties.method = "first")]
    best_r <- apply(r, 1, max)
    rng    <- apply(m, 1, function(x) max(x) - min(x))
    cls <- ifelse(rng < MIN_RANGE, "Static",
                  ifelse(is.na(best_r) | best_r < MIN_COR, "Unassigned", best))
    tibble(class = cls, cor = best_r, range = rng)
  }
  wt_re <- assign_class(zs$WT, traj$WT)
  ko_as <- assign_class(zs$KO, traj$KO)

  proj <- wt_cl |>
    select(feature_id, chr, start, end, wt_class = supercluster_label) |>
    mutate(wt_reassigned = wt_re$class, wt_range = wt_re$range,
           ko_class = ko_as$class, ko_cor = ko_as$cor, ko_range = ko_as$range,
           fate = case_when(ko_class == wt_class   ~ "Retained",
                            ko_class == "Static"     ~ "Static in KO",
                            ko_class == "Unassigned" ~ "Unassigned",
                            TRUE                     ~ "Redistributed"))
  agree <- mean(proj$wt_reassigned == proj$wt_class)
  message(sprintf("  Control: WT re-assignment agrees with workflow labels for %.1f%%", 100 * agree))

  transitions <- proj |>
    count(wt_class, ko_class, name = "n") |>
    group_by(wt_class) |> mutate(frac = n / sum(n)) |> ungroup() |>
    arrange(factor(wt_class, class_levels), desc(n))
  fate_summary <- proj |>
    count(wt_class, fate, name = "n") |>
    group_by(wt_class) |> mutate(frac = n / sum(n)) |> ungroup()
  write_tsv(proj, file.path(outdir, "wt_class_ko_projection.tsv.gz"))
  write_tsv(transitions, file.path(outdir, "wt_to_ko_transitions.tsv"))
  write_tsv(fate_summary, file.path(outdir, "ko_fate_by_wt_class.tsv"))
  write_tsv(tibble(metric = "wt_reassignment_agreement", value = agree),
            file.path(outdir, "projection_control.tsv"))
  print(tidyr::pivot_wider(select(fate_summary, -frac), names_from = fate,
                           values_from = n, values_fill = 0))

  # --- Fig 1: mean trajectory per WT class, WT vs KO (change from Naive) ----
  n_lab <- count(proj, wt_class) |>
    mutate(label = sprintf("%s (n = %s)", wt_class, scales::comma(n)))
  traj_long <- bind_rows(lapply(names(traj), function(g)
    as_tibble(traj[[g]] - traj[[g]][, "D0"]) |>
      mutate(wt_class = wt_cl$supercluster_label, genotype = g))) |>
    tidyr::pivot_longer(starts_with("D"), names_to = "tp", values_to = "delta") |>
    mutate(day = as.integer(sub("D", "", tp))) |>
    group_by(wt_class, genotype, day) |>
    summarise(mean = mean(delta), lo = quantile(delta, 0.25), hi = quantile(delta, 0.75),
              .groups = "drop") |>
    mutate(genotype = factor(genotype, c("WT", "KO")),
           panel = factor(n_lab$label[match(wt_class, n_lab$wt_class)],
                          n_lab$label[match(class_levels, n_lab$wt_class)]))
  p_traj <- ggplot(traj_long, aes(day, mean, colour = genotype, fill = genotype)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = pal_genotype[c("WT", "KO")]) +
    scale_fill_manual(values = pal_genotype[c("WT", "KO")]) +
    scale_x_continuous(breaks = c(0, 3, 5, 8), labels = c("Naive", "D3", "D5", "D8 TE")) +
    labs(x = NULL, y = "Accessibility vs Naive (VST)", colour = NULL, fill = NULL,
         title = sprintf("WT trajectory classes of %s, and the same peaks in ARID1A KO", spec$what),
         subtitle = "Mean and IQR across peaks; Naive WT is the shared baseline") +
    theme_bw(base_size = 9) + theme(legend.position = "top", strip.background = element_blank())
  save_figure(p_traj, "fig1_wt_class_trajectories_wt_vs_ko",
              width = 2.2 * length(class_levels) + 1, height = 3.2, dir = figdir)

  # --- Fig 2: KO fate of each WT class -------------------------------------
  p_fate <- ggplot(mutate(fate_summary, wt_class = factor(wt_class, rev(class_levels)),
                          fate = factor(fate, names(fate_pal))),
                   aes(frac, wt_class, fill = fate)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = fate_pal, drop = FALSE) +
    scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = "Share of WT-class peaks", y = NULL, fill = "In KO",
         title = sprintf("Fate in ARID1A KO: %s", spec$what)) +
    theme_bw(base_size = 9) + theme(legend.position = "top")
  save_figure(p_fate, "fig2_ko_fate_by_wt_class", width = 6,
              height = 0.45 * length(class_levels) + 1.5, dir = figdir)

  # --- Fig 3: alluvial WT class -> KO class --------------------------------
  ko_levels <- c(class_levels, "Static", "Unassigned")
  class_pal <- setNames(c(viridisLite::mako(length(class_levels), begin = 0.15, end = 0.85),
                          "grey70", "grey90"), ko_levels)
  p_alluv <- ggplot(mutate(transitions, wt_class = factor(wt_class, class_levels),
                           ko_class = factor(ko_class, ko_levels)),
                    aes(axis1 = wt_class, axis2 = ko_class, y = n)) +
    geom_alluvium(aes(fill = wt_class), width = 1/4, alpha = 0.75) +
    geom_stratum(width = 1/4, fill = "grey95", colour = "grey30") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6) +
    scale_x_discrete(limits = c("WT class", "KO class"), expand = c(0.12, 0.05)) +
    scale_fill_manual(values = class_pal, guide = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Peaks", title = "WT trajectory class -> class in ARID1A KO") +
    theme_bw(base_size = 9) + theme(panel.grid = element_blank())
  save_figure(p_alluv, "fig3_alluvial_wt_to_ko_class", width = 6, height = 5, dir = figdir)
}

for (nm in names(sets)) project_set(nm, sets[[nm]])
message("\n=== atac_trajectories/03_ko_trajectory_projection.R complete ===")
