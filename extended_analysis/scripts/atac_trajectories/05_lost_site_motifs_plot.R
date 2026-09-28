#!/usr/bin/env Rscript
# =============================================================================
# atac_trajectories/05_lost_site_motifs_plot.R — summarize SEA motif enrichment of lost sites
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Collects every SEA run from atac_trajectories/04_lost_site_motifs.sh and draws:
#   fig4  motifs of each lost-site trajectory class against static peaks
#   fig5  within each WT class, motifs enriched in sites that lose their
#          dynamics in KO (right) or keep them (left), as log2 enrichment
#   fig6  promoter share of lost vs kept sites per WT class. The motifs that
#          separate kept from lost opening sites include GC-rich CpG-promoter
#          motifs (ZBED4, NRF1, ZBTB33), so the annotation is shown alongside
#          to separate a promoter effect from a factor effect.
# Significance: SEA q < 1e-5, as in the timecourse-patterns T cell example.
#
# Inputs:  results/extended_analysis/atac_trajectories/motifs/{comparisons.tsv,sea/,beds/}
#          results/atac/differential/consensus_peaks_annotated.csv (core/02_atacseq_analysis.R)
# Outputs: results/extended_analysis/atac_trajectories/motifs/{sea_all.tsv.gz,promoter_share.tsv}
#          figures/extended_analysis/atac_trajectories/motifs/
# Usage:   Rscript extended_analysis/scripts/atac_trajectories/05_lost_site_motifs_plot.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("extended_analysis/scripts/utils_trajectories.R")

Q_MAX <- 1e-5
TOP_N <- 8

motif_dir <- file.path(traj_dir, "motifs")
figdir    <- file.path(paths$ext_figures, "atac_trajectories/motifs")

cmp <- read_tsv(file.path(motif_dir, "comparisons.tsv"), show_col_types = FALSE)
class_levels <- c("Decreasing", "Transient", "Transient Increasing",
                  "Sustained Increasing", "Late Increasing")

read_sea <- function(comparison, primary) {
  f <- file.path(motif_dir, "sea", paste0(comparison, "__", primary), "sea.tsv")
  if (!file.exists(f)) return(NULL)
  read_tsv(f, comment = "#", show_col_types = FALSE) |>
    filter(!is.na(QVALUE)) |>
    transmute(motif_id = ID, motif = toupper(ALT_ID), tp_pct = `TP%`, fp_pct = `FP%`,
              enr_ratio = ENR_RATIO, qvalue = QVALUE)
}
sea <- cmp |>
  mutate(res = purrr::map2(comparison, primary, read_sea)) |>
  tidyr::unnest(res)
write_tsv(sea, file.path(motif_dir, "sea_all.tsv.gz"))
message(sprintf("  %d SEA rows from %d comparisons", nrow(sea), n_distinct(paste(sea$comparison, sea$primary))))

# Several JASPAR matrices share a TF name; keep the strongest per TF
best_per_tf <- function(d) d |> group_by(across(-c(motif_id, tp_pct, fp_pct, enr_ratio, qvalue))) |>
  slice_min(qvalue, n = 1, with_ties = FALSE) |> ungroup()

# --- Fig 4: lost-site classes vs static peaks ---------------------------------
vs_static <- sea |> filter(comparison == "vs_static", qvalue < Q_MAX) |>
  select(class, motif, motif_id, tp_pct, fp_pct, enr_ratio, qvalue) |> best_per_tf()
top4 <- vs_static |> group_by(class) |> slice_max(enr_ratio, n = TOP_N, with_ties = FALSE) |>
  ungroup() |> pull(motif) |> unique()
d4 <- vs_static |> filter(motif %in% top4) |>
  mutate(class = factor(class, class_levels),
         motif = factor(motif, rev(top4)))
p4 <- ggplot(d4, aes(class, motif, size = -log10(pmax(qvalue, 1e-300)), colour = log2(enr_ratio))) +
  geom_point() +
  scale_colour_viridis_c(option = "mako", direction = -1, name = "log2 enrichment") +
  scale_size_area(max_size = 5, name = "-log10 q") +
  labs(x = NULL, y = NULL,
       title = "Motifs of sites that lose their dynamics in ARID1A KO",
       subtitle = sprintf("Each WT trajectory class vs static peaks; SEA q < %g, top %d per class",
                          Q_MAX, TOP_N)) +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_figure(p4, "fig4_lost_class_motifs_vs_static", width = 5.8,
            height = 0.14 * length(top4) + 1.8, dir = figdir)

# --- Fig 5: lost vs kept within each WT class -----------------------------------
lk <- sea |> filter(comparison %in% c("lost_vs_kept", "kept_vs_lost"), qvalue < Q_MAX) |>
  select(comparison, class, motif, motif_id, tp_pct, fp_pct, enr_ratio, qvalue) |> best_per_tf() |>
  mutate(side = if_else(comparison == "lost_vs_kept", "Lost in KO", "Kept in KO"),
         log2_enr = if_else(side == "Lost in KO", 1, -1) * log2(enr_ratio))
top5 <- lk |> group_by(class, side) |> slice_max(abs(log2_enr), n = 5, with_ties = FALSE) |> ungroup()
d5 <- top5 |>
  mutate(class = factor(class, class_levels),
         label = reorder(paste(motif, class, sep = "___"), log2_enr))  # order within facet
p5 <- ggplot(d5, aes(log2_enr, label, fill = side)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  facet_wrap(~ class, scales = "free_y", nrow = 1) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(values = c(`Lost in KO` = "#2CA02C", `Kept in KO` = "#000000")) +
  labs(x = "log2 enrichment (lost vs kept)", y = NULL, fill = NULL,
       title = "What separates ARID1A-dependent from ARID1A-independent trajectory sites",
       subtitle = sprintf("Within each WT class: sites losing vs keeping their dynamics in KO; SEA q < %g, top 5 per side",
                          Q_MAX)) +
  theme_bw(base_size = 9) + theme(legend.position = "top", strip.background = element_blank())
save_figure(p5, "fig5_lost_vs_kept_motifs", width = 13, height = 3.8, dir = figdir)

# --- Fig 6: promoter share -------------------------------------------------------
anno <- read_csv(file.path(paths$atac, "differential/consensus_peaks_annotated.csv"),
                 show_col_types = FALSE) |>
  transmute(feature_id = peak_id, promoter = grepl("^Promoter", annotation))
set_ids <- function(set) read_tsv(file.path(motif_dir, "beds", paste0(set, ".bed")),
                                  col_names = c("chr", "start", "end", "feature_id"),
                                  show_col_types = FALSE)$feature_id
prom <- bind_rows(
  tibble(class = "Static peaks", side = "Static", feature_id = set_ids("static")),
  bind_rows(lapply(unique(cmp$class[cmp$comparison == "lost_vs_kept"]), function(k) {
    s <- gsub("[^a-z0-9]+", "_", tolower(k))
    bind_rows(tibble(class = k, side = "Lost in KO", feature_id = set_ids(paste0("p1lost_", s))),
              tibble(class = k, side = "Kept in KO", feature_id = set_ids(paste0("p1kept_", s))))
  }))) |>
  left_join(anno, by = "feature_id") |>
  group_by(class, side) |>
  summarise(n = n(), pct_promoter = 100 * mean(promoter, na.rm = TRUE), .groups = "drop")
write_tsv(prom, file.path(motif_dir, "promoter_share.tsv"))
print(as.data.frame(prom))

static_pct <- prom$pct_promoter[prom$side == "Static"]
p6 <- ggplot(filter(prom, side != "Static") |> mutate(class = factor(class, class_levels)),
             aes(class, pct_promoter, fill = side)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_hline(yintercept = static_pct, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c(`Lost in KO` = "#2CA02C", `Kept in KO` = "#000000")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "% of peaks at promoters (ChIPseeker)", fill = NULL,
       title = "Promoter share of sites that lose vs keep their dynamics in KO",
       subtitle = sprintf("Dashed line: static peaks (%.0f%%)", static_pct)) +
  theme_bw(base_size = 9) + theme(legend.position = "top",
                                  axis.text.x = element_text(angle = 30, hjust = 1))
save_figure(p6, "fig6_promoter_share_lost_vs_kept", width = 5.5, height = 3.8, dir = figdir)

message("=== atac_trajectories/05_lost_site_motifs_plot.R complete ===")
