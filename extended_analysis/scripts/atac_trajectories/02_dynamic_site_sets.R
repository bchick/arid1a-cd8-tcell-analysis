#!/usr/bin/env Rscript
# =============================================================================
# atac_trajectories/02_dynamic_site_sets.R — which dynamic sites differ between WT and KO?
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Between the two passes of extended_analysis/scripts/atac_trajectories/01_timecourse_patterns.sh.
#
# Pass 1 called dynamic peaks separately in the WT and KO timecourses (each
# arm: LRT over time padj < 0.01 and range >= 0.5, VST scale). The sites that
# differ are the difference between those two sets:
#   * Lost   = dynamic in WT, not in KO
#   * Gained = dynamic in KO, not in WT
#   * Shared = dynamic in both
# "Not dynamic in KO" can mean flat in KO or short of a threshold, so each lost
# site is split by why it failed in KO: "Flat in KO" (KO range < 0.5, the
# dynamics are gone) or "Sub-threshold in KO" (range >= 0.5 but LRT padj >=
# 0.01, still moving but not called).
#
# Then writes pass-2 inputs: the WT timecourse restricted to the lost sites,
# with no further selection, so timecourse-patterns clusters every lost site
# by its WT trajectory.
#
# Inputs:  results/extended_analysis/atac_trajectories/pass1_wt_ko/results/differential/
# Outputs: results/extended_analysis/atac_trajectories/dynamic_sites/
#          results/extended_analysis/atac_trajectories/pass2_lost_in_ko/{inputs,config.yaml}
#          figures/extended_analysis/atac_trajectories/fig0_dynamic_site_sets.{pdf,png}
# Usage:   Rscript extended_analysis/scripts/atac_trajectories/02_dynamic_site_sets.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("extended_analysis/scripts/utils_trajectories.R")

MIN_RANGE <- 0.5   # timecourse-patterns differential.min_range
FDR       <- 0.01  # timecourse-patterns differential.fdr

setdir <- file.path(traj_dir, "dynamic_sites")
figdir <- file.path(paths$ext_figures, "atac_trajectories")
dir.create(setdir, recursive = TRUE, showWarnings = FALSE)

p1 <- file.path(traj_dir, "pass1_wt_ko/results/differential")
read_arm <- function(a) read_tsv(file.path(p1, paste0(a, "_results.tsv.gz")), show_col_types = FALSE) |>
  select(feature_id, baseMean, padj, range, dynamic)
sites <- full_join(read_arm("WT"), read_arm("KO"), by = "feature_id", suffix = c("_wt", "_ko")) |>
  mutate(set = case_when(dynamic_wt & dynamic_ko ~ "Shared",
                         dynamic_wt              ~ "Lost",
                         dynamic_ko              ~ "Gained",
                         TRUE                    ~ "Static"),
         ko_reason = case_when(set != "Lost" ~ NA_character_,
                               range_ko < MIN_RANGE ~ "Flat in KO",
                               TRUE ~ "Sub-threshold in KO"))

summary_tbl <- sites |> count(set, ko_reason, name = "n")
write_tsv(sites, file.path(setdir, "dynamic_site_sets.tsv.gz"))
write_tsv(summary_tbl, file.path(setdir, "dynamic_site_set_summary.tsv"))
message(sprintf("=== Dynamic: WT %d, KO %d | Shared %d, Lost %d, Gained %d ===",
                sum(sites$dynamic_wt, na.rm = TRUE), sum(sites$dynamic_ko, na.rm = TRUE),
                sum(sites$set == "Shared"), sum(sites$set == "Lost"), sum(sites$set == "Gained")))
print(summary_tbl)

# --- Fig 0: set sizes ---------------------------------------------------------
set_plot <- sites |>
  filter(set != "Static") |>
  mutate(set = factor(set, c("Shared", "Lost", "Gained")),
         part = factor(coalesce(ko_reason, as.character(set)),
                       c("Shared", "Flat in KO", "Sub-threshold in KO", "Gained"))) |>
  count(set, part)
p_sets <- ggplot(set_plot, aes(set, n, fill = part)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::comma(n)), position = position_stack(vjust = 0.5),
            size = 2.6, colour = "white") +
  scale_fill_manual(values = c(Shared = "grey45", `Flat in KO` = "#2CA02C",
                               `Sub-threshold in KO` = "#98D98E", Gained = "#4A90D9")) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Dynamic peaks", fill = NULL,
       title = "Dynamic ATAC peaks in the WT and ARID1A-KO timecourses",
       subtitle = "Lost = dynamic in WT only; Gained = dynamic in KO only") +
  theme_bw(base_size = 9) + theme(legend.position = "right")
save_figure(p_sets, "fig0_dynamic_site_sets", width = 5, height = 3.5, dir = figdir)

# =============================================================================
# Pass-2 inputs: WT timecourse of the lost sites, clustered without selection
# =============================================================================

lost <- sites$feature_id[sites$set == "Lost"]
p2 <- file.path(traj_dir, "pass2_lost_in_ko")
indir <- file.path(p2, "inputs")
dir.create(indir, recursive = TRUE, showWarnings = FALSE)
p2 <- normalizePath(p2); indir <- normalizePath(indir)

count_mat <- read_atac_consensus_counts()
sheet <- atac_trajectory_libraries(colnames(count_mat), genotypes = "WT")
cnt <- count_mat[lost, sheet$sample]
write_tsv(tibble(feature = rownames(cnt), as_tibble(cnt)), file.path(indir, "counts.tsv"))
write_tsv(sheet, file.path(indir, "samplesheet.tsv"))
bed <- attr(count_mat, "bed")
write_tsv(bed[match(lost, bed$name), ], file.path(indir, "features.bed"), col_names = FALSE)
message(sprintf("  Pass 2: %d lost sites x %d WT libraries", nrow(cnt), ncol(cnt)))

# Every lost site already passed both gates in WT; do not re-select. The count
# prefilter is off: these sites passed it in pass 1 over all libraries.
write_tcp_config(file.path(p2, "config.yaml"), indir, file.path(p2, "results"),
                 title = "WT trajectories of ATAC peaks that lose their dynamics in ARID1A KO",
                 extra = list(filter = list(min_mean_counts = 0),
                              differential = list(method = "none", top_n = length(lost),
                                                  min_range = 0)))
