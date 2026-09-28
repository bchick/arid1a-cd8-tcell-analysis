#!/usr/bin/env Rscript
# =============================================================================
# temporal_clustering/01_degpatterns_wt_ko.R — degPatterns temporal clustering, WT vs KO
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Clusters ATAC-seq peaks by temporal dynamics (Naive -> D3 -> D5 -> D8 TE)
# separately for WT and ARID1A-KO with DEGreport::degPatterns, then compares
# which WT patterns are disrupted in KO (cross-tabulation, alluvial, dual
# heatmap, per-cluster retention).
#
# Inputs:  results/atac/differential/checkpoint_sections1to3.RData (main_vsd)
# Outputs: results/extended_analysis/temporal_clustering/wt_ko_degpatterns/ (tables + RData)
#          figures/extended_analysis/temporal_clustering/wt_ko_degpatterns/ (4 figure panels)
# Usage:   Rscript extended_analysis/scripts/temporal_clustering/01_degpatterns_wt_ko.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(DEGreport)
  library(ggalluvial)
  library(RColorBrewer)
  library(cowplot)
})

# Reassign dplyr verbs after Bioc loading
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
count  <- dplyr::count

outdir <- file.path(paths$ext_results, "temporal_clustering/wt_ko_degpatterns")
figdir <- file.path(paths$ext_figures, "temporal_clustering/wt_ko_degpatterns")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

tp_levels <- c("Naive", "D3", "D5", "D8")

# =============================================================================
# name_cluster() — label clusters by temporal profile
# (copied from core/03_atac_temporal_clustering.R)
# =============================================================================

name_cluster <- function(center, tp_names = tp_levels) {
  names(center) <- tp_names
  peak_tp <- which.max(center)
  trough_tp <- which.min(center)
  rng <- max(center) - min(center)

  if (rng < 0.8) return("Constitutive-like")

  # Monotonic patterns
  diffs <- diff(center)
  if (all(diffs > 0.15)) return("Progressive Opening")
  if (all(diffs < -0.15)) return("Progressive Closing")

  # Sharp peak patterns
  if (peak_tp == 1 && center[1] - mean(center[2:4]) > 0.5) return("Naive-specific")
  if (peak_tp == 2 && center[2] - mean(center[c(1,3,4)]) > 0.5) return("D3-specific")
  if (peak_tp == 4 && center[4] - mean(center[1:3]) > 0.5) return("D8-specific")

  # Two-phase patterns
  if (mean(center[1:2]) > 0.2 && mean(center[3:4]) < -0.2) return("Early (Naive+D3)")
  if (mean(center[1:2]) < -0.2 && mean(center[3:4]) > 0.2) return("Late (D5+D8)")
  if (center[1] < -0.3 && all(diffs >= -0.1)) return("Activation-induced")
  if (center[1] > 0.3 && all(diffs <= 0.1)) return("Activation-lost")

  return(paste0("Cluster_", tp_names[peak_tp]))
}

# =============================================================================
# 1. Load data
# =============================================================================

message("=== Loading main_vsd from checkpoint ===")

load(file.path(paths$results, "atac/differential/checkpoint_sections1to3.RData"))
vst_mat <- assay(main_vsd)
cd <- as.data.frame(colData(main_vsd))

message(sprintf("  VST matrix: %d peaks x %d samples", nrow(vst_mat), ncol(vst_mat)))

# Re-source utils.R to restore palettes after load() overwrites
source("scripts/utils.R")
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
count  <- dplyr::count

# =============================================================================
# 2. Identify samples for WT and KO trajectories
# =============================================================================

message("\n=== Identifying WT and KO trajectory samples ===")

# WT trajectory: Naive_WT, D3_WT, D5_WT, D8_WT (all subsets: TE+EEC+MP)
wt_naive <- grep("^Naive_WT_Rep", colnames(vst_mat), value = TRUE)
wt_d3    <- grep("^D3_WT_Rep", colnames(vst_mat), value = TRUE)
wt_d5    <- grep("^D5_WT_Rep", colnames(vst_mat), value = TRUE)
wt_d8_all <- grep("^D8_WT_(TE|EEC|MP)_", colnames(vst_mat), value = TRUE)

# KO trajectory: Naive_WT (shared), D3_KO, D5_KO, D8_KO (all subsets)
ko_d3     <- grep("^D3_KO_Rep", colnames(vst_mat), value = TRUE)
ko_d5     <- grep("^D5_KO_Rep", colnames(vst_mat), value = TRUE)
ko_d8_all <- grep("^D8_KO_(TE|EEC|MP)_", colnames(vst_mat), value = TRUE)

message(sprintf("  WT: Naive=%d, D3=%d, D5=%d, D8_all=%d",
                length(wt_naive), length(wt_d3), length(wt_d5), length(wt_d8_all)))
message(sprintf("  KO: Naive=%d (shared), D3=%d, D5=%d, D8_all=%d",
                length(wt_naive), length(ko_d3), length(ko_d5), length(ko_d8_all)))

# All relevant samples for peak selection
all_relevant <- c(wt_naive, wt_d3, wt_d5, wt_d8_all, ko_d3, ko_d5, ko_d8_all)
message(sprintf("  Total relevant samples: %d", length(all_relevant)))

# =============================================================================
# 3. Select variable peaks using timepoint-averaged MAD
# =============================================================================

message("\n=== Selecting variable peaks ===")

# Average across timepoint-genotype combinations to compute MAD
tp_geno_means <- cbind(
  Naive     = rowMeans(vst_mat[, wt_naive]),
  D3_WT     = rowMeans(vst_mat[, wt_d3]),
  D5_WT     = rowMeans(vst_mat[, wt_d5]),
  D8_WT     = rowMeans(vst_mat[, wt_d8_all]),
  D3_KO     = rowMeans(vst_mat[, ko_d3]),
  D5_KO     = rowMeans(vst_mat[, ko_d5]),
  D8_KO     = rowMeans(vst_mat[, ko_d8_all])
)

peak_mad <- apply(tp_geno_means, 1, mad)

# Take top 5K most variable peaks — degPatterns uses O(n^2) hierarchical
# clustering; 30K takes days. 5K runs in ~minutes and captures the most
# dynamic peaks well. Increase if needed after confirming runtime.
n_top <- 5000
variable_peaks <- names(sort(peak_mad, decreasing = TRUE))[1:min(n_top, sum(peak_mad > 0))]

message(sprintf("  Top %d variable peaks selected (MAD range: %.3f - %.3f)",
                length(variable_peaks),
                min(peak_mad[variable_peaks]), max(peak_mad[variable_peaks])))

# =============================================================================
# 4. Build balanced matrices (N x 8 for each genotype)
# =============================================================================

message("\n=== Building balanced matrices ===")

# --- WT balanced matrix ---
# D8 WT: pseudobulk all subsets (TE+EEC+MP) by experiment
wt_d8_exp1 <- grep("D8_WT_.*_Exp1", wt_d8_all, value = TRUE)
wt_d8_exp2 <- grep("D8_WT_.*_Exp2", wt_d8_all, value = TRUE)

wt_bal <- cbind(
  vst_mat[variable_peaks, wt_naive],
  vst_mat[variable_peaks, wt_d3],
  vst_mat[variable_peaks, wt_d5],
  D8_WT_Exp1 = rowMeans(vst_mat[variable_peaks, wt_d8_exp1]),
  D8_WT_Exp2 = rowMeans(vst_mat[variable_peaks, wt_d8_exp2])
)

wt_meta <- data.frame(
  timepoint = factor(
    c(rep("Naive", 2), rep("D3", 2), rep("D5", 2), rep("D8", 2)),
    levels = tp_levels
  ),
  row.names = colnames(wt_bal)
)

message("  WT balanced matrix:")
message(sprintf("    %d peaks x %d columns", nrow(wt_bal), ncol(wt_bal)))
message(sprintf("    Columns: %s", paste(colnames(wt_bal), collapse = ", ")))
print(table(wt_meta$timepoint))

# --- KO balanced matrix ---
# D8 KO: pseudobulk all subsets (TE+EEC+MP) by experiment
ko_d8_exp1 <- grep("D8_KO_.*_Exp1", ko_d8_all, value = TRUE)
ko_d8_exp2 <- grep("D8_KO_.*_Exp2", ko_d8_all, value = TRUE)

ko_bal <- cbind(
  vst_mat[variable_peaks, wt_naive],  # shared Naive_WT baseline
  vst_mat[variable_peaks, ko_d3],
  vst_mat[variable_peaks, ko_d5],
  D8_KO_Exp1 = rowMeans(vst_mat[variable_peaks, ko_d8_exp1]),
  D8_KO_Exp2 = rowMeans(vst_mat[variable_peaks, ko_d8_exp2])
)

ko_meta <- data.frame(
  timepoint = factor(
    c(rep("Naive", 2), rep("D3", 2), rep("D5", 2), rep("D8", 2)),
    levels = tp_levels
  ),
  row.names = colnames(ko_bal)
)

message("  KO balanced matrix:")
message(sprintf("    %d peaks x %d columns", nrow(ko_bal), ncol(ko_bal)))
message(sprintf("    Columns: %s", paste(colnames(ko_bal), collapse = ", ")))
print(table(ko_meta$timepoint))

# =============================================================================
# 5. Run degPatterns — WT
# =============================================================================

message("\n=== Running degPatterns for WT ===")

old_theme <- theme_get()
theme_set(theme_bw())

pdf(file.path(figdir, "degpatterns_wt_internal.pdf"), width = 10, height = 8)
deg_wt <- tryCatch(
  degPatterns(wt_bal, metadata = wt_meta, time = "timepoint",
              minc = 100, reduce = TRUE, cutoff = 0.7, scale = TRUE),
  error = function(e) {
    message("  WARNING: degPatterns WT failed: ", e$message)
    NULL
  }
)
dev.off()

theme_set(old_theme)

if (is.null(deg_wt)) {
  stop("degPatterns failed for WT — cannot proceed with comparison")
}

# Extract WT clusters
deg_wt_df <- deg_wt$normalized   # $df has only genes+cluster; $normalized has timepoint+value
deg_wt_col <- "cluster"

# Detect gene column name (varies by DEGreport version: "genes" or "id" etc.)
wt_gene_col <- intersect(c("genes", "id", "gene"), colnames(deg_wt_df))[1]
if (is.na(wt_gene_col)) stop("Cannot find gene column in deg_wt$normalized. Columns: ",
                              paste(colnames(deg_wt_df), collapse=", "))
message("  WT normalized columns: ", paste(colnames(deg_wt_df), collapse=", "))
message("  Using gene column: ", wt_gene_col)

wt_clusters <- deg_wt_df %>%
  select(all_of(c(wt_gene_col, deg_wt_col))) %>%
  distinct() %>%
  setNames(c("peak_id", "wt_cluster_num"))

message(sprintf("  WT: %d peaks clustered into %d clusters",
                nrow(wt_clusters), n_distinct(wt_clusters$wt_cluster_num)))

# Name WT clusters
wt_cluster_centers <- deg_wt_df %>%
  group_by(.data[[deg_wt_col]], timepoint) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  setNames(c("cluster_num", "timepoint", "mean_val"))

wt_cluster_nums <- sort(unique(wt_cluster_centers$cluster_num))
wt_name_map <- sapply(wt_cluster_nums, function(cl) {
  vals <- wt_cluster_centers %>%
    filter(cluster_num == cl) %>%
    arrange(factor(timepoint, levels = tp_levels))
  name_cluster(vals$mean_val)
})

# Deduplicate names
name_tab <- table(wt_name_map)
for (dn in names(name_tab[name_tab > 1])) {
  idx <- which(wt_name_map == dn)
  for (j in seq_along(idx)) wt_name_map[idx[j]] <- paste0(dn, " ", j)
}
names(wt_name_map) <- as.character(wt_cluster_nums)

wt_clusters$wt_cluster <- wt_name_map[as.character(wt_clusters$wt_cluster_num)]

message("  WT cluster assignments:")
print(table(wt_clusters$wt_cluster))

# =============================================================================
# 6. Run degPatterns — KO
# =============================================================================

message("\n=== Running degPatterns for KO ===")

theme_set(theme_bw())

pdf(file.path(figdir, "degpatterns_ko_internal.pdf"), width = 10, height = 8)
deg_ko <- tryCatch(
  degPatterns(ko_bal, metadata = ko_meta, time = "timepoint",
              minc = 100, reduce = TRUE, cutoff = 0.7, scale = TRUE),
  error = function(e) {
    message("  WARNING: degPatterns KO failed: ", e$message)
    NULL
  }
)
dev.off()

theme_set(old_theme)

if (is.null(deg_ko)) {
  stop("degPatterns failed for KO — cannot proceed with comparison")
}

# Extract KO clusters
deg_ko_df <- deg_ko$normalized   # $df has only genes+cluster; $normalized has timepoint+value
deg_ko_col <- "cluster"

# Detect gene column name
ko_gene_col <- intersect(c("genes", "id", "gene"), colnames(deg_ko_df))[1]
if (is.na(ko_gene_col)) stop("Cannot find gene column in deg_ko$normalized. Columns: ",
                              paste(colnames(deg_ko_df), collapse=", "))
message("  Using gene column: ", ko_gene_col)

ko_clusters <- deg_ko_df %>%
  select(all_of(c(ko_gene_col, deg_ko_col))) %>%
  distinct() %>%
  setNames(c("peak_id", "ko_cluster_num"))

message(sprintf("  KO: %d peaks clustered into %d clusters",
                nrow(ko_clusters), n_distinct(ko_clusters$ko_cluster_num)))

# Name KO clusters
ko_cluster_centers <- deg_ko_df %>%
  group_by(.data[[deg_ko_col]], timepoint) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  setNames(c("cluster_num", "timepoint", "mean_val"))

ko_cluster_nums <- sort(unique(ko_cluster_centers$cluster_num))
ko_name_map <- sapply(ko_cluster_nums, function(cl) {
  vals <- ko_cluster_centers %>%
    filter(cluster_num == cl) %>%
    arrange(factor(timepoint, levels = tp_levels))
  name_cluster(vals$mean_val)
})

# Deduplicate names
name_tab <- table(ko_name_map)
for (dn in names(name_tab[name_tab > 1])) {
  idx <- which(ko_name_map == dn)
  for (j in seq_along(idx)) ko_name_map[idx[j]] <- paste0(dn, " ", j)
}
names(ko_name_map) <- as.character(ko_cluster_nums)

ko_clusters$ko_cluster <- ko_name_map[as.character(ko_clusters$ko_cluster_num)]

message("  KO cluster assignments:")
print(table(ko_clusters$ko_cluster))

# =============================================================================
# 7. Compare WT vs KO cluster assignments
# =============================================================================

message("\n=== Comparing WT vs KO cluster assignments ===")

comparison <- merge(wt_clusters, ko_clusters, by = "peak_id", all = FALSE)
message(sprintf("  Peaks in both WT and KO: %d", nrow(comparison)))

# Cross-tabulation
crosstab <- table(WT = comparison$wt_cluster, KO = comparison$ko_cluster)
message("\n  Cross-tabulation (WT rows x KO columns):")
print(crosstab)

# Classify: Retained = same cluster name, Redistributed = different
comparison$status <- ifelse(comparison$wt_cluster == comparison$ko_cluster,
                            "Retained", "Redistributed")

status_by_wt <- comparison %>%
  group_by(wt_cluster, status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(wt_cluster) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

message("\n  Disruption summary by WT cluster:")
print(as.data.frame(status_by_wt))

# Save tables
write.csv(comparison[, c("peak_id", "wt_cluster_num", "wt_cluster",
                          "ko_cluster_num", "ko_cluster", "status")],
          file.path(outdir, "peak_cluster_comparison.csv"), row.names = FALSE)

crosstab_df <- as.data.frame.matrix(crosstab)
write.csv(crosstab_df, file.path(outdir, "cluster_crosstab_wt_vs_ko.csv"))

message("  Saved comparison tables")

# =============================================================================
# 8. Figure 1: Side-by-side cluster line plots
# =============================================================================

message("\n=== Figure 1: Side-by-side cluster line plots ===")

# Build profile data for WT
wt_profile <- deg_wt_df %>%
  mutate(cluster_name = wt_name_map[as.character(.data[[deg_wt_col]])]) %>%
  group_by(cluster_name, timepoint) %>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            sd_val = sd(value, na.rm = TRUE),
            n = n_distinct(.data[[wt_gene_col]]),
            .groups = "drop")
wt_profile$timepoint <- factor(wt_profile$timepoint, levels = tp_levels)

# Add peak count label
wt_n <- wt_clusters %>% count(wt_cluster, name = "n_peaks")
wt_profile <- merge(wt_profile, wt_n, by.x = "cluster_name", by.y = "wt_cluster")
wt_profile$label <- sprintf("%s\n(n=%s)", wt_profile$cluster_name,
                             format(wt_profile$n_peaks, big.mark = ","))
wt_profile$genotype <- "WT"

# Build profile data for KO
ko_profile <- deg_ko_df %>%
  mutate(cluster_name = ko_name_map[as.character(.data[[deg_ko_col]])]) %>%
  group_by(cluster_name, timepoint) %>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            sd_val = sd(value, na.rm = TRUE),
            n = n_distinct(.data[[ko_gene_col]]),
            .groups = "drop")
ko_profile$timepoint <- factor(ko_profile$timepoint, levels = tp_levels)

ko_n <- ko_clusters %>% count(ko_cluster, name = "n_peaks")
ko_profile <- merge(ko_profile, ko_n, by.x = "cluster_name", by.y = "ko_cluster")
ko_profile$label <- sprintf("%s\n(n=%s)", ko_profile$cluster_name,
                             format(ko_profile$n_peaks, big.mark = ","))
ko_profile$genotype <- "KO"

# WT panel
p_wt_profiles <- ggplot(wt_profile, aes(x = timepoint, y = mean_val, group = 1)) +
  geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
              alpha = 0.15, fill = "black") +
  geom_line(linewidth = 1.2, color = "black") +
  geom_point(size = 2.5, color = "black") +
  facet_wrap(~ label, scales = "free_y") +
  labs(x = NULL, y = "Scaled accessibility", title = "WT") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 7))

# KO panel
p_ko_profiles <- ggplot(ko_profile, aes(x = timepoint, y = mean_val, group = 1)) +
  geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
              alpha = 0.15, fill = "#2CA02C") +
  geom_line(linewidth = 1.2, color = "#2CA02C") +
  geom_point(size = 2.5, color = "#2CA02C") +
  facet_wrap(~ label, scales = "free_y") +
  labs(x = NULL, y = "Scaled accessibility", title = "ARID1A-KO") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 7))

n_wt_cl <- n_distinct(wt_profile$cluster_name)
n_ko_cl <- n_distinct(ko_profile$cluster_name)
n_cols <- max(n_wt_cl, n_ko_cl)
fig_width <- max(10, n_cols * 2.5)

p_combined <- plot_grid(p_wt_profiles, p_ko_profiles, ncol = 1,
                        labels = c("A", "B"), label_size = 14)
p_titled <- plot_grid(
  ggdraw() + draw_label("degPatterns Temporal Clustering: WT vs ARID1A-KO",
                         fontface = "bold", size = 13),
  p_combined, ncol = 1, rel_heights = c(0.05, 1)
)

save_figure(p_titled, "fig1_cluster_profiles_wt_ko",
            width = fig_width, height = 10, dir = figdir)

# =============================================================================
# 9. Figure 2: Alluvial diagram (WT clusters -> KO clusters)
# =============================================================================

message("\n=== Figure 2: Alluvial diagram ===")

alluvial_data <- comparison %>%
  count(wt_cluster, ko_cluster, name = "n_peaks") %>%
  arrange(wt_cluster, desc(n_peaks))

# Color by WT cluster
wt_cl_names <- sort(unique(alluvial_data$wt_cluster))
n_wt <- length(wt_cl_names)
wt_pal <- if (n_wt <= 8) brewer.pal(max(3, n_wt), "Set2")[1:n_wt] else
  colorRampPalette(brewer.pal(8, "Set2"))(n_wt)
names(wt_pal) <- wt_cl_names

p_alluvial <- ggplot(alluvial_data,
                     aes(axis1 = wt_cluster, axis2 = ko_cluster, y = n_peaks)) +
  geom_alluvium(aes(fill = wt_cluster), width = 1/4, alpha = 0.7) +
  geom_stratum(width = 1/4, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.8) +
  scale_x_discrete(limits = c("WT Cluster", "KO Cluster"), expand = c(0.15, 0.05)) +
  scale_fill_manual(values = wt_pal, name = "WT Cluster") +
  labs(y = "Number of peaks",
       title = "Peak Flow: WT Clusters to KO Clusters",
       subtitle = sprintf("%s peaks tracked across genotypes",
                          format(nrow(comparison), big.mark = ","))) +
  theme_paper +
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())

save_figure(p_alluvial, "fig2_alluvial_wt_to_ko",
            width = 10, height = max(6, n_wt * 0.8), dir = figdir)

# =============================================================================
# 10. Figure 3: Dual heatmap — WT and KO z-scored profiles, same row order
# =============================================================================

message("\n=== Figure 3: Dual heatmap ===")

# Use only peaks present in both WT and KO clustering
shared_peaks <- comparison$peak_id

# Compute per-timepoint averages for WT
wt_avg <- sapply(tp_levels, function(tp) {
  samps <- rownames(wt_meta)[wt_meta$timepoint == tp]
  rowMeans(wt_bal[shared_peaks, samps, drop = FALSE])
})
colnames(wt_avg) <- tp_levels

# Compute per-timepoint averages for KO
ko_avg <- sapply(tp_levels, function(tp) {
  samps <- rownames(ko_meta)[ko_meta$timepoint == tp]
  rowMeans(ko_bal[shared_peaks, samps, drop = FALSE])
})
colnames(ko_avg) <- tp_levels

# Z-score each peak across timepoints
wt_z <- t(scale(t(wt_avg)))
wt_z[is.nan(wt_z)] <- 0
ko_z <- t(scale(t(ko_avg)))
ko_z[is.nan(ko_z)] <- 0

# Order rows by WT cluster assignment
wt_cluster_for_order <- comparison$wt_cluster[match(shared_peaks, comparison$peak_id)]

# Define cluster display order by temporal peak of WT centers
wt_cl_order <- wt_cluster_centers %>%
  arrange(factor(timepoint, levels = tp_levels)) %>%
  group_by(cluster_num) %>%
  summarise(peak_tp = tp_levels[which.max(mean_val)], .groups = "drop") %>%
  mutate(cluster_name = wt_name_map[as.character(cluster_num)]) %>%
  arrange(factor(peak_tp, levels = tp_levels))
cluster_order <- wt_cl_order$cluster_name

# Row order: by cluster, then by peak within cluster
row_order <- order(match(wt_cluster_for_order, cluster_order))

# Cluster colors
n_cl <- length(cluster_order)
cl_pal <- if (n_cl <= 8) brewer.pal(max(3, n_cl), "Set2")[1:n_cl] else
  colorRampPalette(brewer.pal(8, "Set2"))(n_cl)
names(cl_pal) <- cluster_order

# Annotation
row_anno <- rowAnnotation(
  `WT Cluster` = factor(wt_cluster_for_order[row_order], levels = cluster_order),
  col = list(`WT Cluster` = cl_pal),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 9),
  width = unit(5, "mm")
)

# WT heatmap
ht_wt <- Heatmap(
  wt_z[row_order, ],
  name = "z-score\n(WT)",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  left_annotation = row_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 11),
  column_title = "WT",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_split = factor(wt_cluster_for_order[row_order], levels = cluster_order),
  row_gap = unit(1, "mm"),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 8),
  use_raster = TRUE,
  raster_quality = 5,
  border = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(3, "cm")
  )
)

# KO heatmap (same row order)
ht_ko <- Heatmap(
  ko_z[row_order, ],
  name = "z-score\n(KO)",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 11),
  column_title = "ARID1A-KO",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_split = factor(wt_cluster_for_order[row_order], levels = cluster_order),
  row_gap = unit(1, "mm"),
  row_title = NULL,
  use_raster = TRUE,
  raster_quality = 5,
  border = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(3, "cm")
  )
)

ht_combined <- ht_wt + ht_ko

save_heatmap(ht_combined, "fig3_dual_heatmap_wt_ko",
             width = 8, height = 10, dir = figdir)

# =============================================================================
# 11. Figure 4: Disruption barplot
# =============================================================================

message("\n=== Figure 4: Disruption barplot ===")

disruption_data <- comparison %>%
  group_by(wt_cluster) %>%
  summarise(
    total = n(),
    retained = sum(status == "Retained"),
    redistributed = sum(status == "Redistributed"),
    pct_retained = 100 * retained / total,
    pct_redistributed = 100 * redistributed / total,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_redistributed))

# Long format for stacked bar
disruption_long <- disruption_data %>%
  select(wt_cluster, total, pct_retained, pct_redistributed) %>%
  pivot_longer(cols = c(pct_retained, pct_redistributed),
               names_to = "status", values_to = "pct") %>%
  mutate(status = recode(status,
                         "pct_retained" = "Retained",
                         "pct_redistributed" = "Redistributed"))

# Order by disruption
disruption_long$wt_cluster <- factor(disruption_long$wt_cluster,
                                      levels = disruption_data$wt_cluster)

p_disruption <- ggplot(disruption_long,
                       aes(x = wt_cluster, y = pct, fill = status)) +
  geom_col(width = 0.7) +
  geom_text(data = disruption_data,
            aes(x = wt_cluster, y = 102,
                label = sprintf("n=%s", format(total, big.mark = ","))),
            inherit.aes = FALSE, size = 3, hjust = 0) +
  coord_flip(ylim = c(0, 115)) +
  scale_fill_manual(values = c("Retained" = "#4DAF4A", "Redistributed" = "#E41A1C"),
                    name = "Status") +
  labs(x = "WT Cluster", y = "Percentage of peaks (%)",
       title = "Cluster Disruption in ARID1A-KO",
       subtitle = "Fraction of WT cluster peaks that retain same temporal pattern in KO") +
  theme_paper +
  theme(legend.position = "bottom")

save_figure(p_disruption, "fig4_disruption_barplot",
            width = 8, height = max(4, nrow(disruption_data) * 0.6), dir = figdir)

# =============================================================================
# 12. Save workspace
# =============================================================================

message("\n=== Saving workspace ===")

save(
  vst_mat, variable_peaks, tp_levels,
  wt_bal, wt_meta, ko_bal, ko_meta,
  deg_wt, deg_ko,
  wt_clusters, ko_clusters, wt_name_map, ko_name_map,
  comparison, crosstab, status_by_wt, disruption_data,
  wt_z, ko_z,
  file = file.path(outdir, "wt_ko_degpatterns.RData")
)

message(sprintf("Workspace saved: %s", file.path(outdir, "wt_ko_degpatterns.RData")))

# =============================================================================
# Summary
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("SUMMARY")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("Timepoints: %s", paste(tp_levels, collapse = " -> ")))
message(sprintf("Variable peaks (MAD > 0): %d", length(variable_peaks)))
message(sprintf("WT clusters: %d (%d peaks classified)",
                n_distinct(wt_clusters$wt_cluster), nrow(wt_clusters)))
message(sprintf("KO clusters: %d (%d peaks classified)",
                n_distinct(ko_clusters$ko_cluster), nrow(ko_clusters)))
message(sprintf("Peaks in both: %d", nrow(comparison)))
message(sprintf("Retained (same pattern): %d (%.1f%%)",
                sum(comparison$status == "Retained"),
                100 * mean(comparison$status == "Retained")))
message(sprintf("Redistributed: %d (%.1f%%)",
                sum(comparison$status == "Redistributed"),
                100 * mean(comparison$status == "Redistributed")))
message("")
message("WT clusters:")
for (cl in sort(unique(wt_clusters$wt_cluster))) {
  n <- sum(wt_clusters$wt_cluster == cl)
  message(sprintf("  %-25s %6d peaks", cl, n))
}
message("")
message("KO clusters:")
for (cl in sort(unique(ko_clusters$ko_cluster))) {
  n <- sum(ko_clusters$ko_cluster == cl)
  message(sprintf("  %-25s %6d peaks", cl, n))
}
message("")
message(sprintf("Figures: %s", figdir))
message(sprintf("Data:    %s", outdir))
