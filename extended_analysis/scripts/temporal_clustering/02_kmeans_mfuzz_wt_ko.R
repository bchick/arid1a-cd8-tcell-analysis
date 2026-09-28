#!/usr/bin/env Rscript
# =============================================================================
# temporal_clustering/02_kmeans_mfuzz_wt_ko.R — temporal chromatin clustering, k-means + Mfuzz
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Clusters ATAC-seq peaks by temporal dynamics (Naive -> D3 -> D5 -> D8)
# using k-means and Mfuzz (fuzzy c-means), separately for WT and ARID1A-KO.
# KO peaks are projected onto WT cluster centers to measure pattern disruption.
#
# Both methods operate on timepoint-averaged z-scored profiles (4 dimensions).
# 30K top-MAD peaks used — both methods scale to this easily (O(n*k) per iter).
#
# Inputs:  results/atac/differential/checkpoint_sections1to3.RData (main_vsd)
# Outputs: results/extended_analysis/temporal_clustering/kmeans_mfuzz/ (tables + RData)
#          figures/extended_analysis/temporal_clustering/kmeans_mfuzz/
# Usage:   Rscript extended_analysis/scripts/temporal_clustering/02_kmeans_mfuzz_wt_ko.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(Mfuzz)
  library(ggalluvial)
  library(RColorBrewer)
  library(cowplot)
  library(viridisLite)
})

select <- dplyr::select
filter <- dplyr::filter
count  <- dplyr::count

outdir <- file.path(paths$ext_results, "temporal_clustering/kmeans_mfuzz")
figdir <- file.path(paths$ext_figures, "temporal_clustering/kmeans_mfuzz")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

tp_levels <- c("Naive", "D3", "D5", "D8")

# Number of clusters — adjust based on elbow plot output
K_FINAL <- 8
# Number of peaks to cluster
N_TOP <- 30000
# k-means: number of random starts
KMEANS_NSTART <- 25
# k-means: max iterations
KMEANS_ITER <- 500

# =============================================================================
# name_cluster() — label clusters by temporal profile
# =============================================================================

name_cluster <- function(center, tp_names = tp_levels) {
  names(center) <- tp_names
  peak_tp  <- which.max(center)
  rng      <- max(center) - min(center)

  if (rng < 0.8) return("Constitutive-like")

  diffs <- diff(center)
  if (all(diffs >  0.15)) return("Progressive Opening")
  if (all(diffs < -0.15)) return("Progressive Closing")

  if (peak_tp == 1 && center[1] - mean(center[2:4]) > 0.5) return("Naive-specific")
  if (peak_tp == 2 && center[2] - mean(center[c(1,3,4)]) > 0.5) return("D3-specific")
  if (peak_tp == 4 && center[4] - mean(center[1:3]) > 0.5) return("D8-specific")

  if (mean(center[1:2]) > 0.2 && mean(center[3:4]) < -0.2) return("Early (Naive+D3)")
  if (mean(center[1:2]) < -0.2 && mean(center[3:4]) > 0.2) return("Late (D5+D8)")
  if (center[1] < -0.3 && all(diffs >= -0.1)) return("Activation-induced")
  if (center[1] > 0.3  && all(diffs <= 0.1))  return("Activation-lost")

  return(paste0("Cluster_", tp_names[peak_tp]))
}

# =============================================================================
# 1. Load data
# =============================================================================

message("=== Loading main_vsd from checkpoint ===")
load(file.path(paths$results, "atac/differential/checkpoint_sections1to3.RData"))
vst_mat <- assay(main_vsd)
message(sprintf("  VST matrix: %d peaks x %d samples", nrow(vst_mat), ncol(vst_mat)))

source("scripts/utils.R")
select <- dplyr::select
filter <- dplyr::filter

# =============================================================================
# 2. Sample selection
# =============================================================================

message("\n=== Identifying WT and KO trajectory samples ===")

wt_naive  <- grep("^Naive_WT_Rep",          colnames(vst_mat), value = TRUE)
wt_d3     <- grep("^D3_WT_Rep",             colnames(vst_mat), value = TRUE)
wt_d5     <- grep("^D5_WT_Rep",             colnames(vst_mat), value = TRUE)
wt_d8_all <- grep("^D8_WT_(TE|EEC|MP)_",   colnames(vst_mat), value = TRUE)
ko_d3     <- grep("^D3_KO_Rep",             colnames(vst_mat), value = TRUE)
ko_d5     <- grep("^D5_KO_Rep",             colnames(vst_mat), value = TRUE)
ko_d8_all <- grep("^D8_KO_(TE|EEC|MP)_",   colnames(vst_mat), value = TRUE)

message(sprintf("  WT: Naive=%d, D3=%d, D5=%d, D8=%d",
                length(wt_naive), length(wt_d3), length(wt_d5), length(wt_d8_all)))
message(sprintf("  KO: Naive=%d (shared WT), D3=%d, D5=%d, D8=%d",
                length(wt_naive), length(ko_d3), length(ko_d5), length(ko_d8_all)))

# =============================================================================
# 3. Select top-MAD variable peaks
# =============================================================================

message("\n=== Selecting variable peaks ===")

tp_geno_means <- cbind(
  Naive  = rowMeans(vst_mat[, wt_naive]),
  D3_WT  = rowMeans(vst_mat[, wt_d3]),
  D5_WT  = rowMeans(vst_mat[, wt_d5]),
  D8_WT  = rowMeans(vst_mat[, wt_d8_all]),
  D3_KO  = rowMeans(vst_mat[, ko_d3]),
  D5_KO  = rowMeans(vst_mat[, ko_d5]),
  D8_KO  = rowMeans(vst_mat[, ko_d8_all])
)

peak_mad <- apply(tp_geno_means, 1, mad)
n_top <- min(N_TOP, sum(peak_mad > 0))
variable_peaks <- names(sort(peak_mad, decreasing = TRUE))[1:n_top]
message(sprintf("  Selected %d peaks (MAD range: %.3f - %.3f)",
                n_top, min(peak_mad[variable_peaks]), max(peak_mad[variable_peaks])))

# =============================================================================
# 4. Build timepoint-averaged matrices (peaks x 4 timepoints)
# =============================================================================

message("\n=== Building timepoint-averaged matrices ===")

wt_tp <- cbind(
  Naive = rowMeans(vst_mat[variable_peaks, wt_naive]),
  D3    = rowMeans(vst_mat[variable_peaks, wt_d3]),
  D5    = rowMeans(vst_mat[variable_peaks, wt_d5]),
  D8    = rowMeans(vst_mat[variable_peaks, wt_d8_all])
)

ko_tp <- cbind(
  Naive = rowMeans(vst_mat[variable_peaks, wt_naive]),  # shared baseline
  D3    = rowMeans(vst_mat[variable_peaks, ko_d3]),
  D5    = rowMeans(vst_mat[variable_peaks, ko_d5]),
  D8    = rowMeans(vst_mat[variable_peaks, ko_d8_all])
)

# Z-score each peak across the 4 timepoints
wt_z <- t(scale(t(wt_tp)))
ko_z <- t(scale(t(ko_tp)))
wt_z[is.nan(wt_z)] <- 0
ko_z[is.nan(ko_z)] <- 0

message(sprintf("  WT matrix: %d peaks x %d timepoints (z-scored)",
                nrow(wt_z), ncol(wt_z)))

# =============================================================================
# 5. K-means: elbow plot (k = 2:15)
# =============================================================================

message("\n=== K-means elbow plot (k=2:15) ===")

set.seed(42)
k_range  <- 2:15
wss_vals <- sapply(k_range, function(k) {
  km <- kmeans(wt_z, centers = k, nstart = 5, iter.max = 200)
  km$tot.withinss
})

elbow_df <- data.frame(k = k_range, wss = wss_vals)

p_elbow <- ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_point(size = 2.5, color = "black") +
  geom_vline(xintercept = K_FINAL, linetype = "dashed", color = "#D62728", linewidth = 0.7) +
  annotate("text", x = K_FINAL + 0.3, y = max(wss_vals) * 0.95,
           label = sprintf("k=%d (chosen)", K_FINAL), color = "#D62728", hjust = 0, size = 3.5) +
  scale_x_continuous(breaks = k_range) +
  labs(x = "Number of clusters (k)", y = "Total within-cluster SS",
       title = "K-means elbow plot — WT ATAC-seq temporal profiles",
       subtitle = sprintf("%s peaks, 4 timepoints (z-scored)", format(n_top, big.mark = ","))) +
  theme_paper

save_figure(p_elbow, "fig0_kmeans_elbow", width = 7, height = 4, dir = figdir)

# =============================================================================
# 6. K-means final clustering (WT)
# =============================================================================

message(sprintf("\n=== K-means final clustering (k=%d, nstart=%d) ===", K_FINAL, KMEANS_NSTART))

set.seed(42)
km_wt <- kmeans(wt_z, centers = K_FINAL, nstart = KMEANS_NSTART, iter.max = KMEANS_ITER)

message(sprintf("  Cluster sizes: %s",
                paste(sort(km_wt$size, decreasing = TRUE), collapse = ", ")))

# Label clusters by temporal shape
km_centers <- km_wt$centers  # k x 4 matrix
km_names <- sapply(1:K_FINAL, function(k) name_cluster(km_centers[k, ]))

# Deduplicate names
name_tab <- table(km_names)
for (dn in names(name_tab[name_tab > 1])) {
  idx <- which(km_names == dn)
  for (j in seq_along(idx)) km_names[idx[j]] <- paste0(dn, " ", j)
}

message("  Cluster assignments:")
for (k in seq_len(K_FINAL)) {
  message(sprintf("    Cluster %d -> %-30s (%d peaks)", k, km_names[k], km_wt$size[k]))
}

# Data frame: peak -> cluster
km_wt_df <- data.frame(
  peak_id        = variable_peaks,
  km_cluster_num = km_wt$cluster,
  km_cluster     = km_names[km_wt$cluster],
  stringsAsFactors = FALSE
)

# =============================================================================
# 7. K-means: assign KO peaks to nearest WT cluster center
# =============================================================================

message("\n=== Projecting KO onto WT k-means cluster centers ===")

# For each KO peak, find closest WT centroid (Euclidean distance in z-score space)
ko_km_assign <- apply(ko_z, 1, function(x) {
  dists <- apply(km_centers, 1, function(c) sqrt(sum((x - c)^2)))
  which.min(dists)
})

km_ko_df <- data.frame(
  peak_id        = variable_peaks,
  km_cluster_num = ko_km_assign,
  km_cluster     = km_names[ko_km_assign],
  stringsAsFactors = FALSE
)

# Comparison
km_comparison <- data.frame(
  peak_id       = variable_peaks,
  wt_cluster    = km_wt_df$km_cluster,
  ko_cluster    = km_ko_df$km_cluster,
  stringsAsFactors = FALSE
)
km_comparison$status <- ifelse(km_comparison$wt_cluster == km_comparison$ko_cluster,
                                "Retained", "Redistributed")

message(sprintf("  Retained (same cluster): %d / %d (%.1f%%)",
                sum(km_comparison$status == "Retained"), nrow(km_comparison),
                100 * mean(km_comparison$status == "Retained")))

# =============================================================================
# 8. Mfuzz clustering (WT)
# =============================================================================

message(sprintf("\n=== Mfuzz clustering (c=%d) ===", K_FINAL))

# ExpressionSet: peaks x timepoints (rows = genes/peaks, cols = samples/timepoints)
eset_wt <- new("ExpressionSet", exprs = wt_tp)
eset_wt <- filter.std(eset_wt, min.std = 0)
eset_wt <- standardise(eset_wt)  # z-scores each peak across timepoints

# Estimate fuzzifier m
m_est <- mestimate(eset_wt)
message(sprintf("  Estimated fuzzifier m = %.3f", m_est))

set.seed(42)
mf_wt <- mfuzz(eset_wt, c = K_FINAL, m = m_est)

message(sprintf("  Cluster sizes (hard assignment at max membership): %s",
                paste(table(mf_wt$cluster), collapse = ", ")))

# mf_wt$centers: clusters x timepoints (c x 4) — use directly
mf_centers <- mf_wt$centers
mf_names <- sapply(1:K_FINAL, function(k) name_cluster(mf_centers[k, ]))

# Deduplicate
name_tab <- table(mf_names)
for (dn in names(name_tab[name_tab > 1])) {
  idx <- which(mf_names == dn)
  for (j in seq_along(idx)) mf_names[idx[j]] <- paste0(dn, " ", j)
}

message("  Mfuzz cluster assignments:")
for (k in seq_len(K_FINAL)) {
  n_hard <- sum(mf_wt$cluster == k)
  message(sprintf("    Cluster %d -> %-30s (%d peaks)", k, mf_names[k], n_hard))
}

# mf_wt$membership: peaks x clusters (n x c)
mf_membership <- mf_wt$membership
rownames(mf_membership) <- rownames(exprs(eset_wt))

mf_wt_df <- data.frame(
  peak_id         = rownames(exprs(eset_wt)),
  mf_cluster_num  = mf_wt$cluster,
  mf_cluster      = mf_names[mf_wt$cluster],
  max_membership  = apply(mf_membership, 1, max),
  stringsAsFactors = FALSE
)

# =============================================================================
# 9. Mfuzz: project KO onto WT cluster centers
# =============================================================================

message("\n=== Projecting KO onto Mfuzz cluster centers ===")

# mf_wt$centers: clusters x timepoints (c x 4)
# ko_z: peaks x timepoints, z-scored per peak — same standardization as standardise()

compute_membership <- function(x, centers, m) {
  # x: timepoints vector; centers: clusters x timepoints matrix
  # iterate over cluster rows
  dists <- apply(centers, 1, function(c) sqrt(sum((x - c)^2)) + 1e-10)
  exponent <- 2 / (m - 1)
  sapply(seq_along(dists), function(j) 1 / sum((dists[j] / dists)^exponent))
}

# Only project peaks that survived filter.std
ko_z_sub <- ko_z[rownames(mf_membership), , drop = FALSE]

ko_membership <- t(apply(ko_z_sub, 1, function(x)
  compute_membership(x, mf_centers, m_est)))
colnames(ko_membership) <- paste0("C", seq_len(K_FINAL))

ko_mf_cluster <- apply(ko_membership, 1, which.max)

mf_ko_df <- data.frame(
  peak_id        = rownames(ko_z_sub),
  mf_cluster_num = ko_mf_cluster,
  mf_cluster     = mf_names[ko_mf_cluster],
  max_membership = apply(ko_membership, 1, max),
  stringsAsFactors = FALSE
)

# Comparison
mf_comparison <- data.frame(
  peak_id    = mf_wt_df$peak_id,
  wt_cluster = mf_wt_df$mf_cluster,
  ko_cluster = mf_ko_df$mf_cluster[match(mf_wt_df$peak_id, mf_ko_df$peak_id)],
  wt_mem     = mf_wt_df$max_membership,
  stringsAsFactors = FALSE
)
mf_comparison$status <- ifelse(mf_comparison$wt_cluster == mf_comparison$ko_cluster,
                                "Retained", "Redistributed")

message(sprintf("  Retained (same cluster): %d / %d (%.1f%%)",
                sum(mf_comparison$status == "Retained"), nrow(mf_comparison),
                100 * mean(mf_comparison$status == "Retained")))

# =============================================================================
# 10. Cluster color palette
# =============================================================================

# Order clusters by temporal peak for consistent coloring
cluster_order_km <- km_names[order(apply(km_centers, 1, which.max))]
n_cl   <- length(cluster_order_km)
cl_pal <- if (n_cl <= 8) brewer.pal(max(3, n_cl), "Set2")[seq_len(n_cl)] else
  colorRampPalette(brewer.pal(8, "Set2"))(n_cl)
names(cl_pal) <- cluster_order_km

# Match Mfuzz clusters to same palette where names overlap
cluster_order_mf <- mf_names[order(apply(mf_centers, 1, which.max))]
mf_pal <- cl_pal[cluster_order_mf]
mf_pal[is.na(mf_pal)] <- colorRampPalette(brewer.pal(8, "Set2"))(sum(is.na(mf_pal)))
names(mf_pal) <- cluster_order_mf

# =============================================================================
# 11. Figure 1: K-means cluster profiles — WT side-by-side panels
# =============================================================================

message("\n=== Figure 1: K-means cluster profiles ===")

build_profile_df <- function(z_mat, cluster_vec, cluster_names, tp = tp_levels) {
  df <- data.frame(z_mat, peak_id = rownames(z_mat), stringsAsFactors = FALSE)
  colnames(df)[1:4] <- tp
  df$cluster_name <- cluster_names[cluster_vec]
  df_long <- pivot_longer(df, cols = all_of(tp), names_to = "timepoint", values_to = "value")
  df_long$timepoint <- factor(df_long$timepoint, levels = tp)
  df_long
}

km_wt_long <- build_profile_df(wt_z, km_wt$cluster, km_names)
km_wt_sum  <- km_wt_long %>%
  group_by(cluster_name, timepoint) %>%
  summarise(mean_val = mean(value), sd_val = sd(value), n = n_distinct(peak_id), .groups="drop")
km_wt_sum$timepoint <- factor(km_wt_sum$timepoint, levels = tp_levels)

n_peaks_km <- km_wt_df %>% count(km_cluster, name = "n_peaks")
km_wt_sum  <- left_join(km_wt_sum, n_peaks_km, by = c("cluster_name" = "km_cluster"))
km_wt_sum$label <- sprintf("%s\n(n=%s)", km_wt_sum$cluster_name,
                            format(km_wt_sum$n_peaks, big.mark = ","))
km_wt_sum$label <- factor(km_wt_sum$label,
                           levels = unique(km_wt_sum$label[order(
                             match(km_wt_sum$cluster_name, cluster_order_km))]))

p_km_wt <- ggplot(km_wt_sum, aes(x = timepoint, y = mean_val, group = 1)) +
  geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
              alpha = 0.12, fill = "black") +
  geom_line(linewidth = 1.1, color = "black") +
  geom_point(size = 2.2, color = "black") +
  facet_wrap(~ label, scales = "free_y") +
  labs(x = NULL, y = "Scaled accessibility (z-score)",
       title = sprintf("K-means (k=%d): WT temporal chromatin clusters", K_FINAL),
       subtitle = sprintf("%s peaks, Naive → D3 → D5 → D8",
                          format(n_top, big.mark = ","))) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 7.5))

save_figure(p_km_wt, "fig1a_kmeans_wt_profiles",
            width = max(10, K_FINAL * 2.2), height = max(6, ceiling(K_FINAL / 4) * 3),
            dir = figdir)

# KO profiles projected onto WT clusters
km_ko_long <- build_profile_df(ko_z, ko_km_assign, km_names)
km_ko_sum  <- km_ko_long %>%
  group_by(cluster_name, timepoint) %>%
  summarise(mean_val = mean(value), sd_val = sd(value), .groups = "drop")
km_ko_sum$timepoint <- factor(km_ko_sum$timepoint, levels = tp_levels)
km_ko_sum$label <- km_wt_sum$label[match(km_ko_sum$cluster_name, km_wt_sum$cluster_name)]
km_ko_sum$label[is.na(km_ko_sum$label)] <- km_ko_sum$cluster_name[is.na(km_ko_sum$label)]

p_km_ko <- ggplot(km_ko_sum, aes(x = timepoint, y = mean_val, group = 1)) +
  geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
              alpha = 0.12, fill = "#2CA02C") +
  geom_line(linewidth = 1.1, color = "#2CA02C") +
  geom_point(size = 2.2, color = "#2CA02C") +
  facet_wrap(~ cluster_name, scales = "free_y") +
  labs(x = NULL, y = "Scaled accessibility (z-score)",
       title = "K-means: KO peaks projected onto WT cluster centers",
       subtitle = "Mean profile of KO peaks assigned to each WT cluster") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 7.5))

p_km_combined <- plot_grid(p_km_wt, p_km_ko, ncol = 1, labels = c("A", "B"),
                            label_size = 13, rel_heights = c(1, 1))
save_figure(p_km_combined, "fig1b_kmeans_wt_ko_profiles",
            width = max(10, K_FINAL * 2.2), height = max(12, ceiling(K_FINAL / 4) * 6),
            dir = figdir)

# =============================================================================
# 12. Figure 2: Mfuzz cluster profiles — WT
# =============================================================================

message("\n=== Figure 2: Mfuzz cluster profiles ===")

mf_wt_long <- build_profile_df(wt_z[mf_wt_df$peak_id, , drop = FALSE], mf_wt$cluster, mf_names)
mf_wt_sum  <- mf_wt_long %>%
  group_by(cluster_name, timepoint) %>%
  summarise(mean_val = mean(value), sd_val = sd(value), n = n_distinct(peak_id), .groups = "drop")
mf_wt_sum$timepoint <- factor(mf_wt_sum$timepoint, levels = tp_levels)

n_peaks_mf  <- mf_wt_df %>% count(mf_cluster, name = "n_peaks")
mf_wt_sum   <- left_join(mf_wt_sum, n_peaks_mf, by = c("cluster_name" = "mf_cluster"))
mf_wt_sum$label <- sprintf("%s\n(n=%s)", mf_wt_sum$cluster_name,
                            format(mf_wt_sum$n_peaks, big.mark = ","))

p_mf_wt <- ggplot(mf_wt_sum, aes(x = timepoint, y = mean_val, group = 1)) +
  geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
              alpha = 0.12, fill = "black") +
  geom_line(linewidth = 1.1, color = "black") +
  geom_point(size = 2.2, color = "black") +
  facet_wrap(~ label, scales = "free_y") +
  labs(x = NULL, y = "Scaled accessibility (z-score)",
       title = sprintf("Mfuzz (c=%d, m=%.2f): WT temporal chromatin clusters", K_FINAL, m_est),
       subtitle = sprintf("%s peaks, Naive → D3 → D5 → D8",
                          format(n_top, big.mark = ","))) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 7.5))

save_figure(p_mf_wt, "fig2a_mfuzz_wt_profiles",
            width = max(10, K_FINAL * 2.2), height = max(6, ceiling(K_FINAL / 4) * 3),
            dir = figdir)

# Membership score distribution
p_mem <- ggplot(mf_wt_df, aes(x = max_membership)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#D62728") +
  annotate("text", x = 0.52, y = Inf, label = "mem ≥ 0.5", color = "#D62728",
           hjust = 0, vjust = 1.5, size = 3.5) +
  labs(x = "Maximum cluster membership score", y = "Number of peaks",
       title = "Mfuzz: distribution of peak membership scores",
       subtitle = sprintf("%.1f%% peaks with max membership ≥ 0.5",
                          100 * mean(mf_wt_df$max_membership >= 0.5))) +
  theme_paper

save_figure(p_mem, "fig2b_mfuzz_membership_dist", width = 6, height = 4, dir = figdir)

# =============================================================================
# 13. Figure 3: Dual heatmaps (WT vs KO), sorted by WT k-means cluster
# =============================================================================

message("\n=== Figure 3: Dual heatmaps (k-means ordering) ===")

row_cluster <- km_wt_df$km_cluster
names(row_cluster) <- km_wt_df$peak_id
row_order <- order(match(row_cluster, cluster_order_km))

row_anno <- rowAnnotation(
  `WT Cluster` = factor(row_cluster[row_order], levels = cluster_order_km),
  col  = list(`WT Cluster` = cl_pal),
  show_annotation_name  = TRUE,
  annotation_name_gp    = gpar(fontsize = 9),
  width = unit(5, "mm")
)

ht_wt <- Heatmap(
  wt_z[row_order, ],
  name = "z-score\n(WT)",
  col  = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  left_annotation    = row_anno,
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  show_row_names     = FALSE,
  column_names_rot   = 0,
  column_names_centered = TRUE,
  column_names_gp    = gpar(fontsize = 11),
  column_title       = "WT",
  column_title_gp    = gpar(fontsize = 12, fontface = "bold"),
  row_split          = factor(row_cluster[row_order], levels = cluster_order_km),
  row_gap = unit(1, "mm"),
  row_title_rot      = 0,
  row_title_gp       = gpar(fontsize = 8),
  use_raster         = TRUE,
  raster_quality     = 5,
  border             = TRUE,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                               labels_gp = gpar(fontsize = 8),
                               legend_height = unit(3, "cm"))
)

ht_ko <- Heatmap(
  ko_z[row_order, ],
  name = "z-score\n(KO)",
  col  = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  show_row_names     = FALSE,
  column_names_rot   = 0,
  column_names_centered = TRUE,
  column_names_gp    = gpar(fontsize = 11),
  column_title       = "ARID1A-KO",
  column_title_gp    = gpar(fontsize = 12, fontface = "bold"),
  row_split          = factor(row_cluster[row_order], levels = cluster_order_km),
  row_gap  = unit(1, "mm"),
  row_title = NULL,
  use_raster      = TRUE,
  raster_quality  = 5,
  border          = TRUE,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                               labels_gp = gpar(fontsize = 8),
                               legend_height = unit(3, "cm"))
)

save_heatmap(ht_wt + ht_ko, "fig3_dual_heatmap_kmeans_wt_ko",
             width = 8, height = 10, dir = figdir)

# =============================================================================
# 14. Figure 4: Disruption barplots (k-means + Mfuzz side-by-side)
# =============================================================================

message("\n=== Figure 4: Disruption barplots ===")

make_disruption_plot <- function(comparison_df, wt_col, title_label) {
  disruption <- comparison_df %>%
    group_by(wt_cluster) %>%
    summarise(total = n(),
              retained = sum(status == "Retained"),
              pct_retained = 100 * retained / total,
              pct_redistributed = 100 * (1 - retained / total),
              .groups = "drop") %>%
    arrange(desc(pct_redistributed))

  long_df <- disruption %>%
    select(wt_cluster, total, pct_retained, pct_redistributed) %>%
    pivot_longer(cols = c(pct_retained, pct_redistributed),
                 names_to = "status", values_to = "pct") %>%
    mutate(status = recode(status,
                           pct_retained = "Retained",
                           pct_redistributed = "Redistributed"))
  long_df$wt_cluster <- factor(long_df$wt_cluster, levels = disruption$wt_cluster)

  ggplot(long_df, aes(x = wt_cluster, y = pct, fill = status)) +
    geom_col(width = 0.7) +
    geom_text(data = disruption,
              aes(x = wt_cluster, y = 103,
                  label = sprintf("n=%s", format(total, big.mark = ","))),
              inherit.aes = FALSE, size = 2.8, hjust = 0) +
    coord_flip(ylim = c(0, 118)) +
    scale_fill_manual(values = c("Retained" = "#4DAF4A", "Redistributed" = "#E41A1C"),
                      name = "Status in KO") +
    labs(x = "WT Cluster", y = "Peaks (%)",
         title = title_label,
         subtitle = "Fraction retaining WT temporal pattern in ARID1A-KO") +
    theme_paper +
    theme(legend.position = "bottom")
}

p_km_disrupt <- make_disruption_plot(km_comparison, "wt_cluster",
                                     sprintf("K-means (k=%d): cluster disruption in KO", K_FINAL))
p_mf_disrupt <- make_disruption_plot(mf_comparison, "wt_cluster",
                                     sprintf("Mfuzz (c=%d): cluster disruption in KO", K_FINAL))

p_disrupt_combined <- plot_grid(p_km_disrupt, p_mf_disrupt, ncol = 2,
                                 labels = c("K-means", "Mfuzz"), label_size = 11)
save_figure(p_disrupt_combined, "fig4_disruption_barplots",
            width = 14, height = max(5, K_FINAL * 0.65), dir = figdir)

# =============================================================================
# 15. Figure 5: Alluvial — WT k-means to KO assignment
# =============================================================================

message("\n=== Figure 5: Alluvial diagram ===")

alluvial_km <- km_comparison %>%
  count(wt_cluster, ko_cluster, name = "n_peaks") %>%
  arrange(wt_cluster, desc(n_peaks))

p_alluvial <- ggplot(alluvial_km,
                     aes(axis1 = wt_cluster, axis2 = ko_cluster, y = n_peaks)) +
  geom_alluvium(aes(fill = wt_cluster), width = 1/4, alpha = 0.7) +
  geom_stratum(width = 1/4, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6) +
  scale_x_discrete(limits = c("WT Cluster", "KO Assignment"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = cl_pal, name = "WT Cluster") +
  labs(y = "Number of peaks",
       title = "Peak flow: WT k-means clusters → KO cluster assignment",
       subtitle = sprintf("%s peaks", format(n_top, big.mark = ","))) +
  theme_paper +
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())

save_figure(p_alluvial, "fig5_alluvial_kmeans_wt_to_ko",
            width = 10, height = max(6, K_FINAL * 0.9), dir = figdir)

# =============================================================================
# 16. Figure 6: K-means vs Mfuzz cluster comparison (confusion matrix)
# =============================================================================

message("\n=== Figure 6: K-means vs Mfuzz agreement ===")

method_comp <- data.frame(
  peak_id   = km_wt_df$peak_id,
  km_cluster = km_wt_df$km_cluster,
  mf_cluster = mf_wt_df$mf_cluster[match(km_wt_df$peak_id, mf_wt_df$peak_id)]
) %>% filter(!is.na(mf_cluster))

method_ct <- method_comp %>%
  count(km_cluster, mf_cluster, name = "n") %>%
  group_by(km_cluster) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_method_comp <- ggplot(method_ct, aes(x = km_cluster, y = mf_cluster, fill = pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", pct)), size = 3) +
  scale_fill_gradientn(colors = mako(50), name = "% of\nk-means\ncluster",
                        limits = c(0, 100)) +
  labs(x = "K-means cluster", y = "Mfuzz cluster",
       title = "Agreement between k-means and Mfuzz cluster assignments (WT)",
       subtitle = "% of each k-means cluster assigned to each Mfuzz cluster") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8))

save_figure(p_method_comp, "fig6_kmeans_vs_mfuzz_agreement",
            width = max(8, K_FINAL * 1.1), height = max(6, K_FINAL * 0.9), dir = figdir)

# =============================================================================
# 17. Save tables and workspace
# =============================================================================

message("\n=== Saving results ===")

write.csv(km_wt_df,     file.path(outdir, "kmeans_wt_clusters.csv"),     row.names = FALSE)
write.csv(km_comparison, file.path(outdir, "kmeans_wt_ko_comparison.csv"), row.names = FALSE)
write.csv(mf_wt_df,     file.path(outdir, "mfuzz_wt_clusters.csv"),      row.names = FALSE)
write.csv(mf_comparison, file.path(outdir, "mfuzz_wt_ko_comparison.csv"),  row.names = FALSE)

write.csv(as.data.frame(km_centers),
          file.path(outdir, "kmeans_cluster_centers.csv"))
write.csv(as.data.frame(t(mf_wt$centers)),
          file.path(outdir, "mfuzz_cluster_centers.csv"))
write.csv(as.data.frame(mf_membership),
          file.path(outdir, "mfuzz_membership_matrix.csv"))

save(
  variable_peaks, wt_z, ko_z, wt_tp, ko_tp,
  km_wt, km_names, km_wt_df, km_ko_df, km_comparison,
  mf_wt, mf_names, mf_wt_df, mf_ko_df, mf_comparison,
  mf_membership, m_est, K_FINAL, tp_levels, cl_pal, mf_pal,
  file = file.path(outdir, "kmeans_mfuzz_clustering.RData")
)

message(sprintf("  Saved tables and RData to %s", outdir))

# =============================================================================
# Summary
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("SUMMARY")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("Peaks clustered:   %d (top MAD)", n_top))
message(sprintf("Timepoints:        %s", paste(tp_levels, collapse = " -> ")))
message(sprintf("k / c:             %d", K_FINAL))
message(sprintf("Mfuzz m:           %.3f", m_est))
message("")
message("K-MEANS WT clusters:")
for (nm in cluster_order_km) {
  n <- sum(km_wt_df$km_cluster == nm)
  message(sprintf("  %-35s %6d peaks", nm, n))
}
message(sprintf("  KO retained (same cluster): %.1f%%",
                100 * mean(km_comparison$status == "Retained")))
message("")
message("MFUZZ WT clusters:")
for (nm in cluster_order_mf) {
  n <- sum(mf_wt_df$mf_cluster == nm)
  message(sprintf("  %-35s %6d peaks", nm, n))
}
message(sprintf("  KO retained (same cluster): %.1f%%",
                100 * mean(mf_comparison$status == "Retained")))
message("")
message(sprintf("Figures: %s", figdir))
message(sprintf("Data:    %s", outdir))
