#!/usr/bin/env Rscript
# =============================================================================
# core/03_atac_temporal_clustering.R — Multi-method temporal clustering of WT ATAC
# McDonald, Chick et al. 2023 Immunity 56:1303 — core analysis
#
# Clusters WT ATAC-seq consensus peaks by their accessibility dynamics over
# Naive, D3, D5 and D8 (D8 pseudobulked across subsets; 48h excluded because
# it failed QC). Dynamic peaks are selected from DESeq2 VST values, k is
# chosen from several selection metrics, and three methods are compared:
# k-means, Mfuzz (fuzzy c-means) and DEGreport degPatterns.
#
# Inputs:  results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/
#            consensus_peaks.mRp.clN.featureCounts.txt
# Outputs: results/atac/temporal_clustering/ ({kmeans,mfuzz,degpatterns}_clusters.csv,
#            cross_method_comparison.csv, kmeans_*.bed, mfuzz_core_*.bed,
#            temporal_clustering.RData)
#          figures/atacseq/temporal_clustering/
# Usage:   Rscript scripts/core/03_atac_temporal_clustering.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(Mfuzz)
  library(Biobase)
  library(DEGreport)
  library(factoextra)
  library(cluster)
  library(RColorBrewer)
})

# Reassign dplyr verbs after Bioc loading
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename   # masked by S4Vectors
count  <- dplyr::count    # masked by matrixStats
desc   <- dplyr::desc     # masked by IRanges

outdir <- file.path(paths$results, "atac/temporal_clustering")
figdir <- file.path(paths$figures, "atacseq/temporal_clustering")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load featureCounts consensus peak matrix
# =============================================================================

message("=== Loading featureCounts consensus peak matrix ===")

fc_file <- file.path(paths$atac,
  "bowtie2/merged_replicate/macs2/narrow_peak/consensus",
  "consensus_peaks.mRp.clN.featureCounts.txt")

fc_raw <- read.delim(fc_file, comment.char = "#", check.names = FALSE)

peak_info <- fc_raw[, 1:6]
count_mat <- as.matrix(fc_raw[, 7:ncol(fc_raw)])
rownames(count_mat) <- peak_info$Geneid

# Clean column names
colnames(count_mat) <- gsub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(count_mat))

message(sprintf("  %d peaks x %d samples", nrow(count_mat), ncol(count_mat)))

# =============================================================================
# 2. Select WT samples — EXCLUDE all 48h
# =============================================================================

message("=== Selecting WT samples (excluding 48h) ===")

wt_cols <- grep("^(Naive_WT|D3_WT|D5_WT|D8_WT)", colnames(count_mat), value = TRUE)
message(sprintf("  %d WT samples (48h excluded)", length(wt_cols)))

count_wt <- count_mat[, wt_cols]

# Parse sample metadata
parse_sample <- function(s) {
  if (grepl("^Naive", s))
    return(data.frame(timepoint = "Naive", subset = "Bulk", experiment = "Exp1",
                      stringsAsFactors = FALSE))
  if (grepl("^D3", s))
    return(data.frame(timepoint = "D3", subset = "Bulk", experiment = "Exp1",
                      stringsAsFactors = FALSE))
  if (grepl("^D5", s))
    return(data.frame(timepoint = "D5", subset = "Bulk", experiment = "Exp1",
                      stringsAsFactors = FALSE))
  m <- regmatches(s, regexec("D8_WT_(TE|EEC|MP)_(Exp[12])_REP", s))[[1]]
  if (length(m) > 0)
    return(data.frame(timepoint = "D8", subset = m[2], experiment = m[3],
                      stringsAsFactors = FALSE))
  return(data.frame(timepoint = "Unknown", subset = "Unknown", experiment = "Unknown",
                    stringsAsFactors = FALSE))
}

col_data <- do.call(rbind, lapply(wt_cols, parse_sample))
rownames(col_data) <- wt_cols
col_data$timepoint <- factor(col_data$timepoint, levels = c("Naive", "D3", "D5", "D8"))

message("  Samples per timepoint:")
print(table(col_data$timepoint))

# =============================================================================
# 3. DESeq2 VST normalization
# =============================================================================

message("=== DESeq2 VST normalization ===")

keep <- rowSums(count_wt >= 10) >= 2
count_wt_filt <- count_wt[keep, ]
message(sprintf("  Kept %d / %d peaks after filtering (>=10 counts in >=2 samples)",
                nrow(count_wt_filt), nrow(count_wt)))

dds <- DESeqDataSetFromMatrix(
  countData = count_wt_filt,
  colData   = col_data,
  design    = ~ timepoint
)

vst <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_mat <- assay(vst)

message(sprintf("  VST matrix: %d peaks x %d samples", nrow(vst_mat), ncol(vst_mat)))

# =============================================================================
# 4. Average by timepoint (D8 pseudobulked across subsets)
# =============================================================================

message("=== Averaging by timepoint ===")

tp_levels <- c("Naive", "D3", "D5", "D8")
avg_mat <- sapply(tp_levels, function(tp) {
  samps <- rownames(col_data)[col_data$timepoint == tp]
  rowMeans(vst_mat[, samps])
})
colnames(avg_mat) <- tp_levels

for (tp in tp_levels) {
  n <- sum(col_data$timepoint == tp)
  message(sprintf("    %s: %d samples averaged", tp, n))
}

# =============================================================================
# 5. Identify dynamic peaks
# =============================================================================

message("=== Filtering to dynamic peaks ===")

row_var <- apply(avg_mat, 1, var)
row_mad <- apply(avg_mat, 1, mad)

# Constitutive: lowest-variance peaks
n_const <- 1500
const_idx <- order(row_var)[1:n_const]
const_peaks <- rownames(avg_mat)[const_idx]

# Dynamic: everything with MAD > 0 minus constitutive
is_variable <- row_mad > 0
dynamic_peaks <- setdiff(rownames(avg_mat)[is_variable], const_peaks)

message(sprintf("  Constitutive: %d peaks (lowest variance)", length(const_peaks)))
message(sprintf("  Dynamic: %d peaks", length(dynamic_peaks)))

# Z-score matrices
avg_dynamic <- avg_mat[dynamic_peaks, ]
avg_z <- t(scale(t(avg_dynamic)))
avg_z[is.nan(avg_z)] <- 0

avg_const <- avg_mat[const_peaks, ]
const_z <- t(scale(t(avg_const)))
const_z[is.nan(const_z)] <- 0

# =============================================================================
# 6. Optimal k selection
# =============================================================================

message("=== Optimal k selection (k=3 to 10) ===")

k_range <- 3:10

# Subsample for gap statistic speed
set.seed(42)
n_sub <- min(nrow(avg_z), 10000)
sub_idx <- sample(nrow(avg_z), n_sub)
avg_z_sub <- avg_z[sub_idx, ]

# Elbow (WSS)
message("  Computing WSS...")
wss <- sapply(k_range, function(k) {
  km_tmp <- kmeans(avg_z_sub, centers = k, nstart = 25, iter.max = 200)
  km_tmp$tot.withinss
})

# Silhouette
message("  Computing silhouette scores...")
sil_avg <- sapply(k_range, function(k) {
  km_tmp <- kmeans(avg_z_sub, centers = k, nstart = 25, iter.max = 200)
  ss <- silhouette(km_tmp$cluster, dist(avg_z_sub))
  mean(ss[, 3])
})

# Gap statistic
message("  Computing gap statistic (B=50)...")
gap_stat <- clusGap(avg_z_sub, FUNcluster = kmeans, nstart = 25,
                    K.max = max(k_range), B = 50)
gap_tab <- gap_stat$Tab[k_range, ]

# Multi-panel k selection plot
k_df <- data.frame(
  k   = k_range,
  WSS = wss,
  Silhouette = sil_avg,
  Gap = gap_tab[, "gap"],
  Gap_SE = gap_tab[, "SE.sim"]
)

p_wss <- ggplot(k_df, aes(x = k, y = WSS)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_x_continuous(breaks = k_range) +
  labs(x = "k", y = "Total within-cluster SS", title = "Elbow") +
  theme_paper

p_sil <- ggplot(k_df, aes(x = k, y = Silhouette)) +
  geom_line(linewidth = 1, color = "#D62728") +
  geom_point(size = 2.5, color = "#D62728") +
  scale_x_continuous(breaks = k_range) +
  labs(x = "k", y = "Mean silhouette width", title = "Silhouette") +
  theme_paper

p_gap <- ggplot(k_df, aes(x = k, y = Gap)) +
  geom_line(linewidth = 1, color = "#2CA02C") +
  geom_point(size = 2.5, color = "#2CA02C") +
  geom_errorbar(aes(ymin = Gap - Gap_SE, ymax = Gap + Gap_SE), width = 0.2,
                color = "#2CA02C") +
  scale_x_continuous(breaks = k_range) +
  labs(x = "k", y = "Gap statistic", title = "Gap Statistic") +
  theme_paper

p_k_combined <- cowplot::plot_grid(p_wss, p_sil, p_gap, nrow = 1, align = "h")
p_k_titled <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label(
    sprintf("Optimal k Selection (%d dynamic peaks, %d subsampled)", length(dynamic_peaks), n_sub),
    fontface = "bold", size = 12),
  p_k_combined, ncol = 1, rel_heights = c(0.08, 1)
)
save_figure(p_k_titled, "k_selection_metrics", width = 12, height = 4, dir = figdir)

# Determine optimal k
optimal_k_sil <- k_range[which.max(sil_avg)]
optimal_k_gap <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"],
                       method = "Tibs2001SEmax")

message("  k selection summary:")
for (i in seq_along(k_range)) {
  marker <- ""
  if (k_range[i] == optimal_k_sil) marker <- paste0(marker, " <-- silhouette optimum")
  if (k_range[i] == optimal_k_gap) marker <- paste0(marker, " <-- gap optimum")
  message(sprintf("    k=%d: WSS=%.0f  Sil=%.3f  Gap=%.3f (SE=%.3f)%s",
                  k_range[i], wss[i], sil_avg[i],
                  gap_tab[i, "gap"], gap_tab[i, "SE.sim"], marker))
}

K <- optimal_k_sil
message(sprintf("\n  >>> Using k = %d (silhouette optimum) <<<\n", K))

# =============================================================================
# 7. METHOD 1: K-means
# =============================================================================

message(sprintf("=== Method 1: K-means (k=%d) ===", K))

set.seed(42)
km <- kmeans(avg_z, centers = K, nstart = 50, iter.max = 200)

# Name clusters by temporal profile
km_centers <- km$centers
colnames(km_centers) <- tp_levels

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

km_names <- sapply(1:K, function(i) name_cluster(km_centers[i, ]))

# Deduplicate names
name_tab <- table(km_names)
for (dn in names(name_tab[name_tab > 1])) {
  idx <- which(km_names == dn)
  for (j in seq_along(idx)) km_names[idx[j]] <- paste0(dn, " ", j)
}
km_name_map <- setNames(km_names, as.character(1:K))

# Build labels
km_labels <- km_name_map[as.character(km$cluster)]
names(km_labels) <- names(km$cluster)

# Order clusters by temporal peak
peak_order <- order(apply(km_centers, 1, which.max), apply(km_centers, 1, max))
km_display_order <- km_name_map[as.character(peak_order)]

message("  K-means clusters:")
for (i in peak_order) {
  n <- sum(km$cluster == i)
  prof <- paste(sprintf("%.2f", km_centers[i, ]), collapse = "  ")
  message(sprintf("    %s: %d peaks  [%s]", km_name_map[as.character(i)], n, prof))
}

# =============================================================================
# 8. METHOD 2: Mfuzz (fuzzy c-means)
# =============================================================================

message(sprintf("\n=== Method 2: Mfuzz fuzzy c-means (k=%d) ===", K))

# Create ExpressionSet from un-z-scored averages (Mfuzz standardizes internally)
eset <- ExpressionSet(assayData = avg_dynamic)

# Filter zero-variance features (Mfuzz requirement)
eset <- filter.std(eset, min.std = 0)
message(sprintf("  After filter.std: %d peaks", nrow(eset)))

# Standardize (z-score per peak)
eset <- standardise(eset)

# Estimate fuzzifier
m_est <- mestimate(eset)
message(sprintf("  Estimated fuzzifier m = %.2f", m_est))

# Mfuzz Dmin plot for cluster number estimation
message("  Computing Dmin across k=3..10...")
png(file.path(figdir, "mfuzz_Dmin.png"), width = 6, height = 4, units = "in", res = 300)
Dmin_result <- Dmin(eset, m = m_est, crange = k_range, repeats = 3, visu = TRUE)
dev.off()
message("  Saved Dmin plot")

# Run fuzzy c-means
set.seed(42)
mf <- mfuzz(eset, c = K, m = m_est)

membership <- mf$membership
colnames(membership) <- paste0("C", 1:K)

# Name Mfuzz clusters
mf_centers <- mf$centers
colnames(mf_centers) <- tp_levels
mf_names <- sapply(1:K, function(i) name_cluster(mf_centers[i, ]))
name_tab <- table(mf_names)
for (dn in names(name_tab[name_tab > 1])) {
  idx <- which(mf_names == dn)
  for (j in seq_along(idx)) mf_names[idx[j]] <- paste0(dn, " ", j)
}
mf_name_map <- setNames(mf_names, as.character(1:K))

# Core peaks per cluster
core_threshold <- 0.5
core_counts <- sapply(1:K, function(i) sum(membership[, i] > core_threshold))

message("  Mfuzz clusters:")
for (i in 1:K) {
  n <- sum(mf$cluster == i)
  prof <- paste(sprintf("%.2f", mf_centers[i, ]), collapse = "  ")
  message(sprintf("    C%d %s: %d total, %d core  [%s]",
                  i, mf_name_map[as.character(i)], n, core_counts[i], prof))
}

# --- Mfuzz membership heatmap ---

mf_order <- order(mf$cluster, -apply(membership, 1, max))

ht_membership <- Heatmap(
  membership[mf_order, ],
  name = "Membership",
  col = colorRamp2(seq(0, 1, length.out = 100), viridisLite::mako(100)),
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  show_row_names  = FALSE,
  column_names_rot = 45,
  column_labels = paste0("C", 1:K, ": ", mf_name_map),
  column_names_gp = gpar(fontsize = 9),
  column_title     = "Mfuzz Cluster Membership",
  column_title_gp  = gpar(fontsize = 12, fontface = "bold"),
  row_split = factor(paste0("C", mf$cluster[mf_order]), levels = paste0("C", 1:K)),
  row_gap       = unit(1, "mm"),
  row_title_rot = 0,
  row_title_gp  = gpar(fontsize = 9),
  use_raster     = TRUE,
  raster_quality = 5,
  border = TRUE
)

save_heatmap(ht_membership, "mfuzz_membership_heatmap", width = 6, height = 10, dir = figdir)

# --- Mfuzz profile plots (center +/- SD, colored by membership) ---

mf_center_data <- do.call(rbind, lapply(1:K, function(i) {
  idx <- which(mf$cluster == i)
  peak_mat <- exprs(eset)[names(mf$cluster)[idx], , drop = FALSE]
  data.frame(
    cluster  = i,
    name     = mf_name_map[as.character(i)],
    n_total  = length(idx),
    n_core   = core_counts[i],
    timepoint = tp_levels,
    center   = mf_centers[i, ],
    sd       = apply(peak_mat, 2, sd),
    stringsAsFactors = FALSE
  )
}))
mf_center_data$timepoint <- factor(mf_center_data$timepoint, levels = tp_levels)
mf_center_data$label <- sprintf("C%d: %s\n(%d peaks, %d core)",
                                mf_center_data$cluster, mf_center_data$name,
                                mf_center_data$n_total, mf_center_data$n_core)
mf_center_data$label <- factor(mf_center_data$label, levels = unique(mf_center_data$label))

p_mfuzz <- ggplot(mf_center_data, aes(x = timepoint, y = center, group = 1)) +
  geom_ribbon(aes(ymin = center - sd, ymax = center + sd), alpha = 0.2, fill = "steelblue") +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_point(size = 2.5, color = "steelblue") +
  facet_wrap(~ label, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "Standardized accessibility",
       title = "Mfuzz Cluster Profiles (center +/- SD)") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 7))

save_figure(p_mfuzz, "mfuzz_cluster_profiles", width = max(10, K * 2.5), height = 3.5, dir = figdir)

# --- Checkpoint: save k-means + Mfuzz results before degPatterns ---
save(
  count_wt, col_data, vst_mat, avg_mat, avg_z, avg_dynamic, avg_const,
  const_peaks, dynamic_peaks, peak_info, tp_levels, const_z,
  K, optimal_k_sil, optimal_k_gap, sil_avg, wss, gap_stat, k_range,
  km, km_labels, km_name_map, km_centers, km_display_order,
  mf, membership, mf_name_map, mf_centers, m_est, core_threshold, core_counts,
  file = file.path(outdir, "checkpoint_km_mfuzz.RData")
)
message("  Checkpoint saved (k-means + Mfuzz)")

# =============================================================================
# 9. METHOD 3: DEGreport degPatterns
# =============================================================================

message("\n=== Method 3: DEGreport degPatterns ===")

# degPatterns works on individual replicates.
# For scalability, use top N most variable peaks by MAD.
# Pseudobulk D8 by experiment for balanced design (n=2 per timepoint).

n_deg_max <- 5000
peak_mad <- apply(vst_mat[dynamic_peaks, ], 1, mad)
top_var_peaks <- names(sort(peak_mad, decreasing = TRUE))[1:min(n_deg_max, length(peak_mad))]
message(sprintf("  Input: top %d most variable dynamic peaks", length(top_var_peaks)))

# Build balanced matrix: individual reps for Naive/D3/D5, pseudobulked D8 by experiment
early_samps <- c(
  grep("^Naive_WT", rownames(col_data), value = TRUE),
  grep("^D3_WT", rownames(col_data), value = TRUE),
  grep("^D5_WT", rownames(col_data), value = TRUE)
)
d8_exp1 <- rownames(col_data)[col_data$timepoint == "D8" & col_data$experiment == "Exp1"]
d8_exp2 <- rownames(col_data)[col_data$timepoint == "D8" & col_data$experiment == "Exp2"]

bal_mat <- cbind(
  vst_mat[top_var_peaks, early_samps],
  D8_WT_Exp1 = rowMeans(vst_mat[top_var_peaks, d8_exp1]),
  D8_WT_Exp2 = rowMeans(vst_mat[top_var_peaks, d8_exp2])
)

bal_meta <- data.frame(
  timepoint = factor(
    c(rep("Naive", 2), rep("D3", 2), rep("D5", 2), rep("D8", 2)),
    levels = tp_levels
  ),
  row.names = colnames(bal_mat)
)

message("  Balanced design for degPatterns:")
print(table(bal_meta$timepoint))

message("  Running degPatterns...")

# Temporarily reset ggplot theme — degPatterns builds internal plots that
# fail with theme_paper's Arial font on headless servers.
old_theme <- theme_get()
theme_set(theme_bw())

pdf(file.path(figdir, "degpatterns_internal.pdf"), width = 10, height = 8)
deg_result <- tryCatch(
  degPatterns(bal_mat, metadata = bal_meta, time = "timepoint",
              minc = 50, reduce = TRUE, cutoff = 0.7, scale = TRUE),
  error = function(e) {
    message("  WARNING: degPatterns failed: ", e$message)
    message("  Retrying with fewer peaks (3K)...")
    top_3k <- top_var_peaks[1:min(3000, length(top_var_peaks))]
    bal_mat_small <- bal_mat[top_3k, ]
    degPatterns(bal_mat_small, metadata = bal_meta, time = "timepoint",
                minc = 50, reduce = TRUE, cutoff = 0.7, scale = TRUE)
  }
)
dev.off()

# Restore theme
theme_set(old_theme)

# Extract cluster assignments
deg_df <- deg_result$df
deg_col <- if ("merge" %in% colnames(deg_df)) "merge" else "cluster"

deg_peak_clusters <- deg_df %>%
  select(genes, !!sym(deg_col)) %>%
  distinct() %>%
  rename(peak_id = genes, deg_cluster = !!sym(deg_col))

deg_summary <- deg_peak_clusters %>%
  count(deg_cluster, name = "n_peaks") %>%
  arrange(desc(n_peaks))

message("  degPatterns clusters:")
print(as.data.frame(deg_summary))
message(sprintf("  Total peaks clustered: %d / %d input (%.0f%%)",
                nrow(deg_peak_clusters), length(top_var_peaks),
                100 * nrow(deg_peak_clusters) / length(top_var_peaks)))

# Build degPatterns profile plot from cluster data (avoids internal font issues);
# per-timepoint values live in $normalized ($df holds only gene -> cluster)
deg_plot_data <- deg_result$normalized %>%
  group_by(!!sym(deg_col), timepoint) %>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            sd_val   = sd(value, na.rm = TRUE),
            n        = n_distinct(genes),
            .groups  = "drop") %>%
  rename(cluster = !!sym(deg_col))

deg_plot_data$timepoint <- factor(deg_plot_data$timepoint, levels = tp_levels)
n_per_cluster <- deg_peak_clusters %>% count(deg_cluster)
deg_plot_data <- merge(deg_plot_data, n_per_cluster,
                       by.x = "cluster", by.y = "deg_cluster", all.x = TRUE)
deg_plot_data$label <- sprintf("Cluster %s (n=%s)", deg_plot_data$cluster,
                               format(deg_plot_data$n.y, big.mark = ","))

p_deg <- ggplot(deg_plot_data, aes(x = timepoint, y = mean_val, group = 1)) +
  geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
              alpha = 0.15, fill = "#D62728") +
  geom_line(linewidth = 1.2, color = "#D62728") +
  geom_point(size = 2.5, color = "#D62728") +
  facet_wrap(~ label, scales = "free_y") +
  labs(x = NULL, y = "Scaled accessibility",
       title = sprintf("degPatterns Clusters (%d peaks classified)", nrow(deg_peak_clusters))) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 8))

save_figure(p_deg, "degpatterns_clusters", width = 10, height = 8, dir = figdir)

# =============================================================================
# 10. K-means heatmap (main figure)
# =============================================================================

message("\n=== Building K-means heatmap ===")

# Combined labels: constitutive + dynamic
all_labels_km <- c(
  setNames(rep("Constitutive", length(const_peaks)), const_peaks),
  km_labels
)
all_z_km <- rbind(const_z, avg_z[names(km_labels), ])

cluster_display_order <- c("Constitutive", km_display_order)
main_order <- order(match(all_labels_km, cluster_display_order))

# Cluster colors
n_dyn <- length(km_display_order)
dyn_pal <- if (n_dyn <= 8) brewer.pal(max(3, n_dyn), "Set2")[1:n_dyn] else
  colorRampPalette(brewer.pal(8, "Set2"))(n_dyn)
all_colors <- c("Constitutive" = "#1A1A1A", setNames(dyn_pal, km_display_order))

row_anno <- rowAnnotation(
  Cluster = factor(all_labels_km[main_order], levels = cluster_display_order),
  col = list(Cluster = all_colors),
  show_annotation_name = FALSE,
  width = unit(5, "mm")
)

ht_km <- Heatmap(
  all_z_km[main_order, ],
  name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  left_annotation = row_anno,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  show_row_names  = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 12),
  column_title     = "K-means: WT ATAC-seq Temporal Dynamics",
  column_title_gp  = gpar(fontsize = 13, fontface = "bold"),
  row_split = factor(all_labels_km[main_order], levels = cluster_display_order),
  row_gap       = unit(1, "mm"),
  row_title_rot = 0,
  row_title_gp  = gpar(fontsize = 9),
  use_raster     = TRUE,
  raster_quality = 5,
  border = TRUE,
  heatmap_legend_param = list(
    title = "z-score",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(4, "cm")
  )
)

save_heatmap(ht_km, "kmeans_temporal_heatmap", width = 5, height = 10, dir = figdir)

# K-means profile line plots
km_profile <- do.call(rbind, lapply(cluster_display_order, function(cl) {
  idx <- which(all_labels_km == cl)
  if (length(idx) == 0) return(NULL)
  data.frame(
    cluster   = cl,
    n         = length(idx),
    timepoint = tp_levels,
    zscore    = colMeans(all_z_km[idx, , drop = FALSE]),
    sd        = apply(all_z_km[idx, , drop = FALSE], 2, sd),
    stringsAsFactors = FALSE
  )
}))
km_profile$timepoint <- factor(km_profile$timepoint, levels = tp_levels)
km_profile$label <- sprintf("%s (n=%s)", km_profile$cluster,
                            format(km_profile$n, big.mark = ","))
km_profile$label <- factor(km_profile$label, levels = unique(km_profile$label))

p_km <- ggplot(km_profile, aes(x = timepoint, y = zscore, group = 1)) +
  geom_ribbon(aes(ymin = zscore - sd, ymax = zscore + sd), alpha = 0.15) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  facet_wrap(~ label, nrow = 1) +
  labs(x = NULL, y = "Mean z-score",
       title = "K-means Cluster Profiles (mean +/- SD)") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 7))

save_figure(p_km, "kmeans_cluster_profiles", width = max(10, (K + 1) * 2.2), height = 3.5, dir = figdir)

# =============================================================================
# 11. Cross-method comparison
# =============================================================================

message("\n=== Cross-method comparison ===")

# All dynamic peaks have km + mfuzz labels
# degPatterns covers a subset (top variable peaks that passed clustering)
all_dynamic_df <- data.frame(
  peak_id          = dynamic_peaks,
  kmeans_cluster   = km_labels[dynamic_peaks],
  mfuzz_cluster    = mf_name_map[as.character(mf$cluster[dynamic_peaks])],
  mfuzz_max_member = apply(membership[dynamic_peaks, , drop = FALSE], 1, max),
  stringsAsFactors = FALSE
)
all_dynamic_df <- merge(all_dynamic_df, deg_peak_clusters, by = "peak_id", all.x = TRUE)

write.csv(all_dynamic_df, file.path(outdir, "cross_method_comparison.csv"), row.names = FALSE)

# Contingency tables
message("\n  K-means vs Mfuzz:")
print(table(Kmeans = all_dynamic_df$kmeans_cluster,
            Mfuzz  = all_dynamic_df$mfuzz_cluster))

deg_overlap <- all_dynamic_df[!is.na(all_dynamic_df$deg_cluster), ]
message(sprintf("\n  K-means vs degPatterns (%d overlapping peaks):", nrow(deg_overlap)))
print(table(Kmeans = deg_overlap$kmeans_cluster,
            degPat = deg_overlap$deg_cluster))

# Agreement rate: fraction of peaks where km and mfuzz agree on cluster name
# (match by name, since we used the same naming function)
agree <- sum(all_dynamic_df$kmeans_cluster == all_dynamic_df$mfuzz_cluster, na.rm = TRUE)
message(sprintf("\n  K-means / Mfuzz agreement (same cluster name): %d / %d (%.1f%%)",
                agree, nrow(all_dynamic_df), 100 * agree / nrow(all_dynamic_df)))

# =============================================================================
# 12. Save cluster assignments + BED files
# =============================================================================

message("\n=== Saving cluster assignments and BED files ===")

# K-means full table (constitutive + dynamic)
km_out <- data.frame(
  peak_id  = names(all_labels_km),
  chr      = peak_info$Chr[match(names(all_labels_km), peak_info$Geneid)],
  start    = peak_info$Start[match(names(all_labels_km), peak_info$Geneid)],
  end      = peak_info$End[match(names(all_labels_km), peak_info$Geneid)],
  cluster  = all_labels_km,
  stringsAsFactors = FALSE
)
write.csv(km_out, file.path(outdir, "kmeans_clusters.csv"), row.names = FALSE)

# Mfuzz full table with all membership scores
mf_out <- data.frame(
  peak_id    = names(mf$cluster),
  chr        = peak_info$Chr[match(names(mf$cluster), peak_info$Geneid)],
  start      = peak_info$Start[match(names(mf$cluster), peak_info$Geneid)],
  end        = peak_info$End[match(names(mf$cluster), peak_info$Geneid)],
  cluster    = mf_name_map[as.character(mf$cluster)],
  max_membership = apply(membership, 1, max),
  stringsAsFactors = FALSE
)
for (i in 1:K) {
  mf_out[[paste0("membership_C", i)]] <- membership[, i]
}
write.csv(mf_out, file.path(outdir, "mfuzz_clusters.csv"), row.names = FALSE)

# degPatterns table
write.csv(deg_peak_clusters, file.path(outdir, "degpatterns_clusters.csv"), row.names = FALSE)

# BED files — K-means
for (cl in cluster_display_order) {
  rows <- km_out[km_out$cluster == cl, ]
  bed <- data.frame(rows$chr, rows$start, rows$end, rows$peak_id, 0, ".")
  bed_path <- file.path(outdir, paste0("kmeans_", gsub(" ", "_", cl), ".bed"))
  write.table(bed, bed_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message(sprintf("  K-means %s: %d peaks -> %s", cl, nrow(bed), basename(bed_path)))
}

# BED files — Mfuzz core peaks
for (i in 1:K) {
  core_idx <- which(membership[, i] > core_threshold)
  if (length(core_idx) == 0) next
  core_ids <- rownames(membership)[core_idx]
  rows <- peak_info[match(core_ids, peak_info$Geneid), ]
  bed <- data.frame(rows$Chr, rows$Start, rows$End, core_ids,
                    round(membership[core_ids, i] * 1000), ".")
  bed_path <- file.path(outdir, sprintf("mfuzz_core_C%d_%s.bed",
                                         i, gsub(" ", "_", mf_name_map[as.character(i)])))
  write.table(bed, bed_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message(sprintf("  Mfuzz core C%d (%s): %d peaks -> %s",
                  i, mf_name_map[as.character(i)], nrow(bed), basename(bed_path)))
}

# =============================================================================
# 13. Save workspace
# =============================================================================

save(
  count_wt, col_data, vst_mat, avg_mat, avg_z, avg_dynamic, avg_const,
  const_peaks, dynamic_peaks, peak_info, tp_levels,
  K, optimal_k_sil, optimal_k_gap, sil_avg, wss, gap_stat, k_range,
  km, km_labels, km_name_map, km_centers, km_display_order,
  mf, membership, mf_name_map, mf_centers, m_est, core_threshold, core_counts,
  deg_result, deg_peak_clusters, top_var_peaks,
  all_dynamic_df,
  file = file.path(outdir, "temporal_clustering.RData")
)
message(sprintf("\nWorkspace saved: %s", file.path(outdir, "temporal_clustering.RData")))

# =============================================================================
# Summary
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("SUMMARY")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("Timepoints: %s (48h excluded)", paste(tp_levels, collapse = " -> ")))
message(sprintf("WT samples: %d", length(wt_cols)))
message(sprintf("Peaks after filtering: %d", nrow(vst_mat)))
message(sprintf("  Constitutive: %d", length(const_peaks)))
message(sprintf("  Dynamic (clustered): %d", length(dynamic_peaks)))
message(sprintf("Optimal k: %d (silhouette=%.3f)", K, max(sil_avg)))
message("")
message("K-means clusters:")
for (cl in cluster_display_order) {
  n <- sum(all_labels_km == cl)
  message(sprintf("  %-25s %6d peaks (%4.1f%%)", cl, n, 100 * n / length(all_labels_km)))
}
message("")
message(sprintf("Mfuzz clusters (core threshold=%.1f):", core_threshold))
for (i in 1:K) {
  nm <- mf_name_map[as.character(i)]
  message(sprintf("  C%d %-22s %6d total, %5d core", i, nm,
                  sum(mf$cluster == i), core_counts[i]))
}
message("")
message(sprintf("degPatterns: %d clusters, %d/%d peaks clustered",
                nrow(deg_summary), nrow(deg_peak_clusters), length(top_var_peaks)))
message("")
message(sprintf("Figures: %s", figdir))
message(sprintf("Data:    %s", outdir))
