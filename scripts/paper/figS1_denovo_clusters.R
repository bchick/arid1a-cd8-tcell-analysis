#!/usr/bin/env Rscript
# =============================================================================
# paper/figS1_denovo_clusters.R — de novo re-derivation of the Fig 1A OCR clusters
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction (validation only)
#
# The canonical Fig 1 clusters are the published cluster BEDs
# (data/metadata/paper_ocr_clusters/, loaded by fig1_published_clusters.R).
# This script re-derives clusters from the reprocessed nf-core data with the
# paper's stated method, and fig1_published_clusters.R compares the two
# (figures/paper/figS1_cluster_concordance.*).
#
# Paper (STAR Methods): peaks called per condition; differentially accessible
# regions by DESeq2 (2-fold, FDR < 0.05); "k-means clustering of ATAC-seq data
# was performed using DESeq2" normalised signal. Fig 1A clusters: Conserved,
# Naive, Early Activation, Activation, Late Activation.
#
# Implementation:
#   - WT libraries only: Naive, D3 (the paper's "48h/d3" column; 48h ATAC is
#     excluded — REP2 FRiP 0.115), D5, D8 (TE+EEC+MP subsets pooled).
#   - Consensus peaks = nf-core/atacseq merged-replicate consensus.
#   - DESeq2 ~ timepoint; all 6 pairwise WT timepoint contrasts.
#   - Conserved = called in the merged-replicate peak set of every timepoint
#     (boolean matrix) AND not DA in any contrast.
#   - Dynamic = DA (|log2FC| >= 1, padj < 0.05) in >= 1 contrast; k-means
#     (k = K_DYNAMIC) on row z-scores of per-timepoint VST means; clusters are
#     named from their mean profiles and same-named k-means clusters merged.
#   - Dynamic peaks must also be called at >= 1 WT timepoint.
#   - Cluster sizes differ from the paper (HOMER -style dnase peaks, ~43k OCRs)
#     because the MACS2 consensus is larger (~95k after filtering); thresholds
#     are the paper's and were not tuned to match its counts.
#   - Everything else (neither Conserved nor Dynamic) is left unassigned.
#
# Inputs:  consensus featureCounts + boolean matrices (nf-core/atacseq)
# Outputs: results/paper/fig1_denovo_ocr_clusters.csv (peak_id = consensus peak),
#          results/paper/fig1_denovo_cluster_sizes.csv,
#          results/paper/fig1_denovo_vst_timepoint_means.csv,
#          results/paper/fig1_clusters_denovo/*.bed
# Usage:   Rscript scripts/paper/figS1_denovo_clusters.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")

suppressPackageStartupMessages({
  library(DESeq2)
})

K_DYNAMIC <- 6   # over-cluster, then merge by profile name
TIMEPOINTS <- c("Naive", "D3", "D5", "D8")

cons_dir <- file.path(paths$atac, "bowtie2/merged_replicate/macs2/narrow_peak/consensus")
fc_file  <- file.path(cons_dir, "consensus_peaks.mRp.clN.featureCounts.txt")
bool_file <- file.path(cons_dir, "consensus_peaks.mRp.clN.boolean.txt")

# -----------------------------------------------------------------------------
# 1. Counts (WT libraries)
# -----------------------------------------------------------------------------
fc <- read.delim(fc_file, comment.char = "#", check.names = FALSE)
cnt <- as.matrix(fc[, 7:ncol(fc)])
rownames(cnt) <- fc$Geneid
colnames(cnt) <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(cnt))

wt <- grep("^(Naive_WT|D3_WT|D5_WT|D8_WT)_", colnames(cnt), value = TRUE)
coldata <- data.frame(
  row.names = wt,
  timepoint = factor(sub("^(Naive|D3|D5|D8)_.*", "\\1", wt), levels = TIMEPOINTS)
)
message("WT libraries per timepoint:"); print(table(coldata$timepoint))

cnt <- cnt[, wt]
cnt <- cnt[rowSums(cnt >= 10) >= 2, ]
message(sprintf("%d peaks after count filter", nrow(cnt)))

dds <- DESeqDataSetFromMatrix(cnt, coldata, design = ~ timepoint)
dds <- DESeq(dds, parallel = FALSE, quiet = TRUE)

# -----------------------------------------------------------------------------
# 2. Pairwise WT timepoint contrasts
# -----------------------------------------------------------------------------
pairs <- combn(TIMEPOINTS, 2, simplify = FALSE)
is_da <- sapply(pairs, function(p) {
  r <- results(dds, contrast = c("timepoint", p[2], p[1]), alpha = PAPER$atac_padj)
  !is.na(r$padj) & r$padj < PAPER$atac_padj & abs(r$log2FoldChange) >= PAPER$atac_lfc
})
colnames(is_da) <- sapply(pairs, paste, collapse = "_vs_")
message("DA peaks per pairwise contrast:"); print(colSums(is_da))
any_da <- rowSums(is_da) > 0

# -----------------------------------------------------------------------------
# 3. Peak called at every timepoint (merged-replicate boolean)
# -----------------------------------------------------------------------------
bool <- read.delim(bool_file, check.names = FALSE)
b <- bool[match(rownames(cnt), bool$interval_id), ]
called <- cbind(
  Naive = b[["Naive_WT.mRp.clN.bool"]],
  D3    = b[["D3_WT.mRp.clN.bool"]],
  D5    = b[["D5_WT.mRp.clN.bool"]],
  D8    = do.call(pmax, b[grep("^D8_WT_.*\\.bool$", names(b))])
)
called[is.na(called)] <- 0
called_all <- rowSums(called > 0) == 4

# -----------------------------------------------------------------------------
# 4. Per-timepoint VST means; k-means on dynamic peaks
# -----------------------------------------------------------------------------
vst_mat <- assay(vst(dds, blind = TRUE))
tp_mean <- sapply(TIMEPOINTS, function(tp)
  rowMeans(vst_mat[, coldata$timepoint == tp, drop = FALSE]))

dyn <- tp_mean[any_da & rowSums(called > 0) > 0, ]   # DA and called at >= 1 timepoint
z <- t(scale(t(dyn)))

set.seed(42)
km <- kmeans(z, centers = K_DYNAMIC, nstart = 50, iter.max = 500, algorithm = "MacQueen")
prof <- aggregate(as.data.frame(z), list(k = km$cluster), mean)
message("k-means cluster profiles (row z-score):")
print(cbind(prof, n = as.vector(table(km$cluster))), digits = 2)

# Name by where accessibility peaks:
#   Naive            max at Naive
#   Early Activation max at D3 and > 0.5 z above both D5 and D8 (transient)
#   Late Activation  max at D8 and below average at D3
#   Activation       everything else (opened by D3 and largely maintained)
name_profile <- function(p) {
  p <- setNames(unlist(p[TIMEPOINTS]), TIMEPOINTS)
  top <- names(which.max(p))
  if (top == "Naive") return("Naive")
  if (top == "D3" && p["D3"] - max(p[c("D5", "D8")]) > 0.5) return("Early Activation")
  if (top == "D8" && p["D3"] < 0) return("Late Activation")
  "Activation"
}
km_names <- setNames(apply(prof[, TIMEPOINTS], 1, name_profile), prof$k)
message("k-means -> cluster name:"); print(km_names)

cluster <- rep(NA_character_, nrow(cnt))
names(cluster) <- rownames(cnt)
cluster[called_all & !any_da] <- "Conserved"
cluster[rownames(dyn)] <- km_names[as.character(km$cluster)]

# -----------------------------------------------------------------------------
# 5. Write
# -----------------------------------------------------------------------------
peaks <- fc[match(rownames(cnt), fc$Geneid), c("Geneid", "Chr", "Start", "End")]
out <- tibble(peak_id = peaks$Geneid, chr = peaks$Chr,
              start = peaks$Start, end = peaks$End, cluster = cluster) |>
  filter(!is.na(cluster)) |>
  mutate(cluster = factor(cluster, levels = OCR_CLUSTERS)) |>
  arrange(cluster, chr, start)

paper_n <- c(Conserved = 21650, Naive = 3157, `Early Activation` = 7883,
             Activation = 8781, `Late Activation` = 1250)
summ <- out |> dplyr::count(cluster, name = "n_reanalysis") |>
  mutate(n_paper = paper_n[as.character(cluster)])
message("Cluster sizes vs paper:"); print(summ)

write_panel_table(out, "fig1_denovo_ocr_clusters")
write_panel_table(summ, "fig1_denovo_cluster_sizes")
write_panel_table(as_tibble(tp_mean, rownames = "peak_id"), "fig1_denovo_vst_timepoint_means")

bed_dir <- file.path(paths$paper_tab, "fig1_clusters_denovo")
dir.create(bed_dir, showWarnings = FALSE)
for (cl in OCR_CLUSTERS) {
  out |> filter(cluster == cl) |>
    transmute(chr, start = start - 1L, end, name = peak_id, score = 0, strand = ".") |>
    write_tsv(file.path(bed_dir, paste0(gsub(" ", "_", cl), ".bed")), col_names = FALSE)
}
message("Done.")
