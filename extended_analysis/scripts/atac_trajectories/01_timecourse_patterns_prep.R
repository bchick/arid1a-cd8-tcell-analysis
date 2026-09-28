#!/usr/bin/env Rscript
# =============================================================================
# atac_trajectories/01_timecourse_patterns_prep.R — inputs for the WT and KO trajectory run
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Pass 1 of extended_analysis/scripts/atac_trajectories/01_timecourse_patterns.sh:
#   1. this script: two-arm timecourse inputs + config
#   2. timecourse-patterns (github.com/bchick/timecourse-patterns): per arm,
#      DESeq2 LRT over time + range gate selects dynamic peaks, DEGreport
#      degPatterns clusters them into shape-named classes
#   3. atac_trajectories/02_dynamic_site_sets.R: WT-only / KO-only / shared
#      dynamic sites, and pass-2 inputs to cluster the sites lost in KO
#
# Arms: WT and KO, Naive -> D3 -> D5 -> D8 TE. There are no Naive KO
# libraries, so Naive WT is the shared day-0 baseline of both arms. D8 is TE
# only (Exp1 + Exp2, present in both genotypes). All consensus peaks on
# primary chromosomes go in; the workflow selects dynamic peaks per arm
# (LRT padj < 0.01 and range >= 0.5 on the VST scale).
#
# Inputs:  results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt
# Outputs: results/extended_analysis/atac_trajectories/pass1_wt_ko/
#            inputs/{counts.tsv,samplesheet.tsv,features.bed}, config.yaml
# Usage:   Rscript extended_analysis/scripts/atac_trajectories/01_timecourse_patterns_prep.R   (from the repository root;
#          normally called by 01_timecourse_patterns.sh)
# =============================================================================

source("scripts/utils.R")
source("extended_analysis/scripts/utils_trajectories.R")

outdir <- file.path(traj_dir, "pass1_wt_ko")
indir  <- file.path(outdir, "inputs")
dir.create(indir, recursive = TRUE, showWarnings = FALSE)
# Absolute: the workflow runs from its own checkout, not the repository root
outdir <- normalizePath(outdir); indir <- normalizePath(indir)

# =============================================================================
# 1. Counts
# =============================================================================

count_mat <- read_atac_consensus_counts()
std_chr <- paste0("chr", c(1:19, "X", "Y"))
keep <- attr(count_mat, "chr") %in% std_chr

# =============================================================================
# 2. Samplesheet: two arms sharing the Naive WT baseline
# =============================================================================

sheet <- atac_trajectory_libraries(colnames(count_mat))
print(count(sheet, group, time))

cnt <- count_mat[keep, sheet$sample]
write_tsv(tibble(feature = rownames(cnt), as_tibble(cnt)), file.path(indir, "counts.tsv"))
write_tsv(sheet, file.path(indir, "samplesheet.tsv"))
write_tsv(attr(count_mat, "bed")[keep, ], file.path(indir, "features.bed"), col_names = FALSE)
message(sprintf("  %d peaks x %d libraries written to %s", nrow(cnt), ncol(cnt), indir))

# =============================================================================
# 3. Workflow config
# =============================================================================
# Workflow defaults plus the one override the tool's own T cell example uses:
# the second-stage split of "Increasing" needs cross = -0.5 on the 0/3/5/8
# grid. Split labels are per arm, so KO classes are never compared with WT
# classes by name; atac_trajectories/03_ko_trajectory_projection.R projects KO trajectories onto the WT centroids instead.

write_tcp_config(file.path(outdir, "config.yaml"), indir, file.path(outdir, "results"),
                 title = "CD8+ T cell ATAC-seq trajectories: WT and ARID1A KO (Naive to D8 TE)")
