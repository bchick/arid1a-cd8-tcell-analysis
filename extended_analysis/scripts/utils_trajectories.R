#!/usr/bin/env Rscript
# =============================================================================
# utils_trajectories.R — shared helpers for the ATAC trajectory analyses
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# dplyr verb aliases, the trajectory output directory (traj_dir), the nf-core
# consensus featureCounts reader, the two-arm trajectory sample sheet and the
# timecourse-patterns config writer. Source after scripts/utils.R.
#
# Sourced by: atac_trajectories/*.R (01_timecourse_patterns_prep.R through
#             07_glmnet_motif_models.R)
# =============================================================================

# Bioconductor packages mask these; the trajectory scripts use dplyr's
select <- dplyr::select
filter <- dplyr::filter
count  <- dplyr::count

traj_dir <- file.path(paths$ext_results, "atac_trajectories")

#' nf-core consensus featureCounts matrix (replicate-merged peak set), with
#' library names as in the sample sheet (REP1 -> Rep1). Peak chromosomes and a
#' BED4 table (0-based starts) are attached as attributes "chr" and "bed".
read_atac_consensus_counts <- function() {
  fc_file <- file.path(paths$atac, "bowtie2/merged_replicate/macs2/narrow_peak/consensus",
                       "consensus_peaks.mRp.clN.featureCounts.txt")
  raw <- read.delim(fc_file, comment.char = "#", check.names = FALSE)
  m <- as.matrix(raw[, -(1:6)])
  storage.mode(m) <- "integer"
  rownames(m) <- raw$Geneid
  colnames(m) <- gsub("REP(\\d)", "Rep\\1", sub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(m)))
  attr(m, "chr") <- raw$Chr
  attr(m, "bed") <- tibble(chr = raw$Chr, start = raw$Start - 1L, end = raw$End, name = raw$Geneid)
  m
}

#' Trajectory libraries: Naive WT (shared day 0), D3, D5 and D8 TE per
#' genotype. `genotypes` selects the arms; group is "shared" for Naive.
atac_trajectory_libraries <- function(libs, genotypes = c("WT", "KO")) {
  pick <- function(pattern) grep(pattern, libs, value = TRUE)
  arm <- function(s, group, time) tibble(sample = s, group = group, time = time)
  sheet <- bind_rows(
    arm(pick("^Naive_WT_Rep"), "shared", 0),
    bind_rows(lapply(genotypes, function(g) bind_rows(
      arm(pick(sprintf("^D3_%s_Rep", g)), g, 3),
      arm(pick(sprintf("^D5_%s_Rep", g)), g, 5),
      arm(pick(sprintf("^D8_%s_TE_", g)), g, 8))))
  )
  sheet$replicate <- ave(sheet$sample, sheet$group, sheet$time,
                         FUN = function(x) paste0("r", seq_along(x)))
  stopifnot(all(dplyr::count(sheet, group, time)$n >= 2))  # >= 2 replicates everywhere
  sheet
}

#' timecourse-patterns config. `extra` is merged over the defaults used here
#' (the Increasing split needs cross = -0.5 on the 0/3/5/8 grid).
write_tcp_config <- function(file, indir, results_dir, title, extra = list()) {
  cfg <- list(
    input = list(samplesheet = file.path(indir, "samplesheet.tsv"),
                 counts      = file.path(indir, "counts.tsv"),
                 features    = file.path(indir, "features.bed"),
                 baseline_group = "shared", time_unit = "days"),
    output = list(dir = results_dir),
    superclusters = list(split = list(cross = -0.5)),
    report = list(title = title)
  )
  cfg <- utils::modifyList(cfg, extra)
  yaml::write_yaml(cfg, file)
  message("  Config: ", file)
}
