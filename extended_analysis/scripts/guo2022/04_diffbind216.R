#!/usr/bin/env Rscript
# =============================================================================
# guo2022/04_diffbind216.R — reproduce Guo et al. 2022 ED Fig 6a with DiffBind 2.16.0
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Guo 2022 (Nature) Methods, "Analysis of ATAC-seq data":
#   "Differentially accessibility analysis was performed using DiffBind
#    (v.2.16.0) (summit = 250, DESeq2 and other default parameters)."
#   Thresholds (Peak set enrichment analysis): "false positive rate < 0.05 and
#    fold change > 2".
#   Published ED Fig 6a: 11,975 lost / 211 gained, Arid1a KO vs WT, n = 3 v 3.
#
# Our prior result on the same BAMs via nf-core consensus + featureCounts +
# DESeq2 (guo2022/02_atacseq_da.R, after the two-arm batch fix): 50 gained / 65 lost.
# A summit-recentred 501 bp requantification (guo2022/03_summit_requant.R) gave 34 / 37.
#
# This script removes the last software difference by running the actual
# DiffBind 2.16.0 package (Bioconductor 3.11 container), rather than emulating
# its consensus geometry and normalisation inside our own DESeq2 pipeline.
#
# ONE known difference remains and is INTENTIONAL: our BAMs are mm39; Guo used
# mm10. If this run still returns ~tens of peaks, genome
# build is the sole untested variable left.
#
# CRITICAL: only the ARID1A arm. GSE183618 bundles two experiments under one
# "WT" GEO label — WT_REP2/3/4 belong to the MYC arm. Pooling them manufactures
# ~7-9k spurious "lost" peaks. See the arm mapping in guo2022/02_atacseq_da.R
# Section 1.
#
# Requires: Docker, to run inside the DiffBind 2.16.0 container
#   quay.io/biocontainers/bioconductor-diffbind:2.16.0--r40h5f743cb_2
# (the script stops unless DiffBind 2.16.x is loaded).
#
# Inputs:  results/meta_analysis/guo2022/atacseq/bowtie2/merged_library/
#            Guo_g183618_{KO,WT}_REP<N>.mLb.clN.sorted.bam and macs2/narrow_peak/ peaks
# Outputs: results/extended_analysis/guo2022/atacseq/differential/diffbind216/
#            diffbind_samplesheet.csv, dba_{counted,analyzed}.rds,
#            da_KO_vs_WT_diffbind216_{full,significant}.csv,
#            universe_comparison_diffbind216.csv
# Usage:   docker run --rm -u $(id -u):$(id -g) -v "$PWD":/work -w /work \
#            -e ARID1A_PROJECT_DIR=/work -e RENV_ACTIVATE_PROJECT=FALSE \
#            quay.io/biocontainers/bioconductor-diffbind:2.16.0--r40h5f743cb_2 \
#            Rscript extended_analysis/scripts/guo2022/04_diffbind216.R   (from the repository root;
#          `make meta` runs this)
# =============================================================================

suppressPackageStartupMessages(library(DiffBind))

PROJ    <- normalizePath(Sys.getenv("ARID1A_PROJECT_DIR", unset = "."))
BAM_DIR <- file.path(PROJ, "results/meta_analysis/guo2022/atacseq/bowtie2/merged_library")
PK_DIR  <- file.path(BAM_DIR, "macs2/narrow_peak")
OUT_DIR <- file.path(PROJ, "results/extended_analysis/guo2022/atacseq/differential/diffbind216")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# --- Version gate: this MUST be 2.16.0, not the 3.x installed on the host ----
dbv <- as.character(packageVersion("DiffBind"))
log_msg("DiffBind version: ", dbv, " | R: ", R.version.string)
if (!grepl("^2\\.16", dbv)) {
  stop("Expected DiffBind 2.16.x (Guo's version); got ", dbv,
       ". DiffBind 3.x changed normalisation, greylist and counting defaults ",
       "and is NOT their method. Run inside the 2.16.0 container.")
}

# --- Sample sheet: ARID1A arm only, 3 v 3 -----------------------------------
ARID1A_ARM <- list(
  KO = c("KO_REP1", "KO_REP2", "KO_REP3"),
  WT = c("WT_REP1", "WT_REP5", "WT_REP6")   # NOT REP2/3/4 — those are MYC arm
)

samples <- data.frame(
  SampleID   = c(ARID1A_ARM$KO, ARID1A_ARM$WT),
  Condition  = rep(c("KO", "WT"), each = 3),
  Replicate  = rep(1:3, times = 2),
  bamReads   = file.path(BAM_DIR, paste0("Guo_g183618_",
                         c(ARID1A_ARM$KO, ARID1A_ARM$WT), ".mLb.clN.sorted.bam")),
  Peaks      = file.path(PK_DIR, paste0("Guo_g183618_",
                         c(ARID1A_ARM$KO, ARID1A_ARM$WT), ".mLb.clN_peaks.narrowPeak")),
  PeakCaller = "narrow",
  stringsAsFactors = FALSE
)

stopifnot(nrow(samples) == 6, sum(samples$Condition == "WT") == 3)
if (any(grepl("WT_REP[234]", samples$SampleID))) {
  stop("MYC-arm WT sample present — this is the g183618 two-arm batch confound.")
}
write.csv(samples, file.path(OUT_DIR, "diffbind_samplesheet.csv"), row.names = FALSE)
log_msg("Sample sheet OK — ARID1A arm, 3 KO vs 3 WT")

# --- Checkpoints ------------------------------------------------------------
# dba.count() takes ~64 min on these libraries, so reuse a prior run's objects
# when they exist. DIFFBIND_FORCE=1 ignores them and recounts from scratch.
COUNTED_RDS  <- file.path(OUT_DIR, "dba_counted.rds")
ANALYZED_RDS <- file.path(OUT_DIR, "dba_analyzed.rds")
FORCE        <- identical(Sys.getenv("DIFFBIND_FORCE"), "1")

# BAMs + peaks are needed only to (re)count; the bundle ships dba_counted.rds
missing <- c(samples$bamReads, samples$Peaks)[!file.exists(c(samples$bamReads, samples$Peaks))]
if (length(missing) && (FORCE || !any(file.exists(c(COUNTED_RDS, ANALYZED_RDS)))))
  stop("Missing inputs:\n  ", paste(missing, collapse = "\n  "))

dbo      <- NULL
analyzed <- FALSE
if (!FORCE && file.exists(ANALYZED_RDS)) {
  log_msg("Checkpoint: reusing ", basename(ANALYZED_RDS),
          " — skipping dba() / dba.count() / dba.analyze()")
  dbo      <- readRDS(ANALYZED_RDS)
  analyzed <- TRUE
} else if (!FORCE && file.exists(COUNTED_RDS)) {
  log_msg("Checkpoint: reusing ", basename(COUNTED_RDS),
          " — skipping dba() / dba.count()")
  dbo <- readRDS(COUNTED_RDS)
}

# A checkpoint from a different sample set would silently reintroduce the very
# batch confound this script exists to avoid — refuse to use one.
if (!is.null(dbo)) {
  ck <- sort(as.character(dbo$samples$SampleID))
  if (!identical(ck, sort(samples$SampleID))) {
    stop("Checkpoint sample set does not match the ARID1A arm:\n  checkpoint: ",
         paste(ck, collapse = ", "), "\n  expected:   ",
         paste(sort(samples$SampleID), collapse = ", "),
         "\nDelete the .rds files or set DIFFBIND_FORCE=1.")
  }
  log_msg("Checkpoint sample set verified: ", paste(ck, collapse = ", "))
}

if (is.null(dbo)) {
  # --- 1. Consensus peakset (DiffBind defaults: minOverlap = 2) -------------
  log_msg("dba() — building consensus peakset with default parameters")
  dbo <- dba(sampleSheet = samples)
  log_msg("Consensus peakset: ", nrow(dbo$binding), " regions")
  print(dbo)

  # --- 2. Count — THE parameter Guo specify: summits = 250 (=> 501 bp windows)
  # bParallel=FALSE only bounds peak memory on these ultra-deep libraries
  # (300-378M read pairs each); it does not affect the counts or the result.
  log_msg("dba.count(summits = 250) — this is the slow step, serial over 6 deep BAMs")
  t0 <- Sys.time()
  dbo <- dba.count(dbo, summits = 250, bParallel = FALSE)
  log_msg("dba.count done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  log_msg("Counted matrix: ", nrow(dbo$binding), " x ", length(dbo$samples$SampleID))
  print(dbo)

  saveRDS(dbo, COUNTED_RDS)
}

# --- 3. Contrast + DESeq2 analysis, all other parameters default ------------
if (!analyzed) {
  log_msg("dba.contrast(categories = DBA_CONDITION)")
  dbo <- dba.contrast(dbo, categories = DBA_CONDITION, minMembers = 2)
  print(dba.show(dbo, bContrasts = TRUE))

  log_msg("dba.analyze(method = DBA_DESEQ2) — DiffBind 2.16 defaults ",
          "(bFullLibrarySize = TRUE, no LFC shrinkage)")
  dbo <- dba.analyze(dbo, method = DBA_DESEQ2)
  saveRDS(dbo, ANALYZED_RDS)
}

# --- 4. Report every region (th = 1), threshold ourselves -------------------
res <- dba.report(dbo, method = DBA_DESEQ2, th = 1, bNormalized = TRUE,
                  bCounts = TRUE, DataType = DBA_DATA_FRAME)
log_msg("dba.report returned ", nrow(res), " regions")

# Fold sign convention: dba.report gives log2(Group1 / Group2). Resolve which
# group is which from the contrast rather than assuming, then orient so that
# POSITIVE = more accessible in KO ("gained"), NEGATIVE = "lost in KO".
ctr <- dba.show(dbo, bContrasts = TRUE)
# dba.show() names the first group column "Group1" (an earlier version of this
# script read ctr$Group -> NULL -> character(0) -> NA in the tests below).
g1_col <- if ("Group1" %in% names(ctr)) "Group1" else "Group"
g1 <- as.character(ctr[[g1_col]][1]); g2 <- as.character(ctr$Group2[1])
stopifnot(length(g1) == 1L, !is.na(g1), length(g2) == 1L, !is.na(g2))
log_msg("Contrast as built: Group1 = ", g1, " (n=", ctr$Members1[1], "), ",
        "Group2 = ", g2, " (n=", ctr$Members2[1], ")")
if (g1 == "WT" && g2 == "KO") {
  res$Fold <- -res$Fold
  log_msg("Flipped Fold sign so positive = up in KO")
} else if (!(g1 == "KO" && g2 == "WT")) {
  stop("Unexpected contrast groups: ", g1, " vs ", g2)
}

write.csv(res, file.path(OUT_DIR, "da_KO_vs_WT_diffbind216_full.csv"), row.names = FALSE)

# --- 5. Guo's thresholds: FDR < 0.05 AND fold change > 2 --------------------
sig    <- res[res$FDR < 0.05 & abs(res$Fold) > 1, ]
gained <- sum(sig$Fold > 0)
lost   <- sum(sig$Fold < 0)

# FDR-only, for comparison with the guo2022/02_atacseq_da.R reporting
fdr_only <- res[res$FDR < 0.05, ]

summary_tbl <- data.frame(
  method       = "DiffBind 2.16.0 (summit=250, DESeq2, defaults)",
  genome       = "mm39 (Guo used mm10)",
  arm          = "ARID1A arm only (KO_REP1-3 vs WT_REP1/5/6)",
  consensus_regions = nrow(res),
  fdr05_only        = nrow(fdr_only),
  gained_fdr05_fc2  = gained,
  lost_fdr05_fc2    = lost,
  published_gained  = 211,
  published_lost    = 11975,
  stringsAsFactors  = FALSE
)
write.csv(summary_tbl, file.path(OUT_DIR, "universe_comparison_diffbind216.csv"), row.names = FALSE)
write.csv(sig, file.path(OUT_DIR, "da_KO_vs_WT_diffbind216_significant.csv"), row.names = FALSE)

cat("\n", strrep("=", 74), "\n", sep = "")
cat("DiffBind 2.16.0 — Guo ARID1A arm, KO vs WT (mm39 BAMs)\n")
cat(strrep("=", 74), "\n", sep = "")
cat("Consensus regions tested : ", nrow(res), "\n", sep = "")
cat("FDR < 0.05 (any fold)    : ", nrow(fdr_only), "\n", sep = "")
cat("FDR < 0.05 & |FC| > 2    : ", nrow(sig), "\n", sep = "")
cat("  gained in KO           : ", gained, "\n", sep = "")
cat("  lost in KO             : ", lost, "\n", sep = "")
cat("\nGuo published (ED Fig 6a): 211 gained / 11,975 lost\n")
cat("nf-core DESeq2 (02)      : 50 gained / 65 lost\n")
cat("summit501 requant (03)   : 34 gained / 37 lost\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Outputs written to ", OUT_DIR)
log_msg("DONE")
