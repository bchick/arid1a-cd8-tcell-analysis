#!/usr/bin/env Rscript
# =============================================================================
# guo2022/03_summit_requant.R — does the peak universe explain the Guo ED Fig 6a gap?
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Guo report 11,975 lost / 211 gained for Arid1a KO vs WT ATAC (n=3v3,
# FDR<0.05 & FC>2). Our nf-core consensus gives 50 / 65 for the same samples
# (see guo2022/02_atacseq_da.R). Every other explanation has been excluded:
#   - not thresholds  (identical: FDR<0.05, |LFC|>1)
#   - not WT pairing  (fixed in guo2022/02_atacseq_da.R; arm split confirmed in Guo's own peak BEDs)
#   - not power       (MYC arm gives 30,642 peaks at padj<0.05, same pipeline/n)
#   - not QC          (FRiP 0.27-0.37 across all 12 g183618 samples)
#   - not normalization (FRiP KO 0.210 vs WT 0.215; median-of-ratios and
#                        full-library-size agree)
#
# The one remaining difference is the PEAK UNIVERSE. Guo used DiffBind v2.16
# with summit=250 -> fixed 501 bp windows re-centred on summits, over a
# consensus built from the WT+KO samples ONLY. We used nf-core's variable-width
# consensus merged across all 30 Guo samples (median 471 bp, 18.5% >1 kb,
# 114 Mb, 182k peaks). Wide merged peaks dilute a focal signal with unchanged
# flanks, and a 30-sample universe carries a heavier multiple-testing burden.
#
# This script rebuilds the peak universe the DiffBind way and re-runs the same
# DESeq2 contrast, so the ONLY thing that changed is the universe.
#
# Inputs:  results/meta_analysis/guo2022/atacseq/bowtie2/merged_library/
#            <sample>.mLb.clN_summits.bed (macs2/narrow_peak/), <sample>.mLb.clN.sorted.bam
#            (ARID1A arm, 3 KO v 3 WT; if absent, the bundled summit_consensus_501bp.bed
#            and summit_counts.rds are used)
# Outputs: results/extended_analysis/guo2022/atacseq/differential/summit_requant/
#            summit_consensus_501bp.bed, summit_counts.rds,
#            da_KO_vs_WT_summit501.csv, universe_comparison.csv
# Usage:   Rscript extended_analysis/scripts/guo2022/03_summit_requant.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(GenomicRanges)
  library(rtracklayer)
  library(Rsubread)
  library(Rsamtools)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

HALF_WIDTH <- 250L   # DiffBind summit=250 -> 501 bp windows
MIN_REPS   <- 2L     # Guo: "present in at least two replicates of each condition"

atac_dir <- file.path(paths$meta, "atacseq/bowtie2/merged_library")
peak_dir <- file.path(atac_dir, "macs2/narrow_peak")
outdir   <- file.path(paths$ext_results, "guo2022/atacseq/differential/summit_requant")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ARID1A arm only (see the arm mapping in guo2022/02_atacseq_da.R Section 1)
samples <- data.frame(
  sample    = c("Guo_g183618_KO_REP1", "Guo_g183618_KO_REP2", "Guo_g183618_KO_REP3",
                "Guo_g183618_WT_REP1", "Guo_g183618_WT_REP5", "Guo_g183618_WT_REP6"),
  condition = c(rep("KO", 3), rep("WT", 3)),
  stringsAsFactors = FALSE
)
rownames(samples) <- samples$sample

# Sections 1-2 need the MACS2 summits and BAMs (upstream mode). Without them,
# use the consensus and counts they produced, shipped in the data bundle.
summit_files <- file.path(peak_dir, paste0(samples$sample, ".mLb.clN_summits.bed"))
bundled <- file.path(outdir, c("summit_consensus_501bp.bed", "summit_counts.rds"))
raw_ok <- all(file.exists(summit_files))
if (!raw_ok) {
  if (!all(file.exists(bundled)))
    stop("Need the Guo summits/BAMs (upstream mode) or the bundled summit_requant files")
  message("Skipping sections 1-2: summits/BAMs not found; using the bundled consensus + counts")
  consensus <- import(bundled[1], format = "BED")
  names(consensus) <- consensus$name
  cm <- readRDS(bundled[2])
  stopifnot(identical(rownames(cm), names(consensus)))
}

if (raw_ok) {
  # =============================================================================
  # 1. Build a DiffBind-style consensus: summit-centred fixed-width windows
  # =============================================================================
  message("=== Section 1: Building summit-centred consensus (WT+KO only) ===")

  read_summits <- function(s) {
    f <- file.path(peak_dir, paste0(s, ".mLb.clN_summits.bed"))
    stopifnot("summits file not found" = file.exists(f))
    x <- read.table(f, sep = "\t", header = FALSE,
                    col.names = c("chr", "start", "end", "name", "score"))
    GRanges(x$chr, IRanges(x$start + 1, x$end), score = x$score, sample = s)
  }
  sm <- lapply(samples$sample, read_summits)
  names(sm) <- samples$sample
  for (s in samples$sample) message(sprintf("  %-22s %7d summits", s, length(sm[[s]])))

  # Expand each summit to a fixed window, then find windows reproducible across
  # replicates. Reproducibility is assessed per condition (a region kept if it is
  # present in >= MIN_REPS replicates of EITHER condition) so that regions lost in
  # KO are not silently excluded from the test.
  win <- lapply(sm, function(g) resize(g, width = 2L * HALF_WIDTH + 1L, fix = "center"))

  # Pool all windows, merge overlapping, then re-centre each merged region on its
  # strongest contributing summit (DiffBind re-centres on the best summit).
  pooled <- reduce(unlist(GRangesList(win)), ignore.strand = TRUE)
  message(sprintf("  Pooled/merged regions: %d", length(pooled)))

  all_summits <- unlist(GRangesList(sm))
  hits <- findOverlaps(pooled, all_summits, ignore.strand = TRUE)
  best <- data.frame(region = queryHits(hits),
                     score  = all_summits$score[subjectHits(hits)],
                     pos    = start(all_summits)[subjectHits(hits)])
  best <- best[order(best$region, -best$score), ]
  best <- best[!duplicated(best$region), ]

  consensus <- GRanges(
    seqnames(pooled)[best$region],
    IRanges(best$pos - HALF_WIDTH, best$pos + HALF_WIDTH)
  )
  message(sprintf("  Re-centred on best summit: %d regions of %d bp",
                  length(consensus), 2L * HALF_WIDTH + 1L))

  # Reproducibility filter
  sup <- vapply(c("KO", "WT"), function(cond) {
    ss <- samples$sample[samples$condition == cond]
    Reduce(`+`, lapply(ss, function(s) countOverlaps(consensus, win[[s]]) > 0L))
  }, integer(length(consensus)))
  keep <- sup[, "KO"] >= MIN_REPS | sup[, "WT"] >= MIN_REPS
  consensus <- consensus[keep]
  message(sprintf("  Kept %d regions present in >=%d reps of either condition",
                  length(consensus), MIN_REPS))

  # Drop any window that ran off a contig end
  si <- seqlengths(BamFile(file.path(atac_dir,
          paste0(samples$sample[1], ".mLb.clN.sorted.bam"))))
  ok <- start(consensus) >= 1 &
        end(consensus) <= si[as.character(seqnames(consensus))]
  consensus <- consensus[which(ok)]
  names(consensus) <- sprintf("win_%06d", seq_along(consensus))
  message(sprintf("  Final consensus: %d windows (%.1f Mb total)",
                  length(consensus), sum(width(consensus)) / 1e6))

  export(consensus, file.path(outdir, "summit_consensus_501bp.bed"), format = "BED")

  # =============================================================================
  # 2. Re-count reads in those windows (featureCounts, nf-core settings)
  # =============================================================================
  message("=== Section 2: featureCounts over summit windows ===")

  saf <- data.frame(
    GeneID = names(consensus),
    Chr    = as.character(seqnames(consensus)),
    Start  = start(consensus),
    End    = end(consensus),
    Strand = "+",
    stringsAsFactors = FALSE
  )
  bams <- file.path(atac_dir, paste0(samples$sample, ".mLb.clN.sorted.bam"))
  stopifnot("BAM(s) missing" = all(file.exists(bams)))

  fc <- featureCounts(
    files = bams, annot.ext = saf, isPairedEnd = TRUE,
    allowMultiOverlap = TRUE, fracOverlap = 0.2,   # nf-core: -O --fracOverlap 0.2
    strandSpecific = 0, nthreads = 12, verbose = FALSE
  )
  cm <- fc$counts
  colnames(cm) <- samples$sample
  message(sprintf("  Counted %d windows x %d samples", nrow(cm), ncol(cm)))
  saveRDS(cm, file.path(outdir, "summit_counts.rds"))
}  # raw_ok

# =============================================================================
# 3. Same DESeq2 contrast, same thresholds — only the universe changed
# =============================================================================
message("=== Section 3: DESeq2 on the summit universe ===")

cd <- samples
cd$condition <- factor(cd$condition, levels = c("WT", "KO"))
dds <- DESeqDataSetFromMatrix(cm, cd, ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
message(sprintf("  After pre-filter: %d windows", nrow(dds)))
dds <- DESeq(dds)

res_sh <- lfcShrink(dds, coef = "condition_KO_vs_WT", type = "apeglm")
res_un <- results(dds, name = "condition_KO_vs_WT")

tally <- function(r) c(
  gained = sum(r$padj < 0.05 & r$log2FoldChange >  1, na.rm = TRUE),
  lost   = sum(r$padj < 0.05 & r$log2FoldChange < -1, na.rm = TRUE),
  sig    = sum(r$padj < 0.05, na.rm = TRUE)
)
a <- tally(res_sh); b <- tally(res_un)

cat("\n================ SUMMIT-RECENTRED UNIVERSE ================\n")
cat(sprintf("windows tested                 : %d\n", sum(!is.na(res_un$padj))))
cat(sprintf("apeglm-shrunk   gained/lost    : %d / %d   (padj<0.05 total: %d)\n",
            a["gained"], a["lost"], a["sig"]))
cat(sprintf("unshrunken      gained/lost    : %d / %d   (padj<0.05 total: %d)\n",
            b["gained"], b["lost"], b["sig"]))
cat(sprintf("median LFC (unshrunk)          : %+.3f\n",
            median(res_un$log2FoldChange, na.rm = TRUE)))
cat("-----------------------------------------------------------\n")
cat("nf-core variable-width universe: 50 / 65   (179,872 tested)\n")
cat("Guo published (ED Fig 6a)      : 211 / 11,975\n")
cat("===========================================================\n\n")

out <- as.data.frame(res_sh) |>
  tibble::rownames_to_column("window_id") |>
  dplyr::mutate(
    log2FoldChange_unshrunk = res_un$log2FoldChange[match(window_id, rownames(res_un))],
    chr = as.character(seqnames(consensus))[match(window_id, names(consensus))],
    start = start(consensus)[match(window_id, names(consensus))],
    end   = end(consensus)[match(window_id, names(consensus))],
    sig = dplyr::case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Gained",
      padj < 0.05 & log2FoldChange < -1 ~ "Lost",
      TRUE ~ "NS")
  ) |>
  dplyr::arrange(padj)
readr::write_csv(out, file.path(outdir, "da_KO_vs_WT_summit501.csv"))

summ <- data.frame(
  universe = c("nf-core variable-width (all 30 samples)",
               "summit-recentred 501bp (WT+KO only)",
               "Guo published (ED Fig 6a)"),
  n_tested = c(179872, sum(!is.na(res_un$padj)), NA),
  gained   = c(50, a["gained"], 211),
  lost     = c(65, a["lost"], 11975)
)
readr::write_csv(summ, file.path(outdir, "universe_comparison.csv"))
print(summ)
message(sprintf("\nOutputs -> %s", outdir))
