#!/usr/bin/env Rscript
# =============================================================================
# footprinting/export_jaspar_motifs.R — export JASPAR2020 CORE vertebrate PFMs for TOBIAS
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# TOBIAS BINDetect needs motifs in pfm/jaspar/meme format. This exports the same
# motif set used by chromVAR (chromvar/build_chromvar.R; JASPAR2020 CORE
# vertebrates) so the footprinting (tf_tobias_footprints) and chromVAR
# (dose_chromvar) analyses are directly comparable. Run before
# footprinting/run_tobias_footprinting.sh.
#
# Inputs:  JASPAR2020 Bioconductor package (no project files)
# Outputs: results/extended_analysis/footprinting/motifs/jaspar2020_core_vertebrates.jaspar
#          (raw-count JASPAR format: ">ID<TAB>NAME" then A/C/G/T count rows)
# Usage:   Rscript extended_analysis/scripts/footprinting/export_jaspar_motifs.R   (from the repository root)
# =============================================================================
suppressPackageStartupMessages({ library(JASPAR2020); library(TFBSTools) })

outdir <- "results/extended_analysis/footprinting/motifs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, "jaspar2020_core_vertebrates.jaspar")

pfms <- getMatrixSet(JASPAR2020,
                     list(collection = "CORE", tax_group = "vertebrates",
                          matrixtype = "PFM"))
message(sprintf("[export] %d JASPAR2020 CORE vertebrate motifs", length(pfms)))

con <- file(outfile, "w")
for (id in names(pfms)) {
  pf <- pfms[[id]]
  m  <- Matrix(pf)                      # 4 x L integer counts, rows A,C,G,T
  writeLines(sprintf(">%s\t%s", ID(pf), name(pf)), con)
  for (b in c("A", "C", "G", "T")) {
    writeLines(sprintf("%s  [ %s ]", b,
                       paste(format(m[b, ], width = 6), collapse = " ")), con)
  }
}
close(con)
message(sprintf("[export] wrote %s", outfile))
