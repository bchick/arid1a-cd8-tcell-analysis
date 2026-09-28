#!/usr/bin/env Rscript
# =============================================================================
# qc_audit/qc_audit.R — retrospective QC sweep across every differential contrast
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Applies scripts/qc_checks.R to all existing DE/DA results in the project, to
# answer: "is anything else in here confounded the way g183618 was?"
#
# The design_integrity check is the one that matters: it re-derives the group
# structure from the data and asks whether the labels we analysed under are the
# dominant structure. On the g183618 pooled-WT design it FAILS; that failure is
# reproduced here deliberately as a positive control for the check itself.
#
# Inputs:  results/extended_analysis/guo2022/atacseq/differential/checkpoint_sections1to3.RData
#          results/extended_analysis/guo2022/rnaseq/differential/guo_rnaseq_analysis.RData
#          results/atac/differential/atacseq_analysis.RData
#          results/rnaseq/differential/rnaseq_analysis.RData
# Outputs: results/extended_analysis/qc_audit/qc_audit_report.{csv,txt}
# Usage:   Rscript extended_analysis/scripts/qc_audit/qc_audit.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/qc_checks.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

outdir <- file.path(paths$ext_results, "qc_audit")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

all_qc <- list()
note   <- function(...) message(sprintf(...))

add <- function(x) if (!is.null(x) && nrow(x)) all_qc[[length(all_qc) + 1]] <<- x

# =============================================================================
# A. Guo ATAC (guo2022/02_atacseq_da.R) — per-contrast, from the saved checkpoint
# =============================================================================
note("=== A. Guo ATAC contrasts ===")
ck <- file.path(paths$ext_results, "guo2022/atacseq/differential/checkpoint_sections1to3.RData")
if (file.exists(ck)) {
  e <- new.env(); load(ck, envir = e)
  source("scripts/utils.R")   # RData overwrite bug: restore palettes/helpers
  vsd <- e$vsd; meta <- e$meta; da <- e$da_results; cnts <- e$all_counts
  vm <- SummarizedExperiment::assay(vsd)

  specs <- list(
    list(n = "KO_vs_WT",        g = "g183618", a = "KO",    b = "WT"),
    list(n = "MycKO_vs_MycWT",  g = "g183618", a = "MycKO", b = "MycWT"),
    list(n = "vehKO_vs_vehWT",  g = "g198894", a = "vehKO", b = "vehWT"),
    list(n = "inhib_vs_vehWT",  g = "g198894", a = "inhib", b = "vehWT"),
    list(n = "MycHi_vs_MycLo",  g = "g183616", a = "MycHi", b = "MycLo"),
    list(n = "MycHi_vs_Naive",  g = "g183616", a = "MycHi", b = "Naive")
  )
  for (s in specs) {
    keep <- meta$condition %in% c(s$a, s$b)
    if (sum(keep) < 4) next
    add(qc_run_all(mat = vm[, keep, drop = FALSE],
                   groups = as.character(meta$condition)[keep],
                   res = da[[s$n]], counts = cnts[, keep, drop = FALSE],
                   label = paste0("GuoATAC:", s$n)))
    note("  checked %s", s$n)
  }

  # --- positive control for the CHECK: reproduce the known-bad pooled design ---
  bad <- meta$condition %in% c("KO", "WT", "MycWT")
  grp_bad <- ifelse(as.character(meta$condition)[bad] == "KO", "KO", "WT")  # re-pool
  add(qc_run_all(mat = vm[, bad, drop = FALSE], groups = grp_bad,
                 label = "GuoATAC:KO_vs_pooledWT[KNOWN-BAD]"))
  note("  checked KO_vs_pooledWT (known-bad positive control)")
} else {
  note("  checkpoint not found, skipping Guo ATAC")
}

# =============================================================================
# B. Guo RNA (guo2022/01_rnaseq_de.R)
# =============================================================================
note("=== B. Guo RNA contrasts ===")
rd <- file.path(paths$ext_results, "guo2022/rnaseq/differential/guo_rnaseq_analysis.RData")
if (file.exists(rd)) {
  e <- new.env(); load(rd, envir = e); source("scripts/utils.R")
  vobj <- e$vsd %||% e$vst %||% NULL
  meta_r <- e$meta %||% e$coldata %||% NULL
  de <- e$de_results %||% NULL
  if (!is.null(vobj) && !is.null(meta_r) && !is.null(de)) {
    vm <- SummarizedExperiment::assay(vobj)
    gcol <- if ("condition" %in% names(meta_r)) "condition" else names(meta_r)[2]
    pairs <- list(c("KO", "WT"), c("MycKO", "MycWT"), c("inhib", "WT"))
    nms   <- c("KO_vs_WT", "MycKO_vs_MycWT", "inhib_vs_WT")
    for (i in seq_along(pairs)) {
      keep <- as.character(meta_r[[gcol]]) %in% pairs[[i]]
      if (sum(keep) < 4) next
      pc <- if (nms[i] == "KO_vs_WT") list(gene = "Arid1a", dir = "down") else NULL
      add(qc_run_all(mat = vm[, keep, drop = FALSE],
                     groups = as.character(meta_r[[gcol]])[keep],
                     res = de[[nms[i]]], label = paste0("GuoRNA:", nms[i]),
                     positive_control = pc))
      note("  checked %s", nms[i])
    }
  } else note("  RData present but expected objects missing, skipping")
} else {
  note("  Guo RNA RData not found, skipping")
}

# =============================================================================
# C. Main project (McDonald) — ATAC + RNA, if their RData are present
# =============================================================================
note("=== C. Main project contrasts ===")
for (cfg in list(
  list(f = file.path(paths$results, "atac/differential/atacseq_analysis.RData"), tag = "McDonaldATAC"),
  list(f = file.path(paths$results, "rnaseq/differential/rnaseq_analysis.RData"), tag = "McDonaldRNA")
)) {
  if (!file.exists(cfg$f)) { note("  %s: RData not found, skipping", cfg$tag); next }
  e <- new.env(); load(cfg$f, envir = e); source("scripts/utils.R")
  vobj <- e$vsd %||% NULL
  md   <- e$meta %||% e$sample_meta %||% e$coldata %||% NULL
  if (is.null(vobj) || is.null(md)) { note("  %s: objects missing, skipping", cfg$tag); next }
  vm <- SummarizedExperiment::assay(vobj)
  md <- as.data.frame(md)

  # The McDonald design is a timecourse x genotype x cell-subset factorial. The
  # dominant structure is timepoint/subset, NOT genotype — grouping by genotype
  # alone would (correctly) report that genotype isn't the dominant axis, which
  # says nothing about whether the analyses are sound. Use the FULL biological
  # group, which is what the contrasts in 05/06 actually compare.
  g <- if ("group" %in% names(md) && dplyr::n_distinct(md$group) > 1) {
    as.character(md$group)
  } else {
    parts <- intersect(c("genotype", "timepoint", "cell_subset", "treatment_clean"), names(md))
    do.call(paste, c(lapply(parts, function(p) as.character(md[[p]])), sep = "_"))
  }
  add(qc_run_all(mat = vm, groups = g, label = paste0(cfg$tag, ":full_group")))
  note("  checked %s (%d samples, %d groups)", cfg$tag, ncol(vm), dplyr::n_distinct(g))

  # Genotype effect WITHIN each stratum (timepoint x subset) — this is the
  # comparison the analyses actually make, so it is the one worth auditing.
  if (all(c("genotype", "timepoint") %in% names(md))) {
    strat <- if ("cell_subset" %in% names(md)) {
      paste(md$timepoint, md$cell_subset, sep = "_")
    } else as.character(md$timepoint)
    for (s in unique(strat)) {
      k <- strat == s
      if (sum(k) < 4 || dplyr::n_distinct(md$genotype[k]) < 2) next
      add(qc_run_all(mat = vm[, k, drop = FALSE],
                     groups = as.character(md$genotype)[k],
                     label = sprintf("%s:%s[genotype]", cfg$tag, s)))
    }
    note("  checked %s per-stratum genotype effects", cfg$tag)
  }
}

# =============================================================================
# D. Report
# =============================================================================
if (!length(all_qc)) {
  note("No contrasts could be audited.")
  quit(status = 0)
}
report <- dplyr::bind_rows(all_qc)
readr::write_csv(report, file.path(outdir, "qc_audit_report.csv"))

con <- file(file.path(outdir, "qc_audit_report.txt"), open = "wt")
sink(con); sink(con, type = "message")
cat("QC AUDIT — all differential contrasts\n")
cat("=====================================\n")
for (ct in unique(report$contrast)) qc_print(filter(report, contrast == ct), ct)
cat("\n================ SUMMARY ================\n")
summ <- report %>% group_by(contrast) %>%
  summarise(FAIL = sum(status == "FAIL"), WARN = sum(status == "WARN"),
            PASS = sum(status == "PASS"), .groups = "drop") %>%
  arrange(desc(FAIL), desc(WARN))
print(as.data.frame(summ), row.names = FALSE)
sink(type = "message"); sink(); close(con)

for (ct in unique(report$contrast)) qc_print(filter(report, contrast == ct), ct)
cat("\n================ SUMMARY ================\n")
print(as.data.frame(summ), row.names = FALSE)
note("\nReport -> %s", file.path(outdir, "qc_audit_report.{csv,txt}"))
