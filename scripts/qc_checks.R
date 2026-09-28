#!/usr/bin/env Rscript
# =============================================================================
# qc_checks.R — Reusable sanity / QC checks for differential analyses
# McDonald, Chick et al. 2023 Immunity 56:1303 — shared helpers
#
# Motivation: pooling the two experiments deposited under one "WT" label in
# Guo et al. 2022 GSE183618 produced a confident, biologically plausible and
# completely spurious result: a 33:1 accessibility collapse that was actually
# a cross-batch comparison. Every output looked reasonable. These checks are
# designed so that class of error announces itself.
#
# Example:
#   source("scripts/qc_checks.R")
#   qc <- qc_run_all(mat = assay(vsd), groups = meta$condition,
#                    res = da_results$KO_vs_WT, dds = dds_sub,
#                    label = "KO_vs_WT",
#                    positive_control = list(gene = "Arid1a", dir = "down"))
#   qc_print(qc); qc_stop_if_fail(qc)   # optional hard gate
#
# Every check returns a row: check / status (PASS|WARN|FAIL) / value / detail.
# Checks are advisory by default; qc_stop_if_fail() is opt-in.
#
# Inputs:  in-memory matrices / DESeq2 objects passed by the caller
# Outputs: tibbles of check results (nothing is written to disk)
# Sourced by: extended_analysis/scripts/qc_audit/qc_audit.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

qc_row <- function(check, status, value, detail) {
  tibble(check = check, status = status, value = value, detail = detail)
}

# -----------------------------------------------------------------------------
# 1. DESIGN INTEGRITY — the flagship check.
#
# For each sample, find its nearest neighbour by correlation. If a sample's
# closest relative is in a DIFFERENT group, the group labels do not describe
# the dominant structure in the data — usually a batch confound or a swap.
# Also compares mean within-group vs between-group correlation: if between
# >= within, the contrast is measuring batch, not biology.
#
# On the g183618 pooled-WT data this fires immediately (the two WT batches
# correlate with each other LESS than KO correlates with its own-arm WT).
# -----------------------------------------------------------------------------
qc_design_integrity <- function(mat, groups, label = "") {
  groups <- as.character(groups)
  stopifnot(ncol(mat) == length(groups))
  if (ncol(mat) < 4 || length(unique(groups)) < 2) {
    return(qc_row("design_integrity", "SKIP", NA_real_,
                  "needs >=4 samples and >=2 groups"))
  }
  cc <- suppressWarnings(cor(mat, use = "pairwise.complete.obs"))
  diag(cc) <- NA

  # nearest neighbour per sample
  nn      <- colnames(cc)[apply(cc, 2, which.max)]
  nn_grp  <- groups[match(nn, colnames(cc))]
  misfits <- colnames(cc)[nn_grp != groups]

  same <- outer(groups, groups, "==")
  within  <- mean(cc[same & !is.na(cc)], na.rm = TRUE)
  between <- mean(cc[!same & !is.na(cc)], na.rm = TRUE)
  gap <- within - between

  # Severity. A negative gap means the labels do not describe the dominant
  # structure at all (batch or confound is winning) -> FAIL. A healthy positive
  # gap with a couple of stragglers is normal in large factorial designs
  # (and expected where a group is biologically intermediate, e.g. Het in a
  # dose-response) -> WARN, not FAIL. Only when misfits are a large share of
  # samples does a positive gap still indicate a broken contrast.
  frac_misfit <- length(misfits) / ncol(cc)
  status <- if (gap <= 0) "FAIL"
            else if (frac_misfit > 0.25) "FAIL"
            else if (length(misfits) > 0 || gap < 0.005) "WARN"
            else "PASS"
  detail <- sprintf("within=%.4f between=%.4f gap=%+.4f%s",
                    within, between, gap,
                    if (length(misfits))
                      sprintf("; %d/%d (%.0f%%) nearest-neighbour outside own group: %s",
                              length(misfits), ncol(cc), 100 * frac_misfit,
                              paste(misfits, collapse = ", ")) else "")
  qc_row(paste0("design_integrity", if (nzchar(label)) paste0(":", label) else ""),
         status, gap, detail)
}

# -----------------------------------------------------------------------------
# 2. REPLICATE COHESION — flag replicates that don't belong with their group.
# -----------------------------------------------------------------------------
qc_replicate_cohesion <- function(mat, groups, z_thresh = 2) {
  groups <- as.character(groups)
  cc <- suppressWarnings(cor(mat, use = "pairwise.complete.obs")); diag(cc) <- NA
  coh <- vapply(seq_along(groups), function(i) {
    peers <- which(groups == groups[i]); peers <- setdiff(peers, i)
    if (!length(peers)) NA_real_ else mean(cc[i, peers], na.rm = TRUE)
  }, numeric(1))
  names(coh) <- colnames(mat)
  if (all(is.na(coh))) return(qc_row("replicate_cohesion", "SKIP", NA_real_, "no replicates"))
  z <- (coh - mean(coh, na.rm = TRUE)) / (sd(coh, na.rm = TRUE) %||% 1)
  bad <- names(coh)[!is.na(z) & z < -z_thresh]
  status <- if (length(bad)) "WARN" else "PASS"
  qc_row("replicate_cohesion", status, min(coh, na.rm = TRUE),
         sprintf("min within-group r=%.4f (%s)%s",
                 min(coh, na.rm = TRUE), names(which.min(coh)),
                 if (length(bad)) paste0("; outliers: ", paste(bad, collapse = ", ")) else ""))
}

# -----------------------------------------------------------------------------
# 3. POSITIVE CONTROL — the perturbed target must move the expected way.
#   res: data.frame with gene_name/log2FoldChange/padj (or peak-level + gene_name)
# -----------------------------------------------------------------------------
qc_positive_control <- function(res, gene, dir = c("down", "up"), alpha = 0.05) {
  dir <- match.arg(dir)
  gcol <- intersect(c("gene_name", "gene", "symbol", "SYMBOL"), names(res))
  if (!length(gcol)) return(qc_row("positive_control", "SKIP", NA_real_, "no gene column"))
  hit <- res[!is.na(res[[gcol[1]]]) & res[[gcol[1]]] == gene, , drop = FALSE]
  if (!nrow(hit)) return(qc_row(paste0("positive_control:", gene), "SKIP", NA_real_,
                                "gene not present in results"))
  hit <- hit[which.min(hit$padj), ]
  lfc <- hit$log2FoldChange; p <- hit$padj
  right_dir <- (dir == "down" && lfc < 0) || (dir == "up" && lfc > 0)
  status <- if (!right_dir) "FAIL" else if (is.na(p) || p >= alpha) "WARN" else "PASS"
  qc_row(paste0("positive_control:", gene), status, lfc,
         sprintf("LFC=%+.3f padj=%.3g (expected %s)", lfc, p, dir))
}

# -----------------------------------------------------------------------------
# 4. P-VALUE DISTRIBUTION — under a well-calibrated test the null part of the
# histogram is flat. A U-shape or a spike near 1 signals a misspecified model
# (unmodelled batch, wrong dispersion), which is exactly what a confounded
# design produces.
# -----------------------------------------------------------------------------
qc_pvalue_distribution <- function(res) {
  p <- res$pvalue[!is.na(res$pvalue)]
  if (length(p) < 1000) return(qc_row("pvalue_distribution", "SKIP", NA_real_, "too few tests"))
  # density of the upper half relative to uniform expectation
  hi <- mean(p > 0.5) * 2          # ~1.0 if the null half is flat
  top <- mean(p > 0.9) * 10        # ~1.0 if flat at the very top
  status <- if (hi < 0.5 || hi > 1.6 || top > 2) "WARN" else "PASS"
  qc_row("pvalue_distribution", status, hi,
         sprintf("upper-half density=%.2f, top-decile density=%.2f (1.0 = uniform)", hi, top))
}

# -----------------------------------------------------------------------------
# 5. EFFECT SYMMETRY — with median-of-ratios normalisation the median LFC is
# ~0 by construction. A large deviation means normalisation is fighting the
# data (or a real global shift that this normalisation cannot represent).
# -----------------------------------------------------------------------------
qc_effect_symmetry <- function(res, tol = 0.15) {
  lfc <- res$log2FoldChange[!is.na(res$log2FoldChange)]
  if (!length(lfc)) return(qc_row("effect_symmetry", "SKIP", NA_real_, "no LFCs"))
  m <- median(lfc)
  status <- if (abs(m) > tol * 2) "WARN" else "PASS"
  qc_row("effect_symmetry", status, m,
         sprintf("median LFC=%+.3f, %.1f%% negative", m, 100 * mean(lfc < 0)))
}

# -----------------------------------------------------------------------------
# 6. DISPERSION — compare against a reference if one is supplied, else report.
# -----------------------------------------------------------------------------
qc_dispersion <- function(dds, reference = NULL, ratio_warn = 1.8) {
  d <- tryCatch(median(DESeq2::dispersions(dds), na.rm = TRUE), error = function(e) NA_real_)
  if (is.na(d)) return(qc_row("dispersion", "SKIP", NA_real_, "unavailable"))
  status <- "PASS"; det <- sprintf("median dispersion=%.4f", d)
  if (!is.null(reference) && is.finite(reference)) {
    r <- d / reference
    if (r > ratio_warn) status <- "WARN"
    det <- sprintf("%s (%.2fx reference %.4f)", det, r, reference)
  }
  qc_row("dispersion", status, d, det)
}

# -----------------------------------------------------------------------------
# 7. LIBRARY SANITY — depth spread and size-factor range.
# -----------------------------------------------------------------------------
qc_library_sanity <- function(counts, sf = NULL) {
  d <- colSums(counts)
  cv <- sd(d) / mean(d)
  status <- if (cv > 0.75) "WARN" else "PASS"
  det <- sprintf("depth CV=%.2f (range %.1fM-%.1fM)", cv, min(d)/1e6, max(d)/1e6)
  if (!is.null(sf)) {
    rng <- max(sf) / min(sf)
    if (rng > 5) status <- "WARN"
    det <- sprintf("%s; size-factor range=%.2fx", det, rng)
  }
  qc_row("library_sanity", status, cv, det)
}

# -----------------------------------------------------------------------------
# Runner + reporting
# -----------------------------------------------------------------------------
qc_run_all <- function(mat = NULL, groups = NULL, res = NULL, dds = NULL,
                       counts = NULL, label = "", positive_control = NULL,
                       dispersion_reference = NULL) {
  out <- list()
  if (!is.null(mat) && !is.null(groups)) {
    out <- c(out, list(qc_design_integrity(mat, groups),
                       qc_replicate_cohesion(mat, groups)))
  }
  if (!is.null(res)) {
    out <- c(out, list(qc_pvalue_distribution(res), qc_effect_symmetry(res)))
    if (!is.null(positive_control)) {
      out <- c(out, list(qc_positive_control(res, positive_control$gene,
                                             positive_control$dir %||% "down")))
    }
  }
  if (!is.null(dds))    out <- c(out, list(qc_dispersion(dds, dispersion_reference)))
  if (!is.null(counts)) out <- c(out, list(qc_library_sanity(counts,
                          if (!is.null(dds)) tryCatch(DESeq2::sizeFactors(dds),
                                                      error = function(e) NULL) else NULL)))
  bind_rows(out) %>% mutate(contrast = label, .before = 1)
}

qc_print <- function(qc, title = "QC") {
  cat(sprintf("\n---- %s ----\n", title))
  for (i in seq_len(nrow(qc))) {
    mark <- switch(qc$status[i], PASS = "  ok  ", WARN = " WARN ", FAIL = " FAIL ", "  --  ")
    cat(sprintf("[%s] %-28s %s\n", mark, qc$check[i], qc$detail[i]))
  }
  n_f <- sum(qc$status == "FAIL"); n_w <- sum(qc$status == "WARN")
  cat(sprintf("---- %d FAIL, %d WARN, %d PASS ----\n\n",
              n_f, n_w, sum(qc$status == "PASS")))
  invisible(qc)
}

qc_stop_if_fail <- function(qc, allow = character()) {
  bad <- qc %>% filter(status == "FAIL", !check %in% allow)
  if (nrow(bad)) {
    stop(sprintf("QC FAILED (%d check(s)):\n%s", nrow(bad),
                 paste(sprintf("  - %s [%s]: %s", bad$contrast, bad$check, bad$detail),
                       collapse = "\n")), call. = FALSE)
  }
  invisible(TRUE)
}
