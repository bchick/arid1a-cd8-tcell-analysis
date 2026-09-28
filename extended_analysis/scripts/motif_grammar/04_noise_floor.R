#!/usr/bin/env Rscript
# =============================================================================
# motif_grammar/04_noise_floor.R — noise floor for the grammar dAUC, and summary figures
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# WHY THIS EXISTS. 02_grammar_features.R reports a single dAUC point estimate (grammar gain
# over composition) from pooled out-of-fold predictions. "Grammar adds nothing"
# is a NULL claim, and a null claim with no error bar is not evidence — +0.0012
# is only meaningful against the scale of the noise. This script supplies it:
#
#   1. PER-FOLD dAUC — the 5 chromosome-block deltas. Cheap; shows fold spread.
#   2. PAIRED BOOTSTRAP CI — resample peaks, recompute BOTH AUCs on the same
#      resample (paired, so the shared model error cancels), 95% CI on dAUC.
#      Uses the stored out-of-fold predictions, so no refitting.
#   3. CHROMOSOME-CLUSTER BOOTSTRAP — resample whole chromosomes instead of
#      peaks, respecting the dependence that chromosome-held-out CV was designed
#      around. Reported alongside as a robustness check. Note it is not uniformly
#      wider than the peak bootstrap here (narrower for the primary, wider for
#      the control): with ~20 clusters its width tracks between-chromosome dAUC
#      heterogeneity, not peak-level sampling noise. Both are reported; the
#      conclusion does not depend on which is used.
#
# This bootstraps the EVALUATION given fitted models; it does not resample model
# fitting itself (that needs a permutation null, ~60 min/contrast — deferred).
# For a null claim the evaluation CI is the relevant statement: if it spans 0,
# the pooled dAUC is indistinguishable from no gain.
#
# ALSO: composition enrichment, tested two ways.
#   COUNT   mean sites per peak, positive vs negative class.
#   SHARE   each family's fraction of that peak's total sites.
# The share test matters because the control has EVERY family enriched
# (ratios 1.07-1.87). Counts alone cannot separate "specifically ETS-enriched"
# from "globally motif-dense"; share can.
#
# Inputs:  results/extended_analysis/motif_grammar/<contrast>_features.rds (02_grammar_features.R)
# Outputs: results/extended_analysis/motif_grammar/
#            grammar_dauc_noise_floor.csv, grammar_dauc_per_fold.csv,
#            grammar_composition_enrichment.csv
#          figures/extended_analysis/motif_grammar/F1-F5_*.{pdf,png}
# Usage:   Rscript extended_analysis/scripts/motif_grammar/04_noise_floor.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(viridisLite)
})

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate
slice  <- dplyr::slice

set.seed(42)

outdir <- file.path(paths$ext_results, "motif_grammar")
figdir <- file.path(paths$ext_figures, "motif_grammar")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

B <- 2000L   # bootstrap resamples

FAMS <- c("ETS", "RUNX", "T-box", "AP-1/bZIP", "NFkB", "TCF/LEF", "KLF/SP")
fam_tag <- function(f) gsub("[^A-Za-z0-9]", "", f)

CONTRASTS <- list(
  list(tag = "primary_buffered_vs_haplo",
       label = "Primary: buffered vs haploinsufficient",
       pos = "buffered", neg = "haploinsufficient"),
  list(tag = "control_dependent_vs_independent",
       label = "Control: cBAF-dependent vs -independent",
       pos = "cBAF-dependent", neg = "cBAF-independent")
)

auc <- function(truth, sc) {
  r <- rank(sc); n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- 1. noise floor ---------------------------------------------------------
noise_floor <- function(cn) {
  o <- readRDS(file.path(outdir, paste0(cn$tag, "_features.rds")))
  if (is.null(o$preds_comp))
    stop("no stored predictions in ", cn$tag, "_features.rds — rerun ",
         "extended_analysis/scripts/motif_grammar/02_grammar_features.R (it stores out-of-fold predictions)")

  keep <- !is.na(o$preds_comp) & !is.na(o$preds_cg)
  y  <- o$y[keep]
  pc <- o$preds_comp[keep]
  pg <- o$preds_cg[keep]
  ch <- o$meta$chr[keep]
  n  <- length(y)

  d_obs <- auc(y, pg) - auc(y, pc)

  # paired peak-level bootstrap
  boot_peak <- vapply(seq_len(B), function(b) {
    i <- sample.int(n, n, replace = TRUE)
    if (length(unique(y[i])) < 2) return(NA_real_)
    auc(y[i], pg[i]) - auc(y[i], pc[i])
  }, numeric(1))

  # chromosome-cluster bootstrap
  chrs <- unique(ch)
  idx_by_chr <- split(seq_len(n), ch)
  boot_chr <- vapply(seq_len(B), function(b) {
    cs <- sample(chrs, length(chrs), replace = TRUE)
    i <- unlist(idx_by_chr[cs], use.names = FALSE)
    if (length(unique(y[i])) < 2) return(NA_real_)
    auc(y[i], pg[i]) - auc(y[i], pc[i])
  }, numeric(1))

  ci_p <- quantile(boot_peak, c(0.025, 0.975), na.rm = TRUE)
  ci_c <- quantile(boot_chr,  c(0.025, 0.975), na.rm = TRUE)

  per_fold <- tibble(
    contrast = cn$label, fold = seq_along(o$per_fold_comp),
    auc_comp = o$per_fold_comp, auc_cg = o$per_fold_cg,
    d = o$per_fold_cg - o$per_fold_comp)

  list(
    summary = tibble(
      contrast = cn$label, tag = cn$tag, n = n,
      auc_comp = auc(y, pc), auc_cg = auc(y, pg), delta_auc = d_obs,
      ci_lo_peak = ci_p[[1]], ci_hi_peak = ci_p[[2]],
      ci_lo_chr  = ci_c[[1]], ci_hi_chr  = ci_c[[2]],
      spans_zero_peak = ci_p[[1]] <= 0 && ci_p[[2]] >= 0,
      spans_zero_chr  = ci_c[[1]] <= 0 && ci_c[[2]] >= 0,
      fold_d_min = min(per_fold$d, na.rm = TRUE),
      fold_d_max = max(per_fold$d, na.rm = TRUE)),
    per_fold = per_fold,
    boot = tibble(contrast = cn$label,
                  d_peak = boot_peak, d_chr = boot_chr))
}

nf <- lapply(CONTRASTS, noise_floor)
nf_summary <- bind_rows(lapply(nf, `[[`, "summary"))
nf_perfold <- bind_rows(lapply(nf, `[[`, "per_fold"))
nf_boot    <- bind_rows(lapply(nf, `[[`, "boot"))

write_csv(nf_summary, file.path(outdir, "grammar_dauc_noise_floor.csv"))
write_csv(nf_perfold, file.path(outdir, "grammar_dauc_per_fold.csv"))

message("=== dAUC noise floor ===")
print(as.data.frame(nf_summary %>%
  select(contrast, n, auc_comp, auc_cg, delta_auc,
         ci_lo_peak, ci_hi_peak, spans_zero_peak)))

# ---- 2. composition enrichment: count and share -----------------------------
enrich <- function(cn) {
  o <- readRDS(file.path(outdir, paste0(cn$tag, "_features.rds")))
  X <- o$X_comp; y <- o$y
  cols <- paste0("n_", vapply(FAMS, fam_tag, character(1)))
  present <- cols %in% names(X)
  cols <- cols[present]; fams <- FAMS[present]

  M <- as.matrix(X[, cols, drop = FALSE])
  total <- rowSums(M)

  bind_rows(lapply(seq_along(fams), function(j) {
    v <- M[, j]
    m1 <- mean(v[y == 1]); m0 <- mean(v[y == 0])
    # share is only defined for peaks carrying >=1 site of any family
    ok <- total > 0
    s  <- v[ok] / total[ok]; yo <- y[ok]
    s1 <- mean(s[yo == 1]); s0 <- mean(s[yo == 0])

    bs <- vapply(seq_len(500L), function(b) {
      i1 <- sample(which(y == 1), sum(y == 1), replace = TRUE)
      i0 <- sample(which(y == 0), sum(y == 0), replace = TRUE)
      log2(mean(v[i1]) / max(mean(v[i0]), 1e-9))
    }, numeric(1))

    tibble(contrast = cn$label, family = fams[j],
           mean_pos = m1, mean_neg = m0,
           log2_ratio = log2(m1 / max(m0, 1e-9)),
           ci_lo = quantile(bs, 0.025, na.rm = TRUE)[[1]],
           ci_hi = quantile(bs, 0.975, na.rm = TRUE)[[1]],
           p_count = suppressWarnings(wilcox.test(v[y == 1], v[y == 0])$p.value),
           share_pos = s1, share_neg = s0,
           share_diff = s1 - s0,
           p_share = suppressWarnings(wilcox.test(s[yo == 1], s[yo == 0])$p.value))
  }))
}

en <- bind_rows(lapply(CONTRASTS, enrich)) %>%
  group_by(contrast) %>%
  mutate(fdr_count = p.adjust(p_count, "BH"),
         fdr_share = p.adjust(p_share, "BH")) %>%
  ungroup()

write_csv(en, file.path(outdir, "grammar_composition_enrichment.csv"))

message("\n=== composition enrichment (count vs share) ===")
print(as.data.frame(en %>%
  select(contrast, family, log2_ratio, fdr_count, share_pos, share_neg,
         share_diff, fdr_share) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))))

# ---- 3. figures -------------------------------------------------------------
fam_levels <- FAMS
mako7 <- viridisLite::mako(9)[2:8]
names(mako7) <- fam_levels

# Primary leads: it is the new contrast; the control is the validity check.
# Alphabetical faceting would otherwise put the control first everywhere.
ctr_levels <- vapply(CONTRASTS, `[[`, character(1), "label")
as_ctr <- function(x) factor(x, levels = ctr_levels)
nf_summary <- nf_summary %>% mutate(contrast = as_ctr(contrast))
nf_perfold <- nf_perfold %>% mutate(contrast = as_ctr(contrast))
nf_boot    <- nf_boot    %>% mutate(contrast = as_ctr(contrast))
en         <- en         %>% mutate(contrast = as_ctr(contrast))

# F1 — count enrichment by family
p1 <- en %>%
  mutate(family = factor(family, levels = fam_levels),
         sig = ifelse(fdr_count < 0.05, "FDR < 0.05", "n.s.")) %>%
  ggplot(aes(x = log2_ratio, y = family, fill = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_col(aes(alpha = sig), width = 0.7) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                 color = "gray20", linewidth = 0.4) +
  scale_fill_manual(values = mako7, guide = "none") +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "n.s." = 0.35), name = NULL) +
  facet_wrap(~ contrast, scales = "free_x") +
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.08)) +
  theme(panel.spacing.x = unit(1.1, "lines")) +
  labs(title = "Motif family composition, matched sets",
       subtitle = "log2 ratio of mean sites per 300 bp window (positive vs negative class); 95% bootstrap CI",
       x = expression(log[2]~"(mean sites ratio)"), y = NULL)
save_figure(p1, "F1_composition_enrichment", width = 10, height = 4.2, dir = figdir)

# F2 — dAUC with CI against zero
p2 <- nf_summary %>%
  ggplot(aes(x = delta_auc, y = contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = ci_lo_chr, xmax = ci_hi_chr), height = 0,
                 linewidth = 2.2, color = viridisLite::mako(5)[4], alpha = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo_peak, xmax = ci_hi_peak), height = 0.12,
                 linewidth = 0.7, color = viridisLite::mako(5)[2]) +
  geom_point(size = 3.2, color = viridisLite::mako(5)[2]) +
  labs(title = "Does motif syntax predict class beyond motif presence?",
       subtitle = paste0("dAUC = AUC(composition + grammar) - AUC(composition), chromosome-held-out.\n",
                         "Thin bar: paired peak bootstrap 95% CI. Thick pale bar: chromosome-cluster 95% CI. ",
                         B, " resamples."),
       x = expression(Delta~"AUC (grammar gain)"), y = NULL)
save_figure(p2, "F2_dauc_noise_floor", width = 9, height = 3.4, dir = figdir)

# F3 — per-fold paired AUC
p3 <- nf_perfold %>%
  tidyr::pivot_longer(c(auc_comp, auc_cg), names_to = "model", values_to = "auc") %>%
  mutate(model = recode(model, auc_comp = "composition",
                        auc_cg = "composition + grammar")) %>%
  ggplot(aes(x = model, y = auc, group = fold)) +
  geom_line(color = "gray55", linewidth = 0.4) +
  geom_point(aes(color = model), size = 2.4) +
  scale_color_manual(values = unname(viridisLite::mako(5)[c(2, 4)]), guide = "none") +
  facet_wrap(~ contrast) +
  labs(title = "Per-fold AUC, composition vs composition + grammar",
       subtitle = "Each line is one held-out chromosome block (5 blocks)",
       x = NULL, y = "AUC (held-out)")
save_figure(p3, "F3_per_fold_auc", width = 9, height = 4, dir = figdir)

# F4 — density-normalised share
p4 <- en %>%
  mutate(family = factor(family, levels = fam_levels),
         sig = ifelse(fdr_share < 0.05, "FDR < 0.05", "n.s.")) %>%
  ggplot(aes(x = share_diff, y = family, fill = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_col(aes(alpha = sig), width = 0.7) +
  scale_fill_manual(values = mako7, guide = "none") +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "n.s." = 0.35), name = NULL) +
  facet_wrap(~ contrast, scales = "free_x") +
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.08)) +
  theme(panel.spacing.x = unit(1.1, "lines")) +
  labs(title = "Density-normalised composition: each family's share of a window's sites",
       subtitle = "Removes global motif-density differences that inflate every family's count",
       x = "Difference in mean share (positive - negative class)", y = NULL)
save_figure(p4, "F4_composition_share", width = 10, height = 4.2, dir = figdir)

# F5 — bootstrap null distribution of dAUC
p5 <- nf_boot %>%
  ggplot(aes(x = d_peak)) +
  geom_histogram(bins = 60, fill = viridisLite::mako(5)[3], color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(data = nf_summary, aes(xintercept = delta_auc),
             color = "#CC6677", linewidth = 0.8) +
  facet_wrap(~ contrast, scales = "free") +
  labs(title = "Bootstrap distribution of the grammar gain",
       subtitle = "Rose line: observed dAUC. Dashed: zero. Paired peak-level bootstrap.",
       x = expression(Delta~"AUC"), y = "bootstrap resamples")
save_figure(p5, "F5_bootstrap_distribution", width = 10, height = 4, dir = figdir)

message("\n=== motif_grammar/04_noise_floor.R complete. Figures in ", figdir, " ===")
