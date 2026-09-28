#!/usr/bin/env Rscript
# =============================================================================
# motif_grammar/02_grammar_features.R — motif composition vs grammar features and models
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Turns the motif INSTANCES from 01_region_sets.R into (a) composition features and (b) grammar
# features, then asks the only question that matters:
#
#     Does motif SYNTAX predict class beyond motif PRESENCE?
#
# measured as delta AUC of (composition + grammar) over (composition alone),
# with CHROMOSOME-HELD-OUT cross-validation.
#
# Three design decisions that make the delta interpretable:
#
# 1. PWM REDUNDANCY MERGE. JASPAR CORE has 24 ETS motifs; one GGAA site matches
#    many of them, so raw instance counts measure PWM redundancy, not biology.
#    Instances of the same family within MERGE_BP are collapsed to one site
#    (best-scoring representative keeps its position and strand).
#
# 2. GRAMMAR AS FRACTIONS, NOT COUNTS. "pairs within 50 bp" rises with motif
#    count alone, so it would proxy composition. Every grammar feature is a
#    FRACTION of the pairs that exist (or a distance), so it is ~orthogonal to
#    abundance.
#
# 3. COMPOSITION IS IN BOTH MODELS. The baseline is not "no motifs" — it is
#    counts + max scores per family. So the delta isolates syntax.
#
# A near-zero delta is a real, reportable result: it means the dose taxonomy is
# explained by which factors bind, not how their sites are arranged.
#
# Inputs:  results/extended_analysis/motif_grammar/grammar_{primary,control}.rds (01_region_sets.R)
# Outputs: results/extended_analysis/motif_grammar/
#            <contrast>_{coefficients,model_summary}.csv, <contrast>_features.rds
#            (incl. stored CV predictions for 04_noise_floor.R), where <contrast> is
#            primary_buffered_vs_haplo or control_dependent_vs_independent;
#            grammar_model_summary_all.csv
# Usage:   Rscript extended_analysis/scripts/motif_grammar/02_grammar_features.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(glmnet)
  library(Matrix)
})

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate
slice  <- dplyr::slice

set.seed(42)

outdir <- file.path(paths$ext_results, "motif_grammar")
figdir <- file.path(paths$ext_figures, "motif_grammar")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

MERGE_BP  <- 10L    # collapse same-family instances within this distance
CLOSE_BP  <- 50L    # "co-occurring" pair window (composite-element scale)
PHASE_BP  <- 100L   # window for helical-phasing assessment
HELIX     <- 10.5   # bp per DNA helical turn

FAMS <- c("ETS", "RUNX", "T-box", "AP-1/bZIP", "NFkB", "TCF/LEF", "KLF/SP")
fam_tag <- function(f) gsub("[^A-Za-z0-9]", "", f)

# ---- 1. collapse PWM redundancy --------------------------------------------
# Within peak x family, cluster instances whose centres are within MERGE_BP and
# keep the best-scoring one. Without this, ETS "count" is really "how many ETS
# PWMs JASPAR happens to ship".
collapse_sites <- function(pos) {
  pos %>%
    arrange(peak_id, family, pos) %>%
    group_by(peak_id, family) %>%
    mutate(site = cumsum(c(TRUE, diff(pos) > MERGE_BP))) %>%
    group_by(peak_id, family, site) %>%
    slice_max(score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(peak_id, class, family, pos, strand, score)
}

# ---- 2. composition features -----------------------------------------------
comp_features <- function(sites, peaks) {
  n_fam <- sites %>% count(peak_id, family, name = "n") %>%
    pivot_wider(names_from = family, values_from = n, values_fill = 0)
  s_fam <- sites %>% group_by(peak_id, family) %>%
    summarize(s = max(score), .groups = "drop") %>%
    pivot_wider(names_from = family, values_from = s, values_fill = 0)

  for (f in FAMS) {
    if (!f %in% names(n_fam)) n_fam[[f]] <- 0L
    if (!f %in% names(s_fam)) s_fam[[f]] <- 0
  }
  n_fam <- n_fam %>% select(peak_id, all_of(FAMS))
  s_fam <- s_fam %>% select(peak_id, all_of(FAMS))
  names(n_fam)[-1] <- paste0("n_", fam_tag(FAMS))
  names(s_fam)[-1] <- paste0("score_", fam_tag(FAMS))

  out <- tibble(peak_id = peaks) %>%
    left_join(n_fam, by = "peak_id") %>%
    left_join(s_fam, by = "peak_id")
  out[is.na(out)] <- 0
  # log1p the counts: motif counts are heavily right-skewed
  for (f in FAMS) out[[paste0("n_", fam_tag(f))]] <-
      log1p(out[[paste0("n_", fam_tag(f))]])
  out
}

# ---- 3. grammar features ---------------------------------------------------
# For every family pair (A <= B) in each peak:
#   closefrac_A_B  fraction of all A x B pairs that sit within CLOSE_BP
#   mindist_A_B    minimum |centre distance| (WINDOW if the pair never co-occurs)
#   same_A_B       fraction of close pairs on the SAME strand (tandem)
#   conv_A_B       fraction of close pairs that are convergent (-> <-)
#   phase_A_B      fraction of pairs within PHASE_BP on the same helical face
#                  (|cos(2*pi*d/10.5)| high => same side of the helix)
grammar_features <- function(sites, peaks, window) {
  pairs_tbl <- expand.grid(A = FAMS, B = FAMS, stringsAsFactors = FALSE) %>%
    filter(A <= B)

  by_peak <- split(sites, sites$peak_id)

  calc_one <- function(d) {
    out <- list()
    for (k in seq_len(nrow(pairs_tbl))) {
      A <- pairs_tbl$A[k]; B <- pairs_tbl$B[k]
      a <- d[d$family == A, ]; b <- d[d$family == B, ]
      tag <- paste0(fam_tag(A), "_", fam_tag(B))
      if (!nrow(a) || !nrow(b) || (A == B && nrow(a) < 2)) {
        out[[paste0("closefrac_", tag)]] <- 0
        out[[paste0("mindist_",  tag)]] <- window
        out[[paste0("same_",     tag)]] <- 0
        out[[paste0("conv_",     tag)]] <- 0
        out[[paste0("phase_",    tag)]] <- 0
        next
      }
      ii <- expand.grid(i = seq_len(nrow(a)), j = seq_len(nrow(b)))
      if (A == B) ii <- ii[ii$i < ii$j, , drop = FALSE]   # unordered self-pairs
      dd <- b$pos[ii$j] - a$pos[ii$i]
      ad <- abs(dd)
      keep <- ad > 0
      ii <- ii[keep, , drop = FALSE]; dd <- dd[keep]; ad <- ad[keep]
      if (!length(ad)) {
        out[[paste0("closefrac_", tag)]] <- 0
        out[[paste0("mindist_",  tag)]] <- window
        out[[paste0("same_",     tag)]] <- 0
        out[[paste0("conv_",     tag)]] <- 0
        out[[paste0("phase_",    tag)]] <- 0
        next
      }
      sa <- a$strand[ii$i]; sb <- b$strand[ii$j]
      close <- ad <= CLOSE_BP
      near  <- ad <= PHASE_BP

      out[[paste0("closefrac_", tag)]] <- mean(close)
      out[[paste0("mindist_",  tag)]] <- min(ad)
      out[[paste0("same_",     tag)]] <-
        if (any(close)) mean(sa[close] == sb[close]) else 0
      # convergent: upstream motif on +, downstream on -  (-> <-)
      conv <- (sign(dd) > 0 & sa == "+" & sb == "-") |
              (sign(dd) < 0 & sb == "+" & sa == "-")
      out[[paste0("conv_", tag)]] <-
        if (any(close)) mean(conv[close]) else 0
      out[[paste0("phase_", tag)]] <-
        if (any(near)) mean(cos(2 * pi * ad[near] / HELIX) > 0.5) else 0
    }
    as_tibble(out)
  }

  message("  computing grammar features for ", length(by_peak), " peaks…")
  res <- bind_rows(lapply(by_peak, calc_one))
  res$peak_id <- names(by_peak)

  full <- tibble(peak_id = peaks) %>% left_join(res, by = "peak_id")
  # peaks with no motifs at all: neutral defaults
  for (cn in names(full)[-1]) {
    v <- full[[cn]]
    v[is.na(v)] <- if (grepl("^mindist_", cn)) window else 0
    full[[cn]] <- v
  }
  full
}

# ---- 4. model: does grammar beat composition? ------------------------------
# Chromosome-held-out CV. Random CV would leak: nearby/duplicated peaks and
# shared repeat families put near-identical sequences in train and test.
run_contrast <- function(tag, rds, positive) {
  message("=== ", tag, " ===")
  obj <- readRDS(rds)
  gr <- obj$gr; pos <- obj$pos; window <- obj$window

  meta <- tibble(
    peak_id = mcols(gr)$peak_id,
    class   = mcols(gr)$class,
    chr     = as.character(GenomicRanges::seqnames(gr)),
    gc      = mcols(gr)$gc,
    baseMean = mcols(gr)$baseMean
  ) %>% distinct(peak_id, .keep_all = TRUE)

  sites <- collapse_sites(pos)
  message("  instances ", nrow(pos), " -> merged sites ", nrow(sites),
          sprintf(" (%.1f%% collapsed by PWM redundancy)",
                  100 * (1 - nrow(sites) / nrow(pos))))

  X_comp <- comp_features(sites, meta$peak_id)
  X_gram <- grammar_features(sites, meta$peak_id, window)

  stopifnot(identical(X_comp$peak_id, meta$peak_id),
            identical(X_gram$peak_id, meta$peak_id))

  y <- as.integer(meta$class == positive)
  Mc <- as.matrix(X_comp[, -1])
  Mg <- as.matrix(X_gram[, -1])
  Mcg <- cbind(Mc, Mg)

  # drop zero-variance columns (glmnet cannot standardise them)
  Mc  <- Mc[,  apply(Mc,  2, sd) > 0, drop = FALSE]
  Mcg <- Mcg[, apply(Mcg, 2, sd) > 0, drop = FALSE]

  auc <- function(truth, sc) {
    r <- rank(sc); n1 <- sum(truth == 1); n0 <- sum(truth == 0)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }

  chrs <- unique(meta$chr)
  folds <- split(chrs, cut(seq_along(chrs), 5, labels = FALSE))  # 5 chr-blocks

  # preds/fold_id are retained (not just the scalar AUC) so the dAUC can be given
  # a noise floor downstream — a bare point estimate cannot support a null claim.
  eval_model <- function(M, label) {
    preds   <- rep(NA_real_, nrow(M))
    fold_id <- rep(NA_integer_, nrow(M))
    for (i in seq_along(folds)) {
      te <- meta$chr %in% folds[[i]]
      if (!any(te) || !any(!te)) next
      if (length(unique(y[!te])) < 2) next
      fit <- cv.glmnet(M[!te, , drop = FALSE], y[!te], family = "binomial",
                       alpha = 0.5, nfolds = 5, standardize = TRUE)
      preds[te] <- as.numeric(predict(fit, M[te, , drop = FALSE],
                                      s = "lambda.min", type = "link"))
      fold_id[te] <- i
    }
    a <- auc(y[!is.na(preds)], preds[!is.na(preds)])
    per_fold <- vapply(seq_along(folds), function(i) {
      idx <- which(fold_id == i & !is.na(preds))
      if (!length(idx)) return(NA_real_)
      auc(y[idx], preds[idx])
    }, numeric(1))
    message(sprintf("  %-28s AUC = %.4f", label, a))
    list(auc = a, preds = preds, fold_id = fold_id, per_fold = per_fold)
  }

  r_comp <- eval_model(Mc,  "composition only")
  r_cg   <- eval_model(Mcg, "composition + grammar")
  delta <- r_cg$auc - r_comp$auc
  message(sprintf("  %-28s dAUC = %+.4f", "GRAMMAR GAIN", delta))

  # full-data fit for coefficient inspection
  fit_full <- cv.glmnet(Mcg, y, family = "binomial", alpha = 0.5,
                        nfolds = 10, standardize = TRUE)
  co <- as.matrix(coef(fit_full, s = "lambda.min"))
  co_tbl <- tibble(feature = rownames(co), coef = co[, 1]) %>%
    filter(feature != "(Intercept)", coef != 0) %>%
    mutate(type = ifelse(grepl("^n_|^score_", feature), "composition", "grammar")) %>%
    arrange(desc(abs(coef)))

  write_csv(co_tbl, file.path(outdir, paste0(tag, "_coefficients.csv")))

  res <- tibble(contrast = tag, positive_class = positive,
                n = nrow(M <- Mcg), n_pos = sum(y), n_neg = sum(y == 0),
                auc_composition = r_comp$auc, auc_comp_grammar = r_cg$auc,
                delta_auc = delta,
                n_grammar_selected = sum(co_tbl$type == "grammar"),
                n_composition_selected = sum(co_tbl$type == "composition"))
  write_csv(res, file.path(outdir, paste0(tag, "_model_summary.csv")))

  # keep features for downstream inspection / plotting
  saveRDS(list(meta = meta, sites = sites, X_comp = X_comp, X_gram = X_gram,
               y = y, result = res, coefs = co_tbl,
               preds_comp = r_comp$preds, preds_cg = r_cg$preds,
               fold_id = r_comp$fold_id,
               per_fold_comp = r_comp$per_fold, per_fold_cg = r_cg$per_fold),
          file.path(outdir, paste0(tag, "_features.rds")))
  res
}

r1 <- run_contrast("primary_buffered_vs_haplo",
                   file.path(outdir, "grammar_primary.rds"), "buffered")
r2 <- run_contrast("control_dependent_vs_independent",
                   file.path(outdir, "grammar_control.rds"), "cBAF_dependent")

summary_all <- bind_rows(r1, r2)
write_csv(summary_all, file.path(outdir, "grammar_model_summary_all.csv"))
message("=== motif_grammar/02_grammar_features.R complete ===")
print(as.data.frame(summary_all))
