#!/usr/bin/env Rscript
# =============================================================================
# atac_trajectories/07_glmnet_motif_models.R — joint (elastic-net) motif models of ARID1A-dependent sites
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Which motifs predict ARID1A-dependent sites once all motifs are fitted
# together? SEA (atac_trajectories/04-05) and HOMER test each motif on its
# own, so correlated motifs (AP-1 family members, ETS paralogs, GC-rich
# promoter motifs) all come up together. An elastic net fits every motif at
# once and keeps the ones with predictive value of their own. Method follows
# an in-house glmnet motif-enrichment template: binary peak x motif matrix
# from motifmatchr on 200 bp windows, elastic net alpha = 0.5, 10-fold CV, with
# three changes for this data:
#   * Held-out evaluation. In-sample glmnet coefficients looked like signal in
#     the motif-grammar work and were not, so each model is fitted on a
#     stratified 80% split (lambda by cv.glmnet on that split only) and scored
#     by AUC on the other 20%, with bootstrap 95% CIs.
#   * Unpenalized covariates: GC fraction, CpG observed/expected, promoter
#     flag. SEA found the kept side promoter/CpG-rich; without these, glmnet
#     would rediscover "promoter" through GC-rich motifs. A covariates-only
#     model is the baseline, so delta-AUC is what motifs add beyond them.
#   * Stability selection: 100 stratified half-subsamples refitted at the
#     training lambda.1se; a motif is stable at selection frequency >= 0.8.
#
# Models
#   A  multinomial over the 5 WT trajectory classes of lost sites (pass 2)
#   B  binomial, lost vs kept, within each pass-1 WT class
# Motifs: JASPAR2024 CORE vertebrates non-redundant, the same 879 SEA used.
# HOMER known motifs are run on the same region sets for comparison (HOMER
# uses its own motif library and GC-matches its background).
#
# Inputs:  results/extended_analysis/atac_trajectories/motifs/{beds/,sea_all.tsv.gz}
#          (atac_trajectories/04_lost_site_motifs*.{sh,R}, 05_lost_site_motifs_plot.R)
#          results/atac/differential/consensus_peaks_annotated.csv (core/02_atacseq_analysis.R)
#          data/reference/{GRCm39.primary_assembly.genome.fa,jaspar2024_vert_nr.meme}
#            (only to build the match-matrix and HOMER caches; the HOMER
#            comparison also needs findMotifsGenome.pl on PATH)
# Outputs: results/extended_analysis/atac_trajectories/glmnet/
#          figures/extended_analysis/atac_trajectories/glmnet/
# Usage:   Rscript extended_analysis/scripts/atac_trajectories/07_glmnet_motif_models.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("extended_analysis/scripts/utils_trajectories.R")

suppressPackageStartupMessages({
  library(glmnet)
  library(Matrix)
  library(GenomicRanges)
  library(parallel)
})
select <- dplyr::select; filter <- dplyr::filter; count <- dplyr::count

ALPHA   <- 0.5
NFOLDS  <- 10
N_STAB  <- 100
STABLE  <- 0.8
N_BOOT  <- 1000
SEED    <- 42
CORES   <- as.integer(Sys.getenv("THREADS", "16"))
Q_SIG   <- 1e-5   # SEA / HOMER significance, as in 05_lost_site_motifs_plot.R
Q_MISS  <- 0.01   # "not called" by SEA / HOMER: q above this

motif_dir <- file.path(traj_dir, "motifs")
outdir    <- file.path(traj_dir, "glmnet")
homer_dir <- file.path(outdir, "homer")
figdir    <- file.path(paths$ext_figures, "atac_trajectories/glmnet")
dir.create(homer_dir, recursive = TRUE, showWarnings = FALSE)

fa_path   <- genome$fasta
meme_path <- file.path(paths$reference, "jaspar2024_vert_nr.meme")
cache     <- file.path(outdir, "motif_match_cache.rds")

class_levels <- c("Decreasing", "Transient", "Transient Increasing",
                  "Sustained Increasing", "Late Increasing")
slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))

# =============================================================================
# 1. Region sets (200 bp windows from 04_lost_site_motifs_prep.R)
# =============================================================================

read_bed <- function(set) read_tsv(file.path(motif_dir, "beds", paste0(set, ".bed")),
                                   col_names = c("chr", "start", "end", "feature_id"),
                                   show_col_types = FALSE) |> mutate(set = set)
sets <- c("static", paste0("lost_", slug(class_levels)),
          paste0("p1lost_", slug(class_levels)), paste0("p1kept_", slug(class_levels)))
regions <- bind_rows(lapply(sets, read_bed))
uniq <- distinct(regions, feature_id, chr, start, end)
message(sprintf("=== %d region-set rows, %d unique 200 bp windows ===", nrow(regions), nrow(uniq)))

# =============================================================================
# 2. Peak x motif matches + covariates (cached; needs the FASTA to build)
# =============================================================================

if (file.exists(fa_path) && file.exists(meme_path) && !file.exists(cache)) {
  suppressPackageStartupMessages({ library(motifmatchr); library(universalmotif); library(Rsamtools) })
  message("  Building motif match matrix (JASPAR2024, motifmatchr p < 5e-5)")
  mot  <- read_meme(meme_path)
  # count matrices, so motifmatchr builds the log-odds PWMs itself (as for the
  # template's JASPAR2020 PFMatrixList)
  pfms <- do.call(TFBSTools::PFMatrixList, convert_motifs(mot, class = "TFBSTools-PFMatrix"))
  names(pfms) <- vapply(mot, function(m) m@name, "")
  motif_info <- tibble(motif_id = vapply(mot, function(m) m@name, ""),
                       motif    = toupper(vapply(mot, function(m) m@altname, "")))
  gr <- GRanges(uniq$chr, IRanges(uniq$start + 1L, uniq$end), name = uniq$feature_id)
  fa <- FaFile(fa_path)
  mm <- matchMotifs(pfms, gr, genome = fa, p.cutoff = 5e-5)
  hits <- motifMatches(mm); rownames(hits) <- uniq$feature_id; colnames(hits) <- motif_info$motif_id
  seqs <- getSeq(fa, gr)
  nuc  <- Biostrings::letterFrequency(seqs, c("C", "G"))
  cg   <- Biostrings::vcountPattern("CG", seqs)
  L    <- Biostrings::width(seqs)
  anno <- read_csv(file.path(paths$atac, "differential/consensus_peaks_annotated.csv"),
                   show_col_types = FALSE)
  cov <- tibble(feature_id = uniq$feature_id,
                gc = rowSums(nuc) / L,
                cpg_oe = ifelse(nuc[, "C"] * nuc[, "G"] > 0, cg * L / (nuc[, "C"] * nuc[, "G"]), 0),
                promoter = as.numeric(grepl("^Promoter", anno$annotation[match(uniq$feature_id, anno$peak_id)])))
  saveRDS(list(hits = as(hits, "lMatrix"), motif_info = motif_info, cov = cov), cache)
}
if (!file.exists(cache)) {
  message("Skipping glmnet motif models: needs the genome FASTA + JASPAR2024 file, or ", cache)
  quit(save = "no", status = 0)
}
mc <- readRDS(cache)
stopifnot(all(uniq$feature_id %in% rownames(mc$hits)))
message(sprintf("  Match matrix: %d windows x %d motifs, %.1f hits per window",
                nrow(mc$hits), ncol(mc$hits), sum(mc$hits) / nrow(mc$hits)))

motif_name <- setNames(mc$motif_info$motif, mc$motif_info$motif_id)
cov_names  <- c("gc", "cpg_oe", "promoter")

design <- function(ids) {
  m <- as(mc$hits[ids, , drop = FALSE] * 1, "dgCMatrix")
  keep <- colSums(m) >= 20                       # drop motifs almost never matched
  cv <- as.matrix(mc$cov[match(ids, mc$cov$feature_id), cov_names])
  list(x = cbind(Matrix(cv, sparse = TRUE), m[, keep]), n_motif = sum(keep))
}

# =============================================================================
# 3. Held-out AUC, covariates vs covariates + motifs
# =============================================================================

auc <- function(score, pos) {                     # Mann-Whitney AUC
  r <- rank(score); n1 <- sum(pos); n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
strat_split <- function(y, frac, seed) {
  set.seed(seed)
  unlist(lapply(split(seq_along(y), y), function(i) i[sample.int(length(i), round(frac * length(i)))]))
}

fit_model <- function(name, ids, y, family) {
  message(sprintf("\n=== Model %s: %d windows, %s ===", name, length(ids),
                  paste(names(table(y)), table(y), sep = " ", collapse = ", ")))
  d <- design(ids)
  pf <- c(rep(0, length(cov_names)), rep(1, d$n_motif))
  tr <- strat_split(y, 0.8, SEED); te <- setdiff(seq_along(y), tr)
  set.seed(SEED)
  foldid <- integer(length(tr))
  for (lv in levels(y)) { i <- which(y[tr] == lv); foldid[i] <- sample(rep_len(seq_len(NFOLDS), length(i))) }
  doParallel::registerDoParallel(min(CORES, NFOLDS))
  cvf <- cv.glmnet(d$x[tr, ], y[tr], family = family, alpha = ALPHA, foldid = foldid,
                   penalty.factor = pf, type.measure = "deviance", parallel = TRUE)
  base <- glmnet(d$x[tr, cov_names], y[tr], family = family, lambda = 0)

  # lambda.1se is the reported model; lambda.min is a sensitivity check that a
  # zero delta-AUC is not just the more conservative lambda dropping every motif
  p_full <- predict(cvf, d$x[te, ], s = "lambda.1se", type = "response")
  p_min  <- predict(cvf, d$x[te, ], s = "lambda.min", type = "response")
  p_base <- predict(base, d$x[te, cov_names], type = "response")
  if (family == "binomial") {
    p_full <- cbind(p_full); p_min <- cbind(p_min); p_base <- cbind(p_base)
    colnames(p_full) <- colnames(p_min) <- colnames(p_base) <- levels(y)[2]
  } else { p_full <- p_full[, , 1]; p_min <- p_min[, , 1]; p_base <- p_base[, , 1] }
  targets <- colnames(p_full)
  yt <- y[te]

  set.seed(SEED)
  boots <- replicate(N_BOOT, unlist(lapply(split(seq_along(yt), yt), function(i) sample(i, replace = TRUE))),
                     simplify = FALSE)
  auc_tab <- bind_rows(lapply(targets, function(k) {
    pos <- yt == k
    a_f <- auc(p_full[, k], pos); a_b <- auc(p_base[, k], pos); a_m <- auc(p_min[, k], pos)
    bs <- vapply(boots, function(b) c(auc(p_full[b, k], pos[b]), auc(p_base[b, k], pos[b]),
                                      auc(p_min[b, k], pos[b])), numeric(3))
    tibble(model = name, target = k, n_test = length(yt), n_pos = sum(pos),
           auc_cov = a_b, auc_full = a_f, delta = a_f - a_b,
           auc_cov_lo = quantile(bs[2, ], 0.025), auc_cov_hi = quantile(bs[2, ], 0.975),
           auc_full_lo = quantile(bs[1, ], 0.025), auc_full_hi = quantile(bs[1, ], 0.975),
           delta_lo = quantile(bs[1, ] - bs[2, ], 0.025), delta_hi = quantile(bs[1, ] - bs[2, ], 0.975),
           delta_min = a_m - a_b,
           delta_min_lo = quantile(bs[3, ] - bs[2, ], 0.025), delta_min_hi = quantile(bs[3, ] - bs[2, ], 0.975))
  }))
  print(as.data.frame(mutate(auc_tab, across(where(is.numeric), ~ signif(.x, 3)))))

  # --- stability selection at the training lambda.1se ------------------------
  lam <- cvf$lambda.1se
  coefs_one <- function(b) {
    set.seed(SEED + b)
    sub <- unlist(lapply(split(seq_along(y), y), function(i) sample(i, floor(length(i) / 2))))
    f <- glmnet(d$x[sub, ], y[sub], family = family, alpha = ALPHA, lambda = lam, penalty.factor = pf)
    cf <- coef(f)
    if (family == "binomial") cf <- list(cf)
    names(cf) <- targets
    bind_rows(lapply(targets, function(k) {
      v <- as.numeric(cf[[k]][-1, 1]); nm <- rownames(cf[[k]])[-1]
      tibble(target = k, motif_id = nm, coef = v)
    })) |> filter(!motif_id %in% cov_names)
  }
  stab_raw <- bind_rows(mclapply(seq_len(N_STAB), coefs_one, mc.cores = CORES))
  stab <- stab_raw |>
    group_by(target, motif_id) |>
    summarise(sel_freq = mean(coef != 0), median_coef = median(coef[coef != 0]),
              frac_pos = mean(coef[coef != 0] > 0), .groups = "drop") |>
    mutate(model = name, motif = unname(motif_name[motif_id]),
           sign = case_when(is.na(median_coef) ~ NA_character_, median_coef > 0 ~ "+", TRUE ~ "-"))
  # strongest matrix per TF name, as in 05_lost_site_motifs_plot.R / 06_tf_expression_vs_motifs.R
  stab_tf <- stab |> group_by(model, target, motif) |>
    slice_max(order_by = sel_freq, n = 1, with_ties = FALSE) |> ungroup() |>
    mutate(stable = sel_freq >= STABLE)
  message(sprintf("  lambda.1se = %.4g; stable motifs (TF level): %s", lam,
                  paste(stab_tf |> filter(stable) |> count(target) |>
                          mutate(s = paste(target, n)) |> pull(s), collapse = "; ")))
  list(auc = auc_tab, stab = stab_tf, lambda = lam)
}

# Model A: multinomial across lost-site classes (pass 2)
A_reg <- regions |> filter(startsWith(set, "lost_")) |>
  mutate(class = class_levels[match(sub("^lost_", "", set), slug(class_levels))])
res_A <- fit_model("A_lost_classes", A_reg$feature_id, factor(A_reg$class, class_levels), "multinomial")

# Model B: lost vs kept within each pass-1 class
res_B <- lapply(class_levels, function(k) {
  s <- slug(k)
  r <- regions |> filter(set %in% c(paste0("p1lost_", s), paste0("p1kept_", s)))
  y <- factor(ifelse(startsWith(r$set, "p1lost_"), "Lost", "Kept"), c("Kept", "Lost"))
  out <- fit_model(paste0("B_", s), r$feature_id, y, "binomial")
  out$auc$target <- k; out$stab$target <- k
  out
})

auc_all  <- bind_rows(res_A$auc, bind_rows(lapply(res_B, `[[`, "auc")))
stab_all <- bind_rows(res_A$stab, bind_rows(lapply(res_B, `[[`, "stab")))
write_tsv(auc_all, file.path(outdir, "heldout_auc.tsv"))
write_tsv(stab_all, file.path(outdir, "stability_selection.tsv"))

# =============================================================================
# 4. HOMER known motifs on the same windows (cached; needs the FASTA)
# =============================================================================

homer_runs <- bind_rows(
  tibble(model = "A", class = class_levels, primary = paste0("lost_", slug(class_levels)), bg = "static"),
  tibble(model = "B_lost", class = class_levels, primary = paste0("p1lost_", slug(class_levels)),
         bg = paste0("p1kept_", slug(class_levels))),
  tibble(model = "B_kept", class = class_levels, primary = paste0("p1kept_", slug(class_levels)),
         bg = paste0("p1lost_", slug(class_levels))))
homer_ok <- nzchar(Sys.which("findMotifsGenome.pl"))
if (file.exists(fa_path) && homer_ok) {
  bed6 <- function(set) {
    f <- file.path(homer_dir, "beds", paste0(set, ".bed")); dir.create(dirname(f), FALSE, TRUE)
    if (!file.exists(f)) write_tsv(read_bed(set) |> transmute(chr, start, end, feature_id, 0, "+"),
                                   f, col_names = FALSE)
    f
  }
  for (i in seq_len(nrow(homer_runs))) {
    o <- file.path(homer_dir, paste0(homer_runs$model[i], "__", homer_runs$primary[i]))
    if (file.exists(file.path(o, "knownResults.txt"))) next
    message("  HOMER ", basename(o))
    system2("findMotifsGenome.pl", c(bed6(homer_runs$primary[i]), fa_path, o, "-size given", "-nomotif",
                                     "-p", CORES, "-bg", bed6(homer_runs$bg[i])),
            stdout = paste0(o, ".log"), stderr = paste0(o, ".log"))
  }
}
read_homer <- function(model, primary) {
  f <- file.path(homer_dir, paste0(model, "__", primary), "knownResults.txt")
  if (!file.exists(f)) return(NULL)
  h <- read_tsv(f, show_col_types = FALSE, col_types = cols(.default = "c"))
  names(h)[c(1, 5, 7, 9)] <- c("name", "q", "pct_t", "pct_b")
  h |> transmute(motif = toupper(sub("\\(.*$", "", name)),
                 homer_q = as.numeric(q),
                 homer_enr = as.numeric(sub("%", "", pct_t)) / pmax(as.numeric(sub("%", "", pct_b)), 0.01)) |>
    group_by(motif) |> slice_min(homer_q, n = 1, with_ties = FALSE) |> ungroup()
}
homer <- homer_runs |> mutate(res = purrr::map2(model, primary, read_homer)) |> tidyr::unnest(res)
if (!nrow(homer)) message("  No HOMER results (FASTA or HOMER unavailable): comparing with SEA only")

# =============================================================================
# 5. glmnet vs SEA vs HOMER
# =============================================================================

sea <- read_tsv(file.path(motif_dir, "sea_all.tsv.gz"), show_col_types = FALSE) |>
  group_by(comparison, class, motif) |> slice_min(qvalue, n = 1, with_ties = FALSE) |> ungroup()

# Model B: signed on the lost-vs-kept axis (+ = lost)
sea_B <- sea |> filter(comparison %in% c("lost_vs_kept", "kept_vs_lost")) |>
  mutate(side = if_else(comparison == "lost_vs_kept", "+", "-")) |>
  group_by(class, motif) |> slice_min(qvalue, n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(class, motif, sea_side = side, sea_q = qvalue, sea_enr = enr_ratio)
homer_B <- homer |> filter(model %in% c("B_lost", "B_kept")) |>
  mutate(side = if_else(model == "B_lost", "+", "-")) |>
  group_by(class, motif) |> slice_min(homer_q, n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(class, motif, homer_side = side, homer_q, homer_enr)
cmp_B <- stab_all |> filter(startsWith(model, "B_")) |>
  transmute(class = target, motif, glmnet_freq = sel_freq, glmnet_sign = sign, glmnet_stable = stable) |>
  full_join(sea_B, by = c("class", "motif")) |>
  full_join(homer_B, by = c("class", "motif"))

# Model A: glmnet (class vs other lost classes) next to SEA/HOMER (class vs static)
cmp_A <- stab_all |> filter(model == "A_lost_classes") |>
  transmute(class = target, motif, glmnet_freq = sel_freq, glmnet_sign = sign, glmnet_stable = stable) |>
  full_join(sea |> filter(comparison == "vs_static") |> transmute(class, motif, sea_q = qvalue, sea_enr = enr_ratio),
            by = c("class", "motif")) |>
  full_join(homer |> filter(model == "A") |> transmute(class, motif, homer_q, homer_enr),
            by = c("class", "motif"))

flag <- function(d) d |> mutate(
  sea_sig = coalesce(sea_q < Q_SIG, FALSE), homer_sig = coalesce(homer_q < Q_SIG, FALSE),
  sea_miss = coalesce(sea_q >= Q_MISS, TRUE), homer_miss = coalesce(homer_q >= Q_MISS, TRUE),
  agreement = case_when(
    coalesce(glmnet_stable, FALSE) & sea_miss & homer_miss ~ "glmnet only",
    coalesce(glmnet_stable, FALSE) & (sea_sig | homer_sig)  ~ "glmnet + enrichment",
    coalesce(glmnet_stable, FALSE)                         ~ "glmnet + weak enrichment",
    sea_sig & homer_sig                                     ~ "SEA + HOMER, dropped by glmnet",
    sea_sig                                                 ~ "SEA only",
    homer_sig                                               ~ "HOMER only",
    TRUE                                                    ~ "none"))
cmp_A <- flag(cmp_A); cmp_B <- flag(cmp_B)
write_tsv(cmp_A, file.path(outdir, "comparison_modelA_glmnet_sea_homer.tsv"))
write_tsv(cmp_B, file.path(outdir, "comparison_modelB_glmnet_sea_homer.tsv"))
print(as.data.frame(bind_rows(A = count(cmp_A, class, agreement), B = count(cmp_B, class, agreement), .id = "model") |>
                      filter(agreement != "none") |> tidyr::pivot_wider(names_from = agreement, values_from = n, values_fill = 0)))

# =============================================================================
# 6. Figures
# =============================================================================

# --- Fig 1: held-out AUC --------------------------------------------------------
d1 <- auc_all |>
  mutate(panel = if_else(startsWith(model, "A_"), "A: lost-site classes (one vs rest)",
                         "B: lost vs kept within WT class"),
         target = factor(target, class_levels)) |>
  tidyr::pivot_longer(c(auc_cov, auc_full), names_to = "fit", values_to = "auc") |>
  mutate(lo = if_else(fit == "auc_cov", auc_cov_lo, auc_full_lo),
         hi = if_else(fit == "auc_cov", auc_cov_hi, auc_full_hi),
         fit = recode(fit, auc_cov = "GC + CpG + promoter", auc_full = "+ motifs (elastic net)"))
p1 <- ggplot(d1, aes(auc, target, colour = fit)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.5), size = 0.3) +
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(values = c(`GC + CpG + promoter` = "grey55", `+ motifs (elastic net)` = viridisLite::mako(1, begin = 0.35))) +
  labs(x = "Held-out AUC (20% test split, bootstrap 95% CI)", y = NULL, colour = NULL,
       title = "Do motifs predict ARID1A-dependent sites beyond sequence composition?") +
  theme_bw(base_size = 9) + theme(legend.position = "top", strip.background = element_blank())
save_figure(p1, "fig1_heldout_auc", width = 8.5, height = 3.4, dir = figdir)

# --- Fig 2: stable motifs --------------------------------------------------------
d2 <- stab_all |> filter(stable) |>
  mutate(panel = if_else(startsWith(model, "A_"), "A", "B"),
         target = factor(target, class_levels)) |>
  group_by(panel, target) |> slice_max(sel_freq * 100 + abs(coalesce(median_coef, 0)), n = 12, with_ties = FALSE) |>
  ungroup() |>
  mutate(label = reorder(paste(motif, panel, target, sep = "___"), sel_freq))
lab2 <- c(A = "A: this lost-site class vs other lost classes", B = "B: lost (+) vs kept (-) in KO")
p2 <- ggplot(d2, aes(sel_freq, label, fill = sign)) +
  geom_col(width = 0.7) +
  facet_wrap(panel ~ target, scales = "free_y", nrow = 2,
             labeller = labeller(panel = lab2, .multi_line = FALSE)) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(values = c(`+` = "#2CA02C", `-` = "#000000"), name = "Coefficient") +
  coord_cartesian(xlim = c(STABLE, 1)) +
  labs(x = sprintf("Selection frequency (%d half-subsamples)", N_STAB), y = NULL,
       title = "Motifs the elastic net keeps once all motifs and GC/CpG/promoter are in the model") +
  theme_bw(base_size = 8) + theme(strip.background = element_blank(), legend.position = "top")
save_figure(p2, "fig2_stable_motifs", width = 14, height = 7, dir = figdir)

# --- Fig 3: agreement, model B -------------------------------------------------
d3 <- cmp_B |>
  mutate(class = factor(class, class_levels),
         x = if_else(coalesce(sea_side, "+") == "+", 1, -1) * log2(coalesce(sea_enr, 1)),
         y = if_else(coalesce(glmnet_sign, "+") == "+", 1, -1) * coalesce(glmnet_freq, 0),
         homer = case_when(homer_sig & homer_side == "+" ~ "HOMER: lost", homer_sig ~ "HOMER: kept",
                           TRUE ~ "HOMER: n.s."))
lab3 <- d3 |> filter(agreement %in% c("glmnet only", "SEA + HOMER, dropped by glmnet") | abs(y) >= 0.95) |>
  group_by(class) |> slice_max(abs(y) + abs(x), n = 6, with_ties = FALSE) |> ungroup()
p3 <- ggplot(d3, aes(x, y)) +
  geom_hline(yintercept = c(-STABLE, STABLE), linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(aes(colour = homer), size = 0.9, alpha = 0.7) +
  ggrepel::geom_text_repel(data = lab3, aes(label = motif), size = 2.1, max.overlaps = 30, seed = SEED) +
  facet_wrap(~ class, nrow = 1) +
  scale_colour_manual(values = c(`HOMER: lost` = "#2CA02C", `HOMER: kept` = "#000000", `HOMER: n.s.` = "grey80"),
                      name = NULL) +
  labs(x = "SEA log2 enrichment (+ lost, - kept)", y = "glmnet selection frequency (signed)",
       title = "Lost vs kept: joint model (glmnet) vs one-motif-at-a-time tests (SEA, HOMER)") +
  theme_bw(base_size = 8) + theme(legend.position = "top", strip.background = element_blank())
save_figure(p3, "fig3_glmnet_sea_homer_agreement", width = 13, height = 3.8, dir = figdir)

message("=== atac_trajectories/07_glmnet_motif_models.R complete ===")
