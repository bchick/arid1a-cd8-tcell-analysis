#!/usr/bin/env Rscript
# =============================================================================
# motif_grammar/01_region_sets.R — confounder-matched region sets and motif instances
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# First step of the motif-grammar analysis. Builds CONFOUNDER-MATCHED contrast
# sets + fixed-window sequences + motif instance positions (JASPAR2020 CORE
# vertebrates, motifmatchr). Everything downstream (grammar features, glmnet,
# SpaMo) consumes this script's outputs.
#
# Inputs:  results/extended_analysis/het_dose_response/dose_classes_combined.csv
#          results/atac/differential/{consensus_peaks_annotated.csv,da_KO_vs_WT_D8_pseudobulk.csv}
#          data/reference/GRCm39.primary_assembly.genome.fa (genome$fasta; if absent,
#            the bundled grammar_*.rds are used and this script exits early)
# Outputs: results/extended_analysis/motif_grammar/
#            balance_{primary,control}.csv, grammar_{primary,control}.rds,
#            {primary,control}_<class>.fa (per-class FASTA for SpaMo)
# Usage:   Rscript extended_analysis/scripts/motif_grammar/01_region_sets.R   (from the repository root)
# =============================================================================

# Two contrasts:
#   PRIMARY  buffered vs haploinsufficient  — why does ONE Arid1a copy suffice at
#            some cBAF-responsive enhancers but not others? Sequence-encoded
#            dosage sensitivity. Uses the Het samples.
#   CONTROL  cBAF-dependent (lost) vs cBAF-independent (unchanged) — the paper's
#            Fig 5C-style split. Serves as a POSITIVE CONTROL: the pipeline must
#            recover the known ETS enrichment, or it is broken.
#
# WHY MATCHING IS NOT OPTIONAL (measured on this data):
#   buffered vs haploinsufficient : median baseMean 79.4 vs 30.4 (2.6x),
#                                   promoter 31.5% vs 20.2%
#   dependent vs independent      : median baseMean 46.1 vs 21.5 (2.1x),
#                                   promoter 16.1% vs 27.8%
#   Signal strength drives how many motif hits are detectable and promoters have
#   their own motif vocabulary, so an unmatched comparison "discovers" grammar
#   that is really read depth and genomic annotation. Peaks also vary 93-5753 bp,
#   which alone would confound motif counts. We therefore: keep DISTAL only, keep
#   LOST only, fix the window, and match on baseMean x GC deciles.
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Rsamtools)
  library(Biostrings)
  library(motifmatchr)
  library(TFBSTools)
  library(JASPAR2020)
})

# Bioconductor masks dplyr verbs (see core/01_rnaseq_analysis.R)
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

WINDOW   <- 300L   # fixed bp window centred on each peak; peaks are 93-5753 bp
CANON    <- paste0("chr", c(1:19, "X", "Y"))

message("=== 1. Load peak classes + annotation ===")

dose <- read_csv(file.path(paths$ext_results, "het_dose_response/dose_classes_combined.csv"),
          show_col_types = FALSE)
anno <- read_csv(file.path(paths$results,
          "atac/differential/consensus_peaks_annotated.csv"),
          show_col_types = FALSE)
da8  <- read_csv(file.path(paths$results,
          "atac/differential/da_KO_vs_WT_D8_pseudobulk.csv"),
          show_col_types = FALSE)

anno_slim <- anno %>%
  select(peak_id, seqnames, start, end, annotation) %>%
  mutate(region = case_when(
    grepl("^Promoter", annotation)            ~ "Promoter",
    grepl("Intron", annotation)               ~ "Intron",
    grepl("Intergenic|Distal", annotation)    ~ "Intergenic",
    grepl("Exon|UTR|Downstream", annotation)  ~ "Other",
    TRUE                                      ~ "Other"))

message("  dose rows: ", nrow(dose), " | annotated peaks: ", nrow(anno_slim))

# ---- PRIMARY: buffered vs haploinsufficient --------------------------------
# A peak is classed per SUBSET and may differ between subsets. 18,439/19,115
# (96.5%) are consistent; drop the 676 ambiguous rather than pick a winner.
message("=== 2. PRIMARY contrast: buffered vs haploinsufficient ===")

dose_bh <- dose %>% filter(class %in% c("buffered", "haploinsufficient"))

consistent <- dose_bh %>%
  group_by(peak_id) %>%
  summarize(n_class = n_distinct(class), n_lost = sum(direction == "lost"),
            n_sub = n(), class = dplyr::first(class),
            baseMean = mean(baseMean), .groups = "drop") %>%
  filter(n_class == 1)

message("  consistent-class peaks: ", nrow(consistent),
        " (dropped ", n_distinct(dose_bh$peak_id) - nrow(consistent), " ambiguous)")

primary <- consistent %>%
  # LOST only: 'buffered' mixes 15,106 lost with 7,667 gained; gained peaks are a
  # different phenomenon (likely library-size renormalisation drift after global
  # accessibility loss) and must not be pooled with loss.
  filter(n_lost == n_sub) %>%
  left_join(anno_slim, by = "peak_id") %>%
  filter(!is.na(seqnames), seqnames %in% CANON, region != "Promoter") %>%
  select(peak_id, class, baseMean, seqnames, start, end, region)

message("  after lost-only + distal + canonical: ",
        paste(capture.output(print(table(primary$class))), collapse = " "))

# ---- CONTROL: cBAF-dependent vs independent --------------------------------
message("=== 3. CONTROL contrast: cBAF-dependent vs independent ===")

control <- da8 %>%
  mutate(class = case_when(
    padj < 0.05 & log2FoldChange < -1                      ~ "cBAF_dependent",
    !is.na(padj) & padj > 0.5 & abs(log2FoldChange) < 0.25 ~ "cBAF_independent",
    TRUE                                                   ~ NA_character_)) %>%
  filter(!is.na(class)) %>%
  left_join(anno_slim, by = "peak_id") %>%
  filter(!is.na(seqnames), seqnames %in% CANON, region != "Promoter") %>%
  select(peak_id, class, baseMean, seqnames, start, end, region)

message("  after distal + canonical: ",
        paste(capture.output(print(table(control$class))), collapse = " "))

# ---- 4. fixed windows + GC -------------------------------------------------
message("=== 4. Fixed ", WINDOW, " bp windows + GC content ===")

if (!file.exists(genome$fasta)) {
  have <- file.exists(file.path(outdir, c("grammar_primary.rds", "grammar_control.rds")))
  if (all(have)) {
    message("Skipping: genome FASTA not found; using the matched sets from the data bundle")
    quit(save = "no", status = 0)
  }
  stop("Need ", genome$fasta, " (upstream mode) or the bundled grammar_*.rds")
}
fa <- FaFile(genome$fasta)
fa_idx <- scanFaIndex(fa)
chr_len <- setNames(width(fa_idx), as.character(seqnames(fa_idx)))

make_windows <- function(df) {
  ctr <- as.integer((df$start + df$end) / 2)
  half <- WINDOW %/% 2L
  s <- ctr - half; e <- ctr + half - 1L
  keep <- s >= 1 & e <= chr_len[df$seqnames]   # drop windows running off the end
  df <- df[keep, ]; s <- s[keep]; e <- e[keep]
  gr <- GRanges(df$seqnames, IRanges(s, e))
  mcols(gr)$peak_id <- df$peak_id
  mcols(gr)$class   <- df$class
  mcols(gr)$baseMean <- df$baseMean
  gr
}

add_gc <- function(gr) {
  seqs <- getSeq(fa, gr)
  af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
  mcols(gr)$gc <- af[, "G"] + af[, "C"]
  mcols(gr)$n_frac <- 1 - rowSums(af[, c("A", "C", "G", "T")])
  gr
}

gr_primary <- add_gc(make_windows(primary))
gr_control <- add_gc(make_windows(control))

# windows that are mostly N carry no sequence information
gr_primary <- gr_primary[mcols(gr_primary)$n_frac < 0.1]
gr_control <- gr_control[mcols(gr_control)$n_frac < 0.1]

message("  primary windows: ", length(gr_primary),
        " | control windows: ", length(gr_control))

# ---- 5. confounder matching ------------------------------------------------
# Match the minority class 1:1 to the majority on baseMean x GC decile strata.
# Without this, "grammar" differences are really depth and base composition.
message("=== 5. Matching on baseMean x GC deciles ===")

match_sets <- function(gr, minority_class) {
  d <- as.data.frame(mcols(gr))
  d$idx <- seq_len(nrow(d))
  d$bm_dec <- cut(rank(d$baseMean, ties.method = "first"),
                  breaks = 10, labels = FALSE)
  d$gc_dec <- cut(rank(d$gc, ties.method = "first"),
                  breaks = 10, labels = FALSE)
  d$stratum <- paste(d$bm_dec, d$gc_dec, sep = "_")

  min_d <- d[d$class == minority_class, ]
  maj_d <- d[d$class != minority_class, ]

  keep_maj <- unlist(lapply(split(min_d, min_d$stratum), function(chunk) {
    pool <- maj_d$idx[maj_d$stratum == chunk$stratum[1]]
    if (!length(pool)) return(integer(0))
    sample(pool, min(nrow(chunk), length(pool)))
  }))
  # keep only minority rows whose stratum actually found a partner
  matched_strata <- table(d$stratum[d$idx %in% keep_maj])
  keep_min <- unlist(lapply(split(min_d, min_d$stratum), function(chunk) {
    n <- matched_strata[chunk$stratum[1]]
    if (is.na(n) || n == 0) return(integer(0))
    head(chunk$idx, n)
  }))
  gr[sort(c(keep_min, keep_maj))]
}

gr_primary_m <- match_sets(gr_primary, "haploinsufficient")
gr_control_m <- match_sets(gr_control, "cBAF_independent")

report_balance <- function(gr, label) {
  d <- as.data.frame(mcols(gr))
  message("  ", label, ":")
  b <- d %>% group_by(class) %>%
    summarize(n = n(), median_baseMean = round(median(baseMean), 1),
              median_gc = round(median(gc), 3), .groups = "drop")
  for (i in seq_len(nrow(b))) {
    message(sprintf("    %-18s n=%5d  baseMean=%7.1f  GC=%.3f",
                    b$class[i], b$n[i], b$median_baseMean[i], b$median_gc[i]))
  }
  b
}

bal_primary <- report_balance(gr_primary_m, "PRIMARY (matched)")
bal_control <- report_balance(gr_control_m, "CONTROL (matched)")

write_csv(bal_primary, file.path(outdir, "balance_primary.csv"))
write_csv(bal_control, file.path(outdir, "balance_control.csv"))

# ---- 6. motif instances (positions + strand) -------------------------------
# Grammar needs POSITION and STRAND, not the binary hit matrix chromVAR uses.
message("=== 6. Scanning motif instances (positions + strand) ===")

pfm <- getMatrixSet(JASPAR2020, list(collection = "CORE",
                                     tax_group = "vertebrates",
                                     matrixtype = "PWM"))
tf_name <- sapply(pfm, name)

# Reuse the project's family definition verbatim (build_chromvar.R) so
# grammar results speak the same vocabulary as the chromVAR/TOBIAS analyses.
assign_family <- function(tf) {
  u <- toupper(tf)
  dplyr::case_when(
    grepl("ETS|ETV|^ERG|FLI1|GABP|ELK|^ELF|ELF[0-9]|SPI1|SPIB|EHF|FEV", u) ~ "ETS",
    grepl("RUNX", u) ~ "RUNX",
    grepl("TBX|EOMES|^TBR|MGA", u) ~ "T-box",
    grepl("^JUN|^FOS|BATF|^JDP|^ATF3|NFE2|^BACH", u) ~ "AP-1/bZIP",
    grepl("NFKB|^REL$|^RELA|^RELB|^REL[A-B]", u) ~ "NFkB",
    grepl("TCF7|LEF1|TCF7L", u) ~ "TCF/LEF",
    grepl("^KLF|^SP[0-9]", u) ~ "KLF/SP",
    grepl("GATA", u) ~ "GATA",
    grepl("IRF|STAT", u) ~ "IRF/STAT",
    grepl("^EGR", u) ~ "EGR",
    TRUE ~ "other")
}

fam <- assign_family(tf_name)
FAM_KEEP <- c("ETS", "RUNX", "T-box", "AP-1/bZIP", "NFkB", "TCF/LEF", "KLF/SP")
keep_motif <- fam %in% FAM_KEEP
pfm_keep <- pfm[keep_motif]
message("  motifs kept: ", length(pfm_keep), " / ", length(pfm),
        " across families: ", paste(FAM_KEEP, collapse = ", "))

scan_positions <- function(gr, label) {
  message("  scanning ", label, " (", length(gr), " windows)…")
  seqs <- getSeq(fa, gr)
  names(seqs) <- mcols(gr)$peak_id
  mp <- matchMotifs(pfm_keep, seqs, out = "positions", p.cutoff = 1e-4)
  # Structure (verified empirically, not assumed): with a DNAStringSet subject
  # matchMotifs returns a plain list of length n_motifs; each element is a
  # CompressedIRangesList of length n_sequences (one IRanges per window). There
  # are no seqnames — ranges are window-relative — and strand/score live in
  # mcols(), with NO strand() accessor on IRanges.
  motif_names <- tf_name[keep_motif]
  motif_fams  <- fam[keep_motif]
  res <- lapply(seq_along(mp), function(i) {
    g <- mp[[i]]
    nh <- elementNROWS(g)                  # hits per window
    if (!sum(nh)) return(NULL)
    ur <- unlist(g, use.names = FALSE)     # IRanges; mcols: strand, score
    tibble(
      window_i = rep(seq_along(g), nh),    # which window each hit belongs to
      motif    = motif_names[i],
      family   = motif_fams[i],
      pos      = (start(ur) + end(ur)) / 2,  # centre within the window
      strand   = as.character(mcols(ur)$strand),
      score    = mcols(ur)$score
    )
  })
  bind_rows(res) %>%
    mutate(peak_id = mcols(gr)$peak_id[window_i],
           class   = mcols(gr)$class[window_i])
}

pos_primary <- scan_positions(gr_primary_m, "PRIMARY")
pos_control <- scan_positions(gr_control_m, "CONTROL")

message("  primary motif instances: ", nrow(pos_primary))
message("  control motif instances: ", nrow(pos_control))

# ---- 7. save ---------------------------------------------------------------
saveRDS(list(gr = gr_primary_m, pos = pos_primary, window = WINDOW),
        file.path(outdir, "grammar_primary.rds"))
saveRDS(list(gr = gr_control_m, pos = pos_control, window = WINDOW),
        file.path(outdir, "grammar_control.rds"))

# FASTA per class for SpaMo / MEME-suite consumption
write_class_fasta <- function(gr, tag) {
  for (cl in unique(mcols(gr)$class)) {
    sub <- gr[mcols(gr)$class == cl]
    seqs <- getSeq(fa, sub)
    names(seqs) <- mcols(sub)$peak_id
    f <- file.path(outdir, sprintf("%s_%s.fa", tag, gsub("[^A-Za-z0-9]", "", cl)))
    writeXStringSet(seqs, f)
    message("  wrote ", f, " (", length(seqs), " seqs)")
  }
}
write_class_fasta(gr_primary_m, "primary")
write_class_fasta(gr_control_m, "control")

message("=== motif_grammar/01_region_sets.R complete. Outputs in ", outdir, " ===")
