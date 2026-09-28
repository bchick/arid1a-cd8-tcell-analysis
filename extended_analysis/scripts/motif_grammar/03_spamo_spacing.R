#!/usr/bin/env Rscript
# =============================================================================
# motif_grammar/03_spamo_spacing.R — unbiased motif-spacing discovery with MEME-suite SpaMo
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# 01_region_sets.R and 02_grammar_features.R test grammar features WE specified (close-pair fraction, strand
# configuration, helical phasing). SpaMo asks the same question without a
# hypothesis: given a PRIMARY motif, at which exact spacings and orientations is
# a SECONDARY motif enriched, relative to a background of the same sequences?
#
# Run per CLASS, then compare. A spacing that is significant in buffered but not
# in haploinsufficient (or vice versa) is a candidate grammar rule for dosage
# sensitivity.
#
# WHY THIS COMPLEMENTS 02_grammar_features.R: SpaMo is per-sequence-set and reports exact bp
# offsets with strand configuration (same/opposite, upstream/downstream), which
# our binned fraction features deliberately smooth over. The grammar models say whether
# syntax carries information; SpaMo says which syntax.
#
# NOTE: the class FASTAs come from 01_region_sets.R and are already matched on baseMean
# and GC, distal-only and lost-only. SpaMo's own background is derived from the
# same sequences, so the comparison is between like and like.
#
# Requires: MEME suite `spamo` on PATH.
#
# Inputs:  results/extended_analysis/motif_grammar/{primary,control}_<class>.fa (01_region_sets.R)
#          JASPAR2020 CORE vertebrates (Bioconductor package; exported to MEME format)
# Outputs: results/extended_analysis/motif_grammar/
#            motifs/jaspar2020_core_vertebrates.meme, spamo/<class>__<primary>/,
#            spamo_run_manifest.csv, spamo_all_results.csv
# Usage:   Rscript extended_analysis/scripts/motif_grammar/03_spamo_spacing.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(JASPAR2020)
  library(TFBSTools)
  library(universalmotif)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

outdir  <- file.path(paths$ext_results, "motif_grammar")
spamodir <- file.path(outdir, "spamo")
memedir <- file.path(outdir, "motifs")
dir.create(spamodir, recursive = TRUE, showWarnings = FALSE)
dir.create(memedir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. export JASPAR -> MEME ----------------------------------------------
# The existing export (results/extended_analysis/footprinting/motifs/*.jaspar, made by
# extended_analysis/scripts/footprinting/export_jaspar_motifs.R for TOBIAS) is raw-count JASPAR format;
# SpaMo needs MEME. Same motif set either way (JASPAR2020 CORE vertebrates), so
# grammar, chromVAR (dose_chromvar) and TOBIAS (tf_tobias_footprints) all speak the same motif vocabulary.
meme_file <- file.path(memedir, "jaspar2020_core_vertebrates.meme")

pfm <- getMatrixSet(JASPAR2020, list(collection = "CORE",
                                     tax_group = "vertebrates",
                                     matrixtype = "PFM"))
tf_name <- sapply(pfm, name)
tf_id   <- sapply(pfm, ID)

if (!file.exists(meme_file)) {
  message("=== Exporting ", length(pfm), " JASPAR motifs to MEME ===")
  um <- convert_motifs(pfm)
  write_meme(um, meme_file, overwrite = TRUE)
  message("  wrote ", meme_file)
} else {
  message("=== MEME motif file present: ", meme_file, " ===")
}

# ---- 2. pick primaries ------------------------------------------------------
# Reuse the project's family assignment (build_chromvar.R).
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

# One canonical primary per family of interest. ETS first: it is the dose-buffered
# backbone (82% footprint retention at one copy in the TOBIAS dose-response
# footprinting, tf_tobias_dose_response), so "what sits at
# a fixed spacing from ETS" is the mechanistic question.
pick <- function(pattern) {
  i <- grep(pattern, toupper(tf_name))
  if (!length(i)) return(NULL)
  list(id = tf_id[i[1]], name = tf_name[i[1]])
}
primaries <- Filter(Negate(is.null), list(
  ETS   = pick("^ETS1$"),
  RUNX  = pick("^RUNX1$"),
  TBX   = pick("^TBX21$"),
  AP1   = pick("^FOS::JUN$|^BATF$|^JUN$")
))
message("=== Primary motifs ===")
for (nm in names(primaries))
  message(sprintf("  %-5s %s (%s)", nm, primaries[[nm]]$name, primaries[[nm]]$id))

# Fail before launching 16 SpaMo runs if the identifiers we will request are not
# the ones the MEME file actually declares.
meme_ids <- sub("^MOTIF\\s+(\\S+).*$", "\\1",
                grep("^MOTIF ", readLines(meme_file), value = TRUE))
missing <- setdiff(vapply(primaries, function(p) p$name, character(1)), meme_ids)
if (length(missing))
  stop("primary motif(s) absent from ", basename(meme_file), ": ",
       paste(missing, collapse = ", "),
       "\n  MEME declares identifiers like: ", paste(utils::head(meme_ids, 3), collapse = ", "))

# ---- 3. run SpaMo per class -------------------------------------------------
fastas <- list.files(outdir, pattern = "^(primary|control)_.*\\.fa$", full.names = TRUE)
if (!length(fastas)) stop("No class FASTAs found — run extended_analysis/scripts/motif_grammar/01_region_sets.R first")
message("=== Class FASTAs ===")
for (f in fastas) message("  ", basename(f))

run_spamo <- function(fa, prim, tag) {
  od <- file.path(spamodir, sprintf("%s__%s", sub("\\.fa$", "", basename(fa)), tag))
  if (dir.exists(file.path(od)) && file.exists(file.path(od, "spamo.tsv"))) {
    message("  [skip] ", basename(od)); return(od)
  }
  # -primary matches the MEME MOTIF identifier, which write_meme() fills from the
  # universalmotif *name* (TF symbol) and not the JASPAR accession — the accession
  # lands in the alt-name field. Requesting prim$id silently FATALs every run.
  # Safe because JASPAR CORE vertebrates names are unique (746/746).
  cmd <- sprintf(
    "spamo -oc %s -primary %s -margin 150 -range 150 -bin 5 %s %s %s 2>&1",
    shQuote(od), shQuote(prim$name), shQuote(fa), shQuote(meme_file), shQuote(meme_file))
  message("  spamo: ", basename(fa), " primary=", prim$name)
  out <- system(cmd, intern = TRUE)
  if (!file.exists(file.path(od, "spamo.tsv")))
    message("    WARN no spamo.tsv — tail: ", paste(tail(out, 2), collapse = " | "))
  od
}

# A run that executed and found nothing is a RESULT (no spacing preference), not
# an error. A run that never produced spamo.tsv is a genuine failure. Keeping
# zero-row frames preserves the column spec so the two stay distinguishable.
results <- list(); runs <- list()
for (fa in fastas) {
  for (nm in names(primaries)) {
    od  <- run_spamo(fa, primaries[[nm]], nm)
    tsv <- file.path(od, "spamo.tsv")
    ran <- file.exists(tsv)
    n_sig <- NA_integer_
    if (ran) {
      x <- try(read_tsv(tsv, show_col_types = FALSE, comment = "#"), silent = TRUE)
      if (!inherits(x, "try-error")) {
        x$class_set <- sub("\\.fa$", "", basename(fa))
        x$primary_family <- nm
        results[[length(results) + 1]] <- x
        n_sig <- nrow(x)
      }
    }
    runs[[length(runs) + 1]] <- tibble(
      class_set = sub("\\.fa$", "", basename(fa)), primary_family = nm,
      primary_motif = primaries[[nm]]$name, ran = ran, n_significant = n_sig)
  }
}

manifest <- bind_rows(runs)
write_csv(manifest, file.path(outdir, "spamo_run_manifest.csv"))

bad <- manifest %>% dplyr::filter(!ran | is.na(n_significant))
if (nrow(bad))
  stop(nrow(bad), "/", nrow(manifest), " SpaMo runs produced no parseable spamo.tsv — ",
       "check ", spamodir, "\n  first failure: ",
       bad$class_set[1], " / ", bad$primary_motif[1])

all <- bind_rows(results)
write_csv(all, file.path(outdir, "spamo_all_results.csv"))
message(sprintf("=== SpaMo complete: %d/%d runs executed, %d significant spacings ===",
                sum(manifest$ran), nrow(manifest), nrow(all)))
if (!nrow(all)) {
  message("=== NULL RESULT: every run executed cleanly and found no secondary motif ",
          "enriched at any spacing. Consistent with the grammar-model dAUC ~ 0. ===")
} else {
  print(as.data.frame(manifest))
  print(utils::head(as.data.frame(all), 5))
}
