#!/usr/bin/env Rscript
# =============================================================================
# het_dose_response/04_denovo_motif_summary.R — HOMER de novo motif summary across dose classes
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Reads homerResults/motif*.motif files from each HOMER run (driver written by
# 02_feature_enrichment.R). Collapses de novo motifs to TF families via the
# BestGuess tag, and builds:
#   1. family x comparison heatmap, class vs insensitive
#   2. the same for the direct haploinsufficient-vs-buffered contrast
#   3. a table of the top 8 de novo motifs per run with consensus, best-guess,
#      -log10 p, %target, %bg, match score
#
# Inputs:  results/extended_analysis/het_dose_response/homer_motifs/<comparison>/homerResults/
# Outputs: results/extended_analysis/het_dose_response/motif_denovo_top_per_comparison.csv
#          figures/extended_analysis/het_dose_response/
#            motif_denovo_summary_heatmap.{pdf,png}, motif_denovo_haplo_vs_buffered_heatmap.{pdf,png}
# Usage:   Rscript extended_analysis/scripts/het_dose_response/04_denovo_motif_summary.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

homer_root <- file.path(paths$ext_results, "het_dose_response/homer_motifs")
figdir     <- file.path(paths$ext_figures, "het_dose_response")
outdir     <- file.path(paths$ext_results, "het_dose_response")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

comparisons <- list.dirs(homer_root, recursive = FALSE, full.names = FALSE)
message("Found ", length(comparisons), " HOMER runs")

# -----------------------------------------------------------------------------
# 1. Parse all de novo motifs
# -----------------------------------------------------------------------------
# Header line of each motif file:
# >CONSENSUS<TAB>RANK-CONSENSUS,BestGuess:TFNAME(FAMILY)/...(MATCHSCORE)<TAB>logodds<TAB>lnP<TAB>0<TAB>T:n(%),B:n(%),P:1eX
parse_motif_header <- function(path) {
  con <- file(path, "r"); on.exit(close(con))
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0 || !startsWith(line, ">")) return(NULL)
  fields <- strsplit(sub("^>", "", line), "\t", fixed = TRUE)[[1]]
  if (length(fields) < 6) return(NULL)
  consensus <- fields[1]
  name_full <- fields[2]
  lnp       <- suppressWarnings(as.numeric(fields[4]))
  occ       <- fields[6]

  # BestGuess:TFNAME(FAMILY)/.../Homer(0.98)  or  BestGuess:TFNAME/MA1234.2/Jaspar(0.74)
  bg_str <- sub(".*BestGuess:", "", name_full)
  # TF name is everything before the first "(" or "/"
  tf_name <- sub("[/(].*$", "", bg_str)
  # Family tag: text inside first (...) immediately after TF name, if present
  fam_tag <- NA_character_
  m <- regmatches(bg_str, regexpr("^[^/(]+\\(([^)]+)\\)", bg_str))
  if (length(m) && nzchar(m)) {
    fam_tag <- sub("^[^(]+\\(([^)]+)\\).*$", "\\1", m)
  }
  # Match score: last (0.xx) before end of string
  ms <- regmatches(bg_str, regexpr("\\(([0-9.]+)\\)\\s*$", bg_str))
  match_score <- if (length(ms) && nzchar(ms)) {
    as.numeric(sub("[()]", "", sub("\\)\\s*$", "", ms)))
  } else NA_real_

  # Parse T:n(%),B:n(%)
  pct_t <- suppressWarnings(as.numeric(sub("%.*$", "",
                                           sub(".*T:[^,]*\\(", "", occ))))
  pct_b <- suppressWarnings(as.numeric(sub("%.*$", "",
                                           sub(".*B:[^,]*\\(", "", occ))))

  data.frame(
    consensus  = consensus,
    tf_name    = tf_name,
    fam_tag    = fam_tag,
    match      = match_score,
    lnp        = lnp,
    neglogP    = ifelse(is.na(lnp), NA_real_, -lnp / log(10)),
    pct_target = pct_t,
    pct_bg     = pct_b,
    stringsAsFactors = FALSE
  )
}

read_denovo <- function(cmp) {
  hr <- file.path(homer_root, cmp, "homerResults")
  if (!dir.exists(hr)) return(NULL)
  mfiles <- list.files(hr, pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
  if (length(mfiles) == 0) return(NULL)
  # rank = integer in filename
  rank <- as.integer(sub(".*motif([0-9]+)\\.motif$", "\\1", mfiles))
  rows <- lapply(mfiles, parse_motif_header)
  ok   <- !vapply(rows, is.null, logical(1))
  rows <- rows[ok]; rank <- rank[ok]
  out  <- do.call(rbind, rows)
  out$rank       <- rank
  out$comparison <- cmp
  out
}

all_denovo <- bind_rows(lapply(comparisons, read_denovo))
stopifnot(nrow(all_denovo) > 0)
message("Parsed ", nrow(all_denovo), " de novo motifs across ",
        length(unique(all_denovo$comparison)), " comparisons")

# -----------------------------------------------------------------------------
# 2. Assign TF family
#    Primary source: HOMER fam_tag (e.g. "ETS", "Runt", "T-box", "bZIP", "Zf")
#    Secondary: regex on tf_name when tag is missing/ambiguous
# -----------------------------------------------------------------------------
classify_family <- function(tf_name, fam_tag) {
  fam <- rep("other", length(tf_name))
  # Name-first: some family tags are too broad (e.g. "bZIP" includes AP-1 and ATF/CREB)
  # We override with specific TF-name matches where the biology is clear.
  fam[grepl("^(ETS|Fli|ERG|Etv|ETV|GABPA|EWS|Elk|Ets)", tf_name, ignore.case = TRUE)] <- "ETS"
  fam[grepl("^(RUNX|Runx)",                             tf_name, ignore.case = TRUE)] <- "RUNX"
  fam[grepl("^(Tbx|Tbr|Tbet|Eomes|T-box)",              tf_name, ignore.case = TRUE)] <- "T-box"
  fam[grepl("^(Jun|Fos|AP-?1|BATF|Atf|Fra|Bach|Nfe2|Maf)", tf_name, ignore.case = TRUE)] <- "AP-1/bZIP"
  fam[grepl("^(NFkB|NF-kB|Rel|p65|p52|p50)",            tf_name, ignore.case = TRUE)] <- "NFkB"
  fam[grepl("^(KLF|Klf|Sp[0-9])",                       tf_name, ignore.case = TRUE)] <- "KLF/SP"
  fam[grepl("^(TCF|Tcf|LEF|Lef)",                       tf_name, ignore.case = TRUE)] <- "TCF/LEF"
  fam[grepl("^(Egr)",                                   tf_name, ignore.case = TRUE)] <- "EGR"
  fam[grepl("^(Gata)",                                  tf_name, ignore.case = TRUE)] <- "GATA"
  fam[grepl("^(IRF|Stat)",                              tf_name, ignore.case = TRUE)] <- "IRF/STAT"
  # Fallback to HOMER's family tag if still "other"
  use_tag <- fam == "other" & !is.na(fam_tag)
  tag_map <- c("ETS"="ETS","Runt"="RUNX","T-box"="T-box","bZIP"="AP-1/bZIP",
               "RHD"="NFkB","Zf,KLF"="KLF/SP","HMG"="TCF/LEF",
               "Zf"="Zf","Homeobox"="Homeo","NR"="NR","bHLH"="bHLH")
  fam[use_tag] <- ifelse(fam_tag[use_tag] %in% names(tag_map),
                         tag_map[fam_tag[use_tag]], "other")
  fam
}

all_denovo$family <- classify_family(all_denovo$tf_name, all_denovo$fam_tag)

# Quality filter: require best-guess match score >= 0.6 to trust the TF label;
# below that, keep the motif but mark family as "other" (unclassified consensus).
poor_match <- is.na(all_denovo$match) | all_denovo$match < 0.6
all_denovo$family[poor_match] <- "other"

message("Family counts across all runs:")
print(table(all_denovo$family))

# -----------------------------------------------------------------------------
# 3. Top-N de novo motifs per comparison (for CSV)
# -----------------------------------------------------------------------------
top_table <- all_denovo %>%
  group_by(comparison) %>%
  slice_max(order_by = neglogP, n = 8, with_ties = FALSE) %>%
  arrange(comparison, desc(neglogP)) %>%
  select(comparison, rank, consensus, tf_name, family, fam_tag, match,
         neglogP, pct_target, pct_bg) %>%
  ungroup()

write_csv(top_table, file.path(outdir, "motif_denovo_top_per_comparison.csv"))
message("Wrote motif_denovo_top_per_comparison.csv (", nrow(top_table), " rows)")

# -----------------------------------------------------------------------------
# 4. Family-level heatmap — class vs insensitive
# -----------------------------------------------------------------------------
build_family_mat <- function(df, keep_cmps) {
  m <- df %>%
    filter(comparison %in% keep_cmps) %>%
    group_by(comparison, family) %>%
    summarise(neglogP = max(neglogP, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = comparison, values_from = neglogP, values_fill = 0)
  mat <- as.matrix(m[, -1])
  rownames(mat) <- m$family
  mat
}

vs_ins <- grep("_vs_insensitive$", unique(all_denovo$comparison), value = TRUE)
mat_ins <- build_family_mat(all_denovo, vs_ins)

# Drop "other" from the heatmap: it's a mix of unclassified motifs, not a
# single biological family. Keep it in the CSV top-table.
mat_ins <- mat_ins[rownames(mat_ins) != "other", , drop = FALSE]

# Column order: subset in TE, EEC, MP x class in buffered, linear, haplo
col_order <- c(
  "TE_buffered_vs_insensitive", "TE_linear_vs_insensitive", "TE_haplo_vs_insensitive",
  "EEC_buffered_vs_insensitive", "EEC_linear_vs_insensitive", "EEC_haplo_vs_insensitive",
  "MP_buffered_vs_insensitive", "MP_linear_vs_insensitive", "MP_haplo_vs_insensitive"
)
col_order <- intersect(col_order, colnames(mat_ins))
mat_ins   <- mat_ins[, col_order]
col_short <- sub("_vs_insensitive$", "", col_order)
col_short <- sub("_", " ", col_short)

# Row order: by max enrichment across classes, descending
mat_ins <- mat_ins[order(-apply(mat_ins, 1, max)), , drop = FALSE]

# -----------------------------------------------------------------------------
# Scale: buffered ETS reaches ~250, haplo RUNX/Tbox ~30-60, linear <15.
# Use a 2-step colour ramp that gives haplo/linear visible contrast without
# letting buffered ETS dominate the entire palette.
# -----------------------------------------------------------------------------
cap_val <- 100
mat_cap <- pmin(mat_ins, cap_val)

subset_ann <- columnAnnotation(
  subset = sub(" .*$", "", col_short),
  class  = sub("^[^ ]+ ", "", col_short),
  col = list(
    subset = c("TE" = "#CC6677", "EEC" = "#DDCC77", "MP" = "#AA4499"),
    class  = c("buffered" = "#4A90D9", "linear" = "#2CA02C",
               "haplo"    = "#CC3311")
  ),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm")
)

col_fun <- colorRamp2(c(0, 5, 15, 40, 100),
                      c("white", "#FDE0DD", "#FA9FB5", "#C51B8A", "#49006A"))

ht <- Heatmap(
  mat_cap,
  name = "-log10(p)",
  col  = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_labels = col_short,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  top_annotation  = subset_ann,
  column_split = factor(sub(" .*$", "", col_short), levels = c("TE","EEC","MP")),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  border = TRUE,
  rect_gp = gpar(col = "grey92", lwd = 0.3),
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat_ins[i, j]
    if (!is.na(v) && v >= 5) {
      grid.text(sprintf("%.0f", v), x, y,
                gp = gpar(fontsize = 8,
                          col = ifelse(v >= 60, "white", "grey10")))
    }
  },
  heatmap_legend_param = list(
    title = "-log10(p)\n(capped@100)",
    at = c(0, 5, 15, 40, 100),
    legend_height = unit(3, "cm")
  )
)

pdf(file.path(figdir, "motif_denovo_summary_heatmap.pdf"),
    width = 9, height = 0.35 * nrow(mat_ins) + 3)
draw(ht, merge_legends = TRUE)
dev.off()
png(file.path(figdir, "motif_denovo_summary_heatmap.png"),
    width = 9, height = 0.35 * nrow(mat_ins) + 3,
    units = "in", res = 300)
draw(ht, merge_legends = TRUE)
dev.off()
message("Wrote motif_denovo_summary_heatmap.{pdf,png}")

# -----------------------------------------------------------------------------
# 5. Family-level heatmap — haplo vs buffered (direct contrast)
# -----------------------------------------------------------------------------
hvb <- grep("_haplo_vs_buffered$", unique(all_denovo$comparison), value = TRUE)
mat_hvb <- build_family_mat(all_denovo, hvb)
mat_hvb <- mat_hvb[rownames(mat_hvb) != "other", , drop = FALSE]

col_hvb <- intersect(c("TE_haplo_vs_buffered","EEC_haplo_vs_buffered","MP_haplo_vs_buffered"),
                     colnames(mat_hvb))
mat_hvb <- mat_hvb[, col_hvb]
colnames(mat_hvb) <- sub("_haplo_vs_buffered$", "", col_hvb)

mat_hvb <- mat_hvb[order(-apply(mat_hvb, 1, max)), , drop = FALSE]

col_fun2 <- colorRamp2(c(0, 3, 10, 25),
                       c("white", "#FFF7BC", "#FE9929", "#662506"))

ht2 <- Heatmap(
  mat_hvb,
  name = "-log10(p)",
  col  = col_fun2,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 11),
  column_names_rot = 0,
  column_names_centered = TRUE,
  border = TRUE,
  rect_gp = gpar(col = "grey92", lwd = 0.4),
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat_hvb[i, j]
    if (!is.na(v) && v >= 3) {
      grid.text(sprintf("%.0f", v), x, y,
                gp = gpar(fontsize = 9,
                          col = ifelse(v >= 15, "white", "grey10")))
    }
  },
  heatmap_legend_param = list(title = "-log10(p)", legend_height = unit(2.5, "cm"))
)

pdf(file.path(figdir, "motif_denovo_haplo_vs_buffered_heatmap.pdf"),
    width = 5, height = 0.35 * nrow(mat_hvb) + 2.2)
draw(ht2)
dev.off()
png(file.path(figdir, "motif_denovo_haplo_vs_buffered_heatmap.png"),
    width = 5, height = 0.35 * nrow(mat_hvb) + 2.2,
    units = "in", res = 300)
draw(ht2)
dev.off()
message("Wrote motif_denovo_haplo_vs_buffered_heatmap.{pdf,png}")

message("Done: ", Sys.time())
