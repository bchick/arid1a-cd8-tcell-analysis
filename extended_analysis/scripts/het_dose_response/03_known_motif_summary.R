#!/usr/bin/env Rscript
# =============================================================================
# het_dose_response/03_known_motif_summary.R — HOMER known-motif summary across dose classes
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Reads HOMER knownResults.txt for every (subset x class) comparison run by
# run_homer_dose_response.sh (written by 02_feature_enrichment.R) and builds
# a heatmap of -log10(p) for the union of top motifs, highlighting the
# ETS/AP-1 (buffered) vs T-box/RUNX (haploinsufficient) split, plus a
# haploinsufficient-vs-buffered heatmap.
#
# Inputs:  results/extended_analysis/het_dose_response/homer_motifs/<comparison>/knownResults.txt
# Outputs: figures/extended_analysis/het_dose_response/
#            motif_summary_heatmap.{pdf,png}, motif_haplo_vs_buffered_heatmap.{pdf,png}
# Usage:   Rscript extended_analysis/scripts/het_dose_response/03_known_motif_summary.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

homer_root <- file.path(paths$ext_results, "het_dose_response/homer_motifs")
figdir     <- file.path(paths$ext_figures, "het_dose_response")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

comparisons <- list.dirs(homer_root, recursive = FALSE, full.names = FALSE)
message("Found ", length(comparisons), " HOMER runs")

read_known <- function(cmp) {
  f <- file.path(homer_root, cmp, "knownResults.txt")
  if (!file.exists(f)) return(NULL)
  tab <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(tab) <- make.names(colnames(tab))
  tab$motif <- sub("\\(.*$", "", tab$Motif.Name)
  tab$neglogP <- -tab$Log.P.value / log(10)  # natural log -> log10
  tab$comparison <- cmp
  tab
}

all_motifs <- bind_rows(lapply(comparisons, read_known))
stopifnot(nrow(all_motifs) > 0)

# Focus on *_vs_insensitive contrasts (dose class vs flat background)
keep_cmps <- grep("_vs_insensitive$", unique(all_motifs$comparison), value = TRUE)
mat_df <- all_motifs %>% filter(comparison %in% keep_cmps)

# Top 6 motifs per comparison -> union
top_motifs <- mat_df %>%
  group_by(comparison) %>%
  slice_max(order_by = neglogP, n = 6, with_ties = FALSE) %>%
  pull(motif) %>% unique()

message("Union of top motifs: ", length(top_motifs))

# Build matrix: rows = motifs, cols = comparisons
wide <- mat_df %>%
  filter(motif %in% top_motifs) %>%
  group_by(comparison, motif) %>%
  summarise(neglogP = max(neglogP), .groups = "drop") %>%
  pivot_wider(names_from = comparison, values_from = neglogP, values_fill = 0)

mat <- as.matrix(wide[, -1])
rownames(mat) <- wide$motif

# Cap at 250 for viz (some motifs are p~1e-500)
mat_cap <- pmin(mat, 250)

# Column order: subset in TE, EEC, MP x class in buffered, linear, haplo
col_order <- c(
  "TE_buffered_vs_insensitive", "TE_linear_vs_insensitive", "TE_haplo_vs_insensitive",
  "EEC_buffered_vs_insensitive", "EEC_linear_vs_insensitive", "EEC_haplo_vs_insensitive",
  "MP_buffered_vs_insensitive", "MP_linear_vs_insensitive", "MP_haplo_vs_insensitive"
)
col_order <- intersect(col_order, colnames(mat_cap))
mat_cap <- mat_cap[, col_order]

# Relabel columns compactly
col_short <- sub("_vs_insensitive$", "", col_order)
col_short <- sub("_", " ", col_short)

# Row order: cluster within motif family groups
# Manual motif family assignment for annotation
fam <- case_when(
  grepl("^ETS|^Fli|^ERG|^Etv|^ETV|^GABPA|^EWS|^Elk", rownames(mat_cap)) ~ "ETS",
  grepl("^Jun|^Fos|^AP-1|^BATF|^Atf|^Fra",          rownames(mat_cap)) ~ "AP-1/BATF",
  grepl("^RUNX",                                      rownames(mat_cap)) ~ "RUNX",
  grepl("^Tbx|^Tbr|^Tbet|^Eomes",                     rownames(mat_cap)) ~ "T-box",
  grepl("^NFkB|^Rel",                                 rownames(mat_cap)) ~ "NFkB",
  grepl("^KLF|^Klf|^Sp[0-9]",                         rownames(mat_cap)) ~ "KLF/SP",
  grepl("^LEF|^Tcf",                                  rownames(mat_cap)) ~ "TCF/LEF",
  grepl("^Egr",                                       rownames(mat_cap)) ~ "EGR",
  grepl("^Gata",                                      rownames(mat_cap)) ~ "GATA",
  TRUE                                                                   ~ "other"
)

# Order rows: by family, then by max column -log10P descending
row_ord <- order(factor(fam, levels = c("ETS","AP-1/BATF","RUNX","T-box","NFkB","KLF/SP","TCF/LEF","EGR","GATA","other")),
                 -apply(mat_cap, 1, max))
mat_cap <- mat_cap[row_ord, ]
fam     <- fam[row_ord]

row_ann <- rowAnnotation(
  family = fam,
  col = list(family = c(
    "ETS" = "#377EB8", "AP-1/BATF" = "#E41A1C", "RUNX" = "#984EA3",
    "T-box" = "#FF7F00", "NFkB" = "#A65628", "KLF/SP" = "#F781BF",
    "TCF/LEF" = "#4DAF4A", "EGR" = "#999999", "GATA" = "#FFFF33",
    "other" = "#CCCCCC"
  )),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm")
)

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

col_fun <- colorRamp2(c(0, 25, 100, 250),
                      c("white", "#FDE0DD", "#FA9FB5", "#49006A"))

ht <- Heatmap(
  mat_cap,
  name = "-log10(p)",
  col  = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_labels = col_short,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9),
  left_annotation = row_ann,
  top_annotation  = subset_ann,
  row_split = fam,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_title_rot = 0,
  column_split = factor(sub(" .*$", "", col_short), levels = c("TE","EEC","MP")),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  border = TRUE,
  rect_gp = gpar(col = "grey90", lwd = 0.3),
  heatmap_legend_param = list(
    title = "-log10(p)\n(capped@250)",
    at = c(0, 50, 100, 200, 250),
    legend_height = unit(3, "cm")
  )
)

pdf(file.path(figdir, "motif_summary_heatmap.pdf"),
    width = 9, height = 0.18 * nrow(mat_cap) + 2.5)
draw(ht, merge_legends = TRUE)
dev.off()

png(file.path(figdir, "motif_summary_heatmap.png"),
    width = 9, height = 0.18 * nrow(mat_cap) + 2.5,
    units = "in", res = 300)
draw(ht, merge_legends = TRUE)
dev.off()

message("Wrote motif_summary_heatmap.{pdf,png}")

# Also: haplo_vs_buffered direct contrast (what distinguishes dose-sensitive
# from dose-tolerant sites, controlling for ARID1A binding)
keep_hvb <- grep("_haplo_vs_buffered$", unique(all_motifs$comparison), value = TRUE)
hvb_df <- all_motifs %>% filter(comparison %in% keep_hvb)

top_hvb <- hvb_df %>%
  group_by(comparison) %>%
  slice_max(order_by = neglogP, n = 8, with_ties = FALSE) %>%
  pull(motif) %>% unique()

wide_hvb <- hvb_df %>%
  filter(motif %in% top_hvb) %>%
  group_by(comparison, motif) %>%
  summarise(neglogP = max(neglogP), .groups = "drop") %>%
  pivot_wider(names_from = comparison, values_from = neglogP, values_fill = 0)

mat_hvb <- as.matrix(wide_hvb[, -1])
rownames(mat_hvb) <- wide_hvb$motif
col_hvb <- intersect(c("TE_haplo_vs_buffered","EEC_haplo_vs_buffered","MP_haplo_vs_buffered"),
                     colnames(mat_hvb))
mat_hvb <- mat_hvb[, col_hvb]
colnames(mat_hvb) <- sub("_haplo_vs_buffered$", "", col_hvb)

fam_hvb <- case_when(
  grepl("^ETS|^Fli|^ERG|^Etv|^ETV|^GABPA|^EWS", rownames(mat_hvb)) ~ "ETS",
  grepl("^Jun|^Fos|^AP-1|^BATF|^Atf|^Fra",     rownames(mat_hvb)) ~ "AP-1/BATF",
  grepl("^RUNX",                                 rownames(mat_hvb)) ~ "RUNX",
  grepl("^Tbx|^Tbr|^Tbet|^Eomes",                rownames(mat_hvb)) ~ "T-box",
  grepl("^NFkB|^Rel",                            rownames(mat_hvb)) ~ "NFkB",
  grepl("^KLF|^Klf",                             rownames(mat_hvb)) ~ "KLF",
  grepl("^LEF|^Tcf",                             rownames(mat_hvb)) ~ "TCF/LEF",
  grepl("^Egr",                                  rownames(mat_hvb)) ~ "EGR",
  TRUE                                                              ~ "other"
)

row_ord_hvb <- order(factor(fam_hvb,
  levels = c("T-box","NFkB","AP-1/BATF","TCF/LEF","KLF","RUNX","EGR","ETS","other")),
  -apply(mat_hvb, 1, max))
mat_hvb <- mat_hvb[row_ord_hvb, ]
fam_hvb <- fam_hvb[row_ord_hvb]

ht2 <- Heatmap(
  mat_hvb,
  name = "-log10(p)",
  col  = colorRamp2(c(0, 10, 20, 35), c("white","#FFF7BC","#FE9929","#662506")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = fam_hvb,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 0,
  column_names_centered = TRUE,
  border = TRUE,
  rect_gp = gpar(col = "grey90", lwd = 0.4),
  heatmap_legend_param = list(title = "-log10(p)")
)

pdf(file.path(figdir, "motif_haplo_vs_buffered_heatmap.pdf"),
    width = 5.5, height = 0.2 * nrow(mat_hvb) + 2)
draw(ht2)
dev.off()
png(file.path(figdir, "motif_haplo_vs_buffered_heatmap.png"),
    width = 5.5, height = 0.2 * nrow(mat_hvb) + 2,
    units = "in", res = 300)
draw(ht2)
dev.off()

message("Wrote motif_haplo_vs_buffered_heatmap.{pdf,png}")
message("Done: ", Sys.time())
