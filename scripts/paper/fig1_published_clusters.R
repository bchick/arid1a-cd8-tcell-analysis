#!/usr/bin/env Rscript
# =============================================================================
# paper/fig1_published_clusters.R — canonical Fig 1A OCR clusters (published regions)
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# The Fig 1A clusters are the published cluster definitions from the original
# mm10 HOMER analysis, lifted to mm39 (see data/metadata/paper_ocr_clusters/
# README.md). This script:
#   1. writes results/paper/fig1a_ocr_clusters.csv — one row per published
#      region x cluster: peak_id (published HOMER id), chr, start (0-based, BED), end,
#      cluster, plotted (in the 1,500-region Conserved subset drawn in Fig 1A;
#      TRUE for all non-Conserved regions), consensus_peak_id (overlapping
#      nf-core consensus peak(s), ";"-joined; NA if none), consensus_peak_primary
#      (the one with the largest overlap; used to look up DA results);
#   2. writes region BEDs to results/paper/fig1_regions/ for the region-based
#      panels (Conserved = full 21,649 set; Conserved_plotted = drawn subset);
#   3. reports the region -> consensus-peak mapping rate;
#   4. compares published vs de novo clusters. The main mismatch - published
#      Early Activation regions falling in a de novo "Naive" cluster - reflects
#      the absent 48h column in the de novo clustering (48h ATAC failed QC and
#      is excluded from clustering), so regions opened transiently at 48h/d3
#      cannot be separated from Naive-open regions there;
#      uses figS1_denovo_clusters.R output, if run.
#
# Inputs:  data/metadata/paper_ocr_clusters/mm39/
#          nf-core/atacseq consensus peaks (results/atac/.../consensus/)
#          results/paper/fig1_denovo_ocr_clusters.csv (optional)
# Outputs: results/paper/fig1a_ocr_clusters.csv, results/paper/fig1_regions/*.bed
#          results/paper/fig1_cluster_concordance.csv,
#          figures/paper/figS1_cluster_concordance.{pdf,png}
# Usage:   Rscript scripts/paper/fig1_published_clusters.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
})

src <- file.path(paths$metadata, "paper_ocr_clusters/mm39")
files <- c(
  "Conserved"        = "conserved.allSamps.peaks.TotalSet.bed",
  "Naive"            = "naive.specific.sig.bed",
  "Early Activation" = "early.activation.specific.sig.bed",
  "Activation"       = "activation.specific.sig.bed",
  "Late Activation"  = "late.activation.specific.sig.bed"
)
read_bed5 <- function(f) read_tsv(f, col_names = c("chr", "start", "end", "peak_id", "strand"),
                                  col_types = "ciicc")

regions <- imap_dfr(files, ~ read_bed5(file.path(src, .x)) |> mutate(cluster = .y)) |>
  dplyr::select(peak_id, chr, start, end, cluster)
plotted_ids <- read_bed5(file.path(src, "conserved.allSamps.plottedSubSet.peaks.bed"))$peak_id
# One region (Peak_76948) is listed in both Early Activation and Activation in
# the original mm10 files; it is kept in both, as in the published counts.
stopifnot(all(plotted_ids %in% regions$peak_id[regions$cluster == "Conserved"]),
          !anyDuplicated(regions[c("peak_id", "cluster")]))
regions <- regions |>
  mutate(plotted = cluster != "Conserved" | peak_id %in% plotted_ids)

# -----------------------------------------------------------------------------
# Map to nf-core consensus peaks
# -----------------------------------------------------------------------------
cons <- read_tsv(file.path(paths$atac, "bowtie2/merged_replicate/macs2/narrow_peak/consensus",
                           "consensus_peaks.mRp.clN.bed"),
                 col_names = c("chr", "start", "end", "id", "score", "strand"), col_types = "ciicic")
gr_r <- GRanges(regions$chr, IRanges(regions$start + 1L, regions$end))
gr_c <- GRanges(cons$chr, IRanges(cons$start + 1L, cons$end))
hits <- findOverlaps(gr_r, gr_c)
ov <- tibble(q = queryHits(hits), s = subjectHits(hits),
             width = width(pintersect(gr_r[queryHits(hits)], gr_c[subjectHits(hits)]))) |>
  mutate(id = cons$id[s])
map <- ov |> group_by(q) |>
  summarise(consensus_peak_id = paste(id, collapse = ";"),
            consensus_peak_primary = id[which.max(width)], .groups = "drop")

regions <- regions |>
  mutate(q = row_number()) |>
  left_join(map, by = "q") |>
  dplyr::select(-q) |>
  mutate(cluster = factor(cluster, levels = OCR_CLUSTERS)) |>
  arrange(cluster, chr, start)

maprate <- regions |>
  group_by(cluster) |>
  summarise(n_regions = n(), n_plotted = sum(plotted),
            n_mapped = sum(!is.na(consensus_peak_primary)),
            pct_mapped = 100 * n_mapped / n_regions,
            n_multi = sum(grepl(";", consensus_peak_id)), .groups = "drop")
message("Published region -> nf-core consensus peak mapping:"); print(maprate)

write_panel_table(regions, "fig1a_ocr_clusters")
write_panel_table(maprate, "fig1a_consensus_mapping")

# -----------------------------------------------------------------------------
# Region BEDs for region-based panels
# -----------------------------------------------------------------------------
bed_dir <- file.path(paths$paper_tab, "fig1_regions")
dir.create(bed_dir, showWarnings = FALSE)
write_bed <- function(df, name) {
  df |> transmute(chr, start, end, name = peak_id, score = 0, strand = ".") |>
    write_tsv(file.path(bed_dir, paste0(name, ".bed")), col_names = FALSE)
}
for (cl in OCR_CLUSTERS) write_bed(filter(regions, cluster == cl), gsub(" ", "_", cl))
write_bed(filter(regions, cluster == "Conserved", plotted), "Conserved_plotted")
message("Region BEDs in ", bed_dir)

# -----------------------------------------------------------------------------
# Concordance with de novo clusters (validation)
# -----------------------------------------------------------------------------
f_dn <- file.path(paths$paper_tab, "fig1_denovo_ocr_clusters.csv")
if (file.exists(f_dn)) {
  dn <- read_csv(f_dn, show_col_types = FALSE) |>
    dplyr::select(consensus_peak_primary = peak_id, denovo = cluster)
  conc <- regions |>
    left_join(dn, by = "consensus_peak_primary") |>
    mutate(denovo = case_when(
      is.na(consensus_peak_primary) ~ "No consensus peak",
      is.na(denovo) ~ "Unassigned",
      TRUE ~ denovo),
      denovo = factor(denovo, levels = c(OCR_CLUSTERS, "Unassigned", "No consensus peak"))) |>
    dplyr::count(published = cluster, denovo, .drop = FALSE) |>
    group_by(published) |> mutate(pct_of_published = 100 * n / sum(n)) |> ungroup()
  write_panel_table(conc, "fig1_cluster_concordance")

  agree <- conc |> filter(as.character(published) == as.character(denovo))
  message("Published regions recovered in the same de novo cluster:")
  print(agree |> dplyr::select(published, n, pct_of_published))

  m <- conc |> dplyr::select(published, denovo, pct_of_published) |>
    pivot_wider(names_from = denovo, values_from = pct_of_published) |>
    column_to_rownames("published") |> as.matrix()
  n <- conc |> dplyr::select(published, denovo, n) |>
    pivot_wider(names_from = denovo, values_from = n) |>
    column_to_rownames("published") |> as.matrix()
  ht <- Heatmap(m, name = "% of published\ncluster",
                col = colorRamp2(c(0, 50, 100), viridisLite::mako(3, direction = -1)),
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_title = "Published (Fig 1A)", column_title = "De novo (reprocessed data; no 48h ATAC)",
                row_names_side = "left", column_names_rot = 45,
                rect_gp = gpar(col = "white", lwd = 1),
                cell_fun = function(j, i, x, y, w, h, fill)
                  grid.text(sprintf("%d\n%.0f%%", n[i, j], m[i, j]), x, y,
                            gp = gpar(fontsize = 7, col = ifelse(m[i, j] > 50, "white", "black"))))
  save_panel_heatmap(ht, "figS1_cluster_concordance", width = 7, height = 3.6)
} else {
  message("De novo clusters not found; run scripts/paper/figS1_denovo_clusters.R for the concordance panel.")
}
message("Done.")
