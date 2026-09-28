#!/usr/bin/env Rscript
# =============================================================================
# paper/figS_tracks.R — Genome-browser signal tracks: Figures S1A, S2D, S4A, S4B
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
#   S1A  ARID1A CUT&RUN and ATAC-seq (WT time course) at Zeb2, Tbx21, Bhlhe40
#   S2D  ATAC-seq WT / cHet / cKO at Bhlhe40, Tbx21, Zeb2, Batf
#   S4A  ATAC-seq and H3K27ac CUT&RUN at Bhlhe40, Tbx21, Prdm1, Zfp683
#   S4B  ATAC-seq and H3K27me3 CUT&RUN at Tcf7, Cd9, Ccr7, Sell
#
# Signal: nf-core bigWigs (ATAC merged-replicate, scaled to 1M mapped reads;
# CUT&RUN CPM). Signal is averaged in 400 bins per locus; when a condition has
# several bigWigs (e.g. day 8 TE/EEC/MP) they are averaged. Day 8 ATAC uses the
# Exp2 batch only, the one batch with WT, cHet and cKO, so genotypes are
# batch-matched. 48h ATAC is excluded (failed QC). GEO has no day 8 cKO H3K27ac
# or any cKO H3K27me3, so those tracks are WT only.
# Needs the upstream (raw-data) mode: bigWigs are not in the data bundle.
#
# Inputs:  results/atac/bowtie2/merged_replicate/bigwig/*.bigWig
#          results/cutrun/03_peak_calling/03_bed_to_bigwig/*.bigWig
#          data/reference/gencode.vM35.primary_assembly.annotation.gtf
# Outputs: figures/paper/figS{1a,2d,4a,4b}_tracks_<gene>.{pdf,png}
# Usage:   Rscript scripts/paper/figS_tracks.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(patchwork)
})
select <- dplyr::select; filter <- dplyr::filter

atac_bw <- file.path(paths$atac, "bowtie2/merged_replicate/bigwig")
cr_bw   <- file.path(paths$cutrun, "03_peak_calling/03_bed_to_bigwig")

atac <- function(...) file.path(atac_bw, paste0(c(...), ".mRp.clN.bigWig"))
cr   <- function(...) file.path(cr_bw, paste0(c(...), ".bigWig"))
d8_atac <- function(geno) atac(paste0("D8_", geno, "_", c("TE", "EEC", "MP"), "_Exp2"))

if (!require_inputs(c(genome$gtf, list.files(atac_bw, full.names = TRUE)[1],
                      list.files(cr_bw, full.names = TRUE)[1]),
                    "bigWigs / GENCODE GTF")) quit(save = "no", status = 0)

# Track definitions: label, assay (sets the shared y-scale and colour), files
trk <- function(label, assay, files, genotype = "WT")
  list(label = label, assay = assay, files = files, genotype = genotype)

panels <- list(
  figS1a = list(
    genes = c("Zeb2", "Tbx21", "Bhlhe40"),
    tracks = list(
      trk("ATAC Naive", "ATAC", atac("Naive_WT")),
      trk("ATAC d3",    "ATAC", atac("D3_WT")),
      trk("ATAC d5",    "ATAC", atac("D5_WT")),
      trk("ATAC d8",    "ATAC", d8_atac("WT")),
      trk("ARID1A Naive", "ARID1A", cr("ARID1A_Naive_WT_R1")),
      trk("ARID1A 48h",   "ARID1A", cr(paste0("ARID1A_48h_WT_R", 1:2))),
      trk("ARID1A d5",    "ARID1A", cr(paste0("ARID1A_D5_WT_R", 1:2))),
      trk("ARID1A d8",    "ARID1A", cr(paste0("ARID1A_D8_WT_", c("TE", "EEC", "MP"), "_R1"))))),
  figS2d = list(
    genes = c("Bhlhe40", "Tbx21", "Zeb2", "Batf"),
    tracks = list(
      trk("Naive WT", "ATAC", atac("Naive_WT")),
      trk("d3 WT",    "ATAC", atac("D3_WT")),
      trk("d3 cKO",   "ATAC", atac("D3_KO"), "KO"),
      trk("d5 WT",    "ATAC", atac("D5_WT")),
      trk("d5 cKO",   "ATAC", atac("D5_KO"), "KO"),
      trk("d8 WT",    "ATAC", d8_atac("WT")),
      trk("d8 cHet",  "ATAC", d8_atac("Het"), "Het"),
      trk("d8 cKO",   "ATAC", d8_atac("KO"), "KO"))),
  figS4a = list(
    genes = c("Bhlhe40", "Tbx21", "Prdm1", "Zfp683"),
    tracks = list(
      trk("ATAC d5 WT",  "ATAC", atac("D5_WT")),
      trk("ATAC d5 cKO", "ATAC", atac("D5_KO"), "KO"),
      trk("ATAC d8 WT",  "ATAC", d8_atac("WT")),
      trk("ATAC d8 cKO", "ATAC", d8_atac("KO"), "KO"),
      trk("H3K27ac d5 WT",  "H3K27ac", cr("H3K27ac_D5_WT_R1")),
      trk("H3K27ac d5 cKO", "H3K27ac", cr("H3K27ac_D5_KO_R1"), "KO"),
      trk("H3K27ac d8 WT",  "H3K27ac", cr(paste0("H3K27ac_D8_WT_", c("TE", "EEC", "MP"), "_R1"))))),
  figS4b = list(
    genes = c("Tcf7", "Cd9", "Ccr7", "Sell"),
    tracks = list(
      trk("ATAC d8 TE WT",   "ATAC", atac("D8_WT_TE_Exp2")),
      trk("ATAC d8 TE cKO",  "ATAC", atac("D8_KO_TE_Exp2"), "KO"),
      trk("ATAC d8 EEC WT",  "ATAC", atac("D8_WT_EEC_Exp2")),
      trk("ATAC d8 EEC cKO", "ATAC", atac("D8_KO_EEC_Exp2"), "KO"),
      trk("ATAC d8 MP WT",   "ATAC", atac("D8_WT_MP_Exp2")),
      trk("ATAC d8 MP cKO",  "ATAC", atac("D8_KO_MP_Exp2"), "KO"),
      trk("H3K27me3 d8 TE",  "H3K27me3", cr("H3K27me3_D8_WT_TE_R1")),
      trk("H3K27me3 d8 EEC", "H3K27me3", cr("H3K27me3_D8_WT_EEC_R1")),
      trk("H3K27me3 d8 MP",  "H3K27me3", cr("H3K27me3_D8_WT_MP_R1"))))
)

all_files <- unique(unlist(lapply(panels, function(p) lapply(p$tracks, `[[`, "files"))))
stopifnot("Missing bigWigs" = require_inputs(all_files, "track bigWigs"))

# =============================================================================
# Gene models (Ensembl canonical transcript per gene, GENCODE vM35)
# =============================================================================

all_genes <- unique(unlist(lapply(panels, `[[`, "genes")))
gtf_sub <- tempfile(fileext = ".gtf")
system2("grep", c("-E", shQuote(paste0('gene_name "(', paste(all_genes, collapse = "|"), ')";')),
                  genome$gtf), stdout = gtf_sub)
gtf <- import(gtf_sub, format = "gtf")
gtf <- gtf[gtf$gene_name %in% all_genes & gtf$gene_type == "protein_coding"]

# rtracklayer keeps one value of repeated attributes (tag), so read the
# Ensembl_canonical transcript IDs from the GTF text directly
gtf_txt <- readLines(gtf_sub)
canon_ids <- sub('.*transcript_id "([^"]+)".*', "\\1",
                 grep('\ttranscript\t.*tag "Ensembl_canonical"', gtf_txt, value = TRUE))
canonical <- gtf[gtf$type == "transcript" & gtf$transcript_id %in% canon_ids]
stopifnot("No canonical transcript for some genes" = all(all_genes %in% canonical$gene_name))
exons <- gtf[gtf$type == "exon" & gtf$transcript_id %in% canonical$transcript_id]

#' Locus window: gene body plus 30% flank on each side (min 15 kb)
locus_of <- function(g) {
  tx <- canonical[canonical$gene_name == g][1]
  flank <- max(15000, round(0.3 * width(tx)))
  GRanges(seqnames(tx), IRanges(max(1, start(tx) - flank), end(tx) + flank))
}

# =============================================================================
# Signal extraction
# =============================================================================

N_BINS <- 400

bin_signal <- function(files, region) {
  bins <- tile(region, n = N_BINS)[[1]]
  vals <- sapply(files, function(f) {
    cov <- import(f, which = region, as = "RleList")[[as.character(seqnames(region))]]
    as.numeric(Views(cov, ranges(bins)) |> viewMeans())
  })
  tibble(pos = (start(bins) + end(bins)) / 2, signal = rowMeans(matrix(vals, ncol = length(files))))
}

# WT tracks are coloured by assay; cHet / cKO tracks use the genotype palette
assay_cols <- c(ATAC = "#8B1A1A", ARID1A = "#1B7837", H3K27ac = "#762A83", H3K27me3 = "#E08214")

plot_locus <- function(panel_id, gene) {
  p <- panels[[panel_id]]
  region <- locus_of(gene)
  df <- purrr::map_dfr(seq_along(p$tracks), function(i) {
    t <- p$tracks[[i]]
    bin_signal(t$files, region) |>
      mutate(label = t$label, assay = t$assay, genotype = t$genotype, order = i)
  }) |>
    mutate(label = factor(label, levels = vapply(p$tracks, `[[`, "", "label")))

  # One y-scale per assay within the locus
  ymax <- df |> group_by(assay) |> summarise(ymax = max(signal), .groups = "drop")
  fmt <- function(x) trimws(formatC(x, digits = 2, format = "fg"))
  df <- df |> left_join(ymax, by = "assay") |>
    mutate(fill = ifelse(genotype == "WT", unname(assay_cols[assay]),
                         unname(pal_genotype[genotype])),
           label = factor(sprintf("%s [0-%s]", label, fmt(ymax)),
                          levels = unique(sprintf("%s [0-%s]", label, fmt(ymax))[order(order)])))

  tracks <- ggplot(df, aes(pos, signal)) +
    geom_area(aes(fill = fill), outline.type = "upper", linewidth = 0) +
    geom_blank(aes(y = ymax)) +
    scale_fill_identity() +
    facet_grid(label ~ ., scales = "free_y", switch = "y") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = gene) +
    theme_paper +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
          strip.placement = "outside", strip.background = element_blank(),
          panel.spacing = unit(0.5, "mm"), panel.border = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 10))

  ex <- exons[exons$gene_name == gene]
  tx <- canonical[canonical$gene_name == gene][1]
  strand_arrow <- if (as.character(strand(tx)) == "+") "last" else "first"
  model <- ggplot() +
    annotate("segment", x = start(tx), xend = end(tx), y = 0, yend = 0,
             linewidth = 0.4, arrow = arrow(length = unit(1.5, "mm"), ends = strand_arrow)) +
    annotate("rect", xmin = start(ex), xmax = end(ex), ymin = -0.4, ymax = 0.4,
             fill = "black") +
    scale_x_continuous(limits = c(start(region), end(region)), expand = c(0, 0),
                       labels = function(x) sprintf("%.2f Mb", x / 1e6), n.breaks = 3) +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(x = as.character(seqnames(region)), y = NULL) +
    theme_paper +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.border = element_blank(), axis.line.x = element_line(linewidth = 0.3),
          axis.text.x = element_text(size = 6), axis.title.x = element_text(size = 7))

  n <- length(p$tracks)
  out <- tracks / model + plot_layout(heights = c(n, 0.9))
  save_panel(out, sprintf("%s_tracks_%s", panel_id, gene), width = 3.6, height = 0.45 * n + 1.1)
  df |> distinct(label, assay, ymax) |> mutate(panel = panel_id, gene = gene,
                                               locus = as.character(region))
}

scales_used <- purrr::map_dfr(names(panels), function(pid) {
  message("=== ", pid, " ===")
  purrr::map_dfr(panels[[pid]]$genes, function(g) plot_locus(pid, g))
})
write_panel_table(scales_used, "figS_tracks_loci")

message("Signal-track panels complete.")
