#!/usr/bin/env Rscript
# =============================================================================
# core/05_chipseq_analysis.R — T-bet ChIP-seq analysis
# McDonald, Chick et al. 2023 Immunity 56:1303 — core analysis
#
# T-bet ChIP-seq +/- BAF inhibitors at 48h (GSE228546; one replicate each):
# Untreated, IL-12, IL-12 + ACBI1, IL-12 + BRM014.
#   1. T-bet binding sites at baseline and with IL-12
#   2. Effect of the BAF inhibitors on T-bet occupancy (lost / retained /
#      gained peaks relative to IL-12)
#   3. Genomic annotation (ChIPseeker) and a shell script of HOMER motif
#      enrichment commands (written, not run)
#
# Inputs:  results/chipseq/bowtie2/merged_library/macs3/narrow_peak/*.narrowPeak
#          data/reference/gencode.vM35.primary_assembly.annotation.gtf
# Outputs: results/chipseq/differential/ (manifest, overlap BEDs and table,
#            Tbet_IL12_annotated.csv, chipseq_analysis.RData)
#          results/chipseq/run_homer_chipseq.sh
#          figures/chipseq/
# Usage:   Rscript scripts/core/05_chipseq_analysis.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(GenomicFeatures)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})
select <- dplyr::select

outdir <- file.path(paths$results, "chipseq/differential")
figdir <- file.path(paths$figures, "chipseq")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Sample manifest (from actual files on disk)
# =============================================================================

message("=== ChIP-seq sample manifest ===")

bam_dir  <- file.path(paths$chipseq, "bowtie2/merged_library")
peak_dir <- file.path(paths$chipseq, "bowtie2/merged_library/macs3/narrow_peak")

samples <- tribble(
  ~sample_id,                        ~antibody, ~treatment,    ~is_input,
  "ChIP_Untreated_T-bet_REP1",       "Tbet",    "Untreated",   FALSE,
  "ChIP_IL-12_T-bet_REP1",           "Tbet",    "IL12",        FALSE,
  "ChIP_IL-12_ACBI1_T-bet_REP1",     "Tbet",    "IL12_ACBI1",  FALSE,
  "ChIP_IL-12_BRM014_T-bet_REP1",    "Tbet",    "IL12_BRM014", FALSE,
  "ChIP_Untreated_Input_REP1",       "Input",   "Untreated",   TRUE,
  "ChIP_IL-12_Input_REP1",           "Input",   "IL12",        TRUE,
  "ChIP_IL-12_ACBI1_Input_REP1",     "Input",   "IL12_ACBI1",  TRUE,
  "ChIP_IL-12_BRM014_Input_REP1",    "Input",   "IL12_BRM014", TRUE,
) %>%
  mutate(
    bam  = file.path(bam_dir, paste0(sample_id, ".mLb.clN.sorted.bam")),
    peak = file.path(peak_dir, paste0(sample_id, "_peaks.narrowPeak")),
    bam_exists  = file.exists(bam),
    peak_exists = file.exists(peak)
  )

message("Files found:")
samples %>% filter(!is_input) %>%
  select(sample_id, bam_exists, peak_exists) %>%
  { message(capture.output(print(.))); . }

write_csv(samples, file.path(outdir, "chipseq_sample_manifest.csv"))

# =============================================================================
# 2. Load T-bet peaks
# =============================================================================

message("=== Loading T-bet peaks ===")

chip_peaks <- samples %>% filter(!is_input, peak_exists)

peaks_gr <- setNames(
  lapply(chip_peaks$peak, function(p) {
    tryCatch(import(p, format = "narrowPeak"), error = function(e) GRanges())
  }),
  chip_peaks$treatment
)

peak_counts <- sapply(peaks_gr, length)
message("T-bet peak counts:")
for (nm in names(peak_counts)) message(sprintf("  %s: %d peaks", nm, peak_counts[nm]))

# Bar plot: peak counts per condition
p_npeaks <- tibble(
  treatment = factor(names(peak_counts),
                     levels = c("Untreated", "IL12", "IL12_ACBI1", "IL12_BRM014")),
  n_peaks   = as.integer(peak_counts)
) %>%
  ggplot(aes(x = treatment, y = n_peaks, fill = treatment)) +
  geom_col(color = "black", width = 0.6) +
  geom_text(aes(label = scales::comma(n_peaks)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(
    Untreated   = "#AEC7E8",
    IL12        = "#1F77B4",
    IL12_ACBI1  = "#D62728",
    IL12_BRM014 = "#9467BD"
  )) +
  scale_x_discrete(labels = c("Untreated", "IL-12", "IL-12\n+ACBI1", "IL-12\n+BRM014")) +
  labs(title = "T-bet ChIP-seq: peak counts by condition",
       x = "Treatment", y = "Number of peaks") +
  theme_paper + theme(legend.position = "none")

save_figure(p_npeaks, "Tbet_peak_counts", width = 6, height = 5, dir = figdir)

# =============================================================================
# 3. Peak overlap: IL-12 vs inhibitor conditions
# =============================================================================

message("=== T-bet peak overlap across conditions ===")

# Reference: IL-12 T-bet peaks (baseline for inhibitor comparison)
il12_gr <- peaks_gr[["IL12"]]

if (length(il12_gr) > 0) {
  overlap_summary <- list()

  for (cond in c("IL12_ACBI1", "IL12_BRM014")) {
    if (!cond %in% names(peaks_gr) || length(peaks_gr[[cond]]) == 0) next

    cond_gr <- peaks_gr[[cond]]
    shared      <- sum(countOverlaps(il12_gr, cond_gr) > 0)
    il12_only   <- sum(countOverlaps(il12_gr, cond_gr) == 0)
    cond_only   <- sum(countOverlaps(cond_gr, il12_gr) == 0)

    overlap_summary[[cond]] <- tibble(
      comparison  = paste0("IL12_vs_", cond),
      IL12_only   = il12_only,
      Shared      = shared,
      Cond_only   = cond_only,
      IL12_total  = length(il12_gr),
      Cond_total  = length(cond_gr),
      pct_retained = shared / length(il12_gr) * 100
    )

    message(sprintf("  IL-12 vs %s: IL12-only=%d, Shared=%d, %s-only=%d (%.1f%% retained)",
                    cond, il12_only, shared, cond, cond_only,
                    shared / length(il12_gr) * 100))

    # Export BED
    export(il12_gr[countOverlaps(il12_gr, cond_gr) == 0],
           file.path(outdir, paste0("Tbet_IL12_lost_in_", cond, ".bed")))
    export(il12_gr[countOverlaps(il12_gr, cond_gr) > 0],
           file.path(outdir, paste0("Tbet_IL12_retained_in_", cond, ".bed")))
    export(cond_gr[countOverlaps(cond_gr, il12_gr) == 0],
           file.path(outdir, paste0("Tbet_", cond, "_gained.bed")))
  }

  overlap_df <- bind_rows(overlap_summary)
  write_csv(overlap_df, file.path(outdir, "Tbet_inhibitor_peak_overlap.csv"))

  # Bar chart: retained/lost
  overlap_long <- overlap_df %>%
    transmute(
      inhibitor = gsub("IL12_vs_IL12_", "", comparison),
      Retained  = Shared,
      Lost      = IL12_only
    ) %>%
    pivot_longer(-inhibitor, names_to = "category", values_to = "n_peaks") %>%
    mutate(category = factor(category, c("Retained", "Lost")))

  p_inhibitor <- ggplot(overlap_long, aes(x = inhibitor, y = n_peaks, fill = category)) +
    geom_col(position = "stack", color = "black", width = 0.6) +
    scale_fill_manual(values = c(Retained = "#1F77B4", Lost = "#D62728")) +
    labs(title = "T-bet ChIP-seq: effect of BAF inhibitors on peak retention",
         subtitle = "Relative to IL-12 treatment",
         x = "BAF inhibitor", y = "T-bet peaks", fill = "") +
    theme_paper

  save_figure(p_inhibitor, "Tbet_inhibitor_peak_retention", width = 6, height = 5, dir = figdir)
}

# =============================================================================
# 4. Genomic annotation
# =============================================================================

message("=== Annotating T-bet peaks ===")

txdb <- get_txdb()

non_empty <- Filter(function(gr) length(gr) > 0, peaks_gr)
anno_list <- lapply(non_empty, function(gr) {
  annotatePeak(gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
})

if (length(anno_list) > 0) {
  p_anno <- plotAnnoBar(anno_list) + theme_paper +
    ggtitle("T-bet ChIP-seq: genomic annotation by condition") +
    theme(axis.text.y = element_text(size = 8))
  save_figure(p_anno, "Tbet_annotation_bar", width = 9, height = 5, dir = figdir)
}

# Export IL-12 T-bet annotated peak table
if ("IL12" %in% names(anno_list)) {
  il12_anno <- as.data.frame(anno_list[["IL12"]]) %>%
    as_tibble() %>%
    mutate(ensembl_clean = gsub("\\.\\d+$", "", geneId)) %>%
    left_join(
      tryCatch(
        AnnotationDbi::select(org.Mm.eg.db,
          keys = gsub("\\.\\d+$", "", as.data.frame(anno_list[["IL12"]])$geneId),
          keytype = "ENSEMBL", columns = "SYMBOL"),
        error = function(e) data.frame(ENSEMBL = character(), SYMBOL = character())
      ),
      by = c("ensembl_clean" = "ENSEMBL")
    )
  write_csv(il12_anno, file.path(outdir, "Tbet_IL12_annotated.csv"))
}

# =============================================================================
# 5. Write HOMER motif enrichment commands
# =============================================================================

message("=== Writing HOMER motif enrichment commands ===")

homer_cmds <- character()
genome_fa <- genome$fasta

for (cond in names(peaks_gr)) {
  gr <- peaks_gr[[cond]]
  if (length(gr) == 0) next
  bed_file    <- file.path(outdir, paste0("Tbet_", cond, "_peaks.bed"))
  homer_outdir <- file.path(paths$results, "chipseq/motifs", cond)
  export(gr, bed_file)
  homer_cmds <- c(homer_cmds, sprintf(
    "findMotifsGenome.pl %s %s %s -size 200 -mask -p 8",
    bed_file, genome_fa, homer_outdir
  ))
}

writeLines(
  c("#!/bin/bash", "# HOMER motif enrichment for T-bet ChIP-seq", homer_cmds),
  file.path(paths$results, "chipseq/run_homer_chipseq.sh")
)
message("  Wrote HOMER commands to chipseq/run_homer_chipseq.sh")

# =============================================================================
# Save session
# =============================================================================

message("=== Saving R session ===")
save.image(file.path(outdir, "chipseq_analysis.RData"))
message("ChIP-seq analysis complete.")
message("  Results: ", outdir)
message("  Figures: ", figdir)
