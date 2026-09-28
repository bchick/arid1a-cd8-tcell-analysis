#!/usr/bin/env Rscript
# =============================================================================
# chromvar/build_chromvar.R — chromVAR TF-motif deviations by ARID1A dose (D8 WT/Het/KO)
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Precompute step for dose_chromvar (figures/chromvar_dose.R) and
# figures/tobias_dose_response.R. Adapted from the chromVAR section of
# integration/01_multiomic_integration.R, building the SummarizedExperiment
# directly from the nf-core/atacseq featureCounts matrix (no DiffBind object).
# Uses the D8 Exp2 WT/Het/KO × TE/EEC/MP design (matches the het dose-response
# analysis). GC bias + motif matching use the local reference FASTA (seqnames
# match the peaks; avoids the UCSC BSgenome). Motifs: JASPAR2020 CORE
# vertebrates (the JASPAR2024 package download URL is unavailable).
# Needs the raw nf-core/atacseq consensus featureCounts and the mm39 genome FASTA.
# Also writes a family-level trajectory version of dose_chromvar, which
# figures/chromvar_dose.R overwrites with the retention-scaled version.
#
# Inputs:  results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/
#            consensus_peaks.mRp.clN.featureCounts.txt
#          data/reference/GRCm39.primary_assembly.genome.fa (genome$fasta)
# Outputs: results/extended_analysis/chromvar/chromvar_dose.RData  (dev, dev_z, sample_meta)
#          results/extended_analysis/chromvar/chromvar_motif_genotype_means.csv
#          figures/extended_analysis/dose_chromvar.{pdf,png}
# Usage:   Rscript extended_analysis/scripts/chromvar/build_chromvar.R   (from the repository root)
# =============================================================================

t0 <- Sys.time()
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

source("scripts/utils.R")
source("extended_analysis/scripts/utils_figures.R")
suppressPackageStartupMessages({
  library(tidyverse); library(SummarizedExperiment); library(GenomicRanges)
  library(chromVAR); library(motifmatchr); library(TFBSTools); library(JASPAR2020)
  library(Rsamtools); library(Biostrings); library(BiocParallel)
})
register(MulticoreParam(8))

FA  <- genome$fasta
FC  <- file.path(paths$atac,
  "bowtie2/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt")
STD_CHR <- paste0("chr", c(as.character(1:19), "X", "Y"))  # GENCODE mm39 naming (matches FASTA)
outdir <- file.path(paths$ext_results, "chromvar")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. load featureCounts, subset to het dose-response design --------------
log("loading featureCounts matrix")
fc <- read_tsv(FC, comment = "#", show_col_types = FALSE)
samp_cols <- setdiff(colnames(fc), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
samp_ids  <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", samp_cols)
keep <- grepl("^D8_(WT|Het|KO)_(TE|EEC|MP)_Exp2_REP[0-9]+$", samp_ids, ignore.case = TRUE)
sel_cols <- samp_cols[keep]; sel_ids <- samp_ids[keep]

sample_meta <- tibble(sample = sel_ids) %>%
  separate(sample, c("tp", "genotype", "subset", "exp", "rep"), sep = "_", remove = FALSE) %>%
  mutate(genotype = gfac(genotype), subset = sfac(subset),
         dose = c(WT = 2, Het = 1, KO = 0)[as.character(genotype)])
log("samples:", nrow(sample_meta))
print(table(sample_meta$genotype, sample_meta$subset))

# peaks: keep standard chromosomes only (drop scaffolds for clean motif matching)
keep_pk <- fc$Chr %in% STD_CHR
counts  <- as.matrix(fc[keep_pk, sel_cols]); colnames(counts) <- sel_ids
peaks   <- GRanges(fc$Chr[keep_pk], IRanges(fc$Start[keep_pk], fc$End[keep_pk]),
                   peak_id = fc$Geneid[keep_pk])
log("peaks on standard chr:", length(peaks), "of", nrow(fc))
stopifnot("peak/chr filter matched nothing — check seqname naming" = length(peaks) > 1000)

se <- SummarizedExperiment(assays = list(counts = counts), rowRanges = peaks,
                           colData = DataFrame(sample_meta, row.names = sel_ids))

# ---- 2. filter peaks + GC bias from FASTA -----------------------------------
log("filtering peaks")
se <- filterPeaks(se, non_overlapping = TRUE)
log("peaks after filter:", nrow(se))

log("computing GC bias from FASTA")
if (!file.exists(paste0(FA, ".fai"))) indexFa(FA)
seqs <- getSeq(FaFile(FA), rowRanges(se))
gc   <- as.numeric(letterFrequency(seqs, "GC", as.prob = TRUE))
rowData(se)$bias <- gc
rr <- rowRanges(se); rr$bias <- gc; rowRanges(se) <- rr

# ---- 3. JASPAR2020 motifs + matches -----------------------------------------
log("retrieving JASPAR2020 CORE vertebrate PWMs")
motifs <- getMatrixSet(JASPAR2020, list(collection = "CORE", tax_group = "vertebrates",
                                        matrixtype = "PWM"))
log("matching", length(motifs), "motifs across", nrow(se), "peaks")
motif_ix <- matchMotifs(motifs, se, genome = FaFile(FA))
log("matchMotifs done")

# ---- 4. deviations -----------------------------------------------------------
set.seed(42)
log("computing deviations")
dev   <- computeDeviations(object = se, annotations = motif_ix)
dev_z <- deviationScores(dev)                      # motif x sample
tf_name <- name(motifs); names(tf_name) <- TFBSTools::ID(motifs)
rownames(dev_z) <- ifelse(rownames(dev_z) %in% names(tf_name),
                          tf_name[rownames(dev_z)], rownames(dev_z))
save(dev, dev_z, sample_meta, file = file.path(outdir, "chromvar_dose.RData"))
log("deviations saved")

# ---- 5. TF family assignment + genotype means -------------------------------
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

z_long <- as.data.frame(dev_z) %>% rownames_to_column("tf") %>%
  pivot_longer(-tf, names_to = "sample", values_to = "z") %>%
  left_join(sample_meta, by = "sample") %>%
  mutate(family = assign_family(tf))

# family-level deviation per sample, then genotype mean ± SE (pooled subsets)
fam_keep <- c("ETS", "RUNX", "T-box", "AP-1/bZIP", "NFkB", "TCF/LEF", "KLF/SP")
fam_dose <- z_long %>% filter(family %in% fam_keep) %>%
  group_by(family, sample, genotype) %>% summarise(z = mean(z), .groups = "drop") %>%
  group_by(family, genotype) %>%
  summarise(mean = mean(z), se = sd(z) / sqrt(n()), .groups = "drop") %>%
  mutate(family = factor(family, levels = fam_keep))

# per-motif genotype means for the record
z_long %>% group_by(tf, family, genotype) %>% summarise(mean_z = mean(z), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = mean_z) %>%
  write_csv(file.path(outdir, "chromvar_motif_genotype_means.csv"))

# ---- 6. figure: dose_chromvar (family trajectories) -------------------------
p_chromvar <- ggplot(fam_dose, aes(genotype, mean, color = family, group = family)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = family), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_color_manual(values = pal_tf_family, name = "TF family") +
  scale_fill_manual(values = pal_tf_family, guide = "none") +
  scale_x_discrete(labels = c("WT\n(2)", "Het\n(1)", "KO\n(0)")) +
  labs(x = expression(italic("Arid1a")~"copies"),
       y = "Mean chromVAR\ndeviation z-score",
       title = "TF-motif accessibility tracks ARID1A dose",
       subtitle = "chromVAR (JASPAR2020) on D8 OCRs, pooled across TE/EEC/MP") +
  theme_ext
save_ext_figure(p_chromvar, "dose_chromvar", width = 7, height = 5)

log("DONE in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
cat("CHROMVAR_BUILD_OK\n")
