#!/usr/bin/env Rscript
# =============================================================================
# integration/01_multiomic_integration.R — multi-omic integration (ATAC + RNA + CUT&RUN)
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Links D8 KO-vs-WT differential accessibility, differential expression and
# ARID1A CUT&RUN binding at the gene level:
#   - chromatin accessibility vs gene expression correlation (Ensembl-ID join)
#   - ARID1A binding -> accessibility -> expression (Wilcoxon + Fisher tests)
#   - effector TF regulatory summary
#   - chromVAR motif deviations (loaded from chromvar/build_chromvar.R, not recomputed)
#
# Inputs:  results/rnaseq/differential/de_KO_vs_WT_D8_pseudobulk.csv
#          results/atac/differential/{da_KO_vs_WT_D8_pseudobulk.csv,consensus_peaks_annotated.csv}
#          results/cutrun/differential/ARID1A_D5_{WT_only,shared}.bed
#          results/extended_analysis/chromvar/chromvar_dose.RData
# Outputs: results/extended_analysis/integration/ (CSVs, RData)
#          figures/extended_analysis/integration/
# Usage:   Rscript extended_analysis/scripts/integration/01_multiomic_integration.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(DESeq2)
  library(ComplexHeatmap)
  library(ggrepel)
})

# Resolve namespace conflicts: GenomicFeatures pulls in AnnotationDbi, whose S4
# `select` generic masks dplyr::select and has no method for the spec_tbl_df that
# read_csv() returns. Same aliasing as core/01_rnaseq_analysis.R.
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate
slice  <- dplyr::slice

# Output directories
outdir <- file.path(paths$ext_results, "integration")
figdir <- file.path(paths$ext_figures, "integration")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load previous analysis results
# =============================================================================

message("=== Loading analysis results ===")

# RNA-seq DE results
rna_dir <- file.path(paths$results, "rnaseq/differential")
rna_d8 <- read_csv(file.path(rna_dir, "de_KO_vs_WT_D8_pseudobulk.csv"), show_col_types = FALSE)

# ATAC-seq DA results
atac_dir <- file.path(paths$results, "atac/differential")
atac_d8 <- read_csv(file.path(atac_dir, "da_KO_vs_WT_D8_pseudobulk.csv"), show_col_types = FALSE)

# CUT&RUN ARID1A — real per-peak data.
# An earlier version tryCatch'd onto D5_wt_ko_peak_overlap.csv, which is a 6-row
# per-ANTIBODY summary table, not peaks — hence the nonsensical log line
# "ARID1A CUT&RUN: 5 peaks". It was then never used in any computation.
# db_ARID1A_KO_vs_WT.csv does not exist (most antibodies are n=1/genotype, so
# core/04_cutandrun_analysis.R does peak overlap rather than DiffBind). The real
# per-peak products of that overlap are these BEDs + the annotated WT peak table.
cutrun_dir <- file.path(paths$results, "cutrun/differential")

read_peak_bed <- function(f) {
  bed <- read_tsv(f, col_names = c("chr", "start", "end"),
                  col_types = "cii", progress = FALSE)
  GRanges(bed$chr, IRanges(bed$start + 1L, bed$end))  # BED is 0-based half-open
}

arid1a_wt_only <- read_peak_bed(file.path(cutrun_dir, "ARID1A_D5_WT_only.bed"))
arid1a_shared  <- read_peak_bed(file.path(cutrun_dir, "ARID1A_D5_shared.bed"))
# cBAF-dependent ARID1A binding = present in WT, absent in KO
arid1a_dep <- c(arid1a_wt_only, arid1a_shared)

message("  RNA-seq: ", nrow(rna_d8), " genes")
message("  ATAC-seq: ", nrow(atac_d8), " peaks")
message("  ARID1A CUT&RUN: ", length(arid1a_wt_only), " WT-only + ",
        length(arid1a_shared), " shared peaks")

# =============================================================================
# 2. ATAC ↔ RNA correlation (Figure 4E style)
# =============================================================================

message("=== ATAC-RNA correlation ===")

# Load ATAC peak annotations (from ChIPseeker output)
atac_anno <- read_csv(file.path(atac_dir, "consensus_peaks_annotated.csv"),
                       show_col_types = FALSE)

# Link DA peaks to nearest gene, then join with RNA DE.
# Join on Ensembl ID, NOT symbol: ensembl_clean is 100% populated while SYMBOL is
# only 56.7% (peaks whose nearest gene is a Gm-/Rik-/lncRNA locus that
# org.Mm.eg.db has no symbol for — core/02_atacseq_analysis.R already does the
# best mapping available). Measured: the Ensembl join recovers 11,580 genes vs
# 10,233 on symbol (+1,347, +13.2%).
atac_with_gene <- atac_d8 %>%
  left_join(atac_anno %>% select(peak_id, ensembl_clean, SYMBOL),
            by = "peak_id") %>%
  filter(!is.na(ensembl_clean))

# For each gene, summarize ATAC changes (mean fold change of associated peaks)
gene_atac_summary <- atac_with_gene %>%
  group_by(ensembl_clean) %>%
  summarize(
    SYMBOL = dplyr::first(SYMBOL),  # constant within an Ensembl ID
    atac_lfc = mean(log2FoldChange, na.rm = TRUE),
    atac_min_fdr = min(padj, na.rm = TRUE),
    n_peaks = n(),
    n_lost = sum(log2FoldChange < -1 & padj < 0.05, na.rm = TRUE),
    n_gained = sum(log2FoldChange > 1 & padj < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# Join with RNA results on Ensembl ID (strip GENCODE version suffix)
integrated <- gene_atac_summary %>%
  inner_join(rna_d8 %>%
    mutate(ensembl_clean = gsub("\\.\\d+$", "", gene_id)) %>%
    select(ensembl_clean, rna_symbol = gene_name,
           rna_lfc = log2FoldChange, rna_padj = padj),
    by = "ensembl_clean") %>%
  # prefer the ATAC-side symbol, fall back to the GENCODE name from the RNA table
  mutate(SYMBOL = dplyr::coalesce(SYMBOL, rna_symbol)) %>%
  select(-rna_symbol)

message("  integrated: ", nrow(integrated), " genes (ATAC ∩ RNA)")

write_csv(integrated, file.path(outdir, "atac_rna_integrated_d8.csv"))

# Scatter plot: ATAC LFC vs RNA LFC
# Highlight TE and MP signature genes
te_genes <- c("Tbx21", "Zeb2", "Cx3cr1", "Klrg1", "S1pr5", "Gzmb",
              "Id2", "Bhlhe40", "Runx3")
mp_genes <- c("Tcf7", "Id3", "Bach2", "Bcl2", "Il7r", "Sell",
              "Ccr7", "Lef1", "Myb", "Foxo1")

integrated <- integrated %>%
  mutate(signature = case_when(
    SYMBOL %in% te_genes ~ "TE signature",
    SYMBOL %in% mp_genes ~ "MP signature",
    TRUE ~ "Other"
  ))

cor_test <- cor.test(integrated$atac_lfc, integrated$rna_lfc,
                     method = "pearson", use = "complete.obs")

p_scatter <- ggplot(integrated, aes(x = atac_lfc, y = rna_lfc)) +
  geom_point(data = filter(integrated, signature == "Other"),
             color = "grey80", size = 0.5, alpha = 0.5) +
  geom_point(data = filter(integrated, signature == "TE signature"),
             color = pal_subset["TE"], size = 2.5) +
  geom_point(data = filter(integrated, signature == "MP signature"),
             color = pal_subset["MP"], size = 2.5) +
  geom_text_repel(data = filter(integrated, signature != "Other"),
                  aes(label = SYMBOL), size = 3, max.overlaps = 30,
                  fontface = "italic") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.2f\np = %.1e", cor_test$estimate, cor_test$p.value),
           hjust = 1.1, vjust = 1.1, size = 4) +
  scale_color_manual(values = c("TE signature" = pal_subset["TE"],
                                 "MP signature" = pal_subset["MP"],
                                 "Other" = "grey80")) +
  labs(title = "ATAC-seq vs RNA-seq: D8 KO vs WT",
       x = "ATAC log2(FC) (mean per gene)",
       y = "RNA log2(FC)",
       color = "") +
  theme_paper

save_figure(p_scatter, "atac_rna_correlation_d8", width = 8, height = 7, dir = figdir)

# =============================================================================
# 3. chromVAR — motif accessibility deviations
# =============================================================================

# chromVAR is NOT recomputed here — it is built by extended_analysis/scripts/chromvar/build_chromvar.R
# (JASPAR2020, D8 consensus OCRs, WT/Het/KO x TE/EEC/MP) into chromvar_dose.RData.
#
# The previous version gated this entire section on <atac>/diffbind_object.RData,
# a file that has never existed in this project (core/04_cutandrun_analysis.R dropped
# DiffBind because most antibodies are n=1/genotype). So the else branch fired on
# every run: 4 of this script's 9 outputs were silently skipped while the log
# printed a benign-looking "skipping chromVAR". Fail loudly instead.

message("=== chromVAR motif deviation analysis ===")

chromvar_rdata <- file.path(paths$ext_results, "chromvar", "chromvar_dose.RData")
if (!file.exists(chromvar_rdata)) {
  stop("chromvar_dose.RData not found — run extended_analysis/scripts/chromvar/build_chromvar.R first: ",
       chromvar_rdata)
}

suppressPackageStartupMessages(library(chromVAR))

# Load into a private env: a bare load() would clobber palettes/helpers sourced
# from utils.R (known RData-overwrite gotcha).
cv <- new.env()
load(chromvar_rdata, envir = cv)
dev         <- get("dev",         envir = cv)  # chromVARDeviations (746 motifs)
dev_z       <- get("dev_z",       envir = cv)  # 746 x 24 z-scores
sample_meta <- get("sample_meta", envir = cv)  # 24 x 7 (genotype, subset, dose)
stopifnot(identical(colnames(dev_z), sample_meta$sample))

message("  chromVAR: ", nrow(dev_z), " motifs x ", ncol(dev_z), " samples")

# --- top variable motifs ---
# Careful: `dev` is keyed by JASPAR ID (MA0004.1) but build_chromvar.R renamed
# dev_z's rows to TF names (Arnt). computeVariability() inherits dev's IDs, so its
# rownames cannot index dev_z — the `name` column is the bridge.
variability <- computeVariability(dev)
top_motifs <- variability %>%
  as.data.frame() %>%
  rownames_to_column("motif_id") %>%
  arrange(desc(variability)) %>%
  head(30)

stopifnot(all(top_motifs$name %in% rownames(dev_z)))

write_csv(top_motifs, file.path(outdir, "chromvar_top_variable_motifs.csv"))
save(dev, dev_z, sample_meta, variability,
     file = file.path(outdir, "chromvar_results.RData"))

# --- heatmap of top variable motifs ---
# Annotated by Genotype + Subset. The old code used Timepoint, but this matrix is
# D8-only (tp is constant), so subset is the informative second axis.
dev_mat <- dev_z[top_motifs$name, , drop = FALSE]

# as.character() is load-bearing: genotype/subset are FACTORS, and indexing a
# named palette with a factor silently indexes by integer level code, not name —
# pal_subset[unique(subset)] returned Naive/EEC/MP colours for MP/TE/EEC data.
geno_lv <- as.character(unique(sample_meta$genotype))
subs_lv <- as.character(unique(sample_meta$subset))
stopifnot(all(geno_lv %in% names(pal_genotype)), all(subs_lv %in% names(pal_subset)))

ha <- HeatmapAnnotation(
  Genotype = as.character(sample_meta$genotype),
  Subset   = as.character(sample_meta$subset),
  col = list(
    Genotype = pal_genotype[geno_lv],
    Subset   = pal_subset[subs_lv]
  )
)

# col_expression saturates at ±2, but chromVAR deviation z-scores here span
# -23..+28 — using it clips >90% of cells to solid red/blue and the heatmap
# renders as a binary block with no gradient. Rescale the project's diverging
# palette (same colours as col_expression/col_lfc) to this matrix's actual range.
zlim <- as.numeric(quantile(abs(dev_mat), 0.98, na.rm = TRUE))
col_chromvar <- circlize::colorRamp2(
  seq(-zlim, zlim, length.out = 5),
  c("#2166AC", "#4393C3", "#FFFFFF", "#D6604D", "#B2182B")
)

ht_chromvar <- Heatmap(dev_mat,
  name = "Deviation\nz-score",
  col = col_chromvar,
  top_annotation = ha,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "chromVAR: top variable motifs (D8)"
)

save_heatmap(ht_chromvar, "chromvar_top_motifs_heatmap",
             width = 10, height = 10, dir = figdir)

# --- differential motif activity vs WT, per subset ---
# Both KO and Het, using the dose structure (WT=2, Het=1, KO=0 copies) rather
# than a per-timepoint loop.
dev_diff <- lapply(unique(sample_meta$subset), function(sb) {
  lapply(c("KO", "Het"), function(g) {
    wt_idx <- which(sample_meta$genotype == "WT" & sample_meta$subset == sb)
    g_idx  <- which(sample_meta$genotype == g    & sample_meta$subset == sb)
    if (length(wt_idx) < 2 || length(g_idx) < 2) return(NULL)
    diff <- rowMeans(dev_z[, g_idx, drop = FALSE]) -
            rowMeans(dev_z[, wt_idx, drop = FALSE])
    tibble(motif = names(diff), diff_deviation = diff,
           subset = sb, contrast = paste0(g, "_vs_WT"))
  }) %>% compact() %>% bind_rows()
}) %>% compact() %>% bind_rows()

write_csv(dev_diff, file.path(outdir, "chromvar_diff_deviations.csv"))

message("  chromVAR complete: ", nrow(top_motifs), " top variable motifs, ",
        nrow(dev_diff), " motif x contrast rows")

# =============================================================================
# 4. ARID1A binding ↔ accessibility ↔ expression
# =============================================================================
# The file header has always advertised this analysis, but it was never
# implemented: arid1a_res was loaded, row-counted, then referenced only inside an
# exists() guard and never entered a computation. This is the actual three-way link.

message("=== ARID1A binding -> accessibility -> expression ===")

# ATAC consensus OCRs as GRanges (coords live in the annotation table)
atac_gr <- GRanges(
  seqnames = atac_anno$seqnames,
  ranges   = IRanges(atac_anno$start, atac_anno$end),
  peak_id  = atac_anno$peak_id
)

# Which consensus OCRs carry ARID1A in WT?
bound_hits  <- findOverlaps(atac_gr, arid1a_dep)
bound_peaks <- unique(mcols(atac_gr)$peak_id[queryHits(bound_hits)])
message("  ", length(bound_peaks), " / ", length(atac_gr),
        " consensus OCRs overlap a WT ARID1A peak")

# Per-gene ARID1A binding load, keyed the same way as `integrated`
gene_arid1a <- atac_with_gene %>%
  mutate(arid1a_bound = peak_id %in% bound_peaks) %>%
  group_by(ensembl_clean) %>%
  summarize(
    n_arid1a_peaks = sum(arid1a_bound),
    frac_arid1a    = mean(arid1a_bound),
    .groups = "drop"
  )

arid1a_linked <- integrated %>%
  inner_join(gene_arid1a, by = "ensembl_clean") %>%
  mutate(arid1a_target = n_arid1a_peaks > 0)

write_csv(arid1a_linked, file.path(outdir, "arid1a_accessibility_expression.csv"))

# Do ARID1A-bound loci lose more accessibility AND expression in KO?
w_atac <- wilcox.test(atac_lfc ~ arid1a_target, data = arid1a_linked)
w_rna  <- wilcox.test(rna_lfc  ~ arid1a_target, data = arid1a_linked)
message(sprintf("  n ARID1A-target genes: %d / %d",
                sum(arid1a_linked$arid1a_target), nrow(arid1a_linked)))
message(sprintf("  ATAC LFC by ARID1A-target: Wilcoxon p = %.3e", w_atac$p.value))
message(sprintf("  RNA  LFC by ARID1A-target: Wilcoxon p = %.3e", w_rna$p.value))

# --- drift-robust statistics (the ones to quote) ---
# CAVEAT on the mean-LFC metric above: ~18k OCRs collapse in KO, so DESeq2's
# median-of-ratios normalisation redistributes reads and drifts unchanged peaks
# UPWARD — peak-level median LFC is +0.067 while the mean is -0.214, and
# non-significant peaks run 62% positive. Gene-level mean LFC inherits that drift,
# which is why BOTH groups show a positive median atac_lfc despite 18,263 lost vs
# 5,536 gained peaks. With n=11k the Wilcoxon is significant (p~1e-84) on a
# trivial median difference — significance without effect size.
# Counts of SIGNIFICANT peaks are robust to the compositional shift, so the
# lost:gained ratio below is the defensible statistic.
arid1a_counts <- arid1a_linked %>%
  group_by(arid1a_target) %>%
  summarize(genes = n(), lost = sum(n_lost), gained = sum(n_gained),
            peaks = sum(n_peaks), .groups = "drop") %>%
  mutate(lost_gained_ratio = lost / gained,
         # per-peak-surveyed, since ARID1A-bound genes carry more peaks
         # (3.92 vs 2.28) — a confounder for the raw counts
         lost_per_peak   = lost / peaks,
         gained_per_peak = gained / peaks)

write_csv(arid1a_counts, file.path(outdir, "arid1a_target_peak_counts.csv"))

ft <- fisher.test(as.matrix(arid1a_counts[order(arid1a_counts$arid1a_target),
                                          c("lost", "gained")]))
message(sprintf("  lost:gained ratio — not-bound %.2f vs ARID1A-bound %.2f (Fisher p = %.3e, OR = %.2f)",
                arid1a_counts$lost_gained_ratio[!arid1a_counts$arid1a_target],
                arid1a_counts$lost_gained_ratio[arid1a_counts$arid1a_target],
                ft$p.value, ft$estimate))

p_arid1a <- arid1a_linked %>%
  select(SYMBOL, arid1a_target, ATAC = atac_lfc, RNA = rna_lfc) %>%
  pivot_longer(c(ATAC, RNA), names_to = "assay", values_to = "lfc") %>%
  ggplot(aes(x = arid1a_target, y = lfc, fill = arid1a_target)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_violin(scale = "width", alpha = 0.8, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  facet_wrap(~ assay, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "grey75",
                               "TRUE"  = unname(pal_genotype["KO"])),
                    guide = "none") +
  scale_x_discrete(labels = c("FALSE" = "Not bound", "TRUE" = "ARID1A-bound")) +
  labs(title = "ARID1A-bound loci in Arid1a KO: accessibility and expression",
       subtitle = sprintf(paste0("D8 KO vs WT · Wilcoxon ATAC p = %.1e, RNA p = %.1e\n",
                                 "NB: both distributions sit above 0 — DESeq2 renormalisation ",
                                 "drifts unchanged peaks up when ~18k OCRs collapse"),
                          w_atac$p.value, w_rna$p.value),
       x = NULL, y = "log2(FC), KO vs WT") +
  theme_paper

save_figure(p_arid1a, "arid1a_binding_accessibility_expression",
            width = 8, height = 5, dir = figdir)

# Drift-robust companion: significant-peak outcome per peak surveyed
p_arid1a_counts <- arid1a_counts %>%
  select(arid1a_target, Lost = lost_per_peak, Gained = gained_per_peak) %>%
  pivot_longer(c(Lost, Gained), names_to = "direction", values_to = "per_peak") %>%
  ggplot(aes(x = direction, y = per_peak, fill = arid1a_target)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  scale_fill_manual(values = c("FALSE" = "grey75",
                               "TRUE"  = unname(pal_genotype["KO"])),
                    labels = c("FALSE" = "Not bound", "TRUE" = "ARID1A-bound"),
                    name = NULL) +
  labs(title = "ARID1A-bound OCRs preferentially lose accessibility in KO",
       subtitle = sprintf("lost:gained %.2f (ARID1A-bound) vs %.2f (not bound) · Fisher p = %.1e",
                          arid1a_counts$lost_gained_ratio[arid1a_counts$arid1a_target],
                          arid1a_counts$lost_gained_ratio[!arid1a_counts$arid1a_target],
                          ft$p.value),
       x = "Significant DA peaks (padj<0.05, |LFC|>1)",
       y = "Peaks per peak surveyed") +
  theme_paper

save_figure(p_arid1a_counts, "arid1a_target_peak_outcomes",
            width = 7, height = 5, dir = figdir)

# =============================================================================
# 5. Regulatory program summary
# =============================================================================

message("=== Regulatory program: TE vs MP ===")

# Summarize: which TFs lose accessibility AND binding AND expression in KO?
if (exists("integrated")) {
  # TFs with concordant loss across all modalities
  tfs_of_interest <- c("Tbx21", "Batf", "Irf4", "Runx3", "Bhlhe40",
                        "Zeb2", "Id2", "Eomes", "Prdm1", "Zfp683")

  tf_summary <- integrated %>%
    filter(SYMBOL %in% tfs_of_interest) %>%
    select(SYMBOL, atac_lfc, rna_lfc, n_lost, n_gained) %>%
    arrange(rna_lfc)

  write_csv(tf_summary, file.path(outdir, "tf_regulatory_summary.csv"))

  # Dot plot: TF regulatory impact
  p_tf <- ggplot(tf_summary, aes(x = atac_lfc, y = rna_lfc)) +
    geom_point(aes(size = n_lost), color = pal_genotype["KO"]) +
    geom_text_repel(aes(label = SYMBOL), fontface = "italic", size = 4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_size_continuous(range = c(2, 8), name = "# OCRs lost") +
    labs(title = "Effector TF regulatory impact in Arid1a KO",
         x = "Mean ATAC log2(FC) at TF locus",
         y = "RNA log2(FC)") +
    theme_paper

  save_figure(p_tf, "tf_regulatory_impact", width = 8, height = 6, dir = figdir)
}

# =============================================================================
# Save session
# =============================================================================

message("=== Saving R session ===")
save.image(file.path(outdir, "integration_analysis.RData"))
message("Integration analysis complete. Results in: ", outdir)
message("Figures in: ", figdir)
