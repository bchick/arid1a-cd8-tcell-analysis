#!/usr/bin/env Rscript
# =============================================================================
# paper/fig3.R — Figure 3E-I: day 8 TE / EEC / MP RNA-seq, WT vs Arid1a cHet vs cKO
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
#   3E  PCA of all day 8 samples (500 most variable genes)
#   3F  Number of DEGs per subset (cHet and cKO vs WT; >2-fold, padj < 0.05)
#   3G  Heatmap of selected DEGs (gene list and gene clusters from Fig 3G)
#   3H  Volcano plots, cKO vs WT, with TE- and MP-signature genes highlighted
#   3I  GSEA of cKO vs WT against published CD8+ T cell signatures
#
# Signatures (3H, 3I): the paper's MP-vs-TE and activated-vs-naive sets come
# from GSE10239 (Sarkar et al. 2008, KLRG1-high/-int effector, memory and naive
# CD8+ T cells), which is in MSigDB C7/ImmuneSigDB. The legend cites
# "GSE10739", an unrelated human monocyte series, so this is taken as a typo.
# The legend does not say which naive-vs-activated contrast was used; the
# day 4.5 effector contrast is used here, and GSEA against every GSE10239
# contrast is written to results/paper/fig3i_d8_gsea_gse10239_all.csv.
# The d2.5 TbetKO (PRJNA547650), Batf3OE (GSE143504), Runx3KO (GSE81888) and
# d8 BRD4KO (GSE173515) signatures are derived from other raw datasets and are
# not reproduced here; a Hallmark version of 3I is written alongside.
#
# Inputs:  results/rnaseq/star_salmon/salmon.merged.gene_counts.tsv
#          results/rnaseq/differential/de_{KO,Het}_vs_WT_D8_{TE,EEC,MP}.csv
# Outputs: figures/paper/fig3{e..i}_*.{pdf,png}; results/paper/fig3{e..i}_*.csv
# Usage:   Rscript scripts/paper/fig3.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")
source("scripts/paper/rna_helpers.R")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)
  library(ggrepel)
  library(patchwork)
})

subsets   <- c("TE", "EEC", "MP")
genotypes <- c("WT", "Het", "KO")
geno_labels <- c(WT = "WT", Het = "Arid1a cHet", KO = "Arid1a cKO")

de <- list()
for (s in subsets) for (g in c("KO", "Het")) {
  cn <- paste0(g, "_vs_WT_D8_", s)
  de[[cn]] <- read_de(cn)
}

rna <- load_rna_vst()
d8 <- rna$meta |>
  filter(timepoint == "D8") |>
  mutate(cell_subset = factor(as.character(cell_subset), levels = subsets),
         genotype = factor(as.character(genotype), levels = genotypes)) |>
  arrange(cell_subset, genotype, replicate)
vst_d8 <- rna$vst[, d8$sample_name]

# =============================================================================
# 3E — PCA
# =============================================================================

message("=== Fig 3E: day 8 PCA ===")

# 500 most variable genes, as DESeq2::plotPCA (ntop = 500)
top_var <- head(order(matrixStats::rowVars(vst_d8), decreasing = TRUE), 500)
pca <- prcomp(t(vst_d8[top_var, ]), center = TRUE, scale. = FALSE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
pca_df <- d8 |>
  select(sample_name, genotype, cell_subset, replicate) |>
  mutate(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
write_panel_table(pca_df, "fig3e_d8_pca")

p_3e <- ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = cell_subset)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pal_genotype[genotypes], labels = geno_labels,
                     name = "Genotype") +
  scale_shape_manual(values = c(TE = 18, EEC = 16, MP = 17), name = "Subset") +
  labs(x = sprintf("PC1 (%s%% var)", var_pct[1]),
       y = sprintf("PC2 (%s%% var)", var_pct[2]), title = "RNA-seq") +
  theme_paper +
  theme(plot.title = element_text(size = 10))

save_panel(p_3e, "fig3e_d8_pca", width = 5, height = 3.2)

# =============================================================================
# 3F — DEG counts per subset
# =============================================================================

message("=== Fig 3F: day 8 DEG counts ===")

deg_counts <- purrr::imap_dfr(de, count_degs) |>
  tidyr::separate(contrast, c("genotype", NA, NA, NA, "subset"), sep = "_",
                  remove = FALSE)
write_panel_table(deg_counts, "fig3f_d8_deg_counts")
print(deg_counts |> filter(threshold == "legend_2fold") |>
        select(contrast, up, down))

f_df <- deg_counts |>
  filter(threshold == "legend_2fold") |>
  tidyr::pivot_longer(c(up, down), names_to = "direction", values_to = "n") |>
  mutate(n = ifelse(direction == "down", -n, n),
         subset = factor(subset, levels = rev(subsets)),
         genotype = factor(genotype, levels = c("KO", "Het")))

# Het bars are drawn on top of (narrower than) the KO bars, as in the paper
p_3f <- ggplot(f_df, aes(x = n, y = subset, fill = genotype)) +
  geom_col(data = filter(f_df, genotype == "KO"), width = 0.7,
           color = "black", linewidth = 0.3) +
  geom_col(data = filter(f_df, genotype == "Het"), width = 0.4,
           color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = pal_genotype[c("Het", "KO")],
                    labels = geno_labels[c("Het", "KO")], name = NULL) +
  scale_x_continuous(labels = abs) +
  labs(x = "Genes downregulated     Genes upregulated", y = NULL) +
  theme_paper

save_panel(p_3f, "fig3f_d8_deg_counts", width = 4.5, height = 2.4)

# =============================================================================
# 3G — Selected-DEG heatmap
# =============================================================================

message("=== Fig 3G: day 8 gene heatmap ===")

fig3g_genes <- tibble(
  gene_name = c(
    "Klrc2", "Klra1", "Klra4",
    "Id2", "Klrd1", "Klrb1c",
    "Itga4", "Klrc1", "Aqp9", "Klrk1", "Itga1", "Itgax",
    "Tcf3", "Eomes", "S1pr1", "Runx3", "Zfp683", "Tbx21",
    "Ccl5", "Slamf7", "Bhlhe40", "Itgam", "Cx3cr1",
    "Il2rb", "Gzmm", "Ar",
    "Top2a", "Mki67", "Ccl3", "Prdm1", "Gzmk", "Gzmb", "Cdkn2a", "Ccr5",
    "Runx1", "Klra3", "Klra7", "Klf3", "Klre1", "Zeb2", "S1pr5", "Gzma", "Klrg1",
    "Itga6", "Myb", "Trib2", "Klrb1f", "Gzmc", "Irf4", "Maf",
    "Il2ra", "Batf3", "Il21", "Ccr9", "Tox2", "Bhlhe41",
    "Cxcr3", "Tcf7", "Ccr7", "Slamf6", "Id3", "Il7r"),
  cluster = rep(1:9, times = c(3, 3, 6, 11, 3, 8, 9, 13, 6))
)

gm <- rna$gene_map |> distinct(gene_name, .keep_all = TRUE)
g_genes <- fig3g_genes |>
  left_join(gm, by = "gene_name") |>
  filter(!is.na(gene_id), gene_id %in% rownames(vst_d8))
missing <- setdiff(fig3g_genes$gene_name, g_genes$gene_name)
if (length(missing)) message("  Not quantified: ", paste(missing, collapse = ", "))
message(sprintf("  %d / %d genes", nrow(g_genes), nrow(fig3g_genes)))
write_panel_table(g_genes, "fig3g_d8_heatmap_genes")

mat <- vst_d8[g_genes$gene_id, , drop = FALSE]
rownames(mat) <- g_genes$gene_name
mat_z <- pmin(pmax(t(scale(t(mat))), -2), 2)

slice_lv <- as.vector(t(outer(subsets, genotypes, paste, sep = "_")))
col_split <- factor(paste(d8$cell_subset, d8$genotype, sep = "_"), levels = slice_lv)

ht_3g <- Heatmap(
  mat_z,
  name = "Expression\nz-score",
  col = colorRamp2(seq(-2, 2, length.out = 9), mako(9)),
  top_annotation = HeatmapAnnotation(
    Subset = anno_block(
      gp = gpar(fill = rep(unname(pal_subset[subsets]), each = 3), col = NA),
      labels = c("", "TE", "", "", "EEC", "", "", "MP", ""),
      labels_gp = gpar(col = "white", fontsize = 9, fontface = "bold"),
      height = unit(5, "mm")),
    Genotype = anno_block(
      gp = gpar(fill = rep(unname(pal_genotype[genotypes]), 3), col = "black"),
      labels = rep(c("WT", "cHet", "cKO"), 3),
      labels_gp = gpar(col = "white", fontsize = 6, fontface = "bold"),
      height = unit(4, "mm")),
    show_annotation_name = FALSE, gap = unit(0.5, "mm")),
  column_split = col_split, cluster_columns = FALSE, cluster_column_slices = FALSE,
  column_gap = unit(c(0, 0, 2, 0, 0, 2, 0, 0), "mm"), column_title = NULL,
  show_column_names = FALSE,
  row_split = factor(g_genes$cluster, levels = 1:9), cluster_rows = FALSE,
  row_title = NULL, row_gap = unit(1.2, "mm"),
  row_names_gp = gpar(fontsize = 6.5, fontface = "italic"),
  border = TRUE,
  heatmap_legend_param = list(at = c(-2, 0, 2), direction = "horizontal")
)

save_panel_heatmap(ht_3g, "fig3g_d8_gene_heatmap", width = 5.5, height = 9)

# =============================================================================
# Signatures from GSE10239 (MSigDB C7 / ImmuneSigDB)
# =============================================================================

c7 <- msig_sets("C7", "IMMUNESIGDB", pattern = "^GSE10239_")
paper_sets <- list(
  "MPvsTE.Up"              = c7[["GSE10239_KLRG1INT_VS_KLRG1HIGH_EFF_CD8_TCELL_UP"]],
  "MPvsTE.Down"            = c7[["GSE10239_KLRG1INT_VS_KLRG1HIGH_EFF_CD8_TCELL_DN"]],
  "Activated vs Naive.Up"  = c7[["GSE10239_NAIVE_VS_DAY4.5_EFF_CD8_TCELL_DN"]],
  "Activated vs Naive.Down"= c7[["GSE10239_NAIVE_VS_DAY4.5_EFF_CD8_TCELL_UP"]]
)
stopifnot(all(lengths(paper_sets) > 0))
# 3H highlights: GSE10239 KLRG1-high (TE) / KLRG1-int (MP) genes, plus the genes
# the published panels label in red (TE) and blue (MP)
te_labelled <- c("Bhlhe40", "Klrb1c", "Slamf7", "Klrg1", "Havcr2", "Cx3cr1", "Gzma",
                 "Zeb2", "Tbx21", "Prdm1", "Birc5")
mp_labelled <- c("Tcf7", "Sell", "Ltb", "Il7r", "Cxcr3", "Cd27", "Gpr183", "Slamf6",
                 "Id3", "Cxcr5", "Ccr7", "P2rx7", "Pdcd1")
te_sig <- union(te_labelled, setdiff(paper_sets[["MPvsTE.Down"]], mp_labelled))
mp_sig <- union(mp_labelled, setdiff(paper_sets[["MPvsTE.Up"]], te_labelled))

# =============================================================================
# 3H — Volcano plots with TE / MP signature genes
# =============================================================================

message("=== Fig 3H: day 8 volcano plots ===")

# Genes labelled in the published panels
h_labels <- list(
  TE  = c("Bhlhe40", "Klrb1c", "Slamf7", "Klrg1", "Havcr2", "Cx3cr1", "Gzma", "Pdcd1",
          "Tcf7", "Sell", "Ltb", "Il7r", "Cxcr3", "Cd27", "Gpr183", "Slamf6"),
  EEC = c("Cx3cr1", "Bhlhe40", "Zeb2", "Gzma", "Havcr2", "Klrb1c", "Slamf7", "Klrg1",
          "Pdcd1", "Tcf7", "Slamf6", "Ltb", "Id3", "Il7r", "Gpr183", "Cxcr5"),
  MP  = c("Bhlhe40", "Zeb2", "Tbx21", "Cx3cr1", "Slamf7", "Gzma", "P2rx7", "Pdcd1",
          "Prdm1", "Birc5", "Slamf6", "Ccr7")
)
y_cap <- 50   # the published panels cap -log10(padj) at 50

volcano_df <- purrr::map_dfr(subsets, function(s) {
  de[[paste0("KO_vs_WT_D8_", s)]] |>
    filter(!is.na(padj)) |>
    mutate(subset = s,
           signature = case_when(gene_name %in% te_sig ~ "TE signature gene",
                                 gene_name %in% mp_sig ~ "MP signature gene",
                                 TRUE ~ "other"),
           y = pmin(-log10(padj), y_cap),
           label = ifelse(gene_name %in% h_labels[[s]], gene_name, NA))
})
write_panel_table(volcano_df |> filter(signature != "other") |>
                    select(subset, gene_name, log2FoldChange, padj, signature),
                  "fig3h_d8_volcano_signature_genes")

sig_cols <- c("TE signature gene" = unname(pal_subset["TE"]),
              "MP signature gene" = unname(pal_subset["MP"]),
              "other" = "grey80")

make_volcano <- function(s) {
  df <- filter(volcano_df, subset == s)
  ggplot(df, aes(log2FoldChange, y)) +
    geom_point(data = filter(df, signature == "other"), color = "grey80", size = 0.6) +
    geom_point(data = filter(df, signature != "other"), aes(color = signature), size = 1) +
    geom_text_repel(aes(label = label, color = signature), size = 2.4,
                    fontface = "italic", max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.2, na.rm = TRUE, show.legend = FALSE) +
    scale_color_manual(values = sig_cols, breaks = names(sig_cols)[1:2], name = NULL) +
    scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2),
                       oob = scales::squish) +
    scale_y_continuous(limits = c(0, y_cap)) +
    labs(x = "log2FC RNA (WT ← → Arid1a cKO)", y = "-log10(padj)", title = s) +
    theme_paper +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                    color = pal_subset[[s]]))
}

p_3h <- wrap_plots(lapply(subsets, make_volcano), ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")
save_panel(p_3h, "fig3h_d8_volcanos", width = 4, height = 10)

# =============================================================================
# 3I — GSEA against published signatures (+ Hallmark companion)
# =============================================================================

message("=== Fig 3I: day 8 GSEA ===")

gsea_sig <- purrr::map_dfr(subsets, function(s) {
  run_fgsea(de[[paste0("KO_vs_WT_D8_", s)]], paper_sets) |> mutate(subset = s)
})
gsea_h <- purrr::map_dfr(subsets, function(s) {
  run_fgsea(de[[paste0("KO_vs_WT_D8_", s)]], msig_sets("H")) |>
    mutate(subset = s, pathway = hallmark_label(pathway))
})
write_panel_table(gsea_sig |> select(-leadingEdge), "fig3i_d8_gsea_gse10239")
gsea_all_gse10239 <- purrr::map_dfr(subsets, function(s) {
  run_fgsea(de[[paste0("KO_vs_WT_D8_", s)]], c7) |> mutate(subset = s)
})
write_panel_table(gsea_all_gse10239 |> select(-leadingEdge), "fig3i_d8_gsea_gse10239_all")
write_panel_table(gsea_h |> select(-leadingEdge), "fig3i_d8_gsea_hallmark")
print(gsea_sig |> select(subset, pathway, NES, padj))

gsea_dot <- function(df, levels) {
  df |>
    mutate(subset = factor(subset, levels = subsets),
           pathway = factor(pathway, levels = rev(levels)),
           sig = cut(padj, c(-Inf, 0.005, 0.05, Inf),
                     labels = c("padj<0.005", "padj<0.05", "ns"))) |>
    ggplot(aes(subset, pathway)) +
    geom_point(aes(fill = NES, size = sig), shape = 21, stroke = 0.3) +
    scale_fill_gradientn(colors = mako(9), limits = c(-3.2, 3.2),
                         name = "NES\n(Arid1a cKO/WT)") +
    scale_size_manual(values = c("padj<0.005" = 6, "padj<0.05" = 4, ns = 1.2),
                      drop = FALSE, name = NULL) +
    labs(x = NULL, y = NULL) +
    theme_paper +
    theme(axis.text.x = element_text(face = "bold"))
}

p_3i <- gsea_dot(gsea_sig, c("MPvsTE.Up", "Activated vs Naive.Down",
                             "MPvsTE.Down", "Activated vs Naive.Up"))
save_panel(p_3i, "fig3i_d8_gsea_gse10239", width = 4.2, height = 2.6)

top_h <- gsea_h |>
  filter(padj < 0.05) |>
  group_by(pathway) |>
  summarise(m = mean(abs(NES)), .groups = "drop") |>
  slice_max(m, n = 20) |>
  pull(pathway)
p_3i_h <- gsea_dot(filter(gsea_h, pathway %in% top_h), top_h)
save_panel(p_3i_h, "fig3i_d8_gsea_hallmark", width = 5, height = 5.5)

message("Fig 3E-I complete.")
