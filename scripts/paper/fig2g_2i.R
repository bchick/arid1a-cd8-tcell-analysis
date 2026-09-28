#!/usr/bin/env Rscript
# =============================================================================
# paper/fig2g_2i.R — Figure 2G-I: day 3 WT vs Arid1a cKO RNA-seq
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
#   2G  Number of DEGs (>2-fold, padj < 0.05)
#   2H  Heatmap of curated, biologically relevant genes (gene list from Fig 2H)
#   2I  Top Hallmark gene sets (5 down, 4 up) by GSEA (fgsea, 10,000 permutations, padj < 0.01)
#
# Inputs:  results/rnaseq/star_salmon/salmon.merged.gene_counts.tsv
#          results/rnaseq/differential/de_KO_vs_WT_D3.csv (core/01_rnaseq_analysis.R)
# Outputs: figures/paper/fig2{g,h,i}_*.{pdf,png}; results/paper/fig2{g,h,i}_*.csv
# Usage:   Rscript scripts/paper/fig2g_2i.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")
source("scripts/paper/rna_helpers.R")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)
})

de_d3 <- read_de("KO_vs_WT_D3")

# =============================================================================
# 2G — DEG counts
# =============================================================================

message("=== Fig 2G: D3 DEG counts ===")

deg_counts <- count_degs(de_d3, "KO_vs_WT_D3")
write_panel_table(deg_counts, "fig2g_d3_deg_counts")
print(deg_counts)

n_wt <- deg_counts$down[deg_counts$threshold == "legend_2fold"]   # higher in WT
n_ko <- deg_counts$up[deg_counts$threshold == "legend_2fold"]     # higher in KO

g_df <- tibble(
  side  = factor(c("Genes higher in WT", "Genes higher in Arid1a cKO"),
                 levels = c("Genes higher in WT", "Genes higher in Arid1a cKO")),
  n     = c(-n_wt, n_ko),
  genotype = c("WT", "KO")
)
lim <- ceiling(max(abs(g_df$n)) / 250) * 250 + 250

p_2g <- ggplot(g_df, aes(x = n, y = 1, fill = genotype)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_text(aes(label = abs(n), hjust = ifelse(n < 0, 1.15, -0.15)), size = 3) +
  scale_fill_manual(values = pal_genotype[c("WT", "KO")], guide = "none") +
  scale_x_continuous(limits = c(-lim, lim), labels = abs,
                     breaks = seq(-lim + 250, lim - 250, by = 250)) +
  labs(x = "Genes higher in WT          Genes higher in Arid1a cKO", y = NULL,
       title = "Day 3 DEGs (>2-fold, padj < 0.05)") +
  theme_paper +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), plot.title = element_text(size = 10))

save_panel(p_2g, "fig2g_d3_deg_counts", width = 4.5, height = 1.8)

# =============================================================================
# 2H — Curated gene heatmap
# =============================================================================

message("=== Fig 2H: D3 curated gene heatmap ===")

fig2h_genes <- tibble(
  gene_name = c(
    # Transcription factors
    "Bhlhe40", "Eomes", "Ar", "Rxra", "Atf3", "Zeb2", "Id2", "Tcf3",
    "Tbx21", "Batf", "Batf3", "Runx3", "Irf4", "Maf", "Myb", "Tcf7",
    "Klf3", "Foxn3", "Zfp395", "Bhlhe41", "Tox2", "Irf6", "Prdm1",
    "Runx1", "Chd3", "Tnfaip3",
    # Cell cycle
    "Top2a", "Cdkn2d", "Cdkn1b", "Cdkn2a", "Mki67",
    # Cytokine receptors
    "Tnfrsf14", "Il12rb2", "Il2rb", "Il7r", "Il18r1",
    # Chemokine receptors / migration
    "Cx3cr1", "Cxcr3", "S1pr1", "Ccr7", "Ccr5", "Ccr2", "Ccr9",
    # Integrins / adhesion
    "Itga1", "Itga2", "Itgax", "Itgam", "Itga4", "Itga6", "Cd69",
    # Effector molecules
    "Gzma", "Gzmb", "Gzmk", "Gzmm", "Ccl3", "Ccl5",
    # Metabolism / growth
    "Myc", "Srm", "Dusp2", "Tfrc", "Slc7a5"
  ),
  group = rep(c("TFs", "Cell cycle", "Receptors", "Chemokine R",
                "Integrins", "Effector", "Metabolism"),
              times = c(26, 5, 5, 7, 7, 6, 5))
)

rna <- load_rna_vst()
d3 <- rna$meta |>
  filter(timepoint == "D3") |>
  arrange(match(genotype, c("WT", "KO")), replicate)

gm <- rna$gene_map |> distinct(gene_name, .keep_all = TRUE)
h_genes <- fig2h_genes |>
  left_join(gm, by = "gene_name") |>
  filter(!is.na(gene_id), gene_id %in% rownames(rna$vst))
missing <- setdiff(fig2h_genes$gene_name, h_genes$gene_name)
if (length(missing)) message("  Not quantified: ", paste(missing, collapse = ", "))
message(sprintf("  %d / %d genes", nrow(h_genes), nrow(fig2h_genes)))

mat <- t(rna$vst[h_genes$gene_id, d3$sample_name, drop = FALSE])
colnames(mat) <- h_genes$gene_name
mat_z <- pmin(pmax(scale(mat), -2), 2)

sig_genes <- de_d3 |>
  filter(padj < PAPER$rna_padj, abs(log2FoldChange) >= PAPER$rna_lfc) |>
  pull(gene_name)
h_genes <- h_genes |>
  mutate(DEG = ifelse(gene_name %in% sig_genes, "yes", "no"))

write_panel_table(
  h_genes |>
    left_join(de_d3 |> select(gene_id, log2FoldChange, padj), by = "gene_id"),
  "fig2h_d3_heatmap_genes")

row_split <- factor(as.character(d3$genotype), levels = c("WT", "KO"))
col_split <- factor(h_genes$group, levels = unique(fig2h_genes$group))

ht_2h <- Heatmap(
  mat_z,
  name = "Expression\nz-score",
  col = colorRamp2(seq(-2, 2, length.out = 9), mako(9)),
  top_annotation = HeatmapAnnotation(
    `DEG (>2-fold)` = h_genes$DEG,
    col = list(`DEG (>2-fold)` = c(yes = "#3B1F5F", no = "grey80")),
    simple_anno_size = unit(3, "mm"), border = TRUE,
    annotation_name_gp = gpar(fontsize = 8)),
  left_annotation = rowAnnotation(
    Genotype = anno_block(
      gp = gpar(fill = pal_genotype[levels(row_split)], col = "black"),
      labels = c("WT", "Arid1a cKO"),
      labels_gp = gpar(col = "white", fontsize = 8, fontface = "bold"),
      width = unit(6, "mm"))),
  row_split = row_split, row_title = NULL, cluster_rows = FALSE,
  column_split = col_split, cluster_columns = FALSE, cluster_column_slices = FALSE,
  column_title_gp = gpar(fontsize = 7), column_gap = unit(1, "mm"),
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 7, fontface = "italic"),
  border = TRUE,
  heatmap_legend_param = list(at = c(-2, 0, 2))
)

save_panel_heatmap(ht_2h, "fig2h_d3_gene_heatmap", width = 12, height = 3)

# =============================================================================
# 2I — Hallmark GSEA
# =============================================================================

message("=== Fig 2I: D3 Hallmark GSEA ===")

gsea_d3 <- run_fgsea(de_d3, msig_sets("H")) |>
  mutate(pathway = hallmark_label(pathway)) |>
  arrange(padj)

write_panel_table(gsea_d3 |> select(-leadingEdge), "fig2i_d3_hallmark_gsea_all")

# The legend says "Top 8 Hallmark gene sets (adjusted p value < 0.01)" but the
# published panel shows 9 bars: the 5 most negative and 4 most positive NES.
sig_sets <- gsea_d3 |> filter(padj < 0.01)
top_sets <- bind_rows(
  sig_sets |> filter(NES < 0) |> slice_min(NES, n = 5),
  sig_sets |> filter(NES > 0) |> slice_max(NES, n = 4)) |>
  arrange(NES) |>
  mutate(pathway = factor(pathway, levels = pathway),
         log_padj = -log10(padj))

write_panel_table(top_sets |> select(-leadingEdge), "fig2i_d3_hallmark_top9")
print(top_sets |> select(pathway, NES, padj))

p_2i <- ggplot(top_sets, aes(x = NES, y = pathway, fill = log_padj)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "-log10(padj)") +
  labs(x = "NES  (WT ← → Arid1a cKO)", y = NULL, title = "GSEA (Hallmark)") +
  theme_paper +
  theme(axis.text.y = element_text(size = 7), plot.title = element_text(size = 10))

save_panel(p_2i, "fig2i_d3_hallmark_gsea", width = 5, height = 3)

message("Fig 2G-I complete.")
