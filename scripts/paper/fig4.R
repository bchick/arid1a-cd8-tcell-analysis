#!/usr/bin/env Rscript
# =============================================================================
# paper/fig4.R — Figure 4: ARID1A-dependent OCRs in day 8 MP, EEC and TE cells
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
#   4A  ATAC-seq PCA of WT / Het / KO day 8 subsets (+ naive WT)
#   4B  OCRs lost / gained in Arid1a KO and Het vs WT, per subset
#   4C  BED files of subset-specific / shared lost OCRs for the deepTools
#       heatmap (drawn by scripts/paper/fig4_fig5_signal.sh)
#   4D  UpSet plot of OCRs lost in Arid1a KO across subsets
#   4E  ATAC vs RNA log2FC (WT/KO) for OCR-gene pairs, TE/MP signature genes
#
# Thresholds (legend): 2-fold change, adjusted p < 0.05 (BH), for ATAC and RNA.
#
# Inputs:  consensus featureCounts matrix (results/atac/.../consensus/)
#          results/atac/differential/{da_*_vs_WT_D8_*.csv, consensus_peaks_annotated.csv}
#          results/rnaseq/differential/de_KO_vs_WT_D8_*.csv
#          results/rnaseq/star_salmon/salmon.merged.gene_counts.tsv
# Outputs: figures/paper/fig4{a,b,d,e}_*.{pdf,png}; results/paper/fig4*_*.csv
#          results/paper/fig4c_*.bed
# Usage:   Rscript scripts/paper/fig4.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ComplexHeatmap)
  library(ggrepel)
})
select <- dplyr::select
filter <- dplyr::filter

SUBSETS <- c("TE", "EEC", "MP")

consensus_dir <- file.path(paths$atac, "bowtie2/merged_replicate/macs2/narrow_peak/consensus")

# =============================================================================
# 4A. ATAC-seq PCA — day 8 WT / Het / KO subsets plus naive WT
# =============================================================================

message("=== 4A: ATAC PCA ===")

fc <- read.delim(file.path(consensus_dir, "consensus_peaks.mRp.clN.featureCounts.txt"),
                 comment.char = "#", check.names = FALSE)
counts <- as.matrix(fc[, -(1:6)])
rownames(counts) <- fc$Geneid
colnames(counts) <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(counts))

pca_cols <- grep("^(D8_(WT|Het|KO)_(TE|EEC|MP)_|Naive_WT_)", colnames(counts), value = TRUE)
pca_meta <- tibble(sample = pca_cols) |>
  mutate(genotype = ifelse(grepl("^Naive", sample), "WT",
                           sub("^D8_([^_]+)_.*", "\\1", sample)),
         subset   = ifelse(grepl("^Naive", sample), "Naive",
                           sub("^D8_[^_]+_([^_]+)_.*", "\\1", sample)),
         genotype = factor(genotype, levels = c("WT", "Het", "KO")),
         subset   = factor(subset, levels = c("Naive", SUBSETS)))

dds_pca <- DESeqDataSetFromMatrix(counts[, pca_cols],
                                  as.data.frame(pca_meta) |> tibble::column_to_rownames("sample"),
                                  design = ~ 1)
dds_pca <- dds_pca[rowSums(counts(dds_pca) >= 10) >= 2, ]
vsd <- vst(dds_pca, blind = TRUE)
mat <- assay(vsd)
top <- head(order(matrixStats::rowVars(mat), decreasing = TRUE), 5000)
pca <- prcomp(t(mat[top, ]))
pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

pca_df <- pca_meta |> mutate(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
write_panel_table(pca_df, "fig4a_atac_pca")

p4a <- ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = subset)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pal_genotype[c("WT", "Het", "KO")],
                     labels = c("WT", expression(italic(Arid1a)^cHet),
                                expression(italic(Arid1a)^cKO))) +
  scale_shape_manual(values = c(Naive = 1, TE = 18, EEC = 16, MP = 17)) +
  labs(title = "ATAC-seq", color = "Genotype", shape = "Subset",
       x = sprintf("PC1 (%.1f%% var)", pct[1]), y = sprintf("PC2 (%.1f%% var)", pct[2]))
save_panel(p4a, "fig4a_atac_pca", width = 5.5, height = 3.8)

# =============================================================================
# 4B. OCRs lost / gained in KO and Het vs WT per subset
# =============================================================================

message("=== 4B: lost / gained OCR counts ===")

da <- expand_grid(geno = c("KO", "Het"), subset = SUBSETS) |>
  mutate(tab = purrr::map2(geno, subset,
                           ~ read_da(sprintf("%s_vs_WT_D8_%s", .x, .y))))

da_counts <- da |>
  mutate(lost   = purrr::map_int(tab, ~ sum(.x$direction == "lost")),
         gained = purrr::map_int(tab, ~ sum(.x$direction == "gained"))) |>
  select(genotype = geno, subset, lost, gained)
write_panel_table(da_counts, "fig4b_da_counts")
print(da_counts)

p4b <- da_counts |>
  pivot_longer(c(lost, gained), names_to = "direction", values_to = "n") |>
  mutate(value   = ifelse(direction == "lost", -n, n),
         subset  = factor(subset, levels = rev(SUBSETS)),
         genotype = factor(genotype, levels = c("KO", "Het"))) |>
  arrange(genotype) |>                                  # Het drawn on top of KO
  ggplot(aes(x = value, y = subset, fill = genotype)) +
  geom_col(position = "identity", width = 0.7) +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = pal_genotype[c("KO", "Het")],
                    labels = c(KO = expression(italic(Arid1a)^cKO),
                               Het = expression(italic(Arid1a)^cHet))) +
  scale_x_continuous(labels = abs) +
  labs(x = "# OCRs lost  |  # OCRs gained", y = NULL, fill = NULL)
save_panel(p4b, "fig4b_da_counts", width = 5, height = 2.6)

# =============================================================================
# 4C / 4D. Subset-specific vs shared OCRs lost in Arid1a KO
# =============================================================================

message("=== 4C/4D: lost-OCR overlap across subsets ===")

lost_ko <- da |>
  filter(geno == "KO") |>
  mutate(ids = purrr::map(tab, ~ .x$peak_id[.x$direction == "lost"])) |>
  select(subset, ids) |>
  tibble::deframe()

cm <- make_comb_mat(lost_ko[SUBSETS])
upset_df <- tibble(combination = comb_name(cm), size = comb_size(cm),
                   sets = sapply(comb_name(cm), function(code)
                     paste(SUBSETS[strsplit(code, "")[[1]] == "1"], collapse = "+")))
write_panel_table(upset_df, "fig4d_upset_lost_ko")
print(upset_df)

# Subset-specific combinations take the subset colour; shared ones are dark grey
comb_cols <- unname(sapply(comb_name(cm), function(code) {
  bits <- strsplit(code, "")[[1]] == "1"
  if (sum(bits) == 1) pal_subset[SUBSETS[bits]] else "grey20"
}))
ht4d <- UpSet(cm, set_order = SUBSETS, comb_order = order(comb_degree(cm) == 1, -comb_size(cm)),
              comb_col = comb_cols,
              top_annotation = upset_top_annotation(cm, add_numbers = TRUE,
                                                    annotation_name_rot = 90,
                                                    axis_param = list(side = "left"),
                                                    height = unit(4, "cm")),
              left_annotation = upset_left_annotation(cm, gp = gpar(fill = pal_subset[SUBSETS]),
                                                      add_numbers = TRUE),
              right_annotation = NULL,
              column_title = "OCRs lost in Arid1a cKO")
save_panel_heatmap(ht4d, "fig4d_upset_lost_ko", width = 5, height = 3.5)

# 4C groups: lost in one subset only, or lost in all three
saf <- read_tsv(file.path(consensus_dir, "consensus_peaks.mRp.clN.saf"), show_col_types = FALSE)
names(saf)[1:4] <- c("peak_id", "chr", "start", "end")

groups_4c <- list(
  "Down_in_EEC" = setdiff(lost_ko$EEC, union(lost_ko$TE, lost_ko$MP)),
  "Down_in_MP"  = setdiff(lost_ko$MP,  union(lost_ko$TE, lost_ko$EEC)),
  "Down_in_TE"  = setdiff(lost_ko$TE,  union(lost_ko$EEC, lost_ko$MP)),
  "Down_in_all" = Reduce(intersect, lost_ko[SUBSETS])
)
for (g in names(groups_4c)) {
  saf |> filter(peak_id %in% groups_4c[[g]]) |>
    transmute(chr, start = start - 1L, end, peak_id) |>          # SAF is 1-based
    arrange(chr, start) |>
    write_tsv(file.path(paths$paper_tab, paste0("fig4c_", g, ".bed")), col_names = FALSE)
}
write_panel_table(tibble(group = names(groups_4c), n = lengths(groups_4c)), "fig4c_groups")

# =============================================================================
# 4E. ATAC vs RNA log2FC (WT relative to KO) with TE / MP signature genes
# =============================================================================

message("=== 4E: ATAC vs RNA ===")

# TE- and MP-signature genes: DEGs between WT TE and WT MP at day 8
# (2-fold, padj < 0.05) from the same RNA-seq data set.
gc <- read_tsv(file.path(paths$rnaseq, "star_salmon/salmon.merged.gene_counts.tsv"),
               show_col_types = FALSE)
wt_cols <- grep("^D8_WT_(TE|MP)_", names(gc), value = TRUE)
sig_cd <- data.frame(subset = factor(sub("^D8_WT_([^_]+)_.*", "\\1", wt_cols),
                                     levels = c("MP", "TE")), row.names = wt_cols)
dds_sig <- DESeqDataSetFromMatrix(round(as.matrix(gc[, wt_cols])), sig_cd, ~ subset)
rownames(dds_sig) <- gc$gene_id
dds_sig <- DESeq(dds_sig[rowSums(counts(dds_sig) >= 10) >= 2, ], quiet = TRUE)
res_sig <- results(dds_sig, contrast = c("subset", "TE", "MP")) |>
  as.data.frame() |> tibble::rownames_to_column("gene_id") |>
  mutate(ensembl_clean = sub("\\.\\d+$", "", gene_id))
te_sig <- res_sig$ensembl_clean[which(res_sig$padj < PAPER$rna_padj & res_sig$log2FoldChange >=  PAPER$rna_lfc)]
mp_sig <- res_sig$ensembl_clean[which(res_sig$padj < PAPER$rna_padj & res_sig$log2FoldChange <= -PAPER$rna_lfc)]
message(sprintf("  TE-signature genes: %d, MP-signature genes: %d", length(te_sig), length(mp_sig)))
write_panel_table(bind_rows(tibble(signature = "TE", ensembl_clean = te_sig),
                            tibble(signature = "MP", ensembl_clean = mp_sig)) |>
                    left_join(distinct(gc, ensembl_clean = sub("\\.\\d+$", "", gene_id), gene_name),
                              by = "ensembl_clean"),
                  "fig4e_signature_genes")

anno <- read_csv(file.path(paths$atac, "differential/consensus_peaks_annotated.csv"),
                 show_col_types = FALSE) |>
  select(peak_id, ensembl_clean)

label_genes <- c("Bhlhe40", "Klrb1c", "Klrg1", "Tbx21", "Cx3cr1", "Gzma", "Zeb2", "Slamf7",
                 "Tcf7", "Gpr183", "Slamf6", "Cd27", "Pdcd1", "Il7r", "Id3", "Sell", "Ccr7")

e_tabs <- purrr::map(SUBSETS, function(s) {
  rna <- read_csv(file.path(paths$rnaseq, "differential", sprintf("de_KO_vs_WT_D8_%s.csv", s)),
                  show_col_types = FALSE) |>
    mutate(ensembl_clean = sub("\\.\\d+$", "", gene_id)) |>
    filter(!is.na(padj), padj < PAPER$rna_padj, abs(log2FoldChange) >= PAPER$rna_lfc) |>
    select(ensembl_clean, gene_name, rna_lfc = log2FoldChange)
  da |> filter(geno == "KO", subset == s) |> pull(tab) |> purrr::pluck(1) |>
    filter(direction != "ns") |>
    select(peak_id, atac_lfc = log2FoldChange) |>
    inner_join(anno, by = "peak_id") |>
    inner_join(rna, by = "ensembl_clean") |>
    mutate(subset = s,
           # the paper plots WT relative to KO on both axes
           atac_wt_ko = -atac_lfc, rna_wt_ko = -rna_lfc,
           signature = case_when(ensembl_clean %in% te_sig ~ "TE signature gene",
                                 ensembl_clean %in% mp_sig ~ "MP signature gene",
                                 TRUE ~ "Other"))
}) |> bind_rows() |>
  mutate(subset = factor(subset, levels = SUBSETS))
write_panel_table(e_tabs, "fig4e_atac_rna_pairs")

# Quadrant counts: "OCRs (genes)" per signature, as annotated in the paper
quad <- e_tabs |>
  filter(signature != "Other") |>
  mutate(qx = ifelse(atac_wt_ko > 0, "right", "left"), qy = ifelse(rna_wt_ko > 0, "top", "bottom")) |>
  group_by(subset, signature, qx, qy) |>
  summarise(label = sprintf("%d(%d)", n(), n_distinct(ensembl_clean)), .groups = "drop") |>
  mutate(x = ifelse(qx == "right", 9.5, -9.5), hjust = ifelse(qx == "right", 1, 0),
         y = case_when(qy == "top" & signature == "MP signature gene" ~ 9.5,
                       qy == "top" ~ 8.2,
                       signature == "MP signature gene" ~ -8.2,
                       TRUE ~ -9.5))
write_panel_table(select(quad, subset, signature, qx, qy, label), "fig4e_quadrant_counts")

pal_sig <- c("TE signature gene" = "#C0504D", "MP signature gene" = "#2E75B6", Other = "grey75")
clamp <- function(x) pmax(pmin(x, 10), -10)
p4e <- ggplot(e_tabs, aes(clamp(atac_wt_ko), clamp(rna_wt_ko))) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "grey92") +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#E8F3EC") +
  geom_point(data = ~ filter(.x, signature == "Other"), color = pal_sig["Other"], size = 0.5) +
  geom_point(data = ~ filter(.x, signature != "Other"), aes(color = signature), size = 0.8) +
  geom_text_repel(data = ~ filter(.x, gene_name %in% label_genes, signature != "Other") |>
                    distinct(subset, gene_name, .keep_all = TRUE),
                  aes(label = gene_name, color = signature), size = 2.4, fontface = "italic",
                  max.overlaps = 30, show.legend = FALSE) +
  geom_text(data = quad, aes(x = x, y = y, label = label, hjust = hjust, color = signature),
            size = 2.6, show.legend = FALSE, inherit.aes = FALSE) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  facet_wrap(~ subset, ncol = 1) +
  scale_color_manual(values = pal_sig, breaks = names(pal_sig)[1:2]) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
  labs(x = "Log2FC ATAC (WT / Arid1a cKO)", y = "Log2FC RNA (WT / Arid1a cKO)", color = NULL) +
  theme(legend.position = "bottom", strip.text = element_text(color = "black"))
save_panel(p4e, "fig4e_atac_vs_rna", width = 3.6, height = 9)

message("=== fig4.R done ===")
