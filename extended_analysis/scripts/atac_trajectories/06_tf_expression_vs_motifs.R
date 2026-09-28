#!/usr/bin/env Rscript
# =============================================================================
# atac_trajectories/06_tf_expression_vs_motifs.R — does TF mRNA explain the lost accessibility?
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Question: sites that lose their dynamics in ARID1A KO carry particular motifs
# (atac_trajectories/04-05). Is the TF behind each motif itself down-regulated in KO at the
# matching timepoint? If so, lost accessibility could be secondary to lost TF
# expression rather than a direct cBAF requirement at those sites.
#
# Time matching: RNA-seq has KO vs WT at D3 and D8 (TE sorted, like the ATAC
# D8 TE), but no D5 KO libraries and no Naive RNA. So D3 is matched to the
# early classes (Transient, Transient Increasing) and D8 TE to the late ones
# (Sustained, Late Increasing); both timepoints are reported for every class.
#
#   * JASPAR motif -> mouse gene: ALT_ID split on "::" (heterodimers give both
#     genes), matched case-insensitively to the RNA-seq gene symbols.
#   * A TF counts as expressed at a timepoint when its mean normalized count
#     in WT at that timepoint is >= MIN_EXPR.
#   * Motif sets (SEA q < 1e-5): "lost" = enriched in sites losing their
#     dynamics vs sites keeping them (within a WT class); "kept" = the
#     reverse; "program" = enriched in the lost-site class vs static peaks.
#   * Test per class x timepoint: KO vs WT log2FC of the set's TF genes vs
#     every other expressed JASPAR TF gene (two-sided Wilcoxon).
#
# Inputs:  results/extended_analysis/atac_trajectories/motifs/sea_all.tsv.gz (05_lost_site_motifs_plot.R)
#          results/rnaseq/differential/{de_KO_vs_WT_D3,de_KO_vs_WT_D8_TE,normalized_counts}.csv
#          (core/01_rnaseq_analysis.R)
# Outputs: results/extended_analysis/atac_trajectories/tf_expression/
#          figures/extended_analysis/atac_trajectories/tf_expression/
# Usage:   Rscript extended_analysis/scripts/atac_trajectories/06_tf_expression_vs_motifs.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("extended_analysis/scripts/utils_trajectories.R")

Q_MAX    <- 1e-5
MIN_EXPR <- 10

outdir <- file.path(traj_dir, "tf_expression")
figdir <- file.path(paths$ext_figures, "atac_trajectories/tf_expression")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

class_levels <- c("Decreasing", "Transient", "Transient Increasing",
                  "Sustained Increasing", "Late Increasing")
matched_tp <- c(Decreasing = "D8 TE", Transient = "D3", `Transient Increasing` = "D3",
                `Sustained Increasing` = "D8 TE", `Late Increasing` = "D8 TE")

# =============================================================================
# 1. RNA: KO vs WT at D3 and D8 TE, and WT expression over time
# =============================================================================

rna_dir <- file.path(paths$rnaseq, "differential")
de <- bind_rows(
  read_csv(file.path(rna_dir, "de_KO_vs_WT_D3.csv"), show_col_types = FALSE) |> mutate(tp = "D3"),
  read_csv(file.path(rna_dir, "de_KO_vs_WT_D8_TE.csv"), show_col_types = FALSE) |> mutate(tp = "D8 TE")
) |> select(gene_name, tp, log2FoldChange, padj)

nc <- read_csv(file.path(rna_dir, "normalized_counts.csv"), show_col_types = FALSE)
wt_cols <- list(D3 = "^D3_WT_", D5 = "^D5_WT_", `D8 TE` = "^D8_WT_TE_")
wt_expr <- bind_cols(select(nc, gene_name),
                     as_tibble(lapply(wt_cols, function(p) rowMeans(nc[, grep(p, names(nc))])))) |>
  group_by(gene_name) |> summarise(across(everything(), max), .groups = "drop")  # one row per symbol

# =============================================================================
# 2. Motifs -> TF genes
# =============================================================================

sea <- read_tsv(file.path(traj_dir, "motifs/sea_all.tsv.gz"), show_col_types = FALSE)
sym_map <- setNames(wt_expr$gene_name, toupper(wt_expr$gene_name))
motif_genes <- sea |>
  distinct(motif) |>
  mutate(gene = strsplit(sub("\\(VAR\\.[0-9]+\\)$", "", motif), "::", fixed = TRUE)) |>
  tidyr::unnest(gene) |>
  mutate(gene = unname(sym_map[toupper(gene)])) |>
  filter(!is.na(gene))
message(sprintf("  %d of %d JASPAR TF names map to an RNA-seq gene",
                n_distinct(motif_genes$motif), n_distinct(sea$motif)))

sig <- sea |>
  filter(qvalue < Q_MAX) |>
  mutate(set = recode(comparison, lost_vs_kept = "lost", kept_vs_lost = "kept", vs_static = "program")) |>
  group_by(set, class, motif) |> slice_max(enr_ratio, n = 1, with_ties = FALSE) |> ungroup() |>
  inner_join(motif_genes, by = "motif", relationship = "many-to-many") |>
  group_by(set, class, gene) |> summarise(enr_ratio = max(enr_ratio), qvalue = min(qvalue),
                                          motifs = paste(sort(unique(motif)), collapse = ";"),
                                          .groups = "drop")

universe <- motif_genes |> distinct(gene)
tf_rna <- de |> filter(gene_name %in% universe$gene) |>
  left_join(tidyr::pivot_longer(wt_expr, -gene_name, names_to = "tp", values_to = "wt_expr"),
            by = c("gene_name", "tp")) |>
  mutate(expressed = !is.na(wt_expr) & wt_expr >= MIN_EXPR)

# =============================================================================
# 3. Per class x timepoint: are the set's TFs shifted in KO?
# =============================================================================

tests <- tidyr::expand_grid(set = c("lost", "kept", "program"), class = class_levels,
                            tp = c("D3", "D8 TE")) |>
  mutate(res = purrr::pmap(list(set, class, tp), function(s, k, t) {
    genes <- sig$gene[sig$set == s & sig$class == k]
    r <- filter(tf_rna, tp == t, expressed)
    x <- r$log2FoldChange[r$gene_name %in% genes]
    y <- r$log2FoldChange[!r$gene_name %in% genes]
    if (length(x) < 3) return(tibble(n_tf = length(x)))
    tibble(n_tf = length(x), median_lfc = median(x), median_lfc_rest = median(y),
           n_down = sum(r$padj[r$gene_name %in% genes] < 0.05 & x < 0, na.rm = TRUE),
           n_up   = sum(r$padj[r$gene_name %in% genes] < 0.05 & x > 0, na.rm = TRUE),
           p_wilcox = wilcox.test(x, y)$p.value)
  })) |>
  tidyr::unnest(res) |>
  mutate(matched = matched_tp[class] == tp)
write_tsv(tests, file.path(outdir, "tf_set_expression_tests.tsv"))
print(as.data.frame(filter(tests, set != "program") |> mutate(across(where(is.numeric), ~ signif(.x, 3)))))

# Per-TF table for every motif-enriched TF, with its time-matched RNA change
tf_table <- sig |>
  left_join(tf_rna |> select(gene = gene_name, tp, log2FoldChange, padj, wt_expr, expressed),
            by = "gene", relationship = "many-to-many") |>
  mutate(matched = matched_tp[class] == tp) |>
  arrange(set, factor(class, class_levels), desc(enr_ratio), tp)
write_tsv(tf_table, file.path(outdir, "motif_tf_expression.tsv"))

# =============================================================================
# 4. Figures
# =============================================================================

# --- Fig 7: time-matched KO vs WT change of motif TFs, lost vs kept sides -----
d7 <- tf_table |>
  filter(set %in% c("lost", "kept"), matched, expressed) |>
  mutate(class = factor(class, class_levels),
         set = factor(recode(set, lost = "Motif enriched in lost sites",
                             kept = "Motif enriched in kept sites"),
                      c("Motif enriched in lost sites", "Motif enriched in kept sites")))
lab7 <- d7 |> group_by(class, set) |> slice_max(abs(log2FoldChange), n = 3, with_ties = FALSE) |> ungroup()
p7 <- ggplot(d7, aes(set, log2FoldChange, colour = set)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40") +
  geom_jitter(width = 0.15, size = 1.1, alpha = 0.8) +
  ggrepel::geom_text_repel(data = lab7, aes(label = gene), size = 2.3, colour = "black",
                           max.overlaps = 20, seed = 42) +
  facet_wrap(~ class, nrow = 1, labeller = labeller(class = function(k)
    sprintf("%s\nRNA %s", k, matched_tp[k]))) +
  scale_colour_manual(values = c(`Motif enriched in lost sites` = "#2CA02C",
                                 `Motif enriched in kept sites` = "#000000"), guide = "none") +
  scale_x_discrete(labels = c("Lost", "Kept")) +
  labs(x = "Sites whose motif the TF matches", y = "TF mRNA log2FC (KO vs WT)",
       title = "Are the TFs behind lost-site motifs down-regulated in ARID1A KO?",
       subtitle = sprintf("Time-matched RNA-seq; expressed TFs (WT mean normalized count >= %d); SEA q < %g",
                          MIN_EXPR, Q_MAX)) +
  theme_bw(base_size = 9) + theme(strip.background = element_blank())
save_figure(p7, "fig7_motif_tf_expression_lost_vs_kept", width = 12, height = 4, dir = figdir)

# --- Fig 8: key TF families, expression over WT time and KO change -----------
key_tfs <- c("Ets1", "Ets2", "Erg", "Fli1", "Etv1", "Gabpa", "Elf1", "Elf4",
             "Runx1", "Runx3", "Cbfb", "Bcl11b", "Ikzf3",
             "Fos", "Fosb", "Fosl2", "Jun", "Junb", "Jund", "Batf", "Batf3", "Atf3",
             "Tcf7", "Lef1", "Tcf7l2", "Ctcf", "Nrf1", "Tbx21", "Eomes")
d8 <- tf_rna |> filter(gene_name %in% key_tfs) |>
  mutate(gene_name = factor(gene_name, rev(key_tfs)),
         sig = case_when(padj < 0.05 ~ "padj < 0.05", TRUE ~ "n.s."))
p8 <- ggplot(d8, aes(tp, gene_name)) +
  geom_point(aes(size = log10(pmax(wt_expr, 1)), fill = log2FoldChange, colour = sig), shape = 21, stroke = 0.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-2, 2), oob = scales::squish, name = "log2FC\nKO vs WT") +
  scale_colour_manual(values = c(`padj < 0.05` = "black", n.s. = "grey80"), name = NULL) +
  scale_size_area(max_size = 5, name = "log10 WT\nexpression") +
  labs(x = "RNA-seq timepoint", y = NULL,
       title = "Motif-family TFs: KO vs WT mRNA",
       subtitle = "ETS, RUNX, AP-1/BATF, TCF/LEF and controls") +
  theme_bw(base_size = 9)
save_figure(p8, "fig8_key_tf_expression", width = 4.6, height = 6.5, dir = figdir)

message("=== atac_trajectories/06_tf_expression_vs_motifs.R complete ===")
