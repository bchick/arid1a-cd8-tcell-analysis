#!/usr/bin/env Rscript
# =============================================================================
# paper/fig1_fig2_panels.R — Figure 1B/1C/1D/1F, S1B/S1C and Figure 2C/2E/2F
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# Clusters are the published Fig 1A regions (data/metadata/paper_ocr_clusters).
# Region-based panels (1B-1D, 1F, S1C, 2E, 2F) use the regions directly;
# DA-based panels (2C, S1B) map each region to its overlapping nf-core
# consensus peak (consensus_peak_primary; >= 99% mapped) and use our DA tables.
#
# Run after:
#   Rscript scripts/paper/fig1_published_clusters.R  (regions + consensus mapping)
#   bash    scripts/paper/fig1_fig2_deeptools.sh     (profile tables; 1A/2D heatmaps)
#   bash    scripts/paper/fig1d_arid1a_overlap.sh    (1D overlap)
#   bash    scripts/paper/fig1_homer.sh              (1C annotation, 1F/S1C motifs)
# Panels whose inputs are absent are skipped with a message.
#
#   1B  ATAC coverage profiles per cluster (WT Naive/48h/d3/d5/d8). 48h ATAC
#       (merged-replicate bigWig) is shown in the 1A/1B signal panels only, as in
#       the paper: QC Rep1 5.0M read pairs FRiP 0.314, Rep2 4.5M FRiP 0.115,
#       merged FRiP 0.244; excluded from clustering and all statistical tests.
#   1C  genomic annotation per cluster (HOMER annotatePeaks.pl, as in the paper)
#   1D  ARID1A CUT&RUN obs/exp log ratio per cluster (HOMER mergePeaks -matrix)
#   1F  HOMER known-motif enrichment (-log10 p) for the paper's motif panel
#   S1B WT expression of genes annotated to each cluster (d3/d5/d8; the paper's
#       naive RNA-seq is external, GSE152841, and is not included)
#   S1C top 15 known motifs per cluster
#   2C  % of each cluster lost in Arid1a KO (2-fold, FDR < 0.05) at d3/d5/d8
#   2E  ATAC profiles WT vs KO (and d8 Het) per cluster
#   2F  H3K27ac CUT&RUN d5 WT vs KO centred on ATAC peaks
#   1E, 1G use external histone/TF ChIP-seq and are not reproduced.
#
# Inputs:  results/paper/{fig1a_ocr_clusters.csv, deeptools/, homer/}
#          results/atac/differential/{da_*.csv, consensus_peaks_annotated.csv}
#          results/rnaseq/differential/normalized_counts.csv
# Outputs: figures/paper/fig1{b,c,d,f}_*, figS1{b,c}_*, fig2{c,e,f}_*.{pdf,png}
#          results/paper/<panel>.csv
# Usage:   Rscript scripts/paper/fig1_fig2_panels.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")

suppressPackageStartupMessages({
  library(patchwork)
})

clusters <- read_ocr_clusters()
MAT   <- file.path(paths$paper_tab, "deeptools")
HOMER <- file.path(paths$paper_tab, "homer")

pal_cluster <- setNames(viridisLite::mako(6, begin = 0.1, end = 0.85)[1:5], OCR_CLUSTERS)
pal_tp <- c(Naive = "grey60", `48h` = "#E8B04A", d3 = "#D17A22", d5 = "#8B1A1A", d8 = "black")
pal_geno <- c(WT = "black", Het = unname(pal_genotype["Het"]), KO = unname(pal_genotype["KO"]))

theme_panel <- theme_paper + theme(text = element_text(family = "Helvetica"))

#' Parse a deepTools plotProfile --outFileNameData table into long format
read_profile_tab <- function(f) {
  # Rows: sample, region, one value per 20-bp bin (+/-1 kb). The header "bins"
  # row repeats across samples, so bins are taken from each data row.
  body <- readLines(f)[-(1:2)]
  do.call(rbind, lapply(body, function(l) {
    v <- strsplit(l, "\t")[[1]]
    vals <- as.numeric(v[-(1:2)])
    vals <- vals[!is.na(vals)]
    tibble(sample = v[1], region = v[2], bin = seq_along(vals), value = vals)
  })) |>
    mutate(pos_kb = (bin - max(bin) / 2 - 0.5) * 20 / 1000,
           region = factor(region, levels = OCR_CLUSTERS))
}

# =============================================================================
# 1B — ATAC profiles per cluster
# =============================================================================
f <- file.path(MAT, "fig1a_profile.tab")
if (require_inputs(f, "Fig 1B profile table")) {
  pr <- read_profile_tab(f) |>
    filter(grepl("^ATAC_", sample)) |>
    mutate(tp = factor(recode(sub("ATAC_", "", sample), Naive = "Naive", `48h` = "48h",
                              D3 = "d3", D5 = "d5", D8 = "d8"), levels = names(pal_tp)))
  write_panel_table(pr, "fig1b_atac_profiles")
  p <- ggplot(pr, aes(pos_kb, value, colour = tp)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ region, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = pal_tp, name = NULL) +
    labs(x = "Distance from peak center (kb)", y = "ATAC coverage",
         caption = "48h: merged FRiP 0.244 (signal only)") +
    theme_panel
  save_panel(p, "fig1b_atac_profiles", width = 3.2, height = 8)
}

# =============================================================================
# 1C — genomic annotation per cluster
# =============================================================================
annot <- read_csv(file.path(paths$atac, "differential/consensus_peaks_annotated.csv"),
                  show_col_types = FALSE) |>
  dplyr::select(peak_id, annotation, SYMBOL, distanceToTSS) |>
  mutate(category = case_when(
    grepl("^Promoter", annotation) ~ "Promoter",
    grepl("^Distal Intergenic", annotation) ~ "Intergenic",
    grepl("^Intron", annotation) ~ "Intron",
    grepl("^Exon", annotation) ~ "Exon",
    TRUE ~ "Other"),
    category = factor(category, levels = c("Promoter", "Intergenic", "Intron", "Exon", "Other")))

read_homer_annot <- function(cl) {
  f <- file.path(HOMER, "annotation", paste0(gsub(" ", "_", cl), ".annotatePeaks.txt"))
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  read.delim(f, check.names = FALSE) |> as_tibble() |>
    transmute(cluster = cl, homer = sub(" \\(.*", "", Annotation))
}
homer_ann <- bind_rows(lapply(OCR_CLUSTERS, read_homer_annot))
if (nrow(homer_ann) && n_distinct(homer_ann$cluster) == length(OCR_CLUSTERS)) {
  ann_src <- "HOMER annotatePeaks.pl (GENCODE vM35)"
  ann_cl <- homer_ann |>
    mutate(category = case_when(
      homer == "promoter-TSS" ~ "Promoter",
      homer == "Intergenic" ~ "Intergenic",
      homer == "intron" ~ "Intron",
      homer == "exon" ~ "Exon",
      TRUE ~ "Other"))            # TTS, UTRs, non-coding
} else {
  ann_src <- "ChIPseeker (overlapping consensus peak)"
  message("HOMER annotation not found; falling back to ChIPseeker annotation of consensus peaks")
  ann_cl <- clusters |>
    inner_join(annot, by = c("consensus_peak_primary" = "peak_id")) |>
    dplyr::select(cluster, category)
}
ann_cl <- ann_cl |>
  mutate(cluster = factor(cluster, levels = OCR_CLUSTERS),
         category = factor(category, levels = c("Promoter", "Intergenic", "Intron", "Exon", "Other"))) |>
  dplyr::count(cluster, category) |>
  group_by(cluster) |> mutate(pct = 100 * n / sum(n)) |> ungroup() |>
  mutate(source = ann_src)
write_panel_table(ann_cl, "fig1c_annotation")

pal_annot <- c(Promoter = "#2B2B2B", Intergenic = "#1F5A96", Intron = "#F2C12E",
               Exon = "#BFE3CF", Other = "#F7EBC3")
p <- ggplot(ann_cl, aes(x = 1, y = pct, fill = category)) +
  geom_col(width = 1, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(pct >= 8, sprintf("%.1f%%", pct), "")),
            position = position_stack(vjust = 0.5), size = 2.6, colour = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ cluster, ncol = 1) +
  scale_fill_manual(values = pal_annot, name = NULL) +
  theme_void(base_family = "Helvetica") +
  theme(strip.text = element_text(size = 9))
save_panel(p, "fig1c_annotation_pies", width = 2.4, height = 8)

# =============================================================================
# 1D — ARID1A obs/exp (HOMER mergePeaks -matrix log ratio)
# =============================================================================
f <- file.path(HOMER, "overlap/fig1d.logRatio.matrix.txt")
if (require_inputs(f, "Fig 1D HOMER overlap matrix")) {
  m <- read.delim(f, row.names = 1, check.names = FALSE)
  rn <- sub("\\.txt$", "", rownames(m)); cn <- sub("\\.txt$", "", colnames(m))
  dimnames(m) <- list(rn, cn)
  cl_ids <- gsub(" ", "_", OCR_CLUSTERS)
  ar_ids <- paste0("ARID1A_", c("Naive", "48h", "D5", "D8"))
  lr <- as.matrix(m[cl_ids, ar_ids])
  dimnames(lr) <- list(OCR_CLUSTERS, c("Naive", "48h", "d5", "d8"))
  write_panel_table(as_tibble(lr, rownames = "cluster"), "fig1d_arid1a_obs_exp_logratio")
  ht <- Heatmap(lr, name = "log ratio\n(obs/exp)", cluster_rows = FALSE, cluster_columns = FALSE,
                col = colorRamp2(seq(min(0, min(lr)), max(lr), length.out = 5),
                                 c("#FFFFD9", "#C7E9B4", "#41AB5D", "#006D2C", "#00441B")),
                rect_gp = gpar(col = "black", lwd = 0.5),
                column_title = "ARID1A CUT&RUN", column_names_rot = 90,
                row_names_side = "left", border = TRUE)
  save_panel_heatmap(ht, "fig1d_arid1a_overlap", width = 3.4, height = 3)
}

# =============================================================================
# 1F / S1C — HOMER known motifs per cluster
# =============================================================================
read_known <- function(cl) {
  f <- file.path(HOMER, "motifs", gsub(" ", "_", cl), "knownResults.txt")
  if (!file.exists(f)) return(NULL)
  read.delim(f, check.names = FALSE) |>
    as_tibble() |>
    transmute(cluster = cl, motif = `Motif Name`,
              short = sub("/.*", "", motif),
              neglog10p = -`Log P-value` / log(10))
}
known <- bind_rows(lapply(OCR_CLUSTERS, read_known))
if (nrow(known)) {
  known <- known |> mutate(cluster = factor(cluster, levels = OCR_CLUSTERS))
  write_panel_table(known, "fig1f_homer_known_all")

  # Motif columns as shown in paper Fig 1F (family label, HOMER motif prefix)
  f1f <- tribble(
    ~family,    ~short,
    "ETS,Runt", "ETS:RUNX(ETS,Runt)",
    "ETS",      "Fli1(ETS)",
    "ETS",      "ETS1(ETS)",
    "ETS",      "GABPA(ETS)",
    "Zf",       "BORIS(Zf)",
    "Zf",       "CTCF(Zf)",
    "Zf",       "KLF3(Zf)",
    "Zf",       "KLF6(Zf)",
    "bHLH",     "E2A(bHLH)",
    "HMG",      "LEF1(HMG)",
    "HMG",      "Tcf7(HMG)",
    "HMG",      "Tcf3(HMG)",
    "Runt",     "RUNX(Runt)",
    "IRF:bZIP", "IRF:BATF(IRF:bZIP)",
    "bZIP",     "AP-1(bZIP)",
    "bZIP",     "BATF(bZIP)",
    "bZIP",     "Fos(bZIP)",
    "RHD,bZIP", "NFAT:AP1(RHD,bZIP)",
    "NR",       "Nur77(NR)",
    "NR",       "RORg(NR)",
    "NR",       "RORgt(NR)",
    "T-box",    "Eomes(T-box)",
    "T-box",    "Tbx21(T-box)"
  )
  sel <- known |> semi_join(f1f, by = "short") |>
    distinct(cluster, short, .keep_all = TRUE)
  mat <- sel |> dplyr::select(cluster, short, neglog10p) |>
    pivot_wider(names_from = short, values_from = neglog10p) |>
    column_to_rownames("cluster") |> as.matrix()
  mat <- mat[intersect(OCR_CLUSTERS, rownames(mat)), intersect(f1f$short, colnames(mat)), drop = FALSE]
  write_panel_table(as_tibble(mat, rownames = "cluster"), "fig1f_motif_neglog10p")
  fam <- f1f$family[match(colnames(mat), f1f$short)]
  ht <- Heatmap(pmin(mat, 150), name = "-log10 p",
                col = colorRamp2(c(0, 50, 100, 150), c("white", "#9EBCDA", "#8856A7", "#4D004B")),
                cluster_rows = FALSE, cluster_columns = FALSE,
                column_split = factor(fam, levels = unique(fam)), column_title_rot = 45,
                column_title_gp = gpar(fontsize = 8), column_gap = unit(0, "mm"),
                column_labels = sub("\\(.*", "", colnames(mat)),
                column_names_gp = gpar(fontsize = 8), row_names_side = "left",
                border = TRUE, rect_gp = gpar(col = NA))
  save_panel_heatmap(ht, "fig1f_motif_enrichment", width = 8, height = 3.2)

  top15 <- known |> group_by(cluster) |> slice_max(neglog10p, n = 15, with_ties = FALSE) |>
    ungroup() |> mutate(label = paste(sub("/.*", "", motif), cluster, sep = "___"))
  p <- ggplot(top15, aes(neglog10p, reorder(label, neglog10p))) +
    geom_col(fill = "grey30", width = 0.7) +
    facet_wrap(~ cluster, scales = "free", ncol = 5) +
    scale_y_discrete(labels = function(x) sub("___.*", "", x)) +
    labs(x = "-log10 p-value", y = NULL) + theme_panel +
    theme(axis.text.y = element_text(size = 7))
  save_panel(p, "figS1c_top15_motifs", width = 14, height = 3.6)
  write_panel_table(top15 |> dplyr::select(-label), "figS1c_top15_motifs")
}

# =============================================================================
# S1B — expression of genes annotated to each cluster (WT d3/d5/d8)
# =============================================================================
nc <- read_csv(file.path(paths$rnaseq, "differential/normalized_counts.csv"), show_col_types = FALSE)
wt_cols <- list(d3 = grep("^D3_WT", names(nc), value = TRUE),
                d5 = grep("^D5_WT", names(nc), value = TRUE),
                d8 = grep("^D8_WT", names(nc), value = TRUE))
expr <- sapply(wt_cols, function(cc) rowMeans(log2(as.matrix(nc[, cc]) + 1)))
rownames(expr) <- nc$gene_name
expr <- expr[!duplicated(rownames(expr)) & rowMeans(expr) > 2, ]
z <- t(scale(t(expr)))

s1b <- clusters |> inner_join(annot, by = c("consensus_peak_primary" = "peak_id")) |>
  filter(!is.na(SYMBOL), abs(distanceToTSS) <= 50000) |>
  distinct(cluster, SYMBOL) |>
  inner_join(as_tibble(z, rownames = "SYMBOL"), by = "SYMBOL") |>
  pivot_longer(c(d3, d5, d8), names_to = "timepoint", values_to = "z")
write_panel_table(s1b, "figS1b_cluster_gene_expression")
p <- ggplot(s1b, aes(timepoint, z, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  facet_wrap(~ cluster, nrow = 1) +
  scale_fill_manual(values = pal_cluster, guide = "none") +
  coord_cartesian(ylim = c(-1.6, 1.6)) +
  labs(x = NULL, y = "Expression (row z-score)",
       caption = "Genes nearest each OCR (<= 50 kb). Naive RNA (GSE152841) not included.") +
  theme_panel
save_panel(p, "figS1b_cluster_gene_expression", width = 8, height = 2.8)

# =============================================================================
# 2C — % of each cluster lost in Arid1a KO
# =============================================================================
contr <- c(d3 = "KO_vs_WT_D3", d5 = "KO_vs_WT_D5", d8 = "KO_vs_WT_D8_pseudobulk")
# Denominator = all published regions in the cluster (regions without a
# consensus peak or not tested count as not lost); pct_lost_tested uses only
# regions whose consensus peak has a DA result.
lost <- imap_dfr(contr, function(ct, tp) {
  da <- read_da(ct) |> dplyr::select(peak_id, direction)
  clusters |> left_join(da, by = c("consensus_peak_primary" = "peak_id")) |>
    group_by(cluster) |>
    summarise(n = n(), n_tested = sum(!is.na(direction)),
              n_lost = sum(direction == "lost", na.rm = TRUE), .groups = "drop") |>
    mutate(timepoint = tp, pct_lost = 100 * n_lost / n,
           pct_lost_tested = 100 * n_lost / n_tested)
})
write_panel_table(lost, "fig2c_pct_lost_by_cluster")
p <- ggplot(lost, aes(pct_lost, fct_rev(cluster), fill = timepoint)) +
  geom_col(position = position_dodge2(reverse = TRUE), width = 0.8) +
  scale_fill_manual(values = c(d3 = "#D17A22", d5 = "#8B1A1A", d8 = "black"), name = NULL) +
  labs(x = expression("% lost accessibility in "*italic(Arid1a)^cKO), y = "ATAC cluster") +
  theme_panel
save_panel(p, "fig2c_pct_lost", width = 3.8, height = 3.2)

# =============================================================================
# 2E — ATAC profiles WT vs KO per cluster
# =============================================================================
f <- file.path(MAT, "fig2e_profile.tab")
if (require_inputs(f, "Fig 2E profile table")) {
  pr <- read_profile_tab(f) |>
    filter(sample != "Naive_WT") |>
    separate(sample, c("tp", "geno"), sep = "_") |>
    mutate(tp = tolower(tp), geno = factor(geno, levels = names(pal_geno)))
  write_panel_table(pr, "fig2e_atac_profiles")
  p <- ggplot(pr, aes(pos_kb, value, colour = geno)) +
    geom_line(linewidth = 0.6) +
    facet_grid(region ~ tp, scales = "free_y") +
    scale_colour_manual(values = pal_geno, name = NULL,
                        labels = c(WT = "WT", Het = "Arid1a cHet", KO = "Arid1a cKO")) +
    labs(x = "Distance from peak center (kb)", y = "ATAC coverage") +
    theme_panel
  save_panel(p, "fig2e_atac_profiles", width = 6, height = 8)
}

# =============================================================================
# 2F — H3K27ac CUT&RUN d5 WT vs KO on ATAC peaks
# =============================================================================
f  <- file.path(MAT, "fig2f_profile.tab")
fc <- file.path(MAT, "fig2f_profile_by_cluster.tab")
if (require_inputs(c(f, fc), "Fig 2F profile tables")) {
  pr_all <- read_profile_tab(f) |> mutate(region = "All clustered ATAC peaks")
  pr_cl  <- read_profile_tab(fc) |> mutate(region = as.character(region))
  pr <- bind_rows(pr_all, pr_cl) |>
    mutate(geno = factor(sub(".*_", "", sample), levels = c("WT", "KO")),
           region = factor(region, levels = c("All clustered ATAC peaks", OCR_CLUSTERS)))
  write_panel_table(pr, "fig2f_h3k27ac_profiles")
  p <- ggplot(filter(pr, region == "All clustered ATAC peaks"), aes(pos_kb, value, colour = geno)) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = c(WT = "black", KO = unname(pal_genotype["KO"])),
                        labels = c(WT = "WT", KO = "Arid1a cKO"), name = NULL) +
    labs(title = "d5", x = "Distance from ATAC peak (kb)", y = "H3K27ac CUT&RUN coverage") +
    theme_panel
  save_panel(p, "fig2f_h3k27ac_d5", width = 3, height = 2.8)
  p2 <- p
  p2$data <- filter(pr, region != "All clustered ATAC peaks")
  p2 <- p2 +
    facet_wrap(~ region, nrow = 1, scales = "free_y") + labs(title = "d5, by OCR cluster")
  save_panel(p2, "fig2f_h3k27ac_d5_by_cluster", width = 10, height = 2.6)
}

message("Fig 1/2 panels done.")
