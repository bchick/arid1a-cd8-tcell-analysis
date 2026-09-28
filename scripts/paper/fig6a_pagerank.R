#!/usr/bin/env Rscript
# =============================================================================
# paper/fig6a_pagerank.R — Fig 6A + Fig S6A: Taiji-style PageRank TF ranking,
# WT vs Arid1a cKO day-8 MP cells
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# Re-implements the Taiji algorithm (Zhang et al. 2019 NAR 47:e40) in R,
# without chromatin-interaction data:
#   1. Regulatory elements = nf-core ATAC consensus peaks (merged_replicate).
#      Peak -> gene links follow Taiji's no-Hi-C rule, which borrows GREAT's
#      "basal plus extension" domains:
#        promoter: peak overlaps TSS -5 kb/+1 kb (strand-aware)  -> weight 1
#        distal  : any other peak is linked to its nearest gene TSS upstream
#                  and nearest downstream (<= 2 genes, as GREAT/Taiji
#                  effectively do), only if that TSS lies within 50 kb
#                  (tighter than Taiji's 1 Mb cap),
#                  weight = exp(-d / 10 kb), d = |peak midpoint - TSS|.
#      Gene TSSs come from the GENCODE vM35 GTF (gene records).
#   2. TF -> peak: JASPAR2020 CORE vertebrate PWM hits (motifmatchr,
#      p < 5e-5). Dimer motifs (A::B) are credited to each partner; motif
#      variants of one TF are merged (any hit counts).
#   3. Edge weight TF -> gene (per condition c):
#        w = sqrt(expr_c(TF)) * sum_p [ acc_c(p) * link(p, gene) ]
#      over motif-bearing peaks p linked to the gene; acc_c(p) = mean CPM of
#      the peak across the condition's ATAC libraries.
#   4. Node weight (personalisation) = exp(z), z = gene's log2 expression in
#      c standardised across all RNA-seq libraries (as in Taiji).
#   5. Personalised PageRank on the reversed network (gene -> TF), so a TF
#      scores highly when it regulates many highly/specifically expressed
#      genes through accessible elements. TF score = PageRank of its node.
#
# Fig 6A : log2 PageRank ratio (WT/KO) vs log2 mRNA ratio (WT/KO)
# Fig S6A: top 30 TFs per direction with |ratio| > 1.5
#
# Inputs:  consensus_peaks.mRp.clN.featureCounts.txt (merged_replicate)
#          results/rnaseq/differential/{normalized_counts,de_KO_vs_WT_D8_MP}.csv
#          genome FASTA (motif matching), GENCODE vM35 GTF (TSSs)
# Outputs: figures/paper/fig6a_pagerank.{pdf,png}, figS6a_pagerank_ratio.{pdf,png}
#          results/paper/fig6a_pagerank.csv, fig6a_pagerank_replicates.csv
# Usage:   Rscript scripts/paper/fig6a_pagerank.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Rsamtools)
  library(motifmatchr)
  library(TFBSTools)
  library(JASPAR2020)
  library(igraph)
  library(ggrepel)
})

set.seed(42)
PROM_UP <- 5000; PROM_DOWN <- 1000   # Taiji promoter window
DISTAL_MAX <- 50000; DISTAL_DECAY <- 10000
PR_DAMPING <- 0.85
RATIO_FC <- 1.5

fc_file  <- file.path(paths$atac, "bowtie2/merged_replicate/macs2/narrow_peak/consensus",
                      "consensus_peaks.mRp.clN.featureCounts.txt")
rna_file <- file.path(paths$rnaseq, "differential/normalized_counts.csv")
de_file  <- file.path(paths$rnaseq, "differential/de_KO_vs_WT_D8_MP.csv")

# Genome-derived inputs are cached so the panel rebuilds from the data bundle
# without the GTF/FASTA: gene TSSs and the peak x motif match matrix.
tss_cache   <- file.path(paths$paper_tab, "fig6a_gene_tss.tsv.gz")
motif_cache <- file.path(paths$paper_tab, "fig6a_motif_hits.rds")
if (!require_inputs(c(fc_file, rna_file, de_file), "Fig 6A inputs") ||
    !(file.exists(genome$gtf) || file.exists(tss_cache)) ||
    !(file.exists(genome$fasta) || file.exists(motif_cache))) {
  message("Skipping Fig 6A: needs the GTF/FASTA or their caches in results/paper/")
  quit(save = "no", status = 0)
}

# =============================================================================
# 1. Consensus ATAC peaks + condition-specific accessibility (CPM)
# =============================================================================

message("=== ATAC consensus peaks + accessibility (D8 MP WT vs KO) ===")
fc <- as.data.frame(data.table::fread(fc_file, skip = "Geneid"))
cnt <- as.matrix(fc[, -(1:6)])
colnames(cnt) <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(cnt))
rownames(cnt) <- fc$Geneid
cpm <- sweep(cnt, 2, colSums(cnt) / 1e6, "/")

atac_cols <- list(
  WT = grep("^D8_WT_MP_", colnames(cpm), value = TRUE),
  KO = grep("^D8_KO_MP_", colnames(cpm), value = TRUE)
)
message("  WT libs: ", paste(atac_cols$WT, collapse = ", "))
message("  KO libs: ", paste(atac_cols$KO, collapse = ", "))
acc <- sapply(atac_cols, function(cols) rowMeans(cpm[, cols, drop = FALSE]))

std_chr <- paste0("chr", c(1:19, "X", "Y"))
el_gr <- GRanges(fc$Chr, IRanges(fc$Start, fc$End), name = fc$Geneid)
el_gr <- el_gr[as.character(seqnames(el_gr)) %in% std_chr & rowSums(acc[fc$Geneid, ]) > 0]
el_gr <- keepSeqlevels(el_gr, intersect(std_chr, seqlevels(el_gr)))
message(sprintf("  %d peaks on primary chromosomes", length(el_gr)))

# =============================================================================
# 2. Peak -> gene links (Taiji / GREAT basal-plus-extension, no Hi-C)
# =============================================================================

message("=== Peak -> gene links (promoter -5/+1 kb; distal nearest-2 within 50 kb) ===")
rna <- readr::read_csv(rna_file, show_col_types = FALSE)
if (file.exists(genome$gtf)) {
  gtf_genes <- data.table::fread(cmd = paste("awk -F'\\t' '$3==\"gene\"'", shQuote(genome$gtf)),
                                 header = FALSE, sep = "\t", quote = "")
  all_tss <- tibble(chr = gtf_genes$V1, strand = gtf_genes$V7,
                    tss = ifelse(gtf_genes$V7 == "+", gtf_genes$V4, gtf_genes$V5),
                    gene = sub('.*gene_name "([^"]+)".*', "\\1", gtf_genes$V9))
  readr::write_tsv(all_tss, tss_cache)
} else {
  all_tss <- readr::read_tsv(tss_cache, show_col_types = FALSE)
}
tss <- all_tss |>
  filter(chr %in% std_chr, gene %in% rna$gene_name) |>
  distinct(chr, tss, strand, gene)
tss_gr <- GRanges(tss$chr, IRanges(tss$tss, width = 1), strand = tss$strand, gene = tss$gene)
message(sprintf("  %d TSSs for %d genes in the RNA-seq annotation", length(tss_gr), n_distinct(tss$gene)))

# (a) promoter links: peak overlaps TSS -5 kb / +1 kb
prom_gr <- promoters(tss_gr, upstream = PROM_UP, downstream = PROM_DOWN)
ov_p <- findOverlaps(el_gr, prom_gr, ignore.strand = TRUE)
prom_links <- tibble(peak = el_gr$name[queryHits(ov_p)], gene = prom_gr$gene[subjectHits(ov_p)],
                     weight = 1) |> distinct()

# (b) distal links: nearest TSS on each side of the peak midpoint, <= 50 kb
distal_gr <- el_gr[!el_gr$name %in% prom_links$peak]
mid <- resize(distal_gr, width = 1, fix = "center")
tss_us <- unstrand(tss_gr)
side_link <- function(idx) {
  ok <- !is.na(idx)
  d <- abs(start(mid)[ok] - start(tss_us)[idx[ok]])
  tibble(peak = distal_gr$name[ok], gene = tss_us$gene[idx[ok]], d = d)
}
distal_links <- bind_rows(side_link(precede(mid, tss_us)), side_link(follow(mid, tss_us))) |>
  filter(d <= DISTAL_MAX) |>
  group_by(peak, gene) |> summarise(d = min(d), .groups = "drop") |>
  mutate(weight = exp(-d / DISTAL_DECAY)) |>
  select(peak, gene, weight)

links <- bind_rows(prom_links, distal_links) |>
  group_by(peak, gene) |> summarise(weight = max(weight), .groups = "drop")
message(sprintf("  %d promoter links, %d distal links; %d peaks -> %d genes",
                nrow(prom_links), nrow(distal_links), n_distinct(links$peak), n_distinct(links$gene)))
all_el_gr <- el_gr    # motif cache covers every primary-chromosome peak, not just linked ones
el_gr <- el_gr[el_gr$name %in% links$peak]

# =============================================================================
# 3. JASPAR2020 motif hits per element, collapsed to TF genes
# =============================================================================

message("=== Motif matching (JASPAR2020 CORE vertebrates) ===")
motifs <- getMatrixSet(JASPAR2020, list(collection = "CORE", tax_group = "vertebrates",
                                        matrixtype = "PWM"))
if (file.exists(genome$fasta)) {
  fa <- FaFile(genome$fasta)
  mm <- matchMotifs(motifs, all_el_gr, genome = fa, p.cutoff = 5e-5)
  hits <- motifMatches(mm)                    # element x motif (logical, sparse)
  rownames(hits) <- all_el_gr$name
  saveRDS(hits, motif_cache)
  hits <- hits[el_gr$name, , drop = FALSE]
} else {
  hits <- readRDS(motif_cache)
  # Hits are per element, so the cache serves any subset of the elements it covers
  # (it is built on all peaks, so link-rule changes do not invalidate it)
  absent <- setdiff(el_gr$name, rownames(hits))
  if (length(absent))
    stop(length(absent), " elements are not in ", basename(motif_cache),
         "; rebuild it with the genome FASTA (upstream mode)")
  hits <- hits[el_gr$name, , drop = FALSE]
  stopifnot(ncol(hits) == length(motifs))
}
message(sprintf("  %d motifs x %d elements, %d hits", ncol(hits), nrow(hits), sum(hits)))

# Motif name -> TF gene symbol(s) in the RNA-seq annotation
sym_upper <- setNames(rna$gene_name, toupper(rna$gene_name))
motif_tf <- tibble(motif = colnames(hits), motif_name = name(motifs)[colnames(hits)]) |>
  mutate(tf = strsplit(sub("\\(var\\.[0-9]+\\)$", "", motif_name), "::")) |>
  tidyr::unnest(tf) |>
  mutate(tf = unname(sym_upper[toupper(tf)])) |>
  filter(!is.na(tf))
message(sprintf("  %d motifs map to %d TF genes", n_distinct(motif_tf$motif), n_distinct(motif_tf$tf)))

tf_levels <- sort(unique(motif_tf$tf))
M <- Matrix::sparseMatrix(i = match(motif_tf$motif, colnames(hits)),
                          j = match(motif_tf$tf, tf_levels), x = 1,
                          dims = c(ncol(hits), length(tf_levels)))
el_tf <- (hits %*% M) > 0                     # element x TF (any motif variant)
colnames(el_tf) <- tf_levels

# =============================================================================
# 4. Expression: condition means + Taiji node weights
# =============================================================================

expr_mat <- as.matrix(rna[, -(1:2)])
rownames(expr_mat) <- rna$gene_name
expr_mat <- rowsum(expr_mat, rownames(expr_mat))           # collapse duplicate symbols
log_expr <- log2(expr_mat + 1)
mu <- rowMeans(log_expr); sdv <- apply(log_expr, 1, sd); sdv[sdv < 1e-3] <- 1

rna_cols <- list(WT = grep("^D8_WT_MP_", colnames(expr_mat), value = TRUE),
                 KO = grep("^D8_KO_MP_", colnames(expr_mat), value = TRUE))

# =============================================================================
# 5. Network + personalised PageRank
# =============================================================================

# Peak x gene link matrix
links <- filter(links, peak %in% rownames(el_tf), gene %in% rownames(expr_mat))
genes <- sort(unique(links$gene))
A <- Matrix::sparseMatrix(i = match(links$peak, rownames(el_tf)), j = match(links$gene, genes),
                          x = links$weight, dims = c(nrow(el_tf), length(genes)))
el_tf_num <- el_tf * 1

#' TF scores for one condition.
#' @param expr named vector of normalised expression for this sample/condition
#' @param acc_c peak accessibility for this condition (aligned to el_tf rows)
taiji_pagerank <- function(expr, acc_c) {
  # TF -> gene weights: t(el_tf) %*% diag(acc) %*% links
  W <- Matrix::t(el_tf_num) %*% (A * acc_c)                # TF x gene
  tf_expr <- expr[rownames(W)]; tf_expr[is.na(tf_expr)] <- 0
  W <- W * sqrt(tf_expr)
  W <- as(W, "TsparseMatrix")
  edges <- tibble(tf = rownames(W)[W@i + 1], gene = genes[W@j + 1], w = W@x) |>
    filter(w > 0, tf != gene)
  nodes <- union(edges$tf, edges$gene)
  z <- (log2(expr[nodes] + 1) - mu[nodes]) / sdv[nodes]
  z[is.na(z)] <- min(z, na.rm = TRUE)
  z <- pmin(pmax(z, -10), 10)
  # Reversed network: random walk flows from target genes back to their TFs
  g <- graph_from_data_frame(edges |> transmute(from = gene, to = tf, weight = w),
                             directed = TRUE, vertices = nodes)
  pr <- page_rank(g, directed = TRUE, damping = PR_DAMPING,
                  personalized = exp(z)[V(g)$name], weights = E(g)$weight)$vector
  pr[intersect(names(pr), tf_levels)]
}

cond_expr <- function(cols) rowMeans(expr_mat[, cols, drop = FALSE])

message("=== PageRank: condition means ===")
t0 <- Sys.time()
pr <- list(WT = taiji_pagerank(cond_expr(rna_cols$WT), acc[rownames(el_tf), "WT"]),
           KO = taiji_pagerank(cond_expr(rna_cols$KO), acc[rownames(el_tf), "KO"]))
message(sprintf("  done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Per-replicate (RNA replicate i with condition-mean ATAC)
message("=== PageRank: per RNA replicate ===")
rep_pr <- list()
for (cond in names(rna_cols)) for (s in rna_cols[[cond]]) {
  rep_pr[[s]] <- taiji_pagerank(expr_mat[, s], acc[rownames(el_tf), cond])
}

# =============================================================================
# 6. Assemble the table
# =============================================================================

de <- readr::read_csv(de_file, show_col_types = FALSE)
de_tf <- de |> filter(!is.na(gene_name)) |> group_by(gene_name) |>
  slice_max(baseMean, n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(tf = gene_name, expr_log2_WT_over_KO = -log2FoldChange, expr_padj = padj)

tfs <- intersect(names(pr$WT), names(pr$KO))
res <- tibble(tf = tfs, pagerank_WT = pr$WT[tfs], pagerank_KO = pr$KO[tfs]) |>
  mutate(pr_log2_WT_over_KO = log2(pagerank_WT / pagerank_KO),
         tf_expr_WT = cond_expr(rna_cols$WT)[tf],
         tf_expr_KO = cond_expr(rna_cols$KO)[tf]) |>
  left_join(de_tf, by = "tf") |>
  # genes filtered out of DESeq2 (low counts): fall back to the normalised-count ratio
  mutate(expr_log2_WT_over_KO = dplyr::coalesce(expr_log2_WT_over_KO,
                                                log2((tf_expr_WT + 1) / (tf_expr_KO + 1)))) |>
  # Taiji reports TFs that are expressed in the cell state
  filter(pmax(tf_expr_WT, tf_expr_KO) >= 10) |>
  mutate(rank_WT = rank(-pagerank_WT), rank_KO = rank(-pagerank_KO),
         rank_impaired_in_KO = rank(-pr_log2_WT_over_KO)) |>
  arrange(desc(pr_log2_WT_over_KO))

# Replicate stability: ratio from replicate pairs (WT_i / KO_i)
rep_tab <- purrr::map_dfr(seq_along(rna_cols$WT), function(i) {
  if (i > length(rna_cols$KO)) return(NULL)
  w <- rep_pr[[rna_cols$WT[i]]]; k <- rep_pr[[rna_cols$KO[i]]]
  tibble(pair = paste0("Rep", i), tf = res$tf, pr_log2_WT_over_KO = log2(w[res$tf] / k[res$tf]))
})
rep_wide <- tidyr::pivot_wider(rep_tab, names_from = pair, values_from = pr_log2_WT_over_KO,
                               names_prefix = "pr_log2_ratio_")
res <- left_join(res, rep_wide, by = "tf")
rep_cols <- grep("^pr_log2_ratio_Rep", names(res), value = TRUE)
if (length(rep_cols) >= 2) {
  rho <- cor(res[[rep_cols[1]]], res[[rep_cols[2]]], method = "spearman", use = "complete.obs")
  top1 <- res$tf[order(-res[[rep_cols[1]]])][1:30]; top2 <- res$tf[order(-res[[rep_cols[2]]])][1:30]
  message(sprintf("  replicate log2-ratio Spearman rho = %.3f; top-30 impaired overlap = %d/30",
                  rho, length(intersect(top1, top2))))
}

paper_tfs <- c("Zfp683", "Bhlhe40", "Runx3", "Rxra", "Smad3", "Tbx21", "Eomes",
               "Gata3", "Ar", "Runx2", "Arid3b", "Arid2",
               "Irf4", "E2f1", "Egr2", "Batf3", "Myb", "Maf", "Cebpb", "Rorc")
message("=== Paper-highlighted TFs (rank by WT/KO PageRank ratio, of ", nrow(res), ") ===")
print(res |> filter(tf %in% paper_tfs) |>
        select(tf, rank_impaired_in_KO, pr_log2_WT_over_KO, expr_log2_WT_over_KO, all_of(rep_cols)) |>
        as.data.frame(), digits = 3)
missing_tfs <- setdiff(paper_tfs, res$tf)
if (length(missing_tfs))
  message("  Not scored (no JASPAR2020 motif, or not expressed): ", paste(missing_tfs, collapse = ", "))

write_panel_table(res, "fig6a_pagerank")
write_panel_table(rep_tab, "fig6a_pagerank_replicates")

# =============================================================================
# 7. Fig 6A: PageRank ratio vs expression ratio
# =============================================================================

thr <- log2(RATIO_FC)
X_LIM <- 4   # symmetric x-limits covering the bulk; TFs beyond are drawn at the edge
plot_df <- res |> filter(!is.na(expr_log2_WT_over_KO)) |>
  mutate(oob = abs(pr_log2_WT_over_KO) > X_LIM,
         label = ifelse(tf %in% paper_tfs | oob, toupper(tf), NA),
         hl = tf %in% paper_tfs,
         x_plot = scales::oob_squish(pr_log2_WT_over_KO, c(-X_LIM, X_LIM)))
n_oob <- sum(plot_df$oob)

p6a <- ggplot(plot_df, aes(x_plot, expr_log2_WT_over_KO)) +
  annotate("rect", xmin = thr, xmax = Inf, ymin = 0, ymax = Inf, fill = "grey85", alpha = 0.6) +
  annotate("rect", xmin = -Inf, xmax = -thr, ymin = -Inf, ymax = 0,
           fill = pal_genotype[["KO"]], alpha = 0.12) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(-thr, thr), linetype = "dashed") +
  geom_point(data = ~ filter(.x, !hl, !oob), color = "grey65", size = 0.9) +
  geom_point(data = ~ filter(.x, hl, !oob), color = "#CC3311", size = 1.4) +
  # out-of-range TFs: open arrowhead pointing outward, at the axis limit
  geom_point(data = ~ filter(.x, oob, x_plot < 0), shape = "<", size = 3, color = "grey30") +
  geom_point(data = ~ filter(.x, oob, x_plot > 0), shape = ">", size = 3, color = "grey30") +
  geom_label_repel(aes(label = label), size = 2.4, label.size = 0.2, min.segment.length = 0,
                   max.overlaps = Inf, box.padding = 0.3, na.rm = TRUE) +
  scale_x_continuous(limits = c(-X_LIM, X_LIM), oob = scales::oob_squish) +
  labs(x = expression("PageRank score ratio (log"[2]*", WT / "*italic(Arid1a)^cKO*")"),
       y = expression("Expression ratio (log"[2]*", WT / "*italic(Arid1a)^cKO*")"),
       title = "D8 MP: Taiji-style PageRank",
       caption = if (n_oob > 0) sprintf("%d TFs beyond axis limits shown at the edge (WT MP expression \u2248 0)", n_oob)) +
  theme_paper
save_panel(p6a, "fig6a_pagerank", width = 4.5, height = 4.2)

# =============================================================================
# 8. Fig S6A: top 30 TFs per direction, |ratio| > 1.5
# =============================================================================

s6 <- bind_rows(
  res |> filter(pr_log2_WT_over_KO >  thr) |> slice_max(pr_log2_WT_over_KO, n = 30),
  res |> filter(pr_log2_WT_over_KO < -thr) |> slice_min(pr_log2_WT_over_KO, n = 30)
) |>
  mutate(higher_in = ifelse(pr_log2_WT_over_KO > 0, "WT", "KO"),
         tf = forcats::fct_reorder(toupper(tf), pr_log2_WT_over_KO))

ps6 <- ggplot(s6, aes(pr_log2_WT_over_KO, tf, fill = higher_in)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c(WT = "grey40", KO = pal_genotype[["KO"]]),
                    labels = c(WT = "Higher in WT", KO = "Higher in Arid1a cKO"), name = NULL) +
  labs(x = expression("PageRank score ratio (log"[2]*", WT / KO)"), y = NULL,
       title = sprintf("D8 MP: TFs with PageRank ratio > %.1f-fold", RATIO_FC)) +
  theme_paper + theme(axis.text.y = element_text(size = 6), legend.position = "bottom")
save_panel(ps6, "figS6a_pagerank_ratio", width = 4, height = max(4, 0.12 * nrow(s6) + 1.5))

message("=== fig6a_pagerank.R complete ===")
