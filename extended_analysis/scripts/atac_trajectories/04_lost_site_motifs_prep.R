#!/usr/bin/env Rscript
# =============================================================================
# atac_trajectories/04_lost_site_motifs_prep.R — region sets for motif enrichment of lost sites
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Known-motif enrichment follows the timecourse-patterns T cell example
# (examples/tcell_motifs/run_motifs.sh): MEME-suite SEA, JASPAR2024 CORE
# vertebrates non-redundant, every region trimmed to 200 bp around its centre
# so peak width cannot pass for enrichment. Two sets of comparisons:
#
#   vs_static    each WT trajectory class of the lost sites (pass 2 labels)
#                against static peaks, dynamic in neither timecourse: "what
#                does this temporal program carry?"
#   lost_vs_kept within each WT class (pass 1 labels, so both sides share one
#                definition), sites that lose their dynamics in KO against
#                sites that stay dynamic, in both directions: "what makes
#                this program ARID1A-dependent?"
#
# Writes one BED per region set and a manifest of (comparison, primary,
# control) for atac_trajectories/04_lost_site_motifs.sh.
#
# Inputs:  results/extended_analysis/atac_trajectories/{dynamic_sites,pass1_wt_ko,pass2_lost_in_ko}/
# Outputs: results/extended_analysis/atac_trajectories/motifs/{beds/,comparisons.tsv}
# Usage:   Rscript extended_analysis/scripts/atac_trajectories/04_lost_site_motifs_prep.R   (from the repository root;
#          normally called by 04_lost_site_motifs.sh)
# =============================================================================

source("scripts/utils.R")
source("extended_analysis/scripts/utils_trajectories.R")

HALF <- 100L  # 200 bp windows, as in the timecourse-patterns example

motif_dir <- file.path(traj_dir, "motifs")
bed_dir   <- file.path(motif_dir, "beds")
dir.create(bed_dir, recursive = TRUE, showWarnings = FALSE)

sites <- read_tsv(file.path(traj_dir, "dynamic_sites/dynamic_site_sets.tsv.gz"), show_col_types = FALSE)
read_cl <- function(pass) read_tsv(file.path(traj_dir, pass, "results/clusters/WT_clusters.tsv"),
                                   show_col_types = FALSE) |>
  filter(supercluster_label != "Unassigned") |>
  select(feature_id, class = supercluster_label)
pass1 <- read_cl("pass1_wt_ko")
pass2 <- read_cl("pass2_lost_in_ko")

bed <- attr(read_atac_consensus_counts(), "bed")
centre <- function(ids) {
  b <- bed[match(ids, bed$name), ]
  mid <- (b$start + b$end) %/% 2L
  tibble(chr = b$chr, start = pmax(0L, mid - HALF), end = mid + HALF, name = b$name)
}
slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
write_set <- function(name, ids) {
  write_tsv(centre(ids), file.path(bed_dir, paste0(name, ".bed")), col_names = FALSE)
  tibble(set = name, n = length(ids))
}

# --- region sets ---------------------------------------------------------------
sizes <- list(write_set("static", sites$feature_id[sites$set == "Static"]))
comparisons <- list()

for (k in sort(unique(pass2$class))) {
  nm <- paste0("lost_", slug(k))
  sizes[[length(sizes) + 1]] <- write_set(nm, pass2$feature_id[pass2$class == k])
  comparisons[[length(comparisons) + 1]] <-
    tibble(comparison = "vs_static", class = k, primary = nm, control = "static")
}

p1 <- inner_join(pass1, select(sites, feature_id, set), by = "feature_id")
for (k in sort(unique(p1$class))) {
  lost_nm <- paste0("p1lost_", slug(k)); kept_nm <- paste0("p1kept_", slug(k))
  sizes[[length(sizes) + 1]] <- write_set(lost_nm, p1$feature_id[p1$class == k & p1$set == "Lost"])
  sizes[[length(sizes) + 1]] <- write_set(kept_nm, p1$feature_id[p1$class == k & p1$set == "Shared"])
  comparisons[[length(comparisons) + 1]] <- tibble(
    comparison = c("lost_vs_kept", "kept_vs_lost"), class = k,
    primary = c(lost_nm, kept_nm), control = c(kept_nm, lost_nm))
}

sizes <- bind_rows(sizes); comparisons <- bind_rows(comparisons)
write_tsv(sizes, file.path(motif_dir, "region_set_sizes.tsv"))
write_tsv(comparisons, file.path(motif_dir, "comparisons.tsv"))
print(as.data.frame(sizes))
message(sprintf("  %d comparisons -> %s", nrow(comparisons), file.path(motif_dir, "comparisons.tsv")))
