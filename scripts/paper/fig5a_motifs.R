#!/usr/bin/env Rscript
# =============================================================================
# paper/fig5a_motifs.R — Figure 5A: TF motif families enriched in OCRs gained / lost
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# HOMER known motifs in OCRs gained / lost in Arid1a KO day 8 TE, EEC and MP
# cells. Each family is scored by its best motif (largest -log p); bars are
# capped at 500 as in the paper. The day 5 sets are written to the table but,
# as in the paper figure, only d8 subsets are plotted.
#
# Inputs:  results/paper/fig5a_homer/<set>/knownResults.txt (scripts/paper/fig5a_homer.sh)
# Outputs: figures/paper/fig5a_motif_families.{pdf,png}, results/paper/fig5a_motif_families.csv
# Usage:   Rscript scripts/paper/fig5a_motifs.R   (from the repository root)
# =============================================================================

source("scripts/utils.R")
source("scripts/paper/utils_paper.R")
filter <- dplyr::filter

homer_dir <- file.path(paths$paper_tab, "fig5a_homer")
sets <- expand_grid(sample = c("d5", "TE", "EEC", "MP"), direction = c("lost", "gained")) |>
  mutate(file = file.path(homer_dir, paste0(sample, "_", direction), "knownResults.txt"))
if (!require_inputs(sets$file, "HOMER knownResults (run scripts/paper/fig5a_homer.sh)")) quit(save = "no")

# Families shown in the paper, with the HOMER family string(s) that define them.
# Zf is restricted to CTCF / BORIS and bZIP to AP-1-type motifs, matching the
# paper's "Zf(CTCF)" and "bZIP(AP-1)" labels.
families <- tribble(
  ~family,       ~pattern,
  "ETS",         "^[^/]*\\(ETS\\)",
  "Runt",        "^[^/]*\\(Runt\\)",
  "T-box",       "^[^/]*\\(T-box\\)",
  "bZIP(AP-1)",  "^(AP-1|Fra1|Fra2|Fos|Fosl2|Jun-AP1|JunB|Atf3|BATF)\\(bZIP\\)",
  "Zf(CTCF)",    "^(CTCF|BORIS)\\(Zf\\)",
  "ETS:IRF",     "^[^/]*\\(ETS:IRF\\)",
  "Homeobox",    "^[^/]*\\(Homeobox\\)",
  "RHD",         "^[^/]*\\(RHD\\)",
  "IRF",         "^[^/]*\\(IRF\\)",
  "bHLH",        "^[^/]*\\(bHLH\\)"
)

known <- purrr::pmap_dfr(sets, function(sample, direction, file) {
  read_tsv(file, show_col_types = FALSE, name_repair = "minimal") |>
    transmute(sample, direction, motif = `Motif Name`, neglog10p = -`Log P-value` / log(10),
              neglnp = -`Log P-value`)
})

fam_scores <- purrr::pmap_dfr(families, function(family, pattern) {
  known |> filter(grepl(pattern, motif)) |>
    group_by(sample, direction) |>
    slice_max(neglnp, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(family = family)
})
write_panel_table(fam_scores, "fig5a_motif_families")

plot_df <- fam_scores |>
  filter(sample %in% c("TE", "EEC", "MP")) |>
  mutate(score = pmin(neglnp, 500),                   # paper: axis capped at >500
         value = ifelse(direction == "gained", -score, score),
         sample = factor(sample, levels = c("TE", "EEC", "MP")),
         family = factor(family, levels = rev(families$family)))

p5a <- ggplot(plot_df, aes(value, family, fill = sample)) +
  geom_col(position = position_dodge(0.8), width = 0.8, color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = pal_subset[c("TE", "EEC", "MP")]) +
  scale_x_continuous(limits = c(-500, 500), breaks = c(-500, -250, 0, 250, 500),
                     labels = c(">500", "250", "0", "250", ">500")) +
  labs(x = "-log(p value)", y = "TF motif family", fill = NULL,
       title = "OCRs gained in Arid1a cKO  |  OCRs lost in Arid1a cKO") +
  theme(plot.title = element_text(size = 9, face = "plain", hjust = 0.5))
save_panel(p5a, "fig5a_motif_families", width = 4.5, height = 3.8)

message("=== fig5a_motifs.R done ===")
