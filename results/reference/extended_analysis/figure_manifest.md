# Extended-analysis summary figures

Summary figures for the analyses beyond the paper. Each is a vector PDF plus a 300 dpi PNG, written to `figures/extended_analysis/` by `make extensions`; reference renders are in `results/reference/extended_analysis/`. Scripts are in `extended_analysis/scripts/figures/` unless a path is given. Palettes and typography come from `extended_analysis/scripts/utils_figures.R`.

## Overview

| File | Content | Script | Source data |
|------|---------|--------|-------------|
| `overview_study_design` | Assays × timepoints × genotypes | `overview.R` | `data/metadata/master_sample_sheet.tsv` |
| `overview_qc_pca` | RNA and ATAC PCA by genotype and timepoint | `overview.R` | `rnaseq_analysis.RData`, `atacseq_analysis.RData` (vst) |

## cBAF and the effector/memory enhancer landscape

| File | Content | Script | Source data |
|------|---------|--------|-------------|
| `cbaf_ocr_dynamics_heatmap` | WT OCR temporal dynamics (Naive, D3, D5, D8) | `cbaf_landscape.R` | `temporal_clustering.RData` (`all_z_km`) |
| `cbaf_ocr_cluster_profiles` | Temporal cluster profiles | `cbaf_landscape.R` | same |
| `cbaf_arid1a_signal_heatmap` | ARID1A CUT&RUN signal at D5 WT peaks, by timepoint and WT vs KO | `arid1a_signal_heatmap.sh` | CUT&RUN ARID1A bigWigs, D5 WT consensus peaks |
| `cbaf_ko_accessibility_collapse` | D8 KO vs WT volcano; differential OCR counts for KO, Tbx21 KO and Het | `cbaf_landscape.R` | `da_KO_vs_WT_D8_pseudobulk.csv`, `da_summary.csv` |
| `cbaf_rna_signature_heatmap` | Effector/memory signature expression at D8 | `cbaf_landscape.R` | `rnaseq_analysis.RData` |
| `cbaf_rna_de_gsea` | D8 KO vs WT volcano and Hallmark GSEA | `cbaf_landscape.R` | `gsea_KO_vs_WT_D8_pseudobulk.csv` |

## ARID1A dose response (WT vs Het vs KO)

Full write-up: `docs/reports/het_dose_response.md`.

| File | Content | Script | Source data |
|------|---------|--------|-------------|
| `dose_concept` | Dose-response class definitions | `dose_response.R` | — |
| `dose_class_counts_curves` | Class counts per subset and representative dose curves | `dose_response.R` | `class_summary.csv`, `dose_classes_combined.csv` |
| `dose_direction_lost_gained` | Direction of change by class (haploinsufficient ~90% lost, buffered ~70% lost) | `dose_response.R` | `class_summary.csv` |
| `dose_feature_stratification` | Baseline accessibility, peak width and ARID1A occupancy by class | `dose_response.R` | `feature_summary_by_class.csv` |
| `dose_feature_enrichment` | Fisher enrichment of features vs insensitive peaks | `dose_response.R` | `feature_enrichment_vs_insensitive.csv` |
| `dose_motifs_vs_insensitive` | De novo motif families, each class vs insensitive (ETS marks buffered peaks) | `dose_response.R` | `motif_denovo_top_per_comparison.csv` |
| `dose_motifs_haplo_vs_buffered` | Haploinsufficient vs buffered motifs, by subset | `dose_response.R` | same |
| `dose_chromvar` | chromVAR motif deviations across the allelic series | `chromvar/build_chromvar.R` (precompute), `chromvar_dose.R` | `chromvar_dose.RData` |
| `dose_synthesis_model` | Summary schematic | `dose_response.R` | — |

## TF binding, BAF inhibition and footprinting

Footprinting write-up: `docs/reports/tobias_footprinting.md`.

| File | Content | Script | Source data |
|------|---------|--------|-------------|
| `tf_tbet_inhibitor` | T-bet ChIP ± BAF inhibitors: binding expands, mostly at new sites | `tf_inhibitor.R` | `chipseq_analysis.RData`, `Tbet_inhibitor_peak_overlap.csv` |
| `tf_tbet_motifs` | Motif content of T-bet peaks per condition (HOMER known and de novo) | `tbet_motifs.R` | HOMER `knownResults.txt`, `homerMotifs.all.motifs` |
| `tf_arid1a_cooccupancy` | TF CUT&RUN peaks partitioned by ARID1A dependence | `tf_inhibitor.R` | `TF_ARID1A_dependency.csv`, `D5_wt_ko_peak_overlap.csv` |
| `tf_tobias_footprints` | TOBIAS differential footprints, WT vs KO | `tobias_footprints.R` | `footprinting/bindetect/bindetect_results.txt` |
| `tf_tobias_dose_response` | Footprint retention at one *Arid1a* copy, by TF family | `tobias_dose_response.R` | `footprinting/tobias_dose_response_{per_motif,by_family}.csv` |
| `tf_tobias_chromvar_concordance` | TOBIAS vs chromVAR per-motif retention (Spearman ρ = 0.55, n = 232) | `tobias_dose_response.R` | `footprinting/tobias_vs_chromvar_dose_per_motif.csv`, `chromvar_dose.RData` |

Figures for the trajectory, temporal-clustering, motif-grammar, Guo *et al.* 2022 and QC analyses are written to their own subfolders of `figures/extended_analysis/`.
