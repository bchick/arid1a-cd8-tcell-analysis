# =============================================================================
# arid1a-cd8-tcell-analysis — entry points
#
#   make fetch-data    download + verify the processed-data bundle (Zenodo)
#   make figures       every reproducible panel of McDonald, Chick et al. 2023
#   make extensions    analyses beyond the paper (dose response, chromVAR, TOBIAS, ...)
#   make meta          cross-study comparison with Guo et al. 2022, then the QC audit
#   make all           figures + extensions + meta
#   make help          this list (default)
#   make validate      release hygiene + bundle checksums + headline numbers
#
#   make upstream      FASTQ download -> nf-core -> heavy precompute
#                      (needs ~50 TB disk, 400 GB RAM, Docker; see README)
#
# Run from the repository root. Scripts that need raw-data inputs (BAMs,
# bigWigs, genome FASTA) print "Skipping: ..." when those are absent.
# =============================================================================

SHELL  := /bin/bash
RSCRIPT ?= Rscript --no-save --no-restore
C      := scripts/core
U      := scripts/upstream
P      := scripts/paper
E      := extended_analysis/scripts
THREADS ?= 16
# Guo et al. used DiffBind 2.16.0; guo2022/04_diffbind216.R runs it in that exact container
DIFFBIND216 ?= quay.io/biocontainers/bioconductor-diffbind:2.16.0--r40h5f743cb_2

.PHONY: help all fetch-data core figures extensions meta validate upstream \
        upstream-download upstream-nfcore upstream-precompute
.DEFAULT_GOAL := help

help:  ## list the entry points (the header of this Makefile)
	@sed -n '/^# arid1a/,/^# =====/p' $(firstword $(MAKEFILE_LIST)) | sed -e '$$d' -e 's/^# \{0,1\}//'

all: figures extensions meta

fetch-data:
	bash tools/fetch_data.sh

# ---- core downstream analyses (DE/DA tables every figure builds on) ---------
core:
	$(RSCRIPT) $(C)/01_rnaseq_analysis.R
	$(RSCRIPT) $(C)/02_atacseq_analysis.R
	$(RSCRIPT) $(C)/03_atac_temporal_clustering.R
	$(RSCRIPT) $(C)/04_cutandrun_analysis.R
	$(RSCRIPT) $(C)/05_chipseq_analysis.R

# ---- published panels (scripts/paper/) ---------------------------------------
figures: core
	$(RSCRIPT) $(P)/figS1_denovo_clusters.R
	$(RSCRIPT) $(P)/fig1_published_clusters.R
	$(RSCRIPT) $(P)/fig1_fig2_panels.R
	$(RSCRIPT) $(P)/fig2g_2i.R
	$(RSCRIPT) $(P)/fig3.R
	$(RSCRIPT) $(P)/fig4.R
	$(RSCRIPT) $(P)/fig5.R
	$(RSCRIPT) $(P)/fig5a_motifs.R
	$(RSCRIPT) $(P)/fig4_fig5_profiles.R
	$(RSCRIPT) $(P)/fig6a_pagerank.R
	$(RSCRIPT) $(P)/figS_tracks.R

# ---- extensions beyond the paper ---------------------------------------------
extensions: core
	$(RSCRIPT) $(E)/integration/01_multiomic_integration.R
	$(RSCRIPT) $(E)/temporal_clustering/01_degpatterns_wt_ko.R
	$(RSCRIPT) $(E)/temporal_clustering/02_kmeans_mfuzz_wt_ko.R
	THREADS=$(THREADS) RSCRIPT="$(RSCRIPT)" bash $(E)/atac_trajectories/01_timecourse_patterns.sh
	$(RSCRIPT) $(E)/atac_trajectories/03_ko_trajectory_projection.R
	THREADS=$(THREADS) RSCRIPT="$(RSCRIPT)" bash $(E)/atac_trajectories/04_lost_site_motifs.sh
	$(RSCRIPT) $(E)/atac_trajectories/05_lost_site_motifs_plot.R
	$(RSCRIPT) $(E)/atac_trajectories/06_tf_expression_vs_motifs.R
	$(RSCRIPT) $(E)/atac_trajectories/07_glmnet_motif_models.R
	$(RSCRIPT) $(E)/het_dose_response/01_dose_classes.R
	$(RSCRIPT) $(E)/het_dose_response/02_feature_enrichment.R
	$(RSCRIPT) $(E)/het_dose_response/03_known_motif_summary.R
	$(RSCRIPT) $(E)/het_dose_response/04_denovo_motif_summary.R
	$(RSCRIPT) $(E)/motif_grammar/01_region_sets.R
	$(RSCRIPT) $(E)/motif_grammar/02_grammar_features.R
	$(RSCRIPT) $(E)/motif_grammar/03_spamo_spacing.R
	$(RSCRIPT) $(E)/motif_grammar/04_noise_floor.R
	$(RSCRIPT) $(E)/figures/overview.R
	$(RSCRIPT) $(E)/figures/cbaf_landscape.R
	$(RSCRIPT) $(E)/figures/dose_response.R
	$(RSCRIPT) $(E)/figures/chromvar_dose.R
	$(RSCRIPT) $(E)/figures/tf_inhibitor.R
	$(RSCRIPT) $(E)/figures/tbet_motifs.R
	$(RSCRIPT) $(E)/figures/tobias_footprints.R
	$(RSCRIPT) $(E)/figures/tobias_dose_response.R

# ---- Guo et al. 2022 cross-study comparison -----------------------------------
meta: core
	$(RSCRIPT) $(E)/guo2022/01_rnaseq_de.R
	$(RSCRIPT) $(E)/guo2022/02_atacseq_da.R
	$(RSCRIPT) $(E)/guo2022/03_summit_requant.R
	@if command -v docker >/dev/null 2>&1; then \
	  docker run --rm -u $$(id -u):$$(id -g) -v "$(CURDIR)":/work -w /work \
	    -e ARID1A_PROJECT_DIR=/work -e RENV_ACTIVATE_PROJECT=FALSE $(DIFFBIND216) \
	    Rscript --no-save --no-restore $(E)/guo2022/04_diffbind216.R; \
	else echo "Skipping $(E)/guo2022/04_diffbind216.R: needs Docker for the DiffBind 2.16.0 container"; fi
# QC audit of every differential contrast, main project and Guo (needs both)
	$(RSCRIPT) $(E)/qc_audit/qc_audit.R

validate:
	python3 tools/validate_release.py

# ---- upstream: raw data -> processed bundle ----------------------------------
upstream: upstream-download upstream-nfcore upstream-precompute

upstream-download:
	bash $(U)/download_fastqs.sh
	bash $(U)/download_meta_fastqs.sh

upstream-nfcore:
	bash $(U)/run_atacseq.sh
	bash $(U)/run_rnaseq.sh
	bash $(U)/run_cutandrun.sh
	bash $(U)/run_chipseq.sh
	bash $(U)/run_meta_guo_atacseq.sh
	bash $(U)/run_meta_guo_rnaseq.sh

# heavy steps whose compact outputs ship in the bundle
upstream-precompute:
	bash $(E)/footprinting/run_tobias_footprinting.sh
	$(RSCRIPT) $(E)/chromvar/build_chromvar.R
	$(RSCRIPT) $(P)/figS1_denovo_clusters.R
	$(RSCRIPT) $(P)/fig1_published_clusters.R
	bash $(P)/fig1_fig2_deeptools.sh $(THREADS)
	bash $(P)/fig1d_arid1a_overlap.sh $(THREADS)
	bash $(P)/fig1_homer.sh $(THREADS)
	$(RSCRIPT) $(P)/fig4.R
	$(RSCRIPT) $(P)/fig5.R
	bash $(P)/fig4_fig5_signal.sh $(THREADS)
	bash $(P)/fig5a_homer.sh $(THREADS)
