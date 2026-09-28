#!/usr/bin/env bash
# =============================================================================
# upstream/run_rnaseq_kallisto.sh — add Kallisto quantification to nf-core/rnaseq
# McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing
#
# Re-runs nf-core/rnaseq 3.22.2 with --pseudo_aligner kallisto on the study
# RNA-seq (GSE227634; 30 libraries, WT/Het/KO, D3/D5/D8). Launched from the
# repository root, where the original run's .nextflow history lives, so
# -resume reuses the cached STAR + Salmon steps and only Kallisto runs.
# Run after run_rnaseq.sh.
# Requires: Nextflow and Docker (profiles in nextflow/nextflow.config).
#
# Inputs:  nextflow/samplesheets/rnaseq_samplesheet.csv, data/fastq/rnaseq/
# Outputs: results/rnaseq/kallisto/ (log: results/rnaseq/pipeline_info/rnaseq_kallisto_run.log)
# Usage:   bash scripts/upstream/run_rnaseq_kallisto.sh
# =============================================================================

set -euo pipefail

PROJECT_DIR="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
LOG_DIR="${PROJECT_DIR}/results/rnaseq/pipeline_info"
mkdir -p "${LOG_DIR}"

# Launch from project root where .nextflow/history lives (original rnaseq run)
# This enables -resume to reuse cached STAR+Salmon work
cd "${PROJECT_DIR}"

nextflow run nf-core/rnaseq -r 3.22.2 \
    -c nextflow/nextflow.config \
    -profile rnaseq \
    --input nextflow/samplesheets/rnaseq_samplesheet.csv \
    --pseudo_aligner kallisto \
    --save_reference \
    -resume \
    2>&1 | tee "${LOG_DIR}/rnaseq_kallisto_run.log"

echo "Kallisto pipeline finished at $(date)"
