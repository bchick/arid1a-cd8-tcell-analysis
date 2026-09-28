#!/usr/bin/env bash
# =============================================================================
# upstream/run_rnaseq.sh — nf-core/rnaseq 3.22.2 on the study RNA-seq
# McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing
#
# 30 libraries (GSE227634): WT/Het/KO CD8+ T cells, D3/D5/D8; single-end
# 100 bp, reverse-stranded. GRCm39 / GENCODE vM35, STAR + Salmon (star_salmon).
# Requires: Nextflow and Docker (profiles in nextflow/nextflow.config).
#
# Inputs:  nextflow/samplesheets/rnaseq_samplesheet.csv, data/fastq/rnaseq/
# Outputs: results/rnaseq/ (log: results/rnaseq/pipeline_info/rnaseq_run.log)
# Usage:   bash scripts/upstream/run_rnaseq.sh
# =============================================================================

set -euo pipefail

PROJECT_DIR="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
LOG_DIR="${PROJECT_DIR}/results/rnaseq/pipeline_info"
mkdir -p "${LOG_DIR}"

cd "${PROJECT_DIR}"

nextflow run nf-core/rnaseq -r 3.22.2 \
    -c nextflow/nextflow.config \
    -profile rnaseq \
    --input nextflow/samplesheets/rnaseq_samplesheet.csv \
    --save_reference \
    -resume \
    2>&1 | tee "${LOG_DIR}/rnaseq_run.log"

echo "Pipeline finished at $(date)"
