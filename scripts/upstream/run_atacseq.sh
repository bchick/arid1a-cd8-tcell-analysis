#!/usr/bin/env bash
# =============================================================================
# upstream/run_atacseq.sh — nf-core/atacseq 2.1.2 on the study ATAC-seq
# McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing
#
# 62 libraries: 54 ATAC-seq (GSE228171, 75 bp PE) + 8 BAF inhibitor
# (GSE228193, 100 bp PE). GRCm39 / GENCODE vM35, Bowtie2 alignment,
# MACS2 narrow peaks. Resumable (-resume).
# Requires: Nextflow and Docker (profiles in nextflow/nextflow.config).
#
# Inputs:  nextflow/samplesheets/atacseq_samplesheet.csv, data/fastq/atacseq/
# Outputs: results/atac/ (log: results/atac/pipeline_info/atacseq_run.log)
# Usage:   bash scripts/upstream/run_atacseq.sh
# =============================================================================

set -euo pipefail

PROJECT_DIR="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
LOG_DIR="${PROJECT_DIR}/results/atac/pipeline_info"
mkdir -p "${LOG_DIR}"

cd "${PROJECT_DIR}"

nextflow run nf-core/atacseq -r 2.1.2 \
    -c nextflow/nextflow.config \
    -profile atacseq \
    --input nextflow/samplesheets/atacseq_samplesheet.csv \
    --aligner bowtie2 \
    --read_length 75 \
    --narrow_peak \
    --save_reference \
    --mito_name chrM \
    --macs_gsize 2407883318 \
    -resume \
    2>&1 | tee "${LOG_DIR}/atacseq_run.log"

echo "Pipeline finished at $(date)"
