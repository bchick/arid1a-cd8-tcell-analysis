#!/usr/bin/env bash
# =============================================================================
# upstream/run_cutandrun.sh — nf-core/cutandrun 3.2.2 on the CUT&RUN data
# McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing
#
# 28 libraries (GSE228380): ARID1A, H3K27ac, H3K27me3, T-bet, BATF, ETS1, IgG;
# paired-end 42 bp. GRCm39 / GENCODE vM35. Peak callers MACS2 + SEACR. No
# E. coli spike-in, so CPM normalisation.
# Requires: Nextflow and Docker (profiles in nextflow/nextflow.config).
#
# Inputs:  nextflow/samplesheets/cutandrun_samplesheet.csv, data/fastq/cutandrun/
# Outputs: results/cutrun/ (log: results/cutrun/pipeline_info/cutandrun_run.log)
# Usage:   bash scripts/upstream/run_cutandrun.sh
# =============================================================================

set -euo pipefail

PROJECT_DIR="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
LOG_DIR="${PROJECT_DIR}/results/cutrun/pipeline_info"
mkdir -p "${LOG_DIR}"

cd "${PROJECT_DIR}/nextflow/launch_cutandrun"

nextflow run nf-core/cutandrun -r 3.2.2 \
    -c "${PROJECT_DIR}/nextflow/nextflow.config" \
    -profile cutandrun \
    --input "${PROJECT_DIR}/nextflow/samplesheets/cutandrun_samplesheet.csv" \
    --normalisation_mode CPM \
    --peakcaller 'macs2,seacr' \
    --macs_gsize 2494787188 \
    --save_reference \
    -resume \
    2>&1 | tee "${LOG_DIR}/cutandrun_run.log"

echo "Pipeline finished at $(date)"
