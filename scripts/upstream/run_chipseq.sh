#!/usr/bin/env bash
# =============================================================================
# upstream/run_chipseq.sh — nf-core/chipseq 2.1.0 on the T-bet ChIP-seq
# McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing
#
# 8 libraries (GSE228546): T-bet ChIP +/- ACBI1/BRM014 with matched inputs;
# single-end 100 bp. GRCm39 / GENCODE vM35, Bowtie2 (reusing the index built
# by the CUT&RUN run, so run_cutandrun.sh must run first), MACS3 narrow peaks.
# Requires: Nextflow and Docker (profiles in nextflow/nextflow.config).
#
# Inputs:  nextflow/samplesheets/chipseq_samplesheet.csv, data/fastq/chipseq/
#          results/cutrun/00_genome/index/bowtie2/
# Outputs: results/chipseq/ (log: results/chipseq/pipeline_info/chipseq_run.log)
# Usage:   bash scripts/upstream/run_chipseq.sh
# =============================================================================

set -euo pipefail

PROJECT_DIR="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
LOG_DIR="${PROJECT_DIR}/results/chipseq/pipeline_info"
mkdir -p "${LOG_DIR}"

cd "${PROJECT_DIR}/nextflow/launch_chipseq"

nextflow run nf-core/chipseq -r 2.1.0 \
    -c "${PROJECT_DIR}/nextflow/nextflow.config" \
    -profile chipseq \
    --input "${PROJECT_DIR}/nextflow/samplesheets/chipseq_samplesheet.csv" \
    --aligner bowtie2 \
    --bowtie2_index "${PROJECT_DIR}/results/cutrun/00_genome/index/bowtie2/" \
    --read_length 100 \
    --narrow_peak \
    --mito_name chrM \
    --save_reference \
    -resume \
    2>&1 | tee "${LOG_DIR}/chipseq_run.log"

echo "Pipeline finished at $(date)"
