#!/usr/bin/env bash
# =============================================================================
# upstream/run_meta_guo_atacseq.sh — nf-core/atacseq 2.1.2 on Guo et al. 2022
# McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing
#
# Arid1a-relevant ATAC-seq subset of Guo et al. 2022 for the cross-study
# meta-analysis. Same genome and settings as the main run (mm39, GENCODE vM35,
# Bowtie2, MACS2 narrow) so results are directly comparable. Reuses the Bowtie2
# index from the main ATAC run to avoid a slow, memory-hungry rebuild; extra
# resources in nextflow/meta_atac_resources.config.
# Requires: Nextflow and Docker (profiles in nextflow/nextflow.config).
#
# Inputs:  nextflow/samplesheets/meta_guo_atacseq.csv, results/atac/genome/index/bowtie2
# Outputs: results/meta_analysis/guo2022/atacseq/
# Usage:   bash scripts/upstream/run_meta_guo_atacseq.sh
# =============================================================================
set -uo pipefail

PROJECT_DIR="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
OUTDIR="${PROJECT_DIR}/results/meta_analysis/guo2022/atacseq"
LOGDIR="${OUTDIR}/pipeline_info"
LAUNCH="${PROJECT_DIR}/nextflow/launch_meta_guo_atac"
mkdir -p "${LOGDIR}" "${LAUNCH}"
cd "${LAUNCH}" || exit 1

nextflow run nf-core/atacseq -r 2.1.2 \
    -c "${PROJECT_DIR}/nextflow/nextflow.config" \
    -c "${PROJECT_DIR}/nextflow/meta_atac_resources.config" \
    -w "${PROJECT_DIR}/work_meta_guo_atac" \
    --input "${PROJECT_DIR}/nextflow/samplesheets/meta_guo_atacseq.csv" \
    --outdir "${OUTDIR}" \
    --aligner bowtie2 \
    --bowtie2_index "${PROJECT_DIR}/results/atac/genome/index/bowtie2" \
    --narrow_peak \
    --mito_name chrM \
    --macs_gsize 2407883318 \
    --read_length 50 \
    -resume \
    2>&1 | tee "${LOGDIR}/meta_guo_atacseq_run.log"
rc=${PIPESTATUS[0]}
echo "Guo ATAC pipeline finished at $(date) (rc=${rc})"
exit ${rc}
