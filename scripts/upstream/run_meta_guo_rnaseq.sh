#!/usr/bin/env bash
# =============================================================================
# upstream/run_meta_guo_rnaseq.sh — nf-core/rnaseq 3.22.2 on Guo et al. 2022
# McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing
#
# Arid1a-relevant RNA-seq subset of Guo et al. 2022 for the cross-study
# meta-analysis. STAR + Salmon on mm39 / GENCODE vM35, as in the main run.
# Guo RNA-seq is paired-end (ours is single-end); strandedness is inferred per
# sample. Reuses the pre-built STAR index. Do not pass --transcript_fasta
# (transcript ID mismatch seen in the main run); the pipeline derives the
# transcriptome from genome + GTF.
# Requires: Nextflow and Docker (profiles in nextflow/nextflow.config).
#
# Inputs:  nextflow/samplesheets/meta_guo_rnaseq.csv, results/rnaseq/genome/index/star
# Outputs: results/meta_analysis/guo2022/rnaseq/
# Usage:   bash scripts/upstream/run_meta_guo_rnaseq.sh
# =============================================================================
set -uo pipefail

PROJECT_DIR="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
OUTDIR="${PROJECT_DIR}/results/meta_analysis/guo2022/rnaseq"
LOGDIR="${OUTDIR}/pipeline_info"
LAUNCH="${PROJECT_DIR}/nextflow/launch_meta_guo_rna"
mkdir -p "${LOGDIR}" "${LAUNCH}"
cd "${LAUNCH}" || exit 1

nextflow run nf-core/rnaseq -r 3.22.2 \
    -c "${PROJECT_DIR}/nextflow/nextflow.config" \
    -w "${PROJECT_DIR}/work_meta_guo_rna" \
    --input "${PROJECT_DIR}/nextflow/samplesheets/meta_guo_rnaseq.csv" \
    --outdir "${OUTDIR}" \
    --aligner star_salmon \
    --star_index "${PROJECT_DIR}/results/rnaseq/genome/index/star" \
    --save_reference \
    -resume \
    2>&1 | tee "${LOGDIR}/meta_guo_rnaseq_run.log"
rc=${PIPESTATUS[0]}
echo "Guo RNA pipeline finished at $(date) (rc=${rc})"
exit ${rc}
