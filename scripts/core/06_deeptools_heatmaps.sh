#!/usr/bin/env bash
# =============================================================================
# core/06_deeptools_heatmaps.sh — deepTools signal heatmaps
# McDonald, Chick et al. 2023 Immunity 56:1303 — core analysis
#
# computeMatrix + plotHeatmap for ATAC-seq (WT time course at consensus peaks;
# WT vs KO per timepoint), CUT&RUN (WT vs KO at ARID1A-dependent vs
# -independent sites) and T-bet ChIP-seq (+/- BAF inhibitors). Each section is
# skipped if its bigWigs or region files are absent.
# Requires: deepTools, samtools; bigWigs from the nf-core pipelines.
#
# Inputs:  results/atac/bowtie2/merged_replicate/{bigwig,macs2}/
#          results/cutrun/04_reporting/igv/*.bigWig, results/chipseq/
#          results/atac/differential/arid1a_{dependent,independent}_peaks.bed
# Outputs: figures/deeptools/*_heatmap.pdf
#          results/deeptools_matrices/*.gz
# Usage:   bash scripts/core/06_deeptools_heatmaps.sh   (from any directory)
# =============================================================================

set -euo pipefail

PROJECT="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
FIGDIR="${PROJECT}/figures/deeptools"
TMPDIR="${PROJECT}/results/deeptools_matrices"
BLACKLIST="${PROJECT}/data/reference/mm39-blacklist.v2.bed"
ATAC_DIR="${PROJECT}/results/atac/bowtie2/merged_replicate"
CUTRUN_DIR="${PROJECT}/results/cutrun"
CHIPSEQ_DIR="${PROJECT}/results/chipseq"
DA_DIR="${PROJECT}/results/atac/differential"

mkdir -p "${FIGDIR}" "${TMPDIR}"

# Number of processors
NCPU=16

# =============================================================================
# 1. ATAC-seq signal at OCR clusters (Figure 1 style)
# =============================================================================

echo "=== ATAC-seq heatmaps at OCR clusters ==="

# bigWig files from nf-core/atacseq (mergedReplicate)
ATAC_BW_DIR="${ATAC_DIR}/bigwig"

# Compute matrix centered on consensus peaks
if [ -d "${ATAC_BW_DIR}" ]; then
  # Get WT bigWigs for each timepoint
  WT_NAIVE_BW=$(ls ${ATAC_BW_DIR}/Naive_WT*.bigWig 2>/dev/null | head -1)
  WT_D3_BW=$(ls ${ATAC_BW_DIR}/D3_WT*.bigWig 2>/dev/null | head -1)
  WT_D5_BW=$(ls ${ATAC_BW_DIR}/D5_WT*.bigWig 2>/dev/null | head -1)
  WT_D8_BW=$(ls ${ATAC_BW_DIR}/D8_WT_TE*.bigWig 2>/dev/null | head -1)

  CONSENSUS_BED="${PROJECT}/results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.bed"

  if [ -n "${WT_D3_BW}" ]; then
    echo "  Computing matrix for WT ATAC-seq across timepoints..."

    computeMatrix reference-point \
      -S ${WT_NAIVE_BW} ${WT_D3_BW} ${WT_D5_BW} ${WT_D8_BW} \
      -R "${CONSENSUS_BED}" \
      --referencePoint center \
      -a 2000 -b 2000 \
      --binSize 50 \
      --missingDataAsZero \
      --blackListFileName "${BLACKLIST}" \
      -p ${NCPU} \
      -o "${TMPDIR}/atac_wt_timepoints.gz"

    plotHeatmap \
      -m "${TMPDIR}/atac_wt_timepoints.gz" \
      -o "${FIGDIR}/atac_wt_timepoints_heatmap.pdf" \
      --colorMap viridis \
      --samplesLabel "Naive" "D3" "D5" "D8_TE" \
      --heatmapHeight 15 \
      --heatmapWidth 4 \
      --whatToShow "heatmap and colorbar" \
      --zMin 0 \
      --refPointLabel "Peak center"

    echo "  Done: ${FIGDIR}/atac_wt_timepoints_heatmap.pdf"
  fi
fi

# =============================================================================
# 2. ATAC-seq: WT vs KO comparison heatmap (Figure 2 style)
# =============================================================================

echo "=== ATAC-seq WT vs KO heatmaps ==="

for TP in D3 D5 D8; do
  WT_BW=$(ls ${ATAC_BW_DIR}/${TP}_WT*.bigWig 2>/dev/null | head -1)
  KO_BW=$(ls ${ATAC_BW_DIR}/${TP}_KO*.bigWig 2>/dev/null | head -1)

  if [ -n "${WT_BW}" ] && [ -n "${KO_BW}" ]; then
    echo "  ${TP}: WT vs KO..."

    computeMatrix reference-point \
      -S "${WT_BW}" "${KO_BW}" \
      -R "${CONSENSUS_BED}" \
      --referencePoint center \
      -a 2000 -b 2000 \
      --binSize 50 \
      --missingDataAsZero \
      --blackListFileName "${BLACKLIST}" \
      -p ${NCPU} \
      -o "${TMPDIR}/atac_${TP}_wt_vs_ko.gz"

    plotHeatmap \
      -m "${TMPDIR}/atac_${TP}_wt_vs_ko.gz" \
      -o "${FIGDIR}/atac_${TP}_wt_vs_ko_heatmap.pdf" \
      --colorList 'white,#FCA082,#D63B20,#7F0000' \
      --samplesLabel "WT" "KO" \
      --heatmapHeight 15 \
      --heatmapWidth 3 \
      --zMin 0 \
      --refPointLabel "Peak center"
  fi
done

# =============================================================================
# 3. CUT&RUN: ARID1A binding at dependent vs independent sites (Figure 5 style)
# =============================================================================

echo "=== CUT&RUN heatmaps ==="

CUTRUN_BW_DIR="${CUTRUN_DIR}/04_reporting/igv"
ARID1A_DEP="${DA_DIR}/arid1a_dependent_peaks.bed"
ARID1A_INDEP="${DA_DIR}/arid1a_independent_peaks.bed"

if [ -d "${CUTRUN_BW_DIR}" ] && [ -f "${ARID1A_DEP}" ]; then
  for AB in ARID1A Tbet BATF ETS1 H3K27ac; do
    WT_BW=$(ls ${CUTRUN_BW_DIR}/*${AB}*WT*.bigWig 2>/dev/null | head -1)
    KO_BW=$(ls ${CUTRUN_BW_DIR}/*${AB}*KO*.bigWig 2>/dev/null | head -1)

    if [ -n "${WT_BW}" ] && [ -n "${KO_BW}" ]; then
      echo "  ${AB}: WT vs KO at ARID1A-dependent sites..."

      computeMatrix reference-point \
        -S "${WT_BW}" "${KO_BW}" \
        -R "${ARID1A_DEP}" "${ARID1A_INDEP}" \
        --referencePoint center \
        -a 2000 -b 2000 \
        --binSize 50 \
        --missingDataAsZero \
        --blackListFileName "${BLACKLIST}" \
        -p ${NCPU} \
        -o "${TMPDIR}/cutrun_${AB}_wt_vs_ko.gz"

      plotHeatmap \
        -m "${TMPDIR}/cutrun_${AB}_wt_vs_ko.gz" \
        -o "${FIGDIR}/cutrun_${AB}_wt_vs_ko_heatmap.pdf" \
        --colorList 'white,#C7E9C0,#41AB5D,#005A32' \
        --samplesLabel "WT" "KO" \
        --regionsLabel "ARID1A-dependent" "ARID1A-independent" \
        --heatmapHeight 12 \
        --heatmapWidth 3 \
        --zMin 0 \
        --refPointLabel "Peak center"

      echo "  Done: ${FIGDIR}/cutrun_${AB}_wt_vs_ko_heatmap.pdf"
    fi
  done
fi

# =============================================================================
# 4. ChIP-seq: T-bet binding ± BAF inhibitors (Figure 5D style)
# =============================================================================

echo "=== ChIP-seq T-bet heatmaps ==="

CHIPSEQ_BW_DIR="${CHIPSEQ_DIR}/bowtie2/merged_library/bigwig"

if [ -d "${CHIPSEQ_BW_DIR}" ] && [ -f "${ARID1A_DEP}" ]; then
  DMSO_BW=$(ls ${CHIPSEQ_BW_DIR}/*Untreated*T-bet*.bigWig 2>/dev/null | head -1)
  ACBI1_BW=$(ls ${CHIPSEQ_BW_DIR}/*ACBI1*.bigWig 2>/dev/null | head -1)
  BRM014_BW=$(ls ${CHIPSEQ_BW_DIR}/*BRM014*.bigWig 2>/dev/null | head -1)

  if [ -n "${DMSO_BW}" ] && [ -n "${ACBI1_BW}" ]; then
    echo "  T-bet ChIP: DMSO vs ACBI1 vs BRM014..."

    computeMatrix reference-point \
      -S "${DMSO_BW}" "${ACBI1_BW}" "${BRM014_BW}" \
      -R "${ARID1A_DEP}" "${ARID1A_INDEP}" \
      --referencePoint center \
      -a 2000 -b 2000 \
      --binSize 50 \
      --missingDataAsZero \
      --blackListFileName "${BLACKLIST}" \
      -p ${NCPU} \
      -o "${TMPDIR}/chipseq_tbet_inhibitors.gz"

    plotHeatmap \
      -m "${TMPDIR}/chipseq_tbet_inhibitors.gz" \
      -o "${FIGDIR}/chipseq_tbet_inhibitors_heatmap.pdf" \
      --colorList 'white,#FDE0DD,#FA9FB5,#C51B8A' \
      --samplesLabel "DMSO" "ACBI1" "BRM014" \
      --regionsLabel "ARID1A-dependent" "ARID1A-independent" \
      --heatmapHeight 12 \
      --heatmapWidth 3 \
      --zMin 0 \
      --refPointLabel "Peak center"

    echo "  Done: ${FIGDIR}/chipseq_tbet_inhibitors_heatmap.pdf"
  fi
fi

echo "=== deepTools heatmaps complete ==="
