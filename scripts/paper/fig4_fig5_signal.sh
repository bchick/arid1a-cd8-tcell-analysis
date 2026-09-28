#!/usr/bin/env bash
# =============================================================================
# paper/fig4_fig5_signal.sh — deepTools signal heatmaps for Figures 4C, 5C, 5D, 5H, 5J
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# Region sets come from scripts/paper/fig4.R and fig5.R (results/paper/).
# Writes one computeMatrix per panel and the heatmaps. Average profiles
# (4C, 5C, 5E, 5J) are drawn from the same matrices by
# scripts/paper/fig4_fig5_profiles.R.
#
#   4C  ATAC, day 8 TE / EEC / MP x WT / Het / KO, at OCRs lost in one subset or all
#   5C  CUT&RUN ARID1A / BATF / ETS1 / T-bet, d5 WT vs Arid1a KO,
#       at ARID1A-dependent vs -independent OCRs
#   5D  in vitro ATAC and T-bet ChIP (DMSO, DMSO+IL-12, ACBI1+IL-12, BRM014+IL-12)
#       at ACBI1-dependent vs -independent T-bet-bound OCRs
#   5H  ARID1A CUT&RUN, d5 WT vs Tbx21 KO, at the Fig 5G TE / MP OCR groups
#   5J  T-bet CUT&RUN, d5 WT / Arid1a KO +/- T-bet overexpression, at the
#       published Fig 1A Activation / Late Activation OCR clusters (profile only)
#
# Day 8 ATAC uses the Exp2 merged-replicate bigWigs: Exp2 is the batch in which
# WT, Het and KO were all profiled.
# Requires: deepTools >= 3.5; bigWigs from the nf-core pipelines (upstream mode).
#
# Inputs:  results/paper/fig4c_*.bed, results/paper/fig5_beds/*.bed,
#          data/metadata/paper_ocr_clusters/mm39/
#          results/{atac,cutrun,chipseq}/... bigWigs
# Outputs: results/paper/signal/*.gz, figures/paper/<panel>.{pdf,png}
# Usage:   bash scripts/paper/fig4_fig5_signal.sh [threads]   (from the repository root)
# =============================================================================

set -euo pipefail

PROJECT="${ARID1A_PROJECT_DIR:-$(pwd)}"
NCPU="${1:-8}"
TAB="${PROJECT}/results/paper"
BEDS="${TAB}/fig5_beds"
MAT="${TAB}/signal"
FIG="${PROJECT}/figures/paper"
ATAC_BW="${PROJECT}/results/atac/bowtie2/merged_replicate/bigwig"
CR_BW="${PROJECT}/results/cutrun/03_peak_calling/03_bed_to_bigwig"
CHIP_BW="${PROJECT}/results/chipseq/bowtie2/merged_library/bigwig"

mkdir -p "${MAT}" "${FIG}"

for d in "${ATAC_BW}" "${CR_BW}" "${CHIP_BW}"; do
  [ -d "${d}" ] || { echo "Skipping: bigWig directory not found (upstream mode): ${d}" >&2; exit 0; }
done

# matrix <name> <regions...> -- <bigwigs...>
# Sample labels come from the LABELS array and are stored in the matrix header,
# which fig4_fig5_profiles.R reads.
matrix() {
  local name="$1"; shift
  local regions=() scores=()
  while [ "$1" != "--" ]; do regions+=("$1"); shift; done; shift
  scores=("$@")
  if [ -s "${MAT}/${name}.gz" ]; then echo "  ${name}: matrix exists"; return; fi
  echo "  ${name}: computeMatrix (${#regions[@]} region sets x ${#scores[@]} tracks)"
  computeMatrix reference-point --referencePoint center -a 1000 -b 1000 --binSize 20 \
    -R "${regions[@]}" -S "${scores[@]}" --samplesLabel "${LABELS[@]}" \
    --missingDataAsZero --skipZeros -p "${NCPU}" --quiet \
    -o "${MAT}/${name}.gz"
}

heatmap() {  # heatmap <name> <colormap>  (region labels via REGION_LABELS)
  local name="$1" cmap="$2"
  plotHeatmap -m "${MAT}/${name}.gz" -o "${FIG}/${name}.pdf" \
    --colorMap "${cmap}" --whatToShow "heatmap and colorbar" \
    --regionsLabel "${REGION_LABELS[@]}" \
    --sortUsing mean --sortUsingSamples 1 \
    --heatmapHeight 14 --heatmapWidth 2 --zMin 0 \
    --refPointLabel 0 --startLabel "-1" --endLabel "+1 kb" --xAxisLabel "" --dpi 300
  plotHeatmap -m "${MAT}/${name}.gz" -o "${FIG}/${name}.png" \
    --colorMap "${cmap}" --whatToShow "heatmap and colorbar" \
    --regionsLabel "${REGION_LABELS[@]}" \
    --sortUsing mean --sortUsingSamples 1 \
    --heatmapHeight 14 --heatmapWidth 2 --zMin 0 \
    --refPointLabel 0 --startLabel "-1" --endLabel "+1 kb" --xAxisLabel "" --dpi 300
  echo "  ${name}: heatmap written"
}

# -----------------------------------------------------------------------------
# 4C
# -----------------------------------------------------------------------------
echo "=== 4C ==="
S4C=(); L4C=()
for s in TE EEC MP; do
  for g in WT Het KO; do
    S4C+=("${ATAC_BW}/D8_${g}_${s}_Exp2.mRp.clN.bigWig"); L4C+=("${s}_${g}")
  done
done
LABELS=("${L4C[@]}")
matrix fig4c_atac_heatmap \
  "${TAB}/fig4c_Down_in_EEC.bed" "${TAB}/fig4c_Down_in_MP.bed" \
  "${TAB}/fig4c_Down_in_TE.bed"  "${TAB}/fig4c_Down_in_all.bed" -- "${S4C[@]}"
REGION_LABELS=("EEC" "MP" "TE" "All")
heatmap fig4c_atac_heatmap OrRd

# -----------------------------------------------------------------------------
# 5C
# -----------------------------------------------------------------------------
echo "=== 5C ==="
S5C=(); L5C=()
for ab in ARID1A BATF ETS1 Tbet; do
  for g in WT KO; do S5C+=("${CR_BW}/${ab}_D5_${g}_R1.bigWig"); L5C+=("${ab}_${g}"); done
done
LABELS=("${L5C[@]}")
matrix fig5c_cutrun_heatmap "${BEDS}/fig5c_arid1a_dependent.bed" "${BEDS}/fig5c_arid1a_independent.bed" -- "${S5C[@]}"
REGION_LABELS=("Dependent" "Independent")
heatmap fig5c_cutrun_heatmap Greens

# -----------------------------------------------------------------------------
# 5D (heatmaps) — also used for the 5E histograms
# -----------------------------------------------------------------------------
echo "=== 5D ==="
TRT=(Untreated IL-12 IL-12_ACBI1 IL-12_BRM014)
TRT_LAB=(DMSO DMSO+IL-12 ACBI1+IL-12 BRM014+IL-12)
S5D_ATAC=(); S5D_TBET=()
for t in "${TRT[@]}"; do
  S5D_ATAC+=("${ATAC_BW}/ATAC_${t}.mRp.clN.bigWig")
  S5D_TBET+=("${CHIP_BW}/ChIP_${t}_T-bet_REP1.mLb.clN.bigWig")
done
REGION_LABELS=("Dependent" "Independent")
LABELS=("${TRT_LAB[@]}")
matrix fig5d_atac_heatmap "${BEDS}/fig5d_acbi1_dependent.bed" "${BEDS}/fig5d_acbi1_independent.bed" -- "${S5D_ATAC[@]}"
heatmap fig5d_atac_heatmap OrRd
matrix fig5d_tbet_chip_heatmap "${BEDS}/fig5d_acbi1_dependent.bed" "${BEDS}/fig5d_acbi1_independent.bed" -- "${S5D_TBET[@]}"
heatmap fig5d_tbet_chip_heatmap Blues

# -----------------------------------------------------------------------------
# 5H
# -----------------------------------------------------------------------------
echo "=== 5H ==="
S5H=("${CR_BW}/ARID1A_D5_WT_R1.bigWig" "${CR_BW}/ARID1A_D5_TbetKO_R1.bigWig")
REGION_LABELS=("Both" "Tbx21 KO" "Arid1a KO")
LABELS=("WT" "Tbx21 KO")
for s in TE MP; do
  matrix "fig5h_${s}_arid1a_cutrun_heatmap" \
    "${BEDS}/fig5h_${s}_lost_both.bed" "${BEDS}/fig5h_${s}_lost_tbx21ko_only.bed" \
    "${BEDS}/fig5h_${s}_lost_arid1ako_only.bed" -- "${S5H[@]}"
  heatmap "fig5h_${s}_arid1a_cutrun_heatmap" Greens
done

# -----------------------------------------------------------------------------
# 5J (profiles only): the published Fig 1A Activation / Late Activation OCRs
# (original mm10 cluster BEDs lifted to mm39, data/metadata/paper_ocr_clusters/mm39)
# -----------------------------------------------------------------------------
echo "=== 5J ==="
PAPER_CL="${PROJECT}/data/metadata/paper_ocr_clusters/mm39"
if [ -s "${PAPER_CL}/activation.specific.sig.bed" ]; then
  LABELS=(Tbet_WT_EV Tbet_KO_EV Tbet_WT_TbetOE Tbet_KO_TbetOE)
  matrix fig5j_tbet_cutrun \
    "${PAPER_CL}/activation.specific.sig.bed" "${PAPER_CL}/late.activation.specific.sig.bed" -- \
    "${CR_BW}/Tbet_D5_WT_R1.bigWig" "${CR_BW}/Tbet_D5_KO_R1.bigWig" \
    "${CR_BW}/Tbet_D5_WT_TbetOE_R1.bigWig" "${CR_BW}/Tbet_D5_KO_TbetOE_R1.bigWig"
else
  echo "  published cluster BEDs missing: ${PAPER_CL}"
fi

echo "Done. Profiles: Rscript scripts/paper/fig4_fig5_profiles.R"
