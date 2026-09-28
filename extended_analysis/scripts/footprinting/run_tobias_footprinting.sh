#!/usr/bin/env bash
# =============================================================================
# footprinting/run_tobias_footprinting.sh — TOBIAS differential footprinting (D8 WT/Het/KO)
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Footprint-based counterpart of dose_chromvar, independent of HOMER and chromVAR.
#   D8 ATAC, Exp2 batch only (WT=9, Het=9, KO=6 BAMs; TE/EEC/MP subsets pooled)
#   — Het exists only in Exp2, so WT/KO are also restricted to Exp2 to balance
#   batch instead of confounding the WT->Het->KO dose contrast.
# Pipeline: merge BAMs per genotype -> ATACorrect (Tn5 bias) -> ScoreBigwig
#           (footprint) -> BINDetect (differential, JASPAR2020 CORE vertebrates).
# Idempotent: each stage skips if its output already exists. Needs the raw
# nf-core/atacseq BAMs, the mm39 genome FASTA, samtools and TOBIAS on PATH, and
# the motif file from footprinting/export_jaspar_motifs.R.
#
# Inputs:  results/atac/bowtie2/merged_library/D8_{WT,Het,KO}_*_Exp2_*.mLb.clN.sorted.bam
#          results/atac/bowtie2/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.bed
#          data/reference/GRCm39.primary_assembly.genome.fa
#          results/extended_analysis/footprinting/motifs/jaspar2020_core_vertebrates.jaspar
# Outputs: results/extended_analysis/footprinting/{peaks,bam,atacorrect,footprints}/
#          results/extended_analysis/footprinting/bindetect/ (incl. bindetect_results.txt)
# Usage:   bash extended_analysis/scripts/footprinting/run_tobias_footprinting.sh   (from the repository root)
# =============================================================================
set -euo pipefail

PROJ="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
BAMDIR=$PROJ/results/atac/bowtie2/merged_library
GENOME=$PROJ/data/reference/GRCm39.primary_assembly.genome.fa
MOTIFS=$PROJ/results/extended_analysis/footprinting/motifs/jaspar2020_core_vertebrates.jaspar
OUT=$PROJ/results/extended_analysis/footprinting
CORES=16
CONDS=(WT Het KO)

mkdir -p "$OUT"/{peaks,bam,atacorrect,footprints}
log(){ echo "[$(date '+%F %T')] $*"; }

# -- 1. peaks: standard chromosomes only (drop scaffolds) ---------------------
PEAKS=$OUT/peaks/consensus_D8_stdchrom.bed
SRC=$BAMDIR/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.bed
if [[ ! -s $PEAKS ]]; then
  log "filtering consensus peaks to chr1-19,X,Y"
  awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|X|Y)$/ {print $1,$2,$3}' "$SRC" \
    | sort -k1,1 -k2,2n > "$PEAKS"
fi
log "peaks: $(wc -l < "$PEAKS") regions"

# -- 2. merge Exp2 BAMs per genotype (subsets pooled) -------------------------
for g in "${CONDS[@]}"; do
  BAM=$OUT/bam/${g}.bam
  if [[ ! -s $BAM ]]; then
    mapfile -t bams < <(ls "$BAMDIR"/D8_${g}_*_Exp2_*.mLb.clN.sorted.bam)
    log "merging ${#bams[@]} BAMs -> ${g}.bam"
    samtools merge -f -@ "$CORES" "$BAM" "${bams[@]}"
    samtools index -@ "$CORES" "$BAM"
  fi
done

# -- 3. ATACorrect (Tn5 bias correction) --------------------------------------
for g in "${CONDS[@]}"; do
  COR=$OUT/atacorrect/${g}/${g}_corrected.bw
  if [[ ! -s $COR ]]; then
    log "ATACorrect $g"
    TOBIAS ATACorrect --bam "$OUT/bam/${g}.bam" --genome "$GENOME" \
      --peaks "$PEAKS" --prefix "$g" --outdir "$OUT/atacorrect/${g}" \
      --cores "$CORES"
  fi
done

# -- 4. ScoreBigwig (footprint score) -----------------------------------------
for g in "${CONDS[@]}"; do
  FP=$OUT/footprints/${g}_footprints.bw
  if [[ ! -s $FP ]]; then
    log "ScoreBigwig $g"
    TOBIAS ScoreBigwig --signal "$OUT/atacorrect/${g}/${g}_corrected.bw" \
      --regions "$PEAKS" --output "$FP" --score footprint --cores "$CORES"
  fi
done

# -- 5. BINDetect (differential TF footprinting) ------------------------------
if [[ ! -s $OUT/bindetect/bindetect_results.txt ]]; then
  log "BINDetect WT/Het/KO"
  TOBIAS BINDetect --motifs "$MOTIFS" \
    --signals "$OUT/footprints/WT_footprints.bw" \
              "$OUT/footprints/Het_footprints.bw" \
              "$OUT/footprints/KO_footprints.bw" \
    --genome "$GENOME" --peaks "$PEAKS" \
    --cond-names WT Het KO \
    --outdir "$OUT/bindetect" --cores "$CORES"
fi

log "DONE — results in $OUT/bindetect/"
