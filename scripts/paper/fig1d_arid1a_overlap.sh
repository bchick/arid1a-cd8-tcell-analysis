#!/usr/bin/env bash
# =============================================================================
# paper/fig1d_arid1a_overlap.sh — ARID1A peaks per timepoint + HOMER obs/exp (Fig 1D)
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# Paper: overlap heatmaps from HOMER `mergePeaks -matrix` (observed/expected).
#
# nf-core/cutandrun only called ARID1A peaks at D5 (the only timepoint with IgG
# controls). Calling Naive/48h against the D5 WT IgG gave 580 and 6 peaks, so
# every timepoint is called here WITHOUT a control (MACS2 local lambda), with
# otherwise the nf-core/cutandrun MACS2 settings, pooling replicates / D8
# subsets. A FRiP table (reads in D5 WT ARID1A peaks) is written for QC.
# Requires: CUT&RUN BAMs (upstream mode), HOMER, MACS2, samtools; run after
# scripts/paper/fig1_published_clusters.R.
#
# Inputs:  results/cutrun/02_alignment/bowtie2/target/markdup/*.bam
#          results/cutrun/03_peak_calling/04_called_peaks/macs2/ARID1A_D5_WT_R1.macs2_peaks.narrowPeak
#          results/paper/fig1_regions/*.bed
# Outputs: results/paper/homer/arid1a_peaks_nocontrol/, results/paper/homer/overlap/
#          results/paper/fig1d_arid1a_frip.tsv
# Usage:   bash scripts/paper/fig1d_arid1a_overlap.sh [threads]   (from the repository root)
# =============================================================================
set -euo pipefail

THREADS="${1:-16}"
ROOT="${ARID1A_PROJECT_DIR:-$(pwd)}"
cd "$ROOT"

BAM=results/cutrun/02_alignment/bowtie2/target/markdup
NF_PEAKS=results/cutrun/03_peak_calling/04_called_peaks/macs2
CL=results/paper/fig1_regions
OUT=results/paper/homer
PK=$OUT/arid1a_peaks_nocontrol
GSIZE=2494787188
mkdir -p "$PK" "$OUT/overlap"

bams() { for s in "$@"; do echo "$BAM/$s.target.markdup.sorted.bam"; done; }

call() {  # call <name> <samples...>
  local name="$1"; shift
  [[ -s "$PK/${name}_peaks.narrowPeak" ]] && return 0
  macs2 callpeak --nomodel --shift -75 --extsize 150 --keep-dup all -q 0.01 \
    --gsize "$GSIZE" --format BAMPE -t $(bams "$@") \
    --name "$name" --outdir "$PK" 2> "$PK/$name.log"
}
call ARID1A_Naive ARID1A_Naive_WT_R1 &
call ARID1A_48h   ARID1A_48h_WT_R1 ARID1A_48h_WT_R2 &
call ARID1A_D5    ARID1A_D5_WT_R1 ARID1A_D5_WT_R2 &
call ARID1A_D8    ARID1A_D8_WT_TE_R1 ARID1A_D8_WT_EEC_R1 ARID1A_D8_WT_MP_R1 &
wait

# --- QC: fraction of read-1s in D5 WT ARID1A peaks (nf-core MACS2) -----------
REF_BED=$PK/d5_wt_reference_peaks.bed
cut -f1-3 "$NF_PEAKS/ARID1A_D5_WT_R1.macs2_peaks.narrowPeak" > "$REF_BED"
{
  echo -e "library\tread1_total\tread1_in_D5_ARID1A_peaks\tfrip\tn_peaks_nocontrol_call"
  for s in ARID1A_Naive_WT_R1 ARID1A_48h_WT_R1 ARID1A_48h_WT_R2 ARID1A_D5_WT_R1 ARID1A_D5_WT_R2 \
           ARID1A_D8_WT_TE_R1 ARID1A_D8_WT_EEC_R1 ARID1A_D8_WT_MP_R1 IgG_D5_WT_R1 IgG_D5_KO_R1; do
    f=$(bams "$s")
    tot=$(samtools view -@ 4 -c -F 0x904 -f 0x40 "$f")
    inp=$(samtools view -@ 4 -c -F 0x904 -f 0x40 -L "$REF_BED" "$f")
    tp=$(echo "$s" | sed -E 's/ARID1A_([^_]+)_.*/\1/')
    np=NA; [[ -s "$PK/ARID1A_${tp}_peaks.narrowPeak" && $s == ARID1A* ]] && np=$(wc -l < "$PK/ARID1A_${tp}_peaks.narrowPeak")
    awk -v s="$s" -v t="$tot" -v i="$inp" -v n="$np" 'BEGIN{printf "%s\t%d\t%d\t%.4f\t%s\n", s, t, i, i/t, n}'
  done
} > results/paper/fig1d_arid1a_frip.tsv
cat results/paper/fig1d_arid1a_frip.tsv

# --- HOMER mergePeaks -matrix ------------------------------------------------
# HOMER peak files: id chr start end strand
for c in Conserved Naive Early_Activation Activation Late_Activation; do
  awk 'BEGIN{OFS="\t"}{print $4,$1,$2+1,$3,"+"}' "$CL/$c.bed" > "$OUT/overlap/$c.txt"
done
for t in Naive 48h D5 D8; do
  awk -v t=$t 'BEGIN{OFS="\t"}{print t"_"NR,$1,$2+1,$3,"+"}' \
    "$PK/ARID1A_${t}_peaks.narrowPeak" > "$OUT/overlap/ARID1A_$t.txt"
done
( cd "$OUT/overlap" && mergePeaks -d given -gsize "$GSIZE" \
    Conserved.txt Naive.txt Early_Activation.txt Activation.txt Late_Activation.txt \
    ARID1A_Naive.txt ARID1A_48h.txt ARID1A_D5.txt ARID1A_D8.txt \
    -matrix fig1d > /dev/null 2> fig1d.log )
echo "Fig 1D overlap matrices in $OUT/overlap"
