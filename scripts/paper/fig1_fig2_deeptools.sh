#!/usr/bin/env bash
# =============================================================================
# paper/fig1_fig2_deeptools.sh — deepTools signal matrices for Fig 1A/1B and 2D/2E/2F
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# Needs bigWigs from the nf-core runs (upstream mode); run after
# scripts/paper/fig1_published_clusters.R. Regions are the published Fig 1A
# clusters; for Conserved the 1,500-region subset drawn in the paper's heatmap
# is used (Conserved_plotted.bed), for every signal panel. Heatmaps are drawn by
# deepTools; profiles are exported as tables and drawn in R
# (scripts/paper/fig1_fig2_panels.R).
#
# Tracks (paper Fig 1A lists H3K27ac ChIP-seq and H3K4me1 ChIP-seq from
# GSE89036 for Naive; those external tracks are not reproduced. H3K27ac here is
# the deposited D8 WT CUT&RUN, averaged over TE/EEC/MP):
#   ATAC   Naive, 48h, D3, D5, D8 WT. 48h is shown as its own track, as in the
#          paper's "48h/d3" column: merged-replicate bigWig; QC Rep1 5.0M read
#          pairs FRiP 0.314, Rep2 4.5M FRiP 0.115, merged FRiP 0.244. 48h is
#          used in these signal panels only - it is excluded from the de novo
#          clustering and from every differential/statistical test.
#   ARID1A Naive, 48h, D5, D8 WT  (replicates / D8 subsets averaged)
#   H3K27ac D8 WT
# Fig 2F: H3K27ac CUT&RUN D5 WT vs KO.
#
# Inputs:  results/atac/bowtie2/merged_replicate/bigwig/*.bigWig
#          results/cutrun/03_peak_calling/03_bed_to_bigwig/*.bigWig
#          results/paper/fig1_regions/*.bed
# Outputs: figures/paper/fig1a_heatmap.{pdf,png}, figures/paper/fig2d_heatmap.{pdf,png}
#          results/paper/deeptools/ (matrices, *_profile.tab), results/paper/bigwig/
# Usage:   bash scripts/paper/fig1_fig2_deeptools.sh [threads]   (from the repository root)
# =============================================================================
set -euo pipefail

THREADS="${1:-16}"
ROOT="${ARID1A_PROJECT_DIR:-$(pwd)}"
cd "$ROOT"

ATAC=results/atac/bowtie2/merged_replicate/bigwig
CR=results/cutrun/03_peak_calling/03_bed_to_bigwig
CL=results/paper/fig1_regions
BW=results/paper/bigwig
MAT=results/paper/deeptools
FIG=figures/paper
mkdir -p "$BW" "$MAT" "$FIG"

for f in "$CL/Conserved_plotted.bed" "$ATAC/Naive_WT.mRp.clN.bigWig" "$CR/ARID1A_D5_WT_R1.bigWig"; do
  [[ -e "$f" ]] || { echo "Missing $f (run fig1_published_clusters.R / upstream pipelines first)"; exit 1; }
done

REGIONS=("$CL/Conserved_plotted.bed" "$CL/Naive.bed" "$CL/Early_Activation.bed" "$CL/Activation.bed" "$CL/Late_Activation.bed")
LABELS=(Conserved Naive "Early Activation" Activation "Late Activation")

avg() {  # avg <out> <in...>
  local out="$1"; shift
  [[ -s "$out" ]] && bigWigInfo "$out" > /dev/null 2>&1 && return 0
  if [[ $# -eq 1 ]]; then ln -sf "$(realpath "$1")" "$out"; return 0; fi
  bigwigAverage -b "$@" -o "$out" -bs 10 -p "$THREADS"
}

# --- averaged tracks ---------------------------------------------------------
avg "$BW/ATAC_D8_WT.bw"      $ATAC/D8_WT_{TE,EEC,MP}_Exp{1,2}.mRp.clN.bigWig
avg "$BW/ATAC_D8_KO.bw"      $ATAC/D8_KO_{TE,EEC,MP}_Exp{1,2}.mRp.clN.bigWig
avg "$BW/ATAC_D8_Het.bw"     $ATAC/D8_Het_{TE,EEC,MP}_Exp2.mRp.clN.bigWig
avg "$BW/ARID1A_Naive_WT.bw" $CR/ARID1A_Naive_WT_R1.bigWig
avg "$BW/ARID1A_48h_WT.bw"   $CR/ARID1A_48h_WT_R{1,2}.bigWig
avg "$BW/ARID1A_D5_WT.bw"    $CR/ARID1A_D5_WT_R{1,2}.bigWig
avg "$BW/ARID1A_D8_WT.bw"    $CR/ARID1A_D8_WT_{TE,EEC,MP}_R1.bigWig
avg "$BW/H3K27ac_D8_WT.bw"   $CR/H3K27ac_D8_WT_{TE,EEC,MP}_R1.bigWig

cm() {  # cm <name> <bigwigs...>  (reference-point, peak centre +/- 1 kb)
  local name="$1"; shift
  computeMatrix reference-point --referencePoint center -a 1000 -b 1000 -bs 20 \
    -R "${REGIONS[@]}" -S "$@" --skipZeros --missingDataAsZero \
    -p "$THREADS" -o "$MAT/$name.mat.gz" --quiet
}

# --- Fig 1A / 1B -------------------------------------------------------------
cm fig1a \
  $ATAC/Naive_WT.mRp.clN.bigWig $ATAC/48h_WT.mRp.clN.bigWig $ATAC/D3_WT.mRp.clN.bigWig \
  $ATAC/D5_WT.mRp.clN.bigWig "$BW/ATAC_D8_WT.bw" \
  "$BW/ARID1A_Naive_WT.bw" "$BW/ARID1A_48h_WT.bw" "$BW/ARID1A_D5_WT.bw" "$BW/ARID1A_D8_WT.bw" \
  "$BW/H3K27ac_D8_WT.bw"

plotHeatmap -m "$MAT/fig1a.mat.gz" -o "$FIG/fig1a_heatmap.png" --dpi 300 \
  --samplesLabel $'ATAC\nNaive' $'ATAC\n48h' $'ATAC\nd3' $'ATAC\nd5' $'ATAC\nd8' \
                 $'ARID1A\nNaive' $'ARID1A\n48h' $'ARID1A\nd5' $'ARID1A\nd8' $'H3K27ac\nd8' \
  --regionsLabel "${LABELS[@]}" \
  --colorMap Reds Reds Reds Reds Reds Greens Greens Greens Greens magma \
  --sortUsingSamples 1 3 4 5 --sortUsing mean \
  --zMin 0 --zMax 4 4 4 4 4 3 3 3 3 4 \
  --whatToShow "heatmap and colorbar" --heatmapHeight 18 --heatmapWidth 2.2 \
  --refPointLabel 0 --startLabel "-1" --endLabel "+1 kb" --xAxisLabel ""
plotHeatmap -m "$MAT/fig1a.mat.gz" -o "$FIG/fig1a_heatmap.pdf" \
  --samplesLabel $'ATAC\nNaive' $'ATAC\n48h' $'ATAC\nd3' $'ATAC\nd5' $'ATAC\nd8' \
                 $'ARID1A\nNaive' $'ARID1A\n48h' $'ARID1A\nd5' $'ARID1A\nd8' $'H3K27ac\nd8' \
  --regionsLabel "${LABELS[@]}" \
  --colorMap Reds Reds Reds Reds Reds Greens Greens Greens Greens magma \
  --sortUsingSamples 1 3 4 5 --sortUsing mean \
  --zMin 0 --zMax 4 4 4 4 4 3 3 3 3 4 \
  --whatToShow "heatmap and colorbar" --heatmapHeight 18 --heatmapWidth 2.2 \
  --refPointLabel 0 --startLabel "-1" --endLabel "+1 kb" --xAxisLabel ""
plotProfile -m "$MAT/fig1a.mat.gz" -o "$MAT/fig1a_profile_check.png" \
  --outFileNameData "$MAT/fig1a_profile.tab" --regionsLabel "${LABELS[@]}" \
  --samplesLabel ATAC_Naive ATAC_48h ATAC_D3 ATAC_D5 ATAC_D8 ARID1A_Naive ARID1A_48h ARID1A_D5 ARID1A_D8 H3K27ac_D8

# --- Fig 2D / 2E -------------------------------------------------------------
cm fig2d \
  $ATAC/Naive_WT.mRp.clN.bigWig \
  $ATAC/D3_WT.mRp.clN.bigWig $ATAC/D3_KO.mRp.clN.bigWig \
  $ATAC/D5_WT.mRp.clN.bigWig $ATAC/D5_KO.mRp.clN.bigWig \
  "$BW/ATAC_D8_WT.bw" "$BW/ATAC_D8_Het.bw" "$BW/ATAC_D8_KO.bw"

for ext in png pdf; do
  plotHeatmap -m "$MAT/fig2d.mat.gz" -o "$FIG/fig2d_heatmap.$ext" --dpi 300 \
    --samplesLabel "Naive WT" "d3 WT" "d3 KO" "d5 WT" "d5 KO" "d8 WT" "d8 Het" "d8 KO" \
    --regionsLabel "${LABELS[@]}" --colorMap Reds --zMin 0 --zMax 4 \
    --sortUsingSamples 2 4 6 --sortUsing mean \
    --whatToShow "heatmap and colorbar" --heatmapHeight 18 --heatmapWidth 2.2 \
    --refPointLabel 0 --startLabel "-1" --endLabel "+1 kb" --xAxisLabel ""
done
plotProfile -m "$MAT/fig2d.mat.gz" -o "$MAT/fig2d_profile_check.png" \
  --outFileNameData "$MAT/fig2e_profile.tab" --regionsLabel "${LABELS[@]}" \
  --samplesLabel Naive_WT D3_WT D3_KO D5_WT D5_KO D8_WT D8_Het D8_KO

# --- Fig 2F (H3K27ac D5 WT vs KO on all clustered ATAC peaks, and by cluster) -
cm fig2f $CR/H3K27ac_D5_WT_R1.bigWig $CR/H3K27ac_D5_KO_R1.bigWig
plotProfile -m "$MAT/fig2f.mat.gz" -o "$MAT/fig2f_profile_check.png" \
  --outFileNameData "$MAT/fig2f_profile_by_cluster.tab" --regionsLabel "${LABELS[@]}" \
  --samplesLabel H3K27ac_D5_WT H3K27ac_D5_KO
cat "${REGIONS[@]}" > "$MAT/all_clustered_ocrs.bed"
computeMatrix reference-point --referencePoint center -a 1000 -b 1000 -bs 20 \
  -R "$MAT/all_clustered_ocrs.bed" -S $CR/H3K27ac_D5_WT_R1.bigWig $CR/H3K27ac_D5_KO_R1.bigWig \
  --skipZeros --missingDataAsZero -p "$THREADS" -o "$MAT/fig2f_all.mat.gz" --quiet
plotProfile -m "$MAT/fig2f_all.mat.gz" -o "$MAT/fig2f_all_profile_check.png" \
  --outFileNameData "$MAT/fig2f_profile.tab" --regionsLabel "ATAC peaks" \
  --samplesLabel H3K27ac_D5_WT H3K27ac_D5_KO

echo "deepTools matrices in $MAT; heatmaps in $FIG"
