#!/usr/bin/env bash
# =============================================================================
# paper/fig1_homer.sh — HOMER annotation (Fig 1C) and known motifs (Fig 1F, S1C)
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# Runs on the published Fig 1A cluster regions (results/paper/fig1_regions/,
# written by fig1_published_clusters.R; Conserved = the full 21,649 set).
# Paper methods: peaks annotated with HOMER; motifs with findMotifsGenome.pl on
# +/-100 bp of the peak centre against GC-matched random genomic background,
# known HOMER motifs. mm39 has no HOMER genome package, so the GRCm39 FASTA and
# GENCODE vM35 GTF are given explicitly.
# Requires: HOMER.
#
# Inputs:  results/paper/fig1_regions/*.bed
#          data/reference/GRCm39.primary_assembly.genome.fa
#          data/reference/gencode.vM35.primary_assembly.annotation.gtf
# Outputs: results/paper/homer/{annotation,motifs,preparsed}/
# Usage:   bash scripts/paper/fig1_homer.sh [threads]   (from the repository root)
# =============================================================================
set -euo pipefail

THREADS="${1:-16}"
ROOT="${ARID1A_PROJECT_DIR:-$(pwd)}"
cd "$ROOT"

REG=results/paper/fig1_regions
OUT=results/paper/homer
FASTA=data/reference/GRCm39.primary_assembly.genome.fa
GTF=data/reference/gencode.vM35.primary_assembly.annotation.gtf
CLUSTERS=(Conserved Naive Early_Activation Activation Late_Activation)
mkdir -p "$OUT/motifs" "$OUT/annotation" "$OUT/preparsed"

[[ -s "$REG/Conserved.bed" ]] || { echo "Run scripts/paper/fig1_published_clusters.R first"; exit 1; }

# --- Fig 1C: HOMER genomic annotation ---------------------------------------
for c in "${CLUSTERS[@]}"; do
  f="$OUT/annotation/$c.annotatePeaks.txt"
  [[ -s "$f" ]] && continue
  annotatePeaks.pl "$REG/$c.bed" "$FASTA" -gtf "$GTF" -cpu 4 > "$f" 2> "$OUT/annotation/$c.log" &
done
wait

# --- Fig 1F / S1C: known motifs ----------------------------------------------
for c in "${CLUSTERS[@]}"; do
  d="$OUT/motifs/$c"
  [[ -s "$d/knownResults.txt" ]] && continue
  findMotifsGenome.pl "$REG/$c.bed" "$FASTA" "$d" -size 200 -nomotif \
    -p "$THREADS" -preparsedDir "$OUT/preparsed" > "$OUT/motifs/$c.log" 2>&1
done
echo "HOMER outputs in $OUT/{annotation,motifs}"
