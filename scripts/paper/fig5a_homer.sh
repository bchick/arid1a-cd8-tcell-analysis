#!/usr/bin/env bash
# =============================================================================
# paper/fig5a_homer.sh — Figure 5A: HOMER known motifs in OCRs lost / gained in Arid1a KO
# McDonald, Chick et al. 2023 Immunity 56:1303 — paper panel reproduction
#
# Sets: day 5 gained, and day 8 TE / EEC / MP lost and gained.
# STAR Methods: sequences within 100 bp of peak centres (-size 200) compared with
# HOMER known motifs, GC-matched random genomic background (HOMER default).
# Plot with scripts/paper/fig5a_motifs.R.
# Requires: HOMER (findMotifsGenome.pl), soft-masked GRCm39 FASTA.
#
# Inputs:  results/paper/fig5_beds/*.bed (scripts/paper/fig5.R)
#          data/reference/GRCm39.primary_assembly.genome.fa
# Outputs: results/paper/fig5a_homer/<set>/knownResults.txt
# Usage:   bash scripts/paper/fig5a_homer.sh [threads]   (from the repository root)
# =============================================================================

set -euo pipefail

PROJECT="${ARID1A_PROJECT_DIR:-$(pwd)}"
THREADS="${1:-16}"
GENOME="${PROJECT}/data/reference/GRCm39.primary_assembly.genome.fa"
BEDS="${PROJECT}/results/paper/fig5_beds"
OUT="${PROJECT}/results/paper/fig5a_homer"
PREPARSED="${OUT}/preparsed"

[ -f "${GENOME}" ] || { echo "Genome FASTA not found: ${GENOME}" >&2; exit 1; }
[ -d "${BEDS}" ]   || { echo "Run scripts/paper/fig5.R first (${BEDS} missing)" >&2; exit 1; }
mkdir -p "${OUT}" "${PREPARSED}"

# Two runs at a time, THREADS/2 cores each
PER_RUN=$(( THREADS / 2 > 0 ? THREADS / 2 : 1 ))

run_homer() {
  local set="$1"
  local dir="${OUT}/${set}"
  if [ -s "${dir}/knownResults.txt" ]; then
    echo "  ${set}: done (skipping)"; return
  fi
  mkdir -p "${dir}"
  echo "  ${set}: $(wc -l < "${BEDS}/fig5a_${set}.bed") OCRs"
  findMotifsGenome.pl "${BEDS}/fig5a_${set}.bed" "${GENOME}" "${dir}" \
    -size 200 -mask -nomotif -p "${PER_RUN}" -preparsedDir "${PREPARSED}" \
    > "${dir}/homer.log" 2>&1
}

# The first run builds the preparsed background for -size 200; run it alone so
# parallel jobs do not race on writing the same preparsed files.
run_homer d5_lost

SETS=(d5_gained TE_lost TE_gained EEC_lost EEC_gained MP_lost MP_gained)
for ((i = 0; i < ${#SETS[@]}; i += 2)); do
  run_homer "${SETS[i]}" &
  [ -n "${SETS[i+1]:-}" ] && run_homer "${SETS[i+1]}" &
  wait
done

echo "HOMER known-motif results in ${OUT}/*/knownResults.txt"
