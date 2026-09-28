#!/usr/bin/env bash
# =============================================================================
# atac_trajectories/04_lost_site_motifs.sh — known-motif enrichment of sites lost in ARID1A KO
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Runs MEME-suite SEA for every (primary, control) pair written by
# atac_trajectories/04_lost_site_motifs_prep.R, as in the timecourse-patterns T cell example:
# JASPAR2024 CORE vertebrates non-redundant (879 motifs), 200 bp windows,
# fixed seed. atac_trajectories/05_lost_site_motifs_plot.R summarizes the results.
#
# Requires: bedtools, MEME suite 5.5 (sea), curl; mm39 genome FASTA + .fai.
# The JASPAR motif file is downloaded once and checksum-verified.
#
# Inputs:  data/reference/GRCm39.primary_assembly.genome.fa{,.fai}
#          data/reference/jaspar2024_vert_nr.meme (downloaded if absent)
#          results/extended_analysis/atac_trajectories/{dynamic_sites,pass1_wt_ko,pass2_lost_in_ko}/
# Outputs: results/extended_analysis/atac_trajectories/motifs/{beds/,fa/,sea/,comparisons.tsv}
# Usage:   bash extended_analysis/scripts/atac_trajectories/04_lost_site_motifs.sh   (from the repository root)
# =============================================================================
set -euo pipefail

RSCRIPT=${RSCRIPT:-Rscript --no-save --no-restore}
THREADS=${THREADS:-8}
SEED=42
FA=data/reference/GRCm39.primary_assembly.genome.fa
MEME=data/reference/jaspar2024_vert_nr.meme
MEME_URL=https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
MEME_SHA=dd494278d356a4e170908c74d6d8eb746a2c7504d6210cd51017141f21a65b18
OUT=results/extended_analysis/atac_trajectories/motifs

for t in bedtools sea; do command -v $t >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }; done
[[ -f "$FA" && -f "$FA.fai" ]] || { echo "Skipping lost-site motif enrichment: needs $FA and its .fai (upstream mode)"; exit 0; }

if [[ ! -f "$MEME" ]]; then
  curl -L --fail --retry 3 -o "$MEME.part" "$MEME_URL" && mv "$MEME.part" "$MEME"
fi
echo "$MEME_SHA  $MEME" | sha256sum -c --quiet - || { echo "ERROR: $MEME checksum mismatch" >&2; exit 1; }

$RSCRIPT extended_analysis/scripts/atac_trajectories/04_lost_site_motifs_prep.R

# --- sequences, once per region set -----------------------------------------
mkdir -p "$OUT/fa" "$OUT/sea"
for b in "$OUT"/beds/*.bed; do
  n=$(basename "$b" .bed)
  bedtools getfasta -fi "$FA" -bed "$b" -nameOnly -fo "$OUT/fa/$n.fa"
done

# --- one SEA run per comparison ----------------------------------------------
tail -n +2 "$OUT/comparisons.tsv" | cut -f1,3,4 | \
  xargs -P "$THREADS" -L 1 bash -c '
    cmp=$1; p=$2; c=$3; o="'"$OUT"'/sea/${cmp}__${p}"
    [[ -s "'"$OUT"'/fa/$p.fa" && -s "'"$OUT"'/fa/$c.fa" ]] || { echo "  skip $cmp $p (empty set)"; exit 0; }
    sea --p "'"$OUT"'/fa/$p.fa" --n "'"$OUT"'/fa/$c.fa" --m "'"$MEME"'" --o "$o" --seed '"$SEED"' \
      > "$o.log" 2>&1 || { echo "ERROR: sea failed for $cmp $p" >&2; exit 255; }
    printf "  %-13s %-36s %4d motifs\n" "$cmp" "$p" "$(grep -vc "^#\|^RANK\|^$" "$o/sea.tsv")"
  ' _

echo "Done: $OUT/sea/"
