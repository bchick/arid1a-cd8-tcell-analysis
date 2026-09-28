#!/usr/bin/env bash
# =============================================================================
# atac_trajectories/01_timecourse_patterns.sh — WT and ARID1A-KO ATAC trajectory classes
# McDonald, Chick et al. 2023 Immunity 56:1303 — extended analysis
#
# Runs timecourse-patterns (github.com/bchick/timecourse-patterns, MIT) twice:
#   pass 1  WT and KO arms (Naive WT shared): dynamic peaks per arm + classes
#           inputs: atac_trajectories/01_timecourse_patterns_prep.R
#   pass 2  WT trajectories of the sites dynamic in WT but not in KO
#           inputs: atac_trajectories/02_dynamic_site_sets.R (from pass 1's dynamic calls)
# atac_trajectories/03_ko_trajectory_projection.R then projects KO trajectories onto both.
#
# The workflow brings its own pinned R/Bioconductor stack through pixi
# (pixi.lock), separate from this repository's renv library.
# Requires: git, pixi (https://pixi.sh). ~40 min on 8 cores, ~4 GB RAM.
#
# Inputs:  results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt
# Outputs: results/extended_analysis/atac_trajectories/{pass1_wt_ko,dynamic_sites,pass2_lost_in_ko}/
#          (pinned workflow checkout in external/timecourse-patterns/)
# Usage:   bash extended_analysis/scripts/atac_trajectories/01_timecourse_patterns.sh   (from the repository root)
# =============================================================================
set -euo pipefail

TCP_URL=${TCP_URL:-https://github.com/bchick/timecourse-patterns}
TCP_REF=${TCP_REF:-087f02b806138cdba58ec5b4d56db6424df631ec}
TCP_DIR=${TCP_DIR:-external/timecourse-patterns}
THREADS=${THREADS:-8}
RSCRIPT=${RSCRIPT:-Rscript --no-save --no-restore}
E=extended_analysis/scripts
OUT=results/extended_analysis/atac_trajectories

command -v pixi >/dev/null || { echo "ERROR: pixi not found (https://pixi.sh)" >&2; exit 1; }

# --- pinned checkout + environment --------------------------------------------
if [[ ! -d "$TCP_DIR/.git" ]]; then
  git clone --quiet "$TCP_URL" "$TCP_DIR"
fi
if [[ "$(git -C "$TCP_DIR" rev-parse HEAD)" != "$TCP_REF" ]]; then
  git -C "$TCP_DIR" fetch --quiet origin
  git -C "$TCP_DIR" checkout --quiet --detach "$TCP_REF"
fi
echo "timecourse-patterns @ $(git -C "$TCP_DIR" rev-parse --short HEAD)"
# --locked: install exactly what pixi.lock pins, never re-solve
(cd "$TCP_DIR" && pixi install --locked && pixi run postinstall)

run_tcp() {  # $1 = run directory holding config.yaml
  local cfg
  cfg=$(realpath "$1/config.yaml")
  (cd "$TCP_DIR" && pixi run snakemake --configfile "$cfg" -j "$THREADS")
}

# --- pass 1: WT and KO timecourses --------------------------------------------
$RSCRIPT $E/atac_trajectories/01_timecourse_patterns_prep.R
run_tcp "$OUT/pass1_wt_ko"

# --- pass 2: sites dynamic in WT but not in KO --------------------------------
$RSCRIPT $E/atac_trajectories/02_dynamic_site_sets.R
run_tcp "$OUT/pass2_lost_in_ko"

echo "Done: $OUT/pass1_wt_ko/results/report.html, $OUT/pass2_lost_in_ko/results/report.html"
