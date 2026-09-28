#!/usr/bin/env bash
# =============================================================================
# upstream/download_meta_fastqs.sh — Download cross-study FASTQs from ENA
# McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing
#
# Downloads the Guo et al. 2022 and Baxter et al. 2023 FASTQs used in the
# cross-study meta-analysis, filtered by scope tag, with MD5 verification.
# Mirrors scripts/upstream/download_fastqs.sh. By default this is a DRY RUN
# (prints what it would fetch); add --go to download.
#
# Scope tags (default: "arid1a pbaf"):
#   arid1a      Guo Arid1a/WT/DMSO/inhibitor RNA+ATAC  (the core comparison)
#   cmyc        Guo c-Myc arms RNA+ATAC                (effector/memory context)
#   naive       Guo naive ATAC
#   pbaf        Baxter Arid2/Pbrm1 (PBAF) + control bulk RNA+ATAC (contrast)
#   singlecell  Baxter 10x scRNA/scATAC  — not nf-core-processable; excluded by default
#   all-bulk    = arid1a cmyc naive pbaf
# The Guo arid1a scope alone is roughly 364 GB.
#
# Examples:
#   bash scripts/upstream/download_meta_fastqs.sh                      # dry run, core scope
#   bash scripts/upstream/download_meta_fastqs.sh --go arid1a pbaf     # download core
#   bash scripts/upstream/download_meta_fastqs.sh --go all-bulk -p 6   # download all bulk
#
# Inputs:  data/metadata/meta_analysis/ena_download_urls.tsv
#          (from scripts/upstream/make_meta_download_list.py)
# Outputs: data/fastq/meta_analysis/{guo2022,baxter2023}/<assay>/*.fastq.gz
#          data/fastq/meta_analysis/logs/
# Usage:   bash scripts/upstream/download_meta_fastqs.sh [--go] [scope ...] [--parallel N]
# =============================================================================
set -uo pipefail

PROJECT="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
LIST="${PROJECT}/data/metadata/meta_analysis/ena_download_urls.tsv"
LOGDIR="${PROJECT}/data/fastq/meta_analysis/logs"
mkdir -p "$LOGDIR"

GO=0
PARALLEL=4
SCOPES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --go) GO=1; shift ;;
    -p|--parallel) PARALLEL="$2"; shift 2 ;;
    all-bulk) SCOPES+=(arid1a cmyc naive pbaf); shift ;;
    arid1a|cmyc|naive|pbaf|singlecell) SCOPES+=("$1"); shift ;;
    *) echo "unknown arg: $1"; exit 1 ;;
  esac
done
[[ ${#SCOPES[@]} -eq 0 ]] && SCOPES=(arid1a pbaf)
SCOPE_RE=$(IFS='|'; echo "${SCOPES[*]}")

echo "============================================"
echo "Meta-analysis FASTQ download (ENA)"
echo "Scope:    ${SCOPES[*]}"
echo "Mode:     $([[ $GO -eq 1 ]] && echo 'DOWNLOAD' || echo 'DRY RUN (add --go to download)')"
echo "Parallel: ${PARALLEL}"
echo "Started:  $(date)"
echo "============================================"

# Build filtered manifest: url<TAB>dest<TAB>md5 ; rewrite bare ENA host to https
MANIFEST=$(mktemp)
awk -F'\t' -v re="^(${SCOPE_RE})$" 'NR>1 && $1 ~ re {
  url=$6; sub(/^ftp:\/\//,"",url); if (url ~ /^ftp\.sra\.ebi\.ac\.uk/) url="https://" url;
  print url "\t" $8 "\t" $7
}' "$LIST" > "$MANIFEST"

NFILES=$(wc -l < "$MANIFEST")
# estimate size from master sheet
GB=$(awk -F'\t' -v re="^(${SCOPE_RE})$" 'NR>1{ }' "$LIST")
echo "Files to fetch: ${NFILES}"
if [[ $NFILES -eq 0 ]]; then echo "nothing matched scope"; rm -f "$MANIFEST"; exit 0; fi

if [[ $GO -eq 0 ]]; then
  echo ""; echo "--- first 8 files that WOULD be fetched ---"
  head -8 "$MANIFEST" | awk -F'\t' '{print "  "$2}'
  echo "  ... (${NFILES} total)"
  echo ""; echo "DRY RUN — no files downloaded. Re-run with --go to start."
  rm -f "$MANIFEST"; exit 0
fi

fetch_one() {
  local url="$1" dest="$2" md5="$3"
  local abs="${PROJECT}/${dest}"
  mkdir -p "$(dirname "$abs")"
  local verify_ok
  verify_ok() {
    [[ -s "$abs" ]] || return 1
    [[ -z "$md5" ]] && return 0
    echo "${md5}  ${abs}" | md5sum -c --status 2>/dev/null
  }
  if verify_ok; then echo "[skip] $dest"; return 0; fi
  # wget -c must never be combined with -O: on a full-size corrupt file it
  # resumes past the damage instead of replacing it, silently preserving the
  # bad bytes. Always start clean and re-verify.
  local attempt
  for attempt in 1 2 3; do
    rm -f "$abs"
    wget -q --tries=3 --timeout=60 -O "$abs" "$url"
    if verify_ok; then echo "[ok]   $dest (attempt ${attempt})"; return 0; fi
    echo "[retry ${attempt}] checksum mismatch: $dest"
  done
  # Leave nothing corrupt behind: a missing file fails loudly downstream,
  # a corrupt one fails 20 hours into a pipeline.
  rm -f "$abs"
  echo "[MD5 FAIL] $dest — removed after 3 failed attempts"; return 1
}
export -f fetch_one; export PROJECT

cat "$MANIFEST" | xargs -P "$PARALLEL" -d '\n' -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r url dest md5 <<< "{}"; fetch_one "$url" "$dest" "$md5"
' 2>&1 | tee "${LOGDIR}/download_$(date +%Y%m%d_%H%M%S).log"

rm -f "$MANIFEST"
echo "Done: $(date)"
