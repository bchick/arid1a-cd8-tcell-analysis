#!/usr/bin/env bash
# =============================================================================
# fetch_data.sh — Download, verify and unpack the processed-data bundle
#
#   bash tools/fetch_data.sh                  # from Zenodo
#   BUNDLE_TARBALL=/path/to/bundle.tar.gz bash tools/fetch_data.sh   # local copy
#
# The tarball holds repository-relative paths (results/..., data/...), so it is
# unpacked in place at the repository root. Per-file SHA-256 sums are checked
# after unpacking (bundle/SHA256SUMS).
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

ZENODO_RECORD="XXXXXXX"   # set at release
BUNDLE_NAME="arid1a-cd8-tcell-analysis_bundle_v1.0.0.tar.gz"
BUNDLE_URL="${BUNDLE_URL:-https://zenodo.org/records/${ZENODO_RECORD}/files/${BUNDLE_NAME}?download=1}"
TARBALL="${BUNDLE_TARBALL:-bundle/${BUNDLE_NAME}}"

mkdir -p bundle
if [[ ! -f "$TARBALL" ]]; then
  echo "Downloading $BUNDLE_URL"
  curl -L --fail --retry 3 -o "$TARBALL.part" "$BUNDLE_URL"
  mv "$TARBALL.part" "$TARBALL"
fi

echo "Unpacking $TARBALL"
tar -xzf "$TARBALL" -C "$ROOT"

echo "Verifying checksums"
( cd "$ROOT" && sha256sum --quiet -c bundle/SHA256SUMS )
echo "Bundle OK: $(wc -l < bundle/SHA256SUMS) files"
