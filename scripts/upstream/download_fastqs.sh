#!/usr/bin/env bash
# =============================================================================
# upstream/download_fastqs.sh — Download the study FASTQs from ENA
# McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing
#
# Downloads the pre-compressed .fastq.gz files for every sample in the master
# sample sheet (GSE228381 sub-series) from ENA, with MD5 verification.
# Already-downloaded files are skipped; re-run to retry failures.
#   assay:        rnaseq|atacseq|atacseq_inhibitors|cutandrun|chipseq|all (default: all)
#   max_parallel: number of concurrent wget jobs (default: 4)
# Requires: wget, md5sum, python3 + pandas; ample free disk (raw FASTQs are large).
#
# Inputs:  data/metadata/master_sample_sheet.tsv, data/metadata/ena_fastq_urls.tsv
# Outputs: data/fastq/{rnaseq,atacseq,cutandrun,chipseq}/*.fastq.gz
#          data/fastq/download_manifest.txt, data/fastq/logs/
# Usage:   bash scripts/upstream/download_fastqs.sh [assay] [max_parallel]
# =============================================================================

PROJECT="${ARID1A_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
MASTER="${PROJECT}/data/metadata/master_sample_sheet.tsv"
ENA_URLS="${PROJECT}/data/metadata/ena_fastq_urls.tsv"
FASTQ_BASE="${PROJECT}/data/fastq"
LOGDIR="${FASTQ_BASE}/logs"

ASSAY_FILTER="${1:-all}"
MAX_PARALLEL="${2:-4}"

mkdir -p "$LOGDIR" "${FASTQ_BASE}/rnaseq" "${FASTQ_BASE}/atacseq" "${FASTQ_BASE}/cutandrun" "${FASTQ_BASE}/chipseq"

echo "============================================"
echo "FASTQ Download from ENA"
echo "Assay filter: ${ASSAY_FILTER}"
echo "Max parallel: ${MAX_PARALLEL}"
echo "Started: $(date)"
echo "============================================"

# Build download manifest as individual shell commands
python3 -c "
import pandas as pd

master = pd.read_csv('${MASTER}', sep='\t')
ena = pd.read_csv('${ENA_URLS}', sep='\t')
assay_filter = '${ASSAY_FILTER}'

srr_to_assay = {}
for _, r in master.iterrows():
    for srr in str(r['srr']).split(';'):
        srr_to_assay[srr] = r['assay']

assay_dirs = {
    'rnaseq': '${FASTQ_BASE}/rnaseq',
    'atacseq': '${FASTQ_BASE}/atacseq',
    'atacseq_inhibitors': '${FASTQ_BASE}/atacseq',
    'cutandrun': '${FASTQ_BASE}/cutandrun',
    'chipseq': '${FASTQ_BASE}/chipseq',
}

for _, row in ena.iterrows():
    srr = row['run_accession']
    assay = srr_to_assay.get(srr, 'other')
    if assay_filter != 'all' and assay != assay_filter:
        continue
    outdir = assay_dirs.get(assay, '${FASTQ_BASE}/other')
    urls = str(row['fastq_ftp']).split(';')
    md5s = str(row['fastq_md5']).split(';')
    for url, md5 in zip(urls, md5s):
        filename = url.split('/')[-1]
        print(f'ftp://{url}|{outdir}|{filename}|{md5}')
" > "${FASTQ_BASE}/download_manifest.txt"

TOTAL_FILES=$(wc -l < "${FASTQ_BASE}/download_manifest.txt")
echo "Total FASTQ files: ${TOTAL_FILES}"

# Write the per-file download script (self-contained, pipe-delimited args)
cat > "${FASTQ_BASE}/dl_one.sh" << 'DLEOF'
#!/bin/bash
# Args: "url|outdir|filename|md5"
IFS='|' read -r url outdir filename expected_md5 <<< "$1"
logdir="$2"
filepath="${outdir}/${filename}"

if [ -f "$filepath" ]; then
    actual_md5=$(md5sum "$filepath" | cut -d' ' -f1)
    if [ "$actual_md5" = "$expected_md5" ]; then
        echo "[SKIP] ${filename}"
        exit 0
    else
        echo "[REDO] ${filename} (bad MD5)"
        rm -f "$filepath"
    fi
fi

echo "[DOWN] ${filename}"
mkdir -p "$outdir"

if wget -q --tries=3 --timeout=120 --waitretry=10 -O "${filepath}.tmp" "$url" 2>"${logdir}/${filename}.log"; then
    mv "${filepath}.tmp" "$filepath"
    actual_md5=$(md5sum "$filepath" | cut -d' ' -f1)
    if [ "$actual_md5" = "$expected_md5" ]; then
        size=$(du -h "$filepath" | cut -f1)
        echo "[OK]   ${filename} (${size})"
        exit 0
    else
        echo "[FAIL] ${filename} (MD5 mismatch)"
        rm -f "$filepath"
        exit 1
    fi
else
    echo "[FAIL] ${filename} (wget error)"
    rm -f "${filepath}.tmp"
    exit 1
fi
DLEOF
chmod +x "${FASTQ_BASE}/dl_one.sh"

# Count already downloaded
EXISTING=0
while IFS='|' read -r url outdir filename md5; do
    if [ -f "${outdir}/${filename}" ]; then
        EXISTING=$((EXISTING + 1))
    fi
done < "${FASTQ_BASE}/download_manifest.txt"
echo "Already present: ${EXISTING}/${TOTAL_FILES}"
echo "Remaining: $((TOTAL_FILES - EXISTING))"
echo ""

# Run with xargs for parallel execution
cat "${FASTQ_BASE}/download_manifest.txt" | \
    xargs -P "${MAX_PARALLEL}" -I '{}' bash "${FASTQ_BASE}/dl_one.sh" '{}' "${LOGDIR}"

echo ""
echo "============================================"
echo "Finished: $(date)"
echo "============================================"

# Final verification
echo ""
echo "Verifying all downloads..."
MISSING=0
VERIFIED=0
while IFS='|' read -r url outdir filename md5; do
    filepath="${outdir}/${filename}"
    if [ -f "$filepath" ]; then
        actual_md5=$(md5sum "$filepath" | cut -d' ' -f1)
        if [ "$actual_md5" = "$md5" ]; then
            VERIFIED=$((VERIFIED + 1))
        else
            echo "BAD MD5: ${filename}"
            MISSING=$((MISSING + 1))
        fi
    else
        echo "MISSING: ${filename}"
        MISSING=$((MISSING + 1))
    fi
done < "${FASTQ_BASE}/download_manifest.txt"

echo ""
echo "Verified: ${VERIFIED}/${TOTAL_FILES}"
if [ "$MISSING" -gt 0 ]; then
    echo "MISSING/CORRUPT: ${MISSING} -- re-run script to retry"
else
    echo "All ${TOTAL_FILES} files downloaded and verified!"
fi

echo ""
echo "Disk usage:"
for d in rnaseq atacseq cutandrun chipseq; do
    if [ -d "${FASTQ_BASE}/${d}" ]; then
        echo "  ${d}: $(du -sh "${FASTQ_BASE}/${d}" | cut -f1)"
    fi
done
echo "  Total: $(du -sh "${FASTQ_BASE}" | cut -f1)"
