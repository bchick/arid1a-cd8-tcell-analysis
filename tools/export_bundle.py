#!/usr/bin/env python3
"""Build the processed-data bundle from a fully processed working tree.

Usage:
    python tools/export_bundle.py [--manifest tools/bundle_manifest.tsv]
                                  [--out bundle] [--dry-run]

Reads the manifest (glob <TAB> action <TAB> note), collects every matching
file, writes bundle/SHA256SUMS (repository-relative paths) and
bundle/MANIFEST.tsv, and packs everything into
bundle/arid1a-cd8-tcell-analysis_bundle_v<version>.tar.gz, ready for Zenodo.
Patterns that match nothing are reported and make the export fail.
"""
import argparse
import csv
import gzip
import hashlib
import io
import shutil
import sqlite3
import sys
import tarfile
import tempfile
from pathlib import Path

VERSION = "1.0.0"
ABC_FLOOR = 0.01            # keep rows with ABC.Score or powerlaw.Score >= floor
ABC_COLUMNS = ["chr", "start", "end", "name", "class", "activity_base",
               "normalized_atac_enh", "normalized_h3k27ac_enh", "TargetGene",
               "TargetGeneTSS", "TargetGeneIsExpressed", "distance",
               "ABC.Score", "powerlaw.Score"]


TEXT_SUFFIXES = {".tsv", ".csv", ".txt", ".md", ".html", ".json", ".bed", ".narrowPeak"}


def strip_root(data, rel, roots):
    """Make paths in text outputs repository-relative (e.g. SpaMo command lines)."""
    if Path(rel).suffix not in TEXT_SUFFIXES:
        return data
    for r in roots:
        data = data.replace(r + b"/", b"").replace(r, b".")
    return data


def sqlite_strip_root(path, roots):
    """TxDb metadata records the source GTF's absolute path; rewrite it in a copy."""
    with tempfile.TemporaryDirectory() as tmp:
        cp = Path(tmp) / path.name
        shutil.copyfile(path, cp)
        con = sqlite3.connect(cp)
        for r in (r.decode() for r in roots):
            con.execute("UPDATE metadata SET value = replace(replace(value, ?, ''), ?, '.')",
                        (r + "/", r))
        con.commit()
        con.execute("VACUUM")
        con.close()
        return cp.read_bytes()


def read_manifest(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        glob, action, *note = line.split("\t")
        rows.append((glob, action, note[0] if note else ""))
    return rows


def filtered_abc(path):
    """Return gzip bytes of an ABC table reduced to used columns and a score floor."""
    out = io.BytesIO()
    with gzip.open(path, "rt") as fin, gzip.GzipFile(fileobj=out, mode="wb", mtime=0) as gz:
        reader = csv.DictReader(fin, delimiter="\t")
        cols = [c for c in ABC_COLUMNS if c in reader.fieldnames]
        w = io.TextIOWrapper(gz, newline="")
        writer = csv.DictWriter(w, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            score = max(float(row.get("ABC.Score") or 0), float(row.get("powerlaw.Score") or 0))
            if score >= ABC_FLOOR:
                writer.writerow(row)
        w.flush()
        w.detach()
    return out.getvalue()


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", default=".")
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--out", default="bundle")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--strip-prefix", action="append", default=[],
                    help="extra absolute prefix to strip from text files (repeatable)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    manifest = Path(args.manifest) if args.manifest else next(
        p for p in (root / "tools/bundle_manifest.tsv", root / "release/bundle_manifest.tsv") if p.exists())
    out = root / args.out
    out.mkdir(parents=True, exist_ok=True)

    entries, empty = [], []
    for glob, action, note in read_manifest(manifest):
        hits = sorted(p for p in root.glob(glob) if p.is_file())
        if not hits:
            empty.append(glob)
        entries += [(p.relative_to(root).as_posix(), action, note) for p in hits]

    total = sum((root / e[0]).stat().st_size for e in entries)
    print(f"{len(entries)} files, {total / 1e6:.0f} MB before filtering/compression")
    for g in empty:
        print(f"  NO MATCH: {g}")
    if args.dry_run:
        for rel, action, _ in entries:
            print(f"  {action:5s} {(root / rel).stat().st_size / 1e6:8.1f} MB  {rel}")
        sys.exit(1 if empty else 0)
    if empty:
        sys.exit("Fix the manifest: patterns above matched nothing")

    # the project may be reached through a symlinked path; strip both spellings
    roots = sorted({str(root).encode(), str(Path(args.root).absolute()).encode(),
                    *(a.encode() for a in args.strip_prefix)}, key=len, reverse=True)
    tar_path = out / f"arid1a-cd8-tcell-analysis_bundle_v{VERSION}.tar.gz"
    sums, man = [], []
    with tarfile.open(tar_path, "w:gz") as tar:
        for rel, action, note in entries:
            if action == "abc":
                data = filtered_abc(root / rel)
            elif rel.endswith(".sqlite"):
                data = sqlite_strip_root(root / rel, roots)
            else:
                data = strip_root((root / rel).read_bytes(), rel, roots)
            info = tarfile.TarInfo(rel)
            info.size, info.mtime, info.mode = len(data), 0, 0o644
            tar.addfile(info, io.BytesIO(data))
            sums.append(f"{hashlib.sha256(data).hexdigest()}  {rel}")
            man.append(f"{rel}\t{action}\t{len(data)}\t{note}")
            print(f"  {len(data) / 1e6:8.1f} MB  {rel}")
        for name, lines in (("bundle/SHA256SUMS", sums),
                            ("bundle/MANIFEST.tsv", ["path\taction\tbytes\tnote"] + man)):
            data = ("\n".join(lines) + "\n").encode()
            (root / name).write_bytes(data)
            info = tarfile.TarInfo(name)
            info.size, info.mtime, info.mode = len(data), 0, 0o644
            tar.addfile(info, io.BytesIO(data))
    print(f"Wrote {tar_path} ({tar_path.stat().st_size / 1e6:.0f} MB)")


if __name__ == "__main__":
    main()
