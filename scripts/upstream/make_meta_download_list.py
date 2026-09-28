#!/usr/bin/env python3
"""
upstream/make_meta_download_list.py — Per-file ENA download list by scope
McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing

Expands the cross-study master sample sheet into a per-file ENA download list
with scope tags (arid1a, cmyc, naive, pbaf, singlecell), so downloads can be
launched by scope with scripts/upstream/download_meta_fastqs.sh.

Inputs:  data/metadata/meta_analysis/meta_master_sample_sheet.tsv
Outputs: data/metadata/meta_analysis/ena_download_urls.tsv
         (columns: scope_tag study assay datatype run url md5 dest_path)
Usage:   python3 scripts/upstream/make_meta_download_list.py
"""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
MDIR = os.path.join(ROOT, "data", "metadata", "meta_analysis")
FQROOT = "data/fastq/meta_analysis"   # repo-relative destination root

STUDY_DIR = {"Guo2022": "guo2022", "Baxter2023": "baxter2023"}
ASSAY_DIR = {"RNA": "rnaseq", "ATAC": "atacseq", "scRNA": "scrna", "scATAC": "scatac"}


def scope_tag(r):
    if r["datatype"] == "single-cell":
        return "singlecell"
    if r["study"] == "Baxter2023":
        return "pbaf"
    # Guo
    if r["arid1a_relevant"] == "True":
        return "arid1a"
    if r["target"] == "cMyc":
        return "cmyc"
    if r["target"] == "naive":
        return "naive"
    return "other"


def main():
    rows = list(csv.DictReader(open(os.path.join(MDIR, "meta_master_sample_sheet.tsv")), delimiter="\t"))
    out = []
    for r in rows:
        tag = scope_tag(r)
        sdir = STUDY_DIR.get(r["study"], r["study"].lower())
        adir = ASSAY_DIR.get(r["assay"], r["assay"].lower())
        urls = [u for u in r["fastq_ftp"].split(";") if u]
        md5s = r["fastq_md5"].split(";")
        for i, u in enumerate(urls):
            url = u if u.startswith("http") or u.startswith("ftp") else "ftp://" + u
            fname = u.split("/")[-1]
            dest = f"{FQROOT}/{sdir}/{adir}/{fname}"
            md5 = md5s[i] if i < len(md5s) else ""
            out.append([tag, r["study"], r["assay"], r["datatype"], r["run"], url, md5, dest])

    outf = os.path.join(MDIR, "ena_download_urls.tsv")
    with open(outf, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["scope_tag", "study", "assay", "datatype", "run", "url", "md5", "dest_path"])
        w.writerows(out)

    # summary
    import collections
    agg = collections.Counter(r[0] for r in out)
    print(f"wrote {outf}  ({len(out)} files)")
    for tag, n in sorted(agg.items()):
        print(f"   scope={tag:12} files={n}")


if __name__ == "__main__":
    main()
