#!/usr/bin/env python3
"""
upstream/build_meta_sample_sheet.py — Cross-study master sample sheet
McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing

Consolidates the ENA filereport TSVs for the cross-study meta-analysis
(Guo et al. 2022 + Baxter et al. 2023) into one annotated master sample sheet.
Each row is one SRA run, annotated with study / assay / datatype / BAF
complex / target / genotype-condition, plus flags (bulk, arid1a_relevant) so
downstream scope filtering is explicit and reproducible.

Inputs:  data/metadata/meta_analysis/ena_PRJNA*.tsv  (ENA filereport API)
Outputs: data/metadata/meta_analysis/meta_master_sample_sheet.tsv
Usage:   python3 scripts/upstream/build_meta_sample_sheet.py
"""
import csv, re, os, collections

HERE = os.path.dirname(os.path.abspath(__file__))
MDIR = os.path.join(HERE, "..", "..", "data", "metadata", "meta_analysis")

# BioProject -> (study, GSE, nominal assay from GEO sub-series)
PRJ = {
    "PRJNA761410": ("Guo2022",    "GSE183615", "RNA"),
    "PRJNA761411": ("Guo2022",    "GSE183616", "ATAC"),
    "PRJNA761413": ("Guo2022",    "GSE183618", "ATAC"),
    "PRJNA818712": ("Guo2022",    "GSE199184", "RNA"),
    "PRJNA817256": ("Guo2022",    "GSE198894", "ATAC"),
    "PRJNA870147": ("Baxter2023", "GSE211409", "ATAC"),
    "PRJNA870148": ("Baxter2023", "GSE211410", "RNA"),
    "PRJNA960837": ("Baxter2023", "GSE211410", "SC"),   # single-cell multiome
}


def classify(title):
    """Return dict of assay/datatype/complex/target/condition from ENA sample_title."""
    t = title.strip()
    low = t.lower()

    # --- single-cell multiome (Baxter PRJNA960837): "LCMV7: Cl13D30, scRNA" -----
    m = re.search(r"\bsc(rna|atac)\b", low)
    if m:
        assay = "scRNA" if m.group(1) == "rna" else "scATAC"
        cond = re.sub(r",?\s*sc(rna|atac)\s*$", "", t, flags=re.I).strip(" ,")
        return dict(assay=assay, datatype="single-cell", complex="NA",
                    target="timecourse", condition=cond, arid1a=False)

    # --- assay from strategy is passed separately; infer target/condition -------
    # Guo GSE198894 vehicle/inhibitor arms — label DMSO controls consistently
    if "dmso" in low:
        geno = "KO" if "ko" in low else "WT"
        return dict(complex="cBAF", target="Arid1a", condition=f"DMSO_{geno}", arid1a=True)
    # Guo Arid1a arms
    if "arid1a" in low:
        geno = "KO" if re.search(r"arid1a\s*ko", low) else ("WT" if re.search(r"arid1a\s*wt", low) else "?")
        if "dmso" in low:
            cond = f"DMSO_{geno}"
        else:
            cond = geno
        return dict(complex="cBAF", target="Arid1a", condition=cond, arid1a=True)
    # Guo Arid1a pharmacologic inhibitor (BRD-K98645985) — cBAF/Arid1a relevant
    if "brd-k98645985" in low or "brd_k98645985" in low:
        return dict(complex="cBAF", target="Arid1a", condition="inhibitor_BRDK98645985", arid1a=True)
    # Guo c-Myc arms
    if "myc" in low:
        if "high" in low: cond = "Myc_high"
        elif "low" in low: cond = "Myc_low"
        elif re.search(r"myc\s*ko", low): cond = "Myc_KO"
        elif re.search(r"myc\s*wt", low): cond = "Myc_WT"
        else: cond = "Myc"
        return dict(complex="cBAF", target="cMyc", condition=cond, arid1a=False)
    if "naïve" in low or "naive" in low:
        return dict(complex="NA", target="naive", condition="Naive", arid1a=False)
    # Baxter CRISPR sgRNA perturbations
    m = re.search(r"sg([a-z0-9._]+)", low)
    if m:
        gene = m.group(1)
        if gene.startswith("arid2"): comp, tgt = "PBAF", "Arid2"
        elif gene.startswith("pbrm1"): comp, tgt = "PBAF", "Pbrm1"
        elif gene.startswith("ano9"): comp, tgt = "control", "sgControl(Ano9)"
        else: comp, tgt = "other", f"sg{gene}"
        # keep guide id (e.g. e3.1) as condition detail
        guide = re.search(r"sg[a-z0-9]+_?(e[0-9.]+)?", low)
        cond = re.sub(r"\s*\[.*?\]", "", t).strip()
        return dict(complex=comp, target=tgt, condition=cond, arid1a=False)
    # Guo GSE183618 plain "WT ATAC-seq"
    if re.search(r"\bwt\b", low):
        return dict(complex="cBAF", target="WT", condition="WT", arid1a=True)
    return dict(complex="?", target="?", condition=t, arid1a=False)


def rep_from_title(title):
    m = re.search(r"[Rr]eplicate\s*([A-Za-z0-9]+)", title) or re.search(r"[_-][Rr]ep[_-]?(\d+)", title)
    return m.group(1) if m else ""


def main():
    out_rows = []
    for prj, (study, gse, nominal) in PRJ.items():
        f = os.path.join(MDIR, f"ena_{prj}.tsv")
        if not os.path.exists(f):
            print(f"WARN missing {f}")
            continue
        for r in csv.DictReader(open(f), delimiter="\t"):
            if not r.get("run_accession"):
                continue
            strat = r["library_strategy"]
            assay = {"RNA-Seq": "RNA", "ATAC-seq": "ATAC"}.get(strat, strat)
            info = classify(r["sample_title"])
            # single-cell classification overrides assay
            if info.get("datatype") == "single-cell":
                assay = info["assay"]; datatype = "single-cell"
            else:
                datatype = "bulk"
            try:
                gb = round(sum(int(b) for b in r["fastq_bytes"].split(";") if b) / 1e9, 2) if r.get("fastq_bytes") else 0
            except Exception:
                gb = 0
            out_rows.append(dict(
                study=study, gse=gse, bioproject=prj,
                assay=assay, datatype=datatype,
                complex=info.get("complex", "?"), target=info.get("target", "?"),
                condition=info.get("condition", ""), replicate=rep_from_title(r["sample_title"]),
                arid1a_relevant=info.get("arid1a", False),
                run=r["run_accession"], experiment=r["experiment_accession"],
                gsm=r["sample_alias"], layout=r["library_layout"],
                approx_GB=gb, fastq_ftp=r["fastq_ftp"], fastq_md5=r["fastq_md5"],
                sample_title=r["sample_title"],
            ))

    cols = ["study", "gse", "bioproject", "assay", "datatype", "complex", "target",
            "condition", "replicate", "arid1a_relevant", "run", "experiment", "gsm",
            "layout", "approx_GB", "fastq_ftp", "fastq_md5", "sample_title"]
    out = os.path.join(MDIR, "meta_master_sample_sheet.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(out_rows)

    # ---- summary ------------------------------------------------------------
    print(f"wrote {out}  ({len(out_rows)} runs)\n")
    agg = collections.defaultdict(lambda: [0, 0.0])
    for r in out_rows:
        k = (r["study"], r["assay"], r["datatype"])
        agg[k][0] += 1; agg[k][1] += r["approx_GB"]
    print(f"{'study':11} {'assay':7} {'datatype':12} {'runs':>5} {'~GB':>8}")
    for k in sorted(agg):
        n, gb = agg[k]
        print(f"{k[0]:11} {k[1]:7} {k[2]:12} {n:5d} {gb:8.1f}")
    tot_gb = sum(r["approx_GB"] for r in out_rows)
    bulk_gb = sum(r["approx_GB"] for r in out_rows if r["datatype"] == "bulk")
    a1_gb = sum(r["approx_GB"] for r in out_rows if r["arid1a_relevant"])
    print(f"\nTOTAL {len(out_rows)} runs ~{tot_gb:.0f} GB | bulk-only ~{bulk_gb:.0f} GB | "
          f"Arid1a-relevant ~{a1_gb:.0f} GB ({sum(r['arid1a_relevant'] for r in out_rows)} runs)")


if __name__ == "__main__":
    main()
