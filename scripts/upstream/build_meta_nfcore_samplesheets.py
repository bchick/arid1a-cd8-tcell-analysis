#!/usr/bin/env python3
"""
upstream/build_meta_nfcore_samplesheets.py — nf-core samplesheets for Guo et al. 2022
McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing

Builds nf-core samplesheets for the Arid1a-relevant subset of Guo et al. 2022
(the cross-study core) from the cross-study master sample sheet.
  * Guo Arid1a-relevant only (WT/KO/DMSO/inhibitor); c-Myc, naive and Baxter
    samples are excluded.
  * The GSE is encoded in the sample name so the genetic (GSE183615/183618)
    and pharmacologic (GSE199184/198894) experiments stay separable downstream.
  * Runs sharing a GSM are merged (same sample name / same sample+replicate).
  * RNA strandedness = auto (Guo's library prep differs from ours).
  * All Guo libraries used here are paired-end.

Inputs:  data/metadata/meta_analysis/meta_master_sample_sheet.tsv
         data/fastq/meta_analysis/guo2022/  (paths written into the sheets)
Outputs: nextflow/samplesheets/meta_guo_rnaseq.csv   (sample,fastq_1,fastq_2,strandedness)
         nextflow/samplesheets/meta_guo_atacseq.csv  (sample,fastq_1,fastq_2,replicate)
Usage:   python3 scripts/upstream/build_meta_nfcore_samplesheets.py
"""
import os
import csv, os, collections

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
MASTER = os.path.join(ROOT, "data", "metadata", "meta_analysis", "meta_master_sample_sheet.tsv")
FQBASE = os.path.join(os.environ.get("ARID1A_PROJECT_DIR", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "data/fastq/meta_analysis/guo2022")
SSDIR = os.path.join(ROOT, "nextflow", "samplesheets")

COND_SHORT = {"WT": "WT", "KO": "KO", "DMSO_WT": "vehWT", "DMSO_KO": "vehKO",
              "inhibitor_BRDK98645985": "inhib",
              "Myc_WT": "MycWT", "Myc_KO": "MycKO", "Myc_high": "MycHi",
              "Myc_low": "MycLo", "Naive": "Naive"}


def local_fastqs(fastq_ftp, assay_dir):
    files = [u.split("/")[-1] for u in fastq_ftp.split(";") if u]
    r1 = [f for f in files if f.endswith("_1.fastq.gz")]
    r2 = [f for f in files if f.endswith("_2.fastq.gz")]
    f1 = f"{FQBASE}/{assay_dir}/{r1[0]}" if r1 else ""
    f2 = f"{FQBASE}/{assay_dir}/{r2[0]}" if r2 else ""
    return f1, f2


def main():
    rows = [r for r in csv.DictReader(open(MASTER), delimiter="\t")
            if r["study"] == "Guo2022" and r["datatype"] == "bulk"]  # all Guo bulk: Arid1a + c-Myc + naive
    os.makedirs(SSDIR, exist_ok=True)

    # ---------------- RNA (one GSM = one biological sample) ------------------
    rna = [r for r in rows if r["assay"] == "RNA"]
    # order by gse, condition
    rna.sort(key=lambda r: (r["gse"], r["condition"], r["gsm"]))
    ctr = collections.Counter()
    rna_out = []
    for r in rna:
        gse = r["gse"].replace("GSE", "g")
        cond = COND_SHORT.get(r["condition"], r["condition"])
        key = (gse, cond)
        ctr[key] += 1
        sample = f"Guo_{gse}_{cond}_r{ctr[key]}"
        f1, f2 = local_fastqs(r["fastq_ftp"], "rnaseq")
        rna_out.append([sample, f1, f2, "auto"])
    with open(os.path.join(SSDIR, "meta_guo_rnaseq.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
        w.writerows(rna_out)

    # ---------------- ATAC (group = condition; replicate = per GSM) ----------
    atac = [r for r in rows if r["assay"] == "ATAC"]
    # assign replicate index per (gse,cond) group, keyed by GSM
    grp_gsms = collections.defaultdict(list)
    for r in atac:
        gse = r["gse"].replace("GSE", "g"); cond = COND_SHORT.get(r["condition"], r["condition"])
        g = (gse, cond)
        if r["gsm"] not in grp_gsms[g]:
            grp_gsms[g].append(r["gsm"])
    atac.sort(key=lambda r: (r["gse"], r["condition"], r["gsm"], r["run"]))
    atac_out = []
    for r in atac:
        gse = r["gse"].replace("GSE", "g"); cond = COND_SHORT.get(r["condition"], r["condition"])
        g = (gse, cond)
        rep = grp_gsms[g].index(r["gsm"]) + 1
        sample = f"Guo_{gse}_{cond}"
        f1, f2 = local_fastqs(r["fastq_ftp"], "atacseq")
        atac_out.append([sample, f1, f2, rep])
    with open(os.path.join(SSDIR, "meta_guo_atacseq.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["sample", "fastq_1", "fastq_2", "replicate"])
        w.writerows(atac_out)

    # summary
    print(f"RNA  samples: {len(rna_out)} rows")
    for k in sorted(ctr): print(f"   {k[0]}_{k[1]}: {ctr[k]}")
    print(f"\nATAC rows: {len(atac_out)} (runs)  groups:")
    ag = collections.Counter((r[0]) for r in atac_out)
    for k in sorted(ag): print(f"   {k}: {ag[k]} runs / {len(grp_gsms[(k.split('_')[1], k.split('_',2)[2])]) if False else ''}".rstrip())
    for g in sorted(grp_gsms): print(f"   group Guo_{g[0]}_{g[1]}: {len(grp_gsms[g])} replicates")


if __name__ == "__main__":
    main()
