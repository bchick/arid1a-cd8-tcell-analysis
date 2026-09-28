#!/usr/bin/env python3
"""
upstream/build_sample_sheet.py — Master sample sheet and nf-core samplesheets
McDonald, Chick et al. 2023 Immunity 56:1303 — upstream processing

Parses the GEO series matrix files of the five GSE228381 sub-series
(GSE227634 RNA-seq, GSE228171 / GSE228193 ATAC-seq, GSE228380 CUT&RUN,
GSE228546 ChIP-seq) and the SRA run info, and builds one master sample sheet
plus the per-pipeline nf-core samplesheets.

Inputs:  data/metadata/GSE*_series_matrix.txt.gz
         data/metadata/sra_runinfo_PRJNA*.csv
Outputs: data/metadata/master_sample_sheet.tsv
         nextflow/samplesheets/{atacseq,rnaseq,cutandrun,chipseq}_samplesheet.csv
Usage:   python3 scripts/upstream/build_sample_sheet.py
"""

import gzip
import os
import re
import pandas as pd
from pathlib import Path

PROJECT = Path(os.environ.get("ARID1A_PROJECT_DIR", Path(__file__).resolve().parents[2]))
METADATA = PROJECT / "data" / "metadata"
SAMPLESHEETS = PROJECT / "nextflow" / "samplesheets"
SAMPLESHEETS.mkdir(parents=True, exist_ok=True)

# Map sub-series to assay type and BioProject
SUBSERIES = {
    "GSE227634": {"assay": "rnaseq", "bioproject": "PRJNA945945"},
    "GSE228171": {"assay": "atacseq", "bioproject": "PRJNA948430"},
    "GSE228193": {"assay": "atacseq_inhibitors", "bioproject": "PRJNA983073"},
    "GSE228380": {"assay": "cutandrun", "bioproject": "PRJNA949603"},
    "GSE228546": {"assay": "chipseq", "bioproject": "PRJNA950325"},
}


def parse_series_matrix(filepath):
    """Parse a GEO series matrix file and return sample-level metadata."""
    rows = {}
    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!Sample_"):
                key = line.split("\t")[0]
                vals = line.split("\t")[1:]
                vals = [v.strip('"') for v in vals]
                if key not in rows:
                    rows[key] = []
                rows[key].append(vals)
            if line.startswith("!series_matrix_table_begin"):
                break

    n_samples = len(rows["!Sample_geo_accession"][0])
    samples = []

    for i in range(n_samples):
        s = {}
        s["title"] = rows["!Sample_title"][0][i]
        s["gsm"] = rows["!Sample_geo_accession"][0][i]

        # Parse characteristics
        for char_row in rows.get("!Sample_characteristics_ch1", []):
            val = char_row[i]
            if ": " in val:
                k, v = val.split(": ", 1)
                s[f"char_{k.strip()}"] = v.strip()

        # Library strategy
        if "!Sample_library_strategy" in rows:
            s["library_strategy"] = rows["!Sample_library_strategy"][0][i]

        # Instrument
        if "!Sample_instrument_model" in rows:
            s["instrument"] = rows["!Sample_instrument_model"][0][i]

        # SRX from relation
        for rel_row in rows.get("!Sample_relation", []):
            val = rel_row[i]
            if "SRA:" in val:
                srx_match = re.search(r"SRX\d+", val)
                if srx_match:
                    s["srx"] = srx_match.group()

        samples.append(s)

    return samples


def load_sra_runinfo():
    """Load and merge all SRA runinfo CSVs."""
    dfs = []
    for gse, info in SUBSERIES.items():
        csv_path = METADATA / f"sra_runinfo_{info['bioproject']}.csv"
        if csv_path.exists() and csv_path.stat().st_size > 0:
            df = pd.read_csv(csv_path)
            # Filter out empty rows
            df = df[df["Run"].notna() & (df["Run"] != "")]
            df["gse"] = gse
            dfs.append(df)
        else:
            print(f"WARNING: Missing runinfo for {gse} ({info['bioproject']})")
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()


def parse_rnaseq_title(title):
    """Parse RNA-seq sample titles like 'D3_WT_Rep1', 'D8_Het_TE_Rep2'."""
    # Check for cell subset pattern FIRST (more specific)
    m = re.match(r"D(\d+)_([A-Za-z]+)_(TE|EEC|MP)_(Rep\d+)$", title)
    if m:
        return {
            "timepoint": f"D{m.group(1)}",
            "genotype": m.group(2),
            "cell_subset": m.group(3),
            "replicate": m.group(4),
        }
    m = re.match(r"D(\d+)_([A-Za-z]+)_(Rep\d+)$", title)
    if m:
        return {
            "timepoint": f"D{m.group(1)}",
            "genotype": m.group(2),
            "cell_subset": "total",
            "replicate": m.group(3),
        }
    return {"timepoint": "", "genotype": "", "cell_subset": "", "replicate": ""}


def parse_atacseq_title(title):
    """Parse ATAC-seq sample titles with bracket notation."""
    # Extract the short name in brackets, e.g., [D3_WT_Rep1]
    bracket = re.search(r"\[(.+?)\]", title)
    if bracket:
        short = bracket.group(1)
    else:
        short = title

    # Patterns: Naive_WT_Rep1, 48h_WT_Rep1, D3_WT_Rep1, D8_WT_TE_Exp1_Rep1
    # Use [A-Za-z]+ for genotype to avoid matching underscores
    m = re.match(r"(Naive|48h|D\d+)_([A-Za-z]+)_(TE|EEC|MP)_(Exp\d+)_(Rep\d+)$", short)
    if m:
        return {
            "timepoint": m.group(1),
            "genotype": m.group(2),
            "cell_subset": m.group(3),
            "experiment": m.group(4),
            "replicate": m.group(5),
        }
    m = re.match(r"(Naive|48h|D\d+)_([A-Za-z]+)_(TE|EEC|MP)_(Rep\d+)$", short)
    if m:
        return {
            "timepoint": m.group(1),
            "genotype": m.group(2),
            "cell_subset": m.group(3),
            "experiment": "",
            "replicate": m.group(4),
        }
    m = re.match(r"(Naive|48h|D\d+)_([A-Za-z]+)_(Exp\d+)_(Rep\d+)$", short)
    if m:
        return {
            "timepoint": m.group(1),
            "genotype": m.group(2),
            "cell_subset": "total",
            "experiment": m.group(3),
            "replicate": m.group(4),
        }
    m = re.match(r"(Naive|48h|D\d+)_([A-Za-z]+)_(Rep\d+)$", short)
    if m:
        return {
            "timepoint": m.group(1),
            "genotype": m.group(2),
            "cell_subset": "total",
            "experiment": "",
            "replicate": m.group(3),
        }
    return {
        "timepoint": "",
        "genotype": "",
        "cell_subset": "",
        "experiment": "",
        "replicate": "",
    }


def normalize_genotype(raw):
    """Standardize genotype strings."""
    raw_lower = raw.lower().strip()
    if raw_lower in ("wt", "wild type"):
        return "WT"
    if "het" in raw_lower:
        return "Het"
    if "tbet" in raw_lower or "tbx21" in raw_lower:
        return "TbetKO"
    if "ko" in raw_lower or "arid1a ko" in raw_lower:
        return "KO"
    return raw


def build_master_sheet():
    """Build the master sample sheet from series matrices + SRA run info."""
    all_samples = []

    for gse, info in SUBSERIES.items():
        matrix_file = METADATA / f"{gse}_series_matrix.txt.gz"
        if not matrix_file.exists():
            print(f"WARNING: Missing {matrix_file}")
            continue

        samples = parse_series_matrix(matrix_file)
        assay = info["assay"]

        for s in samples:
            row = {
                "gse": gse,
                "gsm": s["gsm"],
                "srx": s.get("srx", ""),
                "title": s["title"],
                "assay": assay,
                "library_strategy": s.get("library_strategy", ""),
                "instrument": s.get("instrument", ""),
            }

            # Parse genotype from characteristics
            genotype_raw = s.get("char_genotype", s.get("char_cell line", ""))
            row["genotype"] = normalize_genotype(genotype_raw)

            # Parse antibody (CUT&RUN and ChIP-seq)
            row["antibody"] = s.get("char_antibody", s.get("char_chip antibody", ""))

            # Parse treatment
            row["treatment"] = s.get("char_treatment", "")

            # Parse timepoint and cell subset from title
            if assay == "rnaseq":
                parsed = parse_rnaseq_title(s["title"])
                row["timepoint"] = parsed["timepoint"]
                row["cell_subset"] = parsed["cell_subset"]
                row["replicate"] = parsed["replicate"]
            elif assay == "atacseq":
                parsed = parse_atacseq_title(s["title"])
                row["timepoint"] = parsed["timepoint"]
                row["cell_subset"] = parsed["cell_subset"]
                row["replicate"] = parsed["replicate"]
                if parsed.get("experiment"):
                    row["experiment"] = parsed["experiment"]
            elif assay == "atacseq_inhibitors":
                # Titles like "WT cells, 48h activated in vitro; Untreated [Control_Rep1]"
                row["timepoint"] = "48h"
                row["cell_subset"] = "total"
                # Extract replicate from bracket
                rep_match = re.search(r"Rep(\d+)", s["title"])
                row["replicate"] = f"Rep{rep_match.group(1)}" if rep_match else "Rep1"
            elif assay == "cutandrun":
                # Parse from title - various formats
                row["antibody"] = s.get("char_antibody", "").split(" (")[0]
                # Try to get timepoint from title
                title = s["title"]
                if "Naive" in title:
                    row["timepoint"] = "Naive"
                elif "48h" in title:
                    row["timepoint"] = "48h"
                elif "d5" in title.lower():
                    row["timepoint"] = "D5"
                elif "d8" in title.lower():
                    row["timepoint"] = "D8"
                # Cell subset
                if " TE " in title:
                    row["cell_subset"] = "TE"
                elif " EEC " in title:
                    row["cell_subset"] = "EEC"
                elif " MP " in title:
                    row["cell_subset"] = "MP"
                else:
                    row["cell_subset"] = "total"
                # Replicate
                rep_match = re.search(r"[Rr]ep\s*(\d+)", title)
                if rep_match:
                    row["replicate"] = f"Rep{rep_match.group(1)}"
                else:
                    row["replicate"] = "Rep1"
                # Tbet-OE
                if "Tbet-OE" in title:
                    row["treatment"] = "Tbet-OE"
            elif assay == "chipseq":
                row["antibody"] = s.get("char_chip antibody", "")
                row["treatment"] = s.get("char_treatment", "")
                row["timepoint"] = "48h"
                row["cell_subset"] = "total"
                # No replicates - single samples per condition
                row["replicate"] = "Rep1"

            # Set defaults for missing fields
            row.setdefault("timepoint", "")
            row.setdefault("cell_subset", "")
            row.setdefault("replicate", "")
            row.setdefault("treatment", "")
            row.setdefault("antibody", "")
            row.setdefault("experiment", "")

            all_samples.append(row)

    df = pd.DataFrame(all_samples)

    # Load SRA run info and merge SRR accessions
    sra = load_sra_runinfo()
    if not sra.empty:
        # Key mapping: SRX (Experiment) -> SRR (Run), LibraryLayout
        srx_to_srr = sra[["Experiment", "Run", "LibraryLayout", "spots", "bases"]].copy()
        srx_to_srr.columns = ["srx", "srr", "library_layout", "spots", "bases"]

        # Some SRX may have multiple SRR (technical replicates / multiple lanes)
        # Group them
        srr_grouped = srx_to_srr.groupby("srx").agg({
            "srr": lambda x: ";".join(sorted(x)),
            "library_layout": "first",
            "spots": "sum",
            "bases": "sum",
        }).reset_index()

        df = df.merge(srr_grouped, on="srx", how="left")
    else:
        df["srr"] = ""
        df["library_layout"] = ""
        df["spots"] = 0
        df["bases"] = 0

    # Helper to sanitize name parts (spaces -> _, remove +)
    def _clean(s):
        return s.replace(" + ", "_").replace(" ", "_").replace("+", "_")

    # Build a clean sample_name
    def make_sample_name(row):
        parts = []
        if row["assay"] == "rnaseq":
            parts = [row["timepoint"], row["genotype"]]
            if row["cell_subset"] != "total":
                parts.append(row["cell_subset"])
            parts.append(row["replicate"])
        elif row["assay"] == "atacseq":
            parts = [row["timepoint"], row["genotype"]]
            if row["cell_subset"] != "total":
                parts.append(row["cell_subset"])
            if row.get("experiment"):
                parts.append(row["experiment"])
            parts.append(row["replicate"])
        elif row["assay"] == "atacseq_inhibitors":
            # Treatment-based naming: ATAC_DMSO_Rep1, ATAC_IL12_Rep1, etc.
            tx = row.get("treatment", "Untreated")
            parts = ["ATAC", _clean(tx), row["replicate"]]
        elif row["assay"] == "cutandrun":
            parts = [row["timepoint"], row["genotype"], row["antibody"].split(" ")[0]]
            if row["cell_subset"] != "total":
                parts.append(row["cell_subset"])
            if row["treatment"]:
                parts.append(_clean(row["treatment"]))
            parts.append(row["replicate"])
        elif row["assay"] == "chipseq":
            tx = _clean(row.get("treatment", "Untreated"))
            ab = row["antibody"] if row["antibody"] and row["antibody"] != "none" else "Input"
            parts = ["ChIP", tx, ab]
        return "_".join([p for p in parts if p])

    df["sample_name"] = df.apply(make_sample_name, axis=1)

    # Reorder columns
    col_order = [
        "sample_name", "gsm", "srx", "srr", "gse", "assay",
        "genotype", "timepoint", "cell_subset", "treatment", "antibody",
        "replicate", "experiment", "library_strategy", "library_layout",
        "instrument", "spots", "bases", "title",
    ]
    col_order = [c for c in col_order if c in df.columns]
    df = df[col_order]

    # Sort
    df = df.sort_values(["assay", "gse", "timepoint", "genotype", "cell_subset", "replicate"]).reset_index(drop=True)

    return df


def write_nfcore_atacseq_samplesheet(master_df):
    """Generate nf-core/atacseq samplesheet CSV.

    nf-core/atacseq format: sample, fastq_1, fastq_2, replicate
    - 'sample' = condition group name (without RepN)
    - 'replicate' = integer replicate number (1, 2, 3...)
    - Rows with same sample + replicate are merged as technical replicates
    """
    atac = master_df[master_df["assay"].isin(["atacseq", "atacseq_inhibitors"])].copy()
    rows = []
    for _, r in atac.iterrows():
        srr_list = str(r.get("srr", "")).split(";") if pd.notna(r.get("srr")) else []
        fastq_dir = PROJECT / "data" / "fastq" / "atacseq"

        # Extract integer replicate from RepN string
        rep_match = re.search(r"Rep(\d+)", str(r.get("replicate", "Rep1")))
        rep_int = int(rep_match.group(1)) if rep_match else 1

        # Build condition name (sample name without RepN)
        # e.g., D3_WT, D8_WT_EEC_Exp1, ATAC_Untreated, ATAC_IL-12_ACBI1
        sample_name = r["sample_name"]
        # Strip trailing _RepN from sample name
        sample_group = re.sub(r"_Rep\d+$", "", sample_name)

        for srr in srr_list:
            if r.get("library_layout") == "PAIRED":
                rows.append({
                    "sample": sample_group,
                    "fastq_1": str(fastq_dir / f"{srr}_1.fastq.gz"),
                    "fastq_2": str(fastq_dir / f"{srr}_2.fastq.gz"),
                    "replicate": rep_int,
                })
            else:
                rows.append({
                    "sample": sample_group,
                    "fastq_1": str(fastq_dir / f"{srr}.fastq.gz"),
                    "fastq_2": "",
                    "replicate": rep_int,
                })
    out = pd.DataFrame(rows)
    out.to_csv(SAMPLESHEETS / "atacseq_samplesheet.csv", index=False)
    return out


def write_nfcore_rnaseq_samplesheet(master_df):
    """Generate nf-core/rnaseq samplesheet CSV."""
    rna = master_df[master_df["assay"] == "rnaseq"].copy()
    rows = []
    for _, r in rna.iterrows():
        srr_list = str(r.get("srr", "")).split(";") if pd.notna(r.get("srr")) else []
        fastq_dir = PROJECT / "data" / "fastq" / "rnaseq"

        for srr in srr_list:
            if r.get("library_layout") == "PAIRED":
                rows.append({
                    "sample": r["sample_name"],
                    "fastq_1": str(fastq_dir / f"{srr}_1.fastq.gz"),
                    "fastq_2": str(fastq_dir / f"{srr}_2.fastq.gz"),
                    "strandedness": "reverse",
                })
            else:
                rows.append({
                    "sample": r["sample_name"],
                    "fastq_1": str(fastq_dir / f"{srr}.fastq.gz"),
                    "fastq_2": "",
                    "strandedness": "reverse",
                })
    out = pd.DataFrame(rows)
    out.to_csv(SAMPLESHEETS / "rnaseq_samplesheet.csv", index=False)
    return out


def write_nfcore_cutandrun_samplesheet(master_df):
    """Generate nf-core/cutandrun samplesheet CSV.

    nf-core/cutandrun format: group, replicate, fastq_1, fastq_2, control
    - 'group' = antibody_condition identifier (e.g., ARID1A_D5_WT, IgG_D5_WT)
    - 'replicate' = integer replicate number
    - 'control' = group name of IgG control (empty for IgG rows themselves)
    - IgG controls matched to targets by condition (timepoint+genotype)
    """
    cr = master_df[master_df["assay"] == "cutandrun"].copy()

    # Build condition key for each sample
    def condition_key(r):
        parts = [r["timepoint"], r["genotype"]]
        if r["cell_subset"] != "total":
            parts.append(r["cell_subset"])
        return "_".join(p for p in parts if p)

    # Build IgG group name lookup: condition -> IgG group name
    # IgG controls exist for D5_WT and D5_KO
    igg_groups = {}
    for _, r in cr.iterrows():
        ab = str(r.get("antibody", "")).strip()
        if ab == "IgG":
            cond = condition_key(r)
            igg_groups[cond] = f"IgG_{cond}"

    rows = []
    for _, r in cr.iterrows():
        srr_list = str(r.get("srr", "")).split(";") if pd.notna(r.get("srr")) else []
        fastq_dir = PROJECT / "data" / "fastq" / "cutandrun"
        ab = str(r.get("antibody", "")).strip()
        cond = condition_key(r)
        is_igg = ab == "IgG"

        # Group name: antibody_condition
        # For Tbet-OE samples, include treatment in condition
        treatment = str(r.get("treatment", "")).strip()
        if treatment and treatment != "nan":
            group_name = f"{ab}_{cond}_{treatment.replace('-', '')}"
        else:
            group_name = f"{ab}_{cond}"

        # Extract integer replicate
        rep_match = re.search(r"Rep(\d+)", str(r.get("replicate", "Rep1")))
        rep_int = int(rep_match.group(1)) if rep_match else 1

        # Control: IgG rows leave control empty; target rows reference IgG group
        if is_igg:
            control = ""
        else:
            # Find matching IgG control by condition
            # First try exact match, then fall back to broader match (same timepoint+genotype)
            base_cond = f"{r['timepoint']}_{r['genotype']}"
            control = igg_groups.get(cond, igg_groups.get(base_cond, ""))

        for srr in srr_list:
            rows.append({
                "group": group_name,
                "replicate": rep_int,
                "fastq_1": str(fastq_dir / f"{srr}_1.fastq.gz"),
                "fastq_2": str(fastq_dir / f"{srr}_2.fastq.gz"),
                "control": control,
            })
    out = pd.DataFrame(rows)
    out.to_csv(SAMPLESHEETS / "cutandrun_samplesheet.csv", index=False)
    return out


def write_nfcore_chipseq_samplesheet(master_df):
    """Generate nf-core/chipseq samplesheet CSV.

    nf-core/chipseq format: sample, fastq_1, fastq_2, replicate, antibody, control, control_replicate
    - Input control rows: antibody, control, control_replicate all empty
    - IP rows: antibody set, control = input sample name, control_replicate = integer
    """
    chip = master_df[master_df["assay"] == "chipseq"].copy()

    # First, build a mapping from treatment -> input sample name
    input_names = {}
    for _, r in chip.iterrows():
        is_input = str(r.get("antibody", "")).lower() == "none"
        if is_input:
            treatment = str(r.get("treatment", "Untreated")).strip()
            input_names[treatment] = r["sample_name"]

    rows = []
    for _, r in chip.iterrows():
        srr_list = str(r.get("srr", "")).split(";") if pd.notna(r.get("srr")) else []
        fastq_dir = PROJECT / "data" / "fastq" / "chipseq"

        is_input = str(r.get("antibody", "")).lower() == "none"
        treatment = str(r.get("treatment", "Untreated")).strip()

        for srr in srr_list:
            if r.get("library_layout") == "PAIRED":
                f1 = str(fastq_dir / f"{srr}_1.fastq.gz")
                f2 = str(fastq_dir / f"{srr}_2.fastq.gz")
            else:
                f1 = str(fastq_dir / f"{srr}.fastq.gz")
                f2 = ""

            if is_input:
                rows.append({
                    "sample": r["sample_name"],
                    "fastq_1": f1,
                    "fastq_2": f2,
                    "replicate": 1,
                    "antibody": "",
                    "control": "",
                    "control_replicate": "",
                })
            else:
                # Look up matching input control by treatment
                control_name = input_names.get(treatment, "")
                rows.append({
                    "sample": r["sample_name"],
                    "fastq_1": f1,
                    "fastq_2": f2,
                    "replicate": 1,
                    "antibody": "Tbet",
                    "control": control_name,
                    "control_replicate": 1,
                })
    out = pd.DataFrame(rows)
    out.to_csv(SAMPLESHEETS / "chipseq_samplesheet.csv", index=False)
    return out


if __name__ == "__main__":
    print("Building master sample sheet...")
    master = build_master_sheet()

    # Save master sheet
    master_path = METADATA / "master_sample_sheet.tsv"
    master.to_csv(master_path, sep="\t", index=False)
    print(f"\nMaster sample sheet: {master_path}")
    print(f"Total samples: {len(master)}")
    print(f"\nSamples per assay:")
    print(master["assay"].value_counts().to_string())
    print(f"\nSamples with SRR accessions: {master['srr'].notna().sum()}")
    print(f"Samples missing SRR: {master['srr'].isna().sum()}")

    # Show library layout
    if "library_layout" in master.columns:
        print(f"\nLibrary layout:")
        print(master.groupby("assay")["library_layout"].first().to_string())

    # Generate nf-core samplesheets
    print("\n--- Generating nf-core samplesheets ---")

    atac_ss = write_nfcore_atacseq_samplesheet(master)
    print(f"ATAC-seq samplesheet: {len(atac_ss)} entries")

    rna_ss = write_nfcore_rnaseq_samplesheet(master)
    print(f"RNA-seq samplesheet: {len(rna_ss)} entries")

    cr_ss = write_nfcore_cutandrun_samplesheet(master)
    print(f"CUT&RUN samplesheet: {len(cr_ss)} entries")

    chip_ss = write_nfcore_chipseq_samplesheet(master)
    print(f"ChIP-seq samplesheet: {len(chip_ss)} entries")

    # Print a summary table
    print("\n\n=== SAMPLE SUMMARY ===")
    for assay in master["assay"].unique():
        sub = master[master["assay"] == assay]
        print(f"\n--- {assay.upper()} ({sub['gse'].iloc[0]}) ---")
        print(f"  Samples: {len(sub)}")
        if "genotype" in sub.columns:
            print(f"  Genotypes: {', '.join(sorted(sub['genotype'].unique()))}")
        if "timepoint" in sub.columns:
            print(f"  Timepoints: {', '.join(sorted(sub['timepoint'].unique()))}")
        if sub["antibody"].any():
            abs_list = sorted(sub[sub["antibody"] != ""]["antibody"].unique())
            if abs_list:
                print(f"  Antibodies: {', '.join(abs_list)}")
        if sub["treatment"].any():
            tx_list = sorted(sub[sub["treatment"] != ""]["treatment"].unique())
            if tx_list:
                print(f"  Treatments: {', '.join(tx_list)}")
