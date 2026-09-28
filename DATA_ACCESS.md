# Data access

This repository contains code, sample sheets and small metadata tables. Sequencing data are in public archives. The processed tables needed to rebuild every figure are in a separate data bundle.

## 1. Processed data bundle (quick start)

| Item | Location |
|------|----------|
| Compact processed tables: peak/gene count matrices, peak sets, motif and footprinting tables, deepTools matrices | Zenodo: `10.5281/zenodo.XXXXXXX` *(DOI assigned at release)* |
| Checksums | `bundle/SHA256SUMS` (verified by `make fetch-data`) |

`make fetch-data` downloads the bundle into `bundle/`, checks each file against its SHA-256 checksum, and unpacks it into the `results/` layout the analysis scripts expect.

## 2. Raw sequencing data: McDonald, Chick et al. 2023

GEO SuperSeries **[GSE228381](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228381)**

| Sub-series | Assay | Samples | Contents |
|------------|-------|--------:|----------|
| [GSE228171](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228171) | ATAC-seq | 54 | WT / Arid1a Het / KO / Tbx21 KO; naive, 48 h, d3, d5, d8 (TE/EEC/MP) |
| [GSE228193](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228193) | ATAC-seq | 8 | In vitro IL-12 ± ACBI1 / BRM014 |
| [GSE227634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227634) | RNA-seq | 30 | WT / Het / KO; d3, d5, d8 subsets |
| [GSE228380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228380) | CUT&RUN | 28 | ARID1A, H3K27ac, H3K27me3, T-bet (± T-bet OE), BATF, ETS1, IgG |
| [GSE228546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228546) | ChIP-seq | 8 | T-bet + input; IL-12 ± ACBI1 / BRM014 |

Run-level metadata: `metadata/master_sample_sheet.tsv`. ENA download URLs with MD5 checksums: `metadata/ena_fastq_urls.tsv`. The download script is `scripts/upstream/download_fastqs.sh`.

## 3. Raw sequencing data: Guo et al. 2022 (cross-study comparison)

Guo A, *et al.* cBAF complex components and MYC cooperate early in CD8+ T cell fate. *Nature* 607(7917), 135–141 (2022). [doi:10.1038/s41586-022-04849-0](https://doi.org/10.1038/s41586-022-04849-0) · PMID 35732731. GEO SuperSeries GSE183619.

| Sub-series | Assay | Contents |
|------------|-------|----------|
| GSE183615 | RNA-seq | Arid1a WT/KO; c-Myc WT/KO |
| GSE199184 | RNA-seq | Vehicle vs Arid1a inhibitor (BRD-K98645985) |
| GSE183618 | ATAC-seq | WT / Arid1a KO / Myc KO. **The "WT" samples come from two separate experiments; never pool them.** |
| GSE198894 | ATAC-seq | DMSO WT / DMSO Arid1a KO / inhibitor |
| GSE183616 | ATAC-seq | Naive / Myc-high / Myc-low |

Run-level metadata: `metadata/meta_analysis/meta_master_sample_sheet.tsv`. The download script is `scripts/upstream/download_meta_fastqs.sh`.

## 4. Reference genome

| File | Source |
|------|--------|
| `GRCm39.primary_assembly.genome.fa` | [GENCODE mouse release M35](https://www.gencodegenes.org/mouse/release_M35.html) |
| `gencode.vM35.primary_assembly.annotation.gtf` | GENCODE mouse release M35 |
| `mm39-blacklist.v2.bed` | [Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist) (lifted to mm39) |

Put these files in `data/reference/`. Only the upstream (raw-data) mode needs them. Three scripts also need the FASTA and skip without it: `het_dose_response/02_feature_enrichment.R`, `motif_grammar/01_region_sets.R` and `chromvar/build_chromvar.R` (all in `extended_analysis/scripts/`).

## 5. Not reproduced from deposited data

These panels in the paper used external data that is not part of GSE228381. `PANEL_MAP.md` lists them with the reason for each.

- Fig 1A H3K4me1 track and naive H3K27ac (GSE89036); Fig 1E enhancer classes (H3K4me1/me3).
- Fig 1G and 5B: public TF ChIP-seq (GSE54191, GSE192390, GSE166718).
- Fig 3I: custom signatures from GSE10739, PRJNA547650, GSE143504, GSE81888 and GSE173515.
- Fig S1B: naive RNA-seq (GSE152841).

## Terms

The code is under the MIT license (`LICENSE`). The raw data are subject to the terms of their GEO/SRA depositions. The processed data bundle is released under CC BY 4.0.
