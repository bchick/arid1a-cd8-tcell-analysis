# Published Fig 1A OCR cluster definitions

These are the ATAC-seq open-chromatin-region (OCR) clusters shown in Figure 1A of
McDonald, Chick et al. 2023, *Immunity* 56:1303–1319
(doi:10.1016/j.immuni.2023.05.005). The authors supplied them from the original
analysis: HOMER `findPeaks -style dnase` peaks on mm10, with clusters defined as
described in the paper's STAR Methods.

They are the **canonical** cluster definitions for the Fig 1 and Fig 2 panel
scripts in `scripts/paper/`. `scripts/paper/figS1_denovo_clusters.R` re-derives
clusters from the reprocessed data only to validate them
(`figures/paper/figS1_cluster_concordance.*`).

## Files

| File | Cluster | mm10 regions | mm39 regions | Lost in liftOver |
|---|---|---:|---:|---:|
| `conserved.allSamps.peaks.TotalSet.bed` | Conserved (all) | 21,650 | 21,649 | 1 |
| `conserved.allSamps.plottedSubSet.peaks.bed` | Conserved (the 1,500 drawn in the Fig 1A heatmap) | 1,500 | 1,500 | 0 |
| `naive.specific.sig.bed` | Naive | 3,157 | 3,157 | 0 |
| `early.activation.specific.sig.bed` | Early Activation | 7,883 | 7,881 | 2 |
| `activation.specific.sig.bed` | Activation | 8,781 | 8,780 | 1 |
| `late.activation.specific.sig.bed` | Late Activation | 1,250 | 1,250 | 0 |

The mm10 counts match the cluster sizes printed in Fig 1A. The files are BED5:
chr, start (0-based), end, HOMER peak id, strand ".". Peak ids are unique within
the full HOMER peak set. One region, `Peak_76948`, is listed in both the Early
Activation and Activation files in the original mm10 data; it is kept in both.

- `mm10/`: the original files as supplied.
- `mm39/`: lifted to GRCm39 with UCSC `liftOver -minMatch=0.95` using the UCSC
  chain `mm10ToMm39.over.chain.gz`. Regions that failed to lift were dropped.

## Use in this repository

`scripts/paper/fig1_published_clusters.R` reads the mm39 files and writes:

- `results/paper/fig1a_ocr_clusters.csv`: one row per region × cluster. It
  includes the overlapping nf-core consensus ATAC peak(s), which are used to
  look up the reprocessed differential-accessibility results (Fig 2C, S1B).
- `results/paper/fig1_regions/*.bed`: region sets for the signal and motif
  panels (Fig 1A–D, 1F, S1C, 2D–F).
- `results/paper/fig1a_consensus_mapping.csv`: the region-to-consensus-peak
  mapping rate per cluster (≥99% for every cluster).
