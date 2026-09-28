# Panel map: McDonald, Chick et al. 2023, *Immunity* 56:1303

This table lists every panel in the paper and supplement against the script that reproduces it, the output it writes, and how closely it matches the published version. Figures go to `figures/paper/<id>.{pdf,png}` and the numbers behind each panel go to `results/paper/`. Reference renders of each figure are in `results/reference/paper/`.

**Status key**

| Mark | Meaning |
|------|---------|
| ✅ | Reproduced. Same analysis and same conclusion; the counts agree closely. |
| 🟡 | Reproduced with a documented difference. The conclusion holds, but counts or thresholds differ (see Notes). |
| ⛔ | Not reproducible from the deposited data (reason given). |
| 🧪 | Wet-lab or flow cytometry panel; not computational. |

Scripts are in `scripts/paper/`. `*.sh` scripts need raw-data inputs (bigWigs, BAMs, genome FASTA, HOMER). Their compact outputs are in the data bundle, so the `*.R` scripts that draw the panels run without them.

## Figure 1: dynamic accessibility and cBAF occupancy

The OCR clusters are the **published cluster definitions**: the original mm10 regions, lifted to mm39. Details are in `data/metadata/paper_ocr_clusters/README.md`. The cluster sizes match the paper exactly: Conserved 21,650 (1,500 drawn in 1A), Naive 3,157, Early Activation 7,883, Activation 8,781, Late Activation 1,250.

| Panel | Content | Script | Status | Notes |
|---|---|---|---|---|
| 1A | ATAC, ARID1A and H3K27ac signal by OCR cluster | `fig1_published_clusters.R`, `fig1_fig2_deeptools.sh` | 🟡 | Tracks from GSE89036 are external and omitted: H3K4me1 and naive H3K27ac. The 48h ATAC track is shown for fidelity, but those libraries failed QC (5.0M/4.5M read pairs, FRiP 0.314/0.115), so it is excluded from all statistics. |
| 1B | ATAC profiles per cluster | `fig1_fig2_panels.R` | ✅ | Naive, 48h, d3, d5, d8 |
| 1C | Genomic annotation per cluster | `fig1_fig2_panels.R` | 🟡 | Intron fractions match (46–52% vs 48–52%). Promoter fractions are higher: GENCODE vM35 has more TSSs than HOMER's mm10 annotation. |
| 1D | ARID1A CUT&RUN overlap, observed/expected | `fig1d_arid1a_overlap.sh`, `fig1_fig2_panels.R` | 🟡 | Uses HOMER `mergePeaks -matrix`, as in the paper, and the pattern matches. ARID1A peaks were called without a control, because only the D5 libraries have a matched IgG. The 48h column rests on 58 peaks. |
| 1E | Enhancer classes (H3K4me1/me3, H3K27me3) | — | ⛔ | Needs external histone data (GSE89036) |
| 1F | HOMER known motifs per cluster | `fig1_homer.sh`, `fig1_fig2_panels.R` | ✅ | All 23 motifs shown in the paper were recovered, with the same family pattern |
| 1G | Public TF ChIP overlap | — | ⛔ | Needs external ChIP data (GSE54191, GSE192390, GSE166718) |
| S1A | Signal tracks: Zeb2, Tbx21, Bhlhe40 | `figS_tracks.R` | ✅ | |
| S1B | Expression of cluster-annotated genes | `fig1_fig2_panels.R` | 🟡 | The naive RNA-seq is external (GSE152841), so only d3/d5/d8 are shown |
| S1C | Top 15 motifs per cluster | `fig1_fig2_panels.R` | ✅ | |
| (new) | Published vs de novo clusters | `figS1_denovo_clusters.R`, `fig1_published_clusters.R` | — | Validation only. Share of each published cluster recovered in the same de novo cluster: Naive 90%, Late Activation 84%, Conserved 77%, Activation 62%, Early Activation 46%. Most of the Early Activation disagreement comes from the missing 48h column. |

## Figure 2: ARID1A and activation-induced enhancers

| Panel | Content | Script | Status | Notes |
|---|---|---|---|---|
| 2A, 2B | Cell numbers; CTV dilution | — | 🧪 | |
| 2C | % of each cluster lost in KO at d3/d5/d8 | `fig1_fig2_panels.R` | 🟡 | d3 agrees closely (Activation 36% vs ~34%). At d8, more loss is found than in the paper, e.g. Activation 54% vs ~32%. |
| 2D, 2E | ATAC heatmaps and profiles, WT/Het/KO | `fig1_fig2_deeptools.sh`, `fig1_fig2_panels.R` | ✅ | D8 uses the Exp2 batch (the only batch with WT, Het and KO) |
| 2F | H3K27ac, WT vs KO | `fig1_fig2_panels.R` | 🟡 | D5 WT vs KO is shown (−17% vs ~−18% in the paper). |
| 2G | D3 DEG counts | `fig2g_2i.R` | 🟡 | See **DEG thresholds** below |
| 2H | D3 curated heatmap | `fig2g_2i.R` | ✅ | All 61 genes |
| 2I | D3 Hallmark GSEA | `fig2g_2i.R` | ✅ | The same 9 gene sets in the same order |

## Figure 3: dose-dependent effector gene expression

| Panel | Content | Script | Status | Notes |
|---|---|---|---|---|
| 3A–D | Flow cytometry | — | 🧪 | |
| 3E | RNA PCA | `fig3.R` | ✅ | Same layout; PC1/PC2 = 38.7/14.7% (paper 49.4/18.8%) |
| 3F | D8 DEG counts per subset | `fig3.R` | 🟡 | See **DEG thresholds** |
| 3G | Curated heatmap | `fig3.R` | ✅ | 62/62 genes, in the paper's 9 gene clusters |
| 3H | Volcanoes with TE/MP signature genes | `fig3.R` | ✅ | |
| 3I | GSEA against published signatures | `fig3.R` | 🟡 | Only the signatures that can be rebuilt are used. MP-vs-TE matches the paper in TE and EEC. The legend's "GSE10739" appears to be a typo: that accession is an unrelated human monocyte series, so GSE10239 (Sarkar 2008) is used. The TbetKO, Batf3OE, Runx3KO and BRD4KO signatures need external data (⛔). |

## Figure 4: ARID1A-dependent OCRs at day 8

| Panel | Content | Script | Status | Notes |
|---|---|---|---|---|
| 4A | ATAC PCA | `fig4.R` | ✅ | |
| 4B | Lost/gained OCRs per subset | `fig4.R` | ✅ | KO lost TE 10,804 / EEC 14,365 / MP 15,547 |
| 4C | Heatmaps by lost-OCR category | `fig4_fig5_signal.sh`, `fig4_fig5_profiles.R` | ✅ | WT > Het > KO in every subset |
| 4D | UpSet of lost OCRs | `fig4.R` | ✅ | Lost in all three subsets: 8,416 (paper ~7,700) |
| 4E | ATAC log2FC vs RNA log2FC | `fig4.R` | ✅ | TE and MP signatures defined as WT TE vs WT MP DEGs |

## Figure 5: cBAF and T-bet targeting

| Panel | Content | Script | Status | Notes |
|---|---|---|---|---|
| 5A | Motif families in lost/gained OCRs | `fig5a_homer.sh`, `fig5a_motifs.R` | ✅ | ETS, Runt and T-box at the >500 cap in lost OCRs |
| 5B | In vitro TF ChIP overlap | — | ⛔ | Needs external ChIP data (GSE192390) |
| 5C | Annotation and CUT&RUN at ARID1A-dependent OCRs | `fig5.R`, `fig4_fig5_signal.sh` | 🟡 | 5,453 dependent / 12,700 independent (paper 9,965 / 21,375). About 30% are dependent in both analyses; the absolute counts are lower because MACS2 calls fewer ARID1A peaks than SEACR. ARID1A, BATF, ETS1 and T-bet binding all collapse in KO. |
| 5D–F | In vitro ATAC and T-bet ChIP ± ACBI1/BRM014 | `fig5.R`, `fig4_fig5_signal.sh`, `fig4_fig5_profiles.R` | 🟡 | See **Inhibitor normalization** below |
| 5G | Arid1a KO vs Tbx21 KO lost-OCR overlap | `fig5.R` | 🟡 | Most OCRs lost in Tbx21 KO are also lost in Arid1a KO, as in the paper. More Tbx21-KO-only OCRs are found. |
| 5H | ARID1A CUT&RUN in Tbx21 KO | `fig5.R`, `fig4_fig5_signal.sh` | ✅ | |
| 5I | T-bet overexpression and %TE | — | 🧪 | |
| 5J | T-bet CUT&RUN ± T-bet overexpression | `fig4_fig5_signal.sh`, `fig4_fig5_profiles.R` | ✅ | Overexpression raises T-bet 2.8× in WT and has no effect in KO |
| S5A | Tbx21 KO DA counts | `fig5.R` | ✅ | |

## Figure 6: ARID1A and Trm formation

| Panel | Content | Script | Status | Notes |
|---|---|---|---|---|
| 6A | PageRank TF scores vs mRNA, D8 MP WT vs KO | `fig6a_pagerank.R` | 🟡 | Taiji-style personalised PageRank, reimplemented in R. Peaks are linked to genes by promoter overlap plus a ±50 kb distance decay (GREAT-style, as Taiji does without Hi-C). The paper's hits rank near the top of 251 TFs: BHLHE40 2, EOMES 9, RXRA 7, SMAD3 13, TBX21 15, RUNX3 23; the KO-high TFs rank 235–247. ZFP683 has no JASPAR2020 motif and is not scored. |
| 6B–I | Memory, Trm and recall | — | 🧪 | |
| S6A | Top PageRank ratios | `fig6a_pagerank.R` | ✅ | |
| S2A–C, S2D | Deletion and titres; tracks | `figS_tracks.R` (S2D) | 🧪 / ✅ | |
| S3 | Flow cytometry | — | 🧪 | |
| S4A, S4B | Signal tracks | `figS_tracks.R` | 🟡 | H3K27ac at D8 and H3K27me3 tracks are WT only |
| S5B, S6B–H | Flow cytometry | — | 🧪 | |

## Cross-cutting differences

**DEG thresholds (2G, 3F).** The legends say ">2-fold", but the Methods specify |log2FC| ≥ 0.585 (1.5-fold) with HOMER/DESeq2 on unshrunken fold changes. Both thresholds are reported here:

| D3, WT-high / KO-high | Paper | 2-fold (legend) | 1.5-fold (Methods) |
|---|---|---|---|
| DEGs | 534 / 595 | 376 / 516 | 1,105 / 1,283 |

The direction and relative sizes agree in every comparison. The per-subset tables in `results/paper/` carry both thresholds.

**Inhibitor normalization (5D–F).** The inhibitor-treated libraries have lower FRiP: 0.134/0.143 with IL-12, 0.117/0.113 with ACBI1, 0.088/0.117 with BRM014. DESeq2's median-of-ratios normalization absorbs this global loss. The primary analysis instead normalizes to total counted reads, which is closer to the paper's HOMER total-tag normalization. Both are reported (`fig5_norm_sensitivity`):

| Lost vs DMSO + IL-12 (2-fold, padj < 0.05) | ACBI1 | BRM014 |
|---|---|---|
| Total-read normalization (primary) | 2,296 | 1,963 |
| Median-of-ratios normalization | 283 | 127 |
| Paper (ACBI1-dependent, T-bet-bound OCRs) | 7,487 | — |

T-bet ChIP reproduces the key result. IL-12 raises T-bet binding at ACBI1-dependent OCRs (1.63 → 3.36), and ACBI1 blocks it (1.24). Neither normalization reaches the paper's 7,487. The global FRiP loss is consistent with BAF disruption, but it is confounded with library quality, and each condition has only two replicates.

**Genome and pipelines.** The paper used mm10 with STAR and HOMER. This reanalysis uses mm39 (GENCODE vM35), nf-core pipelines, MACS2/SEACR and DESeq2 with apeglm. The full comparison is in `docs/methods.md` §10.
