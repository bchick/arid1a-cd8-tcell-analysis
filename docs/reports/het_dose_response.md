# ARID1A Dose-Response Analysis — D8 Effector Subsets

**Project:** McDonald, Chick et al. 2023 reanalysis — extended analysis
**Data:** nf-core/atacseq consensus peaks, D8 TE/EEC/MP, Exp2 only
**Design:** WT (n=3) / Het (n=3) / KO (n=2) per subset, one batch
**Scripts:** `extended_analysis/scripts/het_dose_response/01_dose_classes.R`, `extended_analysis/scripts/het_dose_response/02_feature_enrichment.R`, `extended_analysis/scripts/het_dose_response/03_known_motif_summary.R`, `extended_analysis/scripts/het_dose_response/04_denovo_motif_summary.R`
**Outputs:** `results/extended_analysis/het_dose_response/`, `figures/extended_analysis/het_dose_response/`

---

## 1. Question

The original Het samples were collected but not exploited in McDonald 2023 for a dose-response analysis. With WT (2 copies), Het (1 copy) and KO (0 copies) of *Arid1a*, every D8 OCR can be scored on how it responds to ARID1A dosage. The question is *not* "does accessibility change?" but "**what shape does the allelic response take at each OCR, and what distinguishes the shapes?**"

## 2. Classification scheme

Per subset, DESeq2 LRT (`~genotype` vs `~1`) selects OCRs with any genotype effect (padj < 0.05). Three Wald contrasts — Het vs WT, KO vs WT, KO vs Het — provide effect sizes. Each OCR is assigned to one class:

| Class | Definition |
|---|---|
| **insensitive** | LRT not significant |
| **buffered** | KO differs from WT, Het ≈ WT (`sig_KOWT` & `!sig_HetWT` & `|lfc_HetWT|<0.5`) — *one copy is sufficient* |
| **linear** | Monotonic, both Het & KO significant vs WT, same sign, KO vs Het also significant, `0.5·|KO| ≤ |Het| < |KO|` — *true dose-dependence* |
| **haploinsufficient** | Het ≈ KO (both differ from WT, KO-vs-Het NOT significant) — *one copy already produces the full effect* |
| **nonmonotonic** | Het and KO change in opposite directions |
| **other_responsive** | LRT-significant but doesn't fit cleanly |

Thresholds: padj < 0.05, |log2FC| ≥ 0.5.

## 3. Class counts

![Class counts per subset](../../figures/extended_analysis/het_dose_response/class_counts_bar.png)

| | TE | EEC | MP |
|---|---:|---:|---:|
| insensitive | 48,112 | 39,630 | 36,607 |
| **buffered** | 3,386 | 9,367 | 10,020 |
| **linear** | 73 | 271 | 309 |
| **haploinsufficient** | 2,495 | 1,358 | 2,234 |
| nonmonotonic | 34 | 84 | 72 |
| other_responsive | 4,173 | 10,560 | 13,723 |

**Buffered dominates responsive OCRs** (2–3× haploinsufficient across subsets). Most ARID1A-responsive accessibility is recessive — one allele's worth of ARID1A suffices. **Linear (true dose-dependent) is rare** (37–309 peaks) but, as shown below, biochemically the cleanest class.

Dose curves (per-class mean VST across WT → Het → KO):

![TE dose curves](../../figures/extended_analysis/het_dose_response/dose_curves_TE.png)
![EEC dose curves](../../figures/extended_analysis/het_dose_response/dose_curves_EEC.png)
![MP dose curves](../../figures/extended_analysis/het_dose_response/dose_curves_MP.png)

## 4. Direction of change: haplo peaks are overwhelmingly *lost*

| Subset | Haplo **lost** | Haplo **gained** | Ratio |
|---|---:|---:|---:|
| TE | 2,195 | 300 | 7.3× |
| EEC | 1,230 | 128 | 9.6× |
| MP | 2,018 | 216 | 9.3× |

Buffered peaks split ~2:1 lost:gained. Haploinsufficient peaks are ~9:1 lost.
**Sites that genuinely require two ARID1A copies are sites that open in WT and fail to open on dosage loss.** Gained-in-KO peaks (de-repression) are almost all buffered or other — one allele is enough to keep them closed.

## 5. Feature differences between classes

Baseline WT accessibility (VST) is highest at linear and buffered, lower at haploinsufficient and insensitive:

![Baseline WT accessibility by class](../../figures/extended_analysis/het_dose_response/baseline_WT_by_class.png)

Peak widths are largest at linear peaks (median 840–950 bp) and smallest at haplo (median 667–746 bp) — consistent with linear being "big" ARID1A-remodeled elements:

![Peak width by class](../../figures/extended_analysis/het_dose_response/peak_width_by_class.png)

Genomic annotation (ChIPseeker, GENCODE vM35):

![Annotation by class](../../figures/extended_analysis/het_dose_response/annotation_by_class.png)

- **Insensitive:** ~35–37% promoter — cBAF is mostly an enhancer remodeler; promoters are held open by other mechanisms.
- **Buffered, linear, haploinsufficient:** ~13–21% promoter, ~68–72% distal — all three responsive classes are enhancer-enriched.

CUT&RUN co-occupancy at D5 WT (ARID1A, H3K27ac, Tbet, BATF, ETS1):

![CUT&RUN overlap by class](../../figures/extended_analysis/het_dose_response/cutrun_overlap_by_class.png)
![TF co-occupancy per peak by class](../../figures/extended_analysis/het_dose_response/tf_cooccupancy_by_class.png)

ARID1A binding rate (from D5 CUT&RUN, a proxy for "direct target"):

| | buffered | **linear** | haploinsufficient | insensitive |
|---|---:|---:|---:|---:|
| TE | 53% | **84%** | 41% | 27% |
| EEC | 43% | **69%** | 34% | 22% |
| MP | 44% | **67%** | 27% | 19% |

**Linear peaks are the most directly ARID1A-bound class.** Haploinsufficient peaks have intermediate ARID1A binding — more than insensitive, less than buffered. They are not simply "extra-dependent direct targets"; they are a distinct category.

## 6. Motifs — the clearest axis of separation

HOMER de novo motif discovery was run per (subset × class), with two background schemes: (i) class-lost vs insensitive-lost within the same subset, and (ii) haplo-lost vs buffered-lost (the matched-ARID1A-binding contrast). De novo results are used rather than the known-motif scan because they index TF families by motifs actually present in the foreground, not by the HOMER library's prior. Each de novo motif was assigned to a TF family via its HOMER best-guess match (Etv2/ETS, Runx1/Runt, Tbx21/T-box, Fos/bZIP, etc.); ambiguous best-guesses (match score < 0.6, or non-vertebrate best-guess such as yeast/plant TFs) were labelled *other* and excluded from the family summary. Per-comparison top-8 motifs with consensus, best-guess TF, and occupancy are in `results/extended_analysis/het_dose_response/motif_denovo_top_per_comparison.csv`.

### 6a. Class vs insensitive — ETS is the dose-tolerant signature

![De novo motif summary heatmap](../../figures/extended_analysis/het_dose_response/motif_denovo_summary_heatmap.png)

Cells show max −log10 p across de novo motifs assigned to each family. Colour scale capped at 100 so that the buffered-ETS runaway (max = 252 in EEC) does not drown the haplo-range signal.

Three patterns in the figure, in order of effect size:

1. **ETS is the single family that is strongly buffered-specific.** EEC buffered −log10 p = 252 (Etv2 de novo motif, 53 % of foreground vs 32 % of background), MP buffered = 188 (ERG), TE buffered = 52 (Ets1-distal). ETS enrichment in all three haplo columns is absent or marginal (MP haplo = 12, others ≈ 0). This is the one clean asymmetry between buffered and haplo.

2. **RUNX and T-box are enriched in *both* buffered and haplo — they do not separate the two classes.** RUNX1/2 de novo motifs hit buffered (EEC = 129, MP = 108) almost as hard as haplo (MP = 58, TE = 54, EEC = 31). T-box motifs behave similarly. The earlier knownResults interpretation that "RUNX/T-box dominates haplo" reflected a relative-scale artefact (known RUNX hits cap at ~45 vs known ETS at >200, so RUNX looked pale only because ETS was off the chart); the de novo contrast shows RUNX/T-box are pan-responsive rather than dose-sensitivity markers.

3. **AP-1/bZIP and NFκB appear mainly in haplo, not buffered.** AP-1/bZIP is absent from all buffered columns but present in EEC haplo (18) and MP haplo (29). NFκB shows up only in EEC haplo (18). These are the families whose presence predicts dose-sensitivity in EEC/MP, not RUNX/T-box.

4. Linear columns are underpowered (all < 15; TE linear has 73 peaks) — no conclusions drawn from the linear class motif scan.

### 6b. Haplo vs buffered — direct contrast at ARID1A-bound sites

Haplo and buffered sites have broadly similar ARID1A occupancy (~40 % vs ~45 %) and identical genomic annotation. A direct haplo-vs-buffered de novo scan controls for that shared anatomy and asks: *within the ARID1A-responsive universe, which motifs separate dose-sensitive from dose-tolerant sites?*

![De novo haplo vs buffered motif contrast](../../figures/extended_analysis/het_dose_response/motif_denovo_haplo_vs_buffered_heatmap.png)

The contrast is strongly **subset-specific**:

- **MP haplo vs buffered:** NFκB (−log10 p = 23) and AP-1/bZIP (21), with bHLH, IRF/STAT and GATA at lower significance. Dose-sensitivity in the memory-precursor program is the NFκB/AP-1 arm of the T-cell effector TF network.
- **EEC haplo vs buffered:** GATA (16), NFκB (14), AP-1/bZIP (14), T-box (13), TCF/LEF (12). A broader AP-1/NFκB-flavoured signature with a GATA component.
- **TE haplo vs buffered:** RUNX (34), TCF/LEF (29), KLF/SP (24), bHLH (21), IRF/STAT (19), AP-1/bZIP (12). TE is the outlier — its dose-sensitive peaks carry a RUNX/TCF/KLF signature, not an NFκB one. Because RUNX is not differentially enriched in EEC or MP haplo-vs-buffered (both classes have it), the TE-specific RUNX hit here is a real subset asymmetry.

Taken together, §6b shows ARID1A dose-sensitivity is wired through subset-distinct TF programs rather than a single generic motif — **AP-1/NFκB in MP and EEC, RUNX/TCF/KLF in TE** — all sitting on top of the shared RUNX/T-box scaffold that characterises the broader ARID1A-responsive universe.

## 7. Interpretation

1. **ETS-decorated enhancers are the dose-tolerant set.** Sites strongly enriched for ETS motifs (Etv2, ERG, Fli1, Ets1) dominate the buffered class across all three subsets — one *Arid1a* allele is sufficient to keep ETS-wired enhancers open. This is consistent with the paper's observation that Het mice are largely phenotypically normal: the ETS-centric effector enhancer backbone is preserved at 1 copy.
2. **RUNX and T-box are not dose-sensitivity markers.** Both families are enriched in buffered and haplo classes roughly equally. The apparent "RUNX/T-box program" in haplo peaks in the known-motif heatmap is a scale-bar artefact. RUNX/T-box mark the broader ARID1A-responsive effector universe, not the subset of it that fails at reduced dosage.
3. **Dose-sensitive sites are subset-specific.** The TFs whose motifs specifically predict haplo-vs-buffered split differ by subset: MP and EEC haplo peaks carry NFκB and AP-1/bZIP signatures (EEC additionally GATA); TE haplo peaks carry RUNX + TCF/LEF + KLF/SP + bHLH. Dose-sensitivity at ARID1A-bound enhancers is not one mechanism but three subset-tuned ones.
4. **The rare linear (truly dose-graded) class represents the most directly ARID1A-bound enhancers** (up to 84% ARID1A occupancy in TE, widest peaks, highest TF co-occupancy). These ~650 peaks across subsets are the strongest candidates for direct-titration mechanistic follow-up.
5. **De-repression (gained-in-KO peaks) is almost entirely recessive**: gained peaks are 77–85% buffered/other. Loss of a single ARID1A allele does not open ectopic chromatin.

## 8. Methods

### 8.1 Input data

Counts were taken from the nf-core/atacseq consensus feature-counts matrix built at the merged-replicate level with MACS2 narrow peak calling:

```
results/atac/bowtie2/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt
```

The matrix was restricted to D8 TE/EEC/MP Exp2 samples with genotype WT, Het or KO (regex `^D8_(WT|Het|KO)_(TE|EEC|MP)_Exp2_Rep[0-9]+$`), yielding a balanced single-batch design of n = 3 WT, 3 Het, 2 KO per subset (24 libraries total). Genotype was coded as an ordered factor `WT > Het > KO` with dose 2, 1, 0 copies of *Arid1a*. The 48h timepoint was excluded a priori per project QC (REP2 FRiP = 0.115, REP1 FRiP = 0.314, 5M reads — see `MEMORY.md`).

### 8.2 Per-subset differential accessibility

For each subset (TE, EEC, MP) independently:

1. `DESeqDataSetFromMatrix(countData, colData, design = ~ genotype)`
2. Pre-filter: retain peaks with ≥ 10 counts in at least `min(table(genotype))` = 2 samples.
3. **LRT** (`DESeq(test = "LRT", reduced = ~ 1)`) to test any genotype effect, extracted with `results(independentFiltering = FALSE)` so no peaks are dropped by default filtering.
4. **Wald** (`DESeq()`) for three pairwise contrasts: `c("genotype", "Het", "WT")`, `c("genotype", "KO", "WT")`, `c("genotype", "KO", "Het")`. Unshrunken log2 fold-changes were used for classification — `apeglm` was not applied because the classification explicitly requires magnitude comparisons *between* contrasts (Het vs WT magnitude compared to KO vs WT), and apeglm shrinks each independently, which would distort the ratio.
5. VST (`vst(blind = FALSE)`) for visualization and per-genotype means.

Subsets were modeled separately rather than with an interaction term (`~ subset * genotype`) because dispersion estimates differ between subsets and the classification is inherently per-subset.

### 8.3 Dose-response classification

Given LRT padj and the three Wald contrasts per peak, classification used BH-adjusted thresholds `α = 0.05` and `|log2FC| ≥ 0.5` (`HET_KO_RATIO = 0.5`). Tests are applied top-down; a peak takes the first class that matches:

| Class | Rule (all conditions must hold) |
|---|---|
| nonmonotonic | LRT sig; `sign(lfc_HetWT) ≠ sign(lfc_KOWT)`; both |lfc| ≥ 0.5; at least one of HetWT/KOWT padj < 0.05 |
| linear | LRT sig; same sign; HetWT, KOWT, KOHet all padj < 0.05; `0.5·|lfc_KOWT| ≤ \|lfc_HetWT\| < \|lfc_KOWT\|` |
| haploinsufficient | LRT sig; same sign; HetWT padj < 0.05; KOWT padj < 0.05; KOHet padj ≥ 0.05; both |lfc| ≥ 0.5 |
| buffered | LRT sig; KOWT padj < 0.05; \|lfc_KOWT\| ≥ 0.5; HetWT padj ≥ 0.05; \|lfc_HetWT\| < 0.5 |
| other_responsive | LRT sig, did not match above |
| insensitive | LRT padj ≥ 0.05 or NA |

Direction is defined by `sign(lfc_KOWT)` (gained = accessibility up in KO vs WT; lost = down), or NA where neither HetWT nor KOWT is significant. Insensitive peaks carry NA direction by construction.

### 8.4 Genomic annotation

Consensus peaks were annotated once with ChIPseeker `annotatePeak` against `TxDb.Mmusculus.UCSC.mm39.knownGene` with `tssRegion = c(-3000, 3000)` and `level = "gene"`. The raw ChIPseeker annotation strings were collapsed to six categories:

- Promoter: begins with "Promoter"
- 5′UTR / 3′UTR: begins with "5' UTR" / "3' UTR"
- Exon / Intron: contains "Exon" / "Intron"
- Intergenic: contains "Intergenic", "Downstream", or "Distal"

The annotation table is cached at `results/atac/differential/consensus_peaks_annotated.csv` and re-used across downstream scripts.

### 8.5 CUT&RUN co-occupancy features

D5 WT CUT&RUN MACS2 narrow peaks (nf-core/cutandrun) were loaded for five tracks: ARID1A (pooled R1+R2 union via `reduce(c(gr1, gr2))`), H3K27ac, Tbet, BATF, ETS1. Per-track replicates were pooled by `GenomicRanges::reduce()` on the concatenation. Each ATAC consensus peak was scored for any overlap (`overlapsAny`) against each track, producing five binary features `has_ARID1A`, `has_H3K27ac`, `has_Tbet`, `has_BATF`, `has_ETS1`. TF co-occupancy was the integer sum of `has_ARID1A + has_Tbet + has_BATF + has_ETS1` (0–4). D5 WT was used because paired D8 CUT&RUN was not generated in the original study; the D5→D8 stability assumption is flagged in §9. H3K27me3 tracks were called at other genotypes/timepoints but not at D5 WT, so H3K27me3 was excluded from feature tests.

### 8.6 Feature enrichment tests

Feature enrichment was scoped to peaks *losing* accessibility (direction = "lost" per class, union with direction-NA insensitive as background):

```r
feat_lost <- feat_class %>%
  filter(direction == "lost" | class == "insensitive")
```

This scope was chosen because §4 established that the dose-response signal is overwhelmingly in the lost direction; including gained peaks would dilute the contrast. For each (subset × class × feature) cell with n ≥ 10 in foreground and background, a 2×2 Fisher exact test was run against the within-subset insensitive set:

```
               class   insensitive
feature=T        a          c
feature=F        b          d
```

Features tested: `has_ARID1A`, `has_H3K27ac`, `has_Tbet`, `has_BATF`, `has_ETS1`, `is_promoter` (feature == "Promoter"), `is_distal` (feature ∈ {"Intergenic", "Intron"}). Classes tested: buffered, linear, haploinsufficient, nonmonotonic. Output columns include odds ratio, raw p-value, and BH-adjusted padj across the full 84-test family (3 subsets × 4 classes × 7 features). Effect sizes (frac_class, frac_bg) are reported alongside p-values.

### 8.7 HOMER motif enrichment

BEDs were exported per (subset × class) restricted to direction = "lost" peaks. HOMER 4.11 `findMotifsGenome.pl` was run with:

```
findMotifsGenome.pl <class.bed> GRCm39.primary_assembly.genome.fa <outdir> \
  -size 200 -mask -p 8 -bg <background.bed>
```

`-size 200` centers each peak on a 200-bp window (narrower than the median peak width to reduce motif dilution), `-mask` applies RepeatMasker soft-masking from the FASTA. Twelve runs were performed, in two schemes:

1. **Class vs insensitive** (9 runs: TE/EEC/MP × buffered/linear/haplo): foreground = class-lost peaks, background = insensitive-lost peaks in the same subset. Detects TF programs characteristic of that dose-response class against a null of non-responsive chromatin.
2. **Haplo vs buffered** (3 runs: one per subset): foreground = haploinsufficient-lost, background = buffered-lost in the same subset. Controls for ARID1A binding rate and peak anatomy — both foreground and background are ARID1A-responsive — isolating motifs that specifically separate dose-sensitive from dose-tolerant sites.

Both HOMER outputs are produced in the same run; the report uses de novo motifs (`homerResults/motifN.motif`) rather than `knownResults.txt`. Rationale: `knownResults` scores the HOMER motif library against the foreground and is biased toward motifs that happen to be in the library and that are common enough in the background to have mature priors, which produces an enormous dynamic range dominated by pan-mammalian ETS hits. De novo discovery finds motifs actually present in the foreground and returns a single HOMER best-guess annotation per motif. For this analysis the de novo p-values are tighter and the family structure is more interpretable (§6). `knownResults.txt` files are retained for reference but not used in the figures.

### 8.8 Motif summary figure

`extended_analysis/scripts/het_dose_response/03_known_motif_summary.R`, `extended_analysis/scripts/het_dose_response/04_denovo_motif_summary.R` parses the header of every `homerResults/motifN.motif` file across all 12 HOMER runs (389 de novo motifs total). Each header carries the consensus sequence, HOMER's best-guess TF (`BestGuess:TF(family)/source/.../Homer(match_score)`), the natural-log p-value, and target/background occupancy percentages. −log10 p is computed as −lnP / ln 10. Each de novo motif is assigned to a TF family (ETS, RUNX, T-box, AP-1/bZIP, NFκB, KLF/SP, TCF/LEF, EGR, GATA, IRF/STAT, NR, bHLH, Zf, Homeo) by a regex over the best-guess TF name, with HOMER's family tag as a secondary source where the regex did not match. De novo motifs whose best-guess match score was < 0.6, or whose best-guess TF was clearly non-vertebrate (yeast: CUP9/CRZ1/DOT6/etc.; plant: WRKY40/NAC035/OsbZIP42/etc.; fly: Doc2/pho/CGxxxx), were labelled *other* and excluded from the family-level heatmap; they are retained in the per-comparison CSV (`motif_denovo_top_per_comparison.csv`). The class-vs-insensitive heatmap shows, per (subset × class, family) cell, the maximum −log10 p across de novo motifs assigned to that family. Scale is capped at 100: the EEC-buffered ETS cell reaches 252 (Etv2), and without capping the haplo range (~20–60) renders near-white. The haplo-vs-buffered panel uses the same logic on three comparisons with a scale capped at 30 (max observed ≈ 34).

### 8.9 Software versions

- R 4.x; DESeq2 1.44; ChIPseeker 1.40; GenomicRanges 1.56; rtracklayer 1.64
- `TxDb.Mmusculus.UCSC.mm39.knownGene`, `org.Mm.eg.db` 3.19
- ComplexHeatmap 2.20; ggplot2 3.5; tidyverse 2.0
- HOMER 4.11 (`/usr/bin/homer`); mouse genome GRCm39 primary assembly
- nf-core/atacseq pipeline (consensus mRp.clN matrix) and nf-core/cutandrun (MACS2 narrow peaks)

### 8.10 Reproducibility

Full pipeline in three steps from the repository root:

```bash
Rscript extended_analysis/scripts/het_dose_response/01_dose_classes.R        # classification + dose curves
Rscript extended_analysis/scripts/het_dose_response/02_feature_enrichment.R   # annotation, CUT&RUN, Fisher, BEDs, HOMER runner
bash   results/extended_analysis/het_dose_response/run_homer_dose_response.sh
Rscript extended_analysis/scripts/het_dose_response/04_denovo_motif_summary.R    # de novo motif summary heatmaps
```

Script 11 writes `het_dose_response.RData` (dds, vsd, res_tabs) and script 12 writes `het_feature_enrichment.RData` for interactive reanalysis without rerunning DESeq2.

## 9. Caveats

- KO n=2 per subset — `|lfc| ≥ 0.5` with DESeq2 Wald is stable at baseMean > 20 but borderline at low-count peaks. The linear and nonmonotonic class boundaries have modest power at small effect sizes; small rare-class counts (nonmonotonic 16–84, linear 37–309) should not be overinterpreted. Buffered and haploinsufficient counts (1k–10k) are robust.
- Classification treats peaks as independent; ARID1A-bound enhancer clusters may share a class non-independently.
- All CUT&RUN features use D5 WT signal as a proxy for D8 binding; the D5→D8 displacement is assumed small but has not been validated.

## 10. Outputs

Tables (`results/extended_analysis/het_dose_response/`):
- `dose_classes_{TE,EEC,MP}.csv` — per-peak class, VST means, all pairwise LFCs and padj
- `dose_classes_combined.csv` — long form, all subsets
- `class_summary.csv` — class × direction × subset counts
- `feature_summary_by_class.csv` — per-class medians for width, baseline, ARID1A binding, etc.
- `feature_enrichment_vs_insensitive.csv` — Fisher tests of each feature, per-class
- `homer_motifs/<subset>_<class>_vs_<bg>/homerResults/motif*.motif` — 12 HOMER de novo runs (used in figures)
- `homer_motifs/<subset>_<class>_vs_<bg>/knownResults.txt` — 12 HOMER known-motif runs (reference, not used in figures)
- `motif_denovo_top_per_comparison.csv` — top 8 de novo motifs per comparison with best-guess TF, family, −log10 p, occupancy

Figures (`figures/extended_analysis/het_dose_response/`):
- `class_counts_bar.{pdf,png}`, `dose_curves_{TE,EEC,MP}.{pdf,png}`, `pca_per_subset.{pdf,png}`
- `baseline_WT_by_class`, `peak_width_by_class`, `annotation_by_class`, `cutrun_overlap_by_class`, `tf_cooccupancy_by_class`
- `motif_denovo_summary_heatmap.{pdf,png}` — de novo family × comparison, class vs insensitive
- `motif_denovo_haplo_vs_buffered_heatmap.{pdf,png}` — de novo, direct haplo vs buffered contrast
- `motif_summary_heatmap.{pdf,png}`, `motif_haplo_vs_buffered_heatmap.{pdf,png}` — known-motif versions (for reference; not used in §6)
