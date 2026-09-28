# Computational Methods

Reanalysis of McDonald, Chick *et al.* (2023) "Canonical BAF complex activity shapes the enhancer
landscape that licenses CD8⁺ T cell effector and memory fates", *Immunity* 56:1303–1319
(GEO SuperSeries GSE228381), with a cross-study comparison against Guo *et al.* (2022),
*Nature* 607:135–141 (GEO SuperSeries GSE183619).

All paths below are relative to the repository root. Manuscript code is in `scripts/` (`core/`, `paper/`, `upstream/`); analyses beyond the paper are in `extended_analysis/scripts/<analysis>/`.

---

## 1. Overview

### 1.1 Data

All raw sequencing data deposited under GSE228381 were reprocessed from FASTQ:

| Sub-series | Assay | Samples | Layout | Design |
|---|---|---:|---|---|
| GSE227634 | RNA-seq | 30 | SE | WT/Het/KO P14 CD8⁺ T cells; D3 (WT, KO), D5 (WT), D8 TE/EEC/MP (WT, Het, KO) |
| GSE228171 | ATAC-seq | 54 | PE 75 bp | WT/Het/KO/Tbet-KO; Naive, 48 h, D3, D5, D8 TE/EEC/MP (two D8 experiments, Exp1 and Exp2) |
| GSE228193 | ATAC-seq | 8 | PE 100 bp | WT, 48 h: untreated, IL-12, IL-12 + ACBI1, IL-12 + BRM014 (n = 2 each) |
| GSE228380 | CUT&RUN | 28 | PE | ARID1A, H3K27ac, H3K27me3, T-bet, BATF, ETS1, IgG; Naive, 48 h, D5, D8 |
| GSE228546 | ChIP-seq | 8 | SE | T-bet and input; untreated, IL-12, IL-12 + ACBI1, IL-12 + BRM014 (n = 1 each) |

Sample metadata were parsed from the five GEO series-matrix files and joined to SRA run
accessions to produce a master sample sheet (`data/metadata/master_sample_sheet.tsv`, 128 samples;
`scripts/upstream/build_sample_sheet.py`), from which all nf-core samplesheets in `nextflow/samplesheets/`
were generated. In the ATAC-seq D8 data, Het libraries exist only in Exp2 and Tbet-KO libraries
only in Exp1; this batch structure is handled explicitly in every D8 analysis (§5.2).

### 1.2 Reference genome and annotation

- Genome: GRCm39 (mm39) primary assembly, GENCODE release M35
  (`GRCm39.primary_assembly.genome.fa`); UCSC-style `chr` sequence names.
- Gene annotation: GENCODE vM35 primary-assembly GTF.
- Blacklist: the ENCODE mm10 blacklist v2 (Amemiya *et al.* 2019) was lifted to mm39 with UCSC
  liftOver (`mm10ToMm39.over.chain.gz`), giving `data/reference/mm39-blacklist.v2.bed`. No mm39
  blacklist was supplied to the nf-core pipelines, so blacklist filtering was not applied at the
  alignment/peak-calling stage; the lifted blacklist was applied to deepTools signal matrices.
- Peak annotation used either a TxDb built from the GENCODE vM35 GTF
  (`GenomicFeatures::makeTxDbFromGFF`; `core/02_atacseq_analysis.R`, `integration/`, `guo2022/`) or
  `TxDb.Mmusculus.UCSC.mm39.knownGene` (`het_dose_response/`), with `org.Mm.eg.db` for identifier
  mapping. Ensembl version suffixes were stripped before mapping.

### 1.3 Execution environment

Primary processing used nf-core pipelines under Nextflow 25.10.0 with Docker containers
(`nextflow/nextflow.config`; local executor). Pipeline-internal tool versions are those of the
pinned pipeline release containers and are recorded in each pipeline's `pipeline_info/` output.
Downstream analyses were run in R 4.4.2 and Python 3 on the host.

---

## 2. Data acquisition

FASTQ files were downloaded from the European Nucleotide Archive (ENA) FTP mirror rather than
via the SRA toolkit, using the ENA file reports for the five BioProjects (PRJNA945945,
PRJNA948430, PRJNA983073, PRJNA949603, PRJNA950325). URLs and ENA-supplied MD5 checksums are
stored in `data/metadata/ena_fastq_urls.tsv`; `scripts/upstream/download_fastqs.sh` downloads in parallel
and verifies every file against its MD5. All 218 FASTQ files (128 samples) passed verification.

Guo *et al.* 2022 data (§8) were obtained the same way (`scripts/upstream/download_meta_fastqs.sh`), with
each file MD5-verified, retried up to three times, and deleted if verification ultimately failed.
Before pipeline launch, all Guo FASTQs were additionally re-verified by MD5 and tested with
`gzip -t`; two files that failed were re-downloaded and
re-verified.

---

## 3. Primary processing

Launch scripts with the exact command lines are `scripts/upstream/run_{atacseq,rnaseq,cutandrun,chipseq}.sh`.
Only parameters that differ from pipeline defaults are listed.

### 3.1 ATAC-seq — nf-core/atacseq 2.1.2

All 62 ATAC-seq libraries (GSE228171 and GSE228193) were processed in a single run so that they
share one consensus peak set. Non-default parameters: `--aligner bowtie2`, `--narrow_peak`,
`--mito_name chrM`, `--macs_gsize 2407883318`, `--read_length 75`.

The pipeline performed adapter/quality trimming (Trim Galore 0.6.7 / cutadapt 3.4), alignment with
Bowtie2 2.4.4, duplicate marking (Picard 3.0.0) and filtering of the merged-library BAMs to remove
duplicates, mitochondrial reads, unmapped/multi-mapping reads (MAPQ < 1), improperly paired or
orphan reads, fragments > 2 kb, reads with > 4 mismatches and soft-clipped reads (pipeline
defaults). Peaks were called with MACS2 2.2.7.1 (`--nomodel --keep-dup all`, narrow mode,
default q-value threshold) on each library and on replicate-merged BAMs. The merged-replicate
(`mRp`) consensus peak set (union of peaks across condition groups) was quantified on every
filtered library BAM with featureCounts (Subread 2.0.1; `-F SAF -O --fracOverlap 0.2`). This
matrix, `bowtie2/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt`,
is the input to all ATAC-seq count-based analyses. Library QC (FRiP, ataqv 1.3.1, fragment-size
distributions) was taken from the pipeline MultiQC report.

### 3.2 RNA-seq — nf-core/rnaseq 3.22.2

Thirty single-end libraries were processed with `--aligner star_salmon` (STAR 2.7.11b alignment,
Salmon 1.10.3 quantification of the STAR transcriptome BAMs), declared as reverse-stranded. Trim
Galore 0.6.10 was used for trimming. The transcriptome FASTA was generated by the pipeline from
the genome FASTA and GTF (the GENCODE transcript FASTA was not supplied because its transcript
identifiers do not match the GTF-derived transcriptome). Gene-level counts
(`star_salmon/salmon.merged.gene_counts.tsv`) and TPMs were used downstream. A second run with
`--pseudo_aligner kallisto` (kallisto 0.51.1; `scripts/upstream/run_rnaseq_kallisto.sh`) was performed as
an alternative quantification and is not used in any reported analysis.

### 3.3 CUT&RUN — nf-core/cutandrun 3.2.2

Twenty-eight paired-end libraries were processed with `--normalisation_mode CPM`,
`--peakcaller macs2,seacr` and `--macs_gsize 2494787188`. The original experiment did not include
an *E. coli* spike-in, so spike-in normalisation was not used; bigWig tracks are CPM-normalised.

Reads were trimmed (Trim Galore 0.6.6), aligned with Bowtie2 2.4.4
(`--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 --minins 10 --maxins 700
--dovetail`) and filtered to properly paired reads with MAPQ ≥ 20. Duplicates were marked but
retained in target libraries and removed from IgG libraries (pipeline defaults). Each target
library was paired with the IgG library from the same timepoint/genotype (samplesheet `control`
column; IgG was available for D5 WT and D5 KO). Peaks were called against IgG with:

- **MACS2** 2.2.7.1: `--format BAMPE --nomodel --shift -75 --extsize 150 --keep-dup all -q 0.01`;
- **SEACR** 1.3: stringent mode against the IgG bedGraph with `non` normalisation (pipeline default).

MACS2 narrowPeak calls were used for downstream overlap analyses. Libraries without a matched IgG
(ARID1A at Naive, 48 h and D8; H3K27ac and H3K27me3 at D8; ARID1A in D5 Tbet-KO) were aligned and
converted to signal tracks but not peak-called; they are analysed as signal only.

### 3.4 ChIP-seq — nf-core/chipseq 2.1.0

Eight single-end libraries (four T-bet ChIPs, each with its own input) were processed with
`--aligner bowtie2` (reusing the Bowtie2 index built by the CUT&RUN run), `--narrow_peak`,
`--mito_name chrM` and `--read_length 100` (effective genome size estimated by the pipeline with
khmer). Alignment used Bowtie2 2.5.2 and peaks were called against the matched input with MACS3
3.0.1 at default thresholds. Each condition has a single replicate.

---

## 4. Quality control and sample exclusions

### 4.1 Library QC

FastQC, alignment, duplication, FRiP, fragment-size and library-complexity metrics were reviewed
from the MultiQC reports of each pipeline. Sample-level structure was examined by PCA and
Pearson correlation of variance-stabilised counts (DESeq2 `vst`).

### 4.2 Exclusions

- **All 48 h ATAC-seq libraries from GSE228171 were excluded** from every reported analysis.
  REP2 failed QC (FRiP = 0.115) and REP1 was borderline (FRiP = 0.314 with ~5 M reads), leaving no
  usable replicate pair. The libraries were retained only in joint peak calling and in the
  consensus count matrix; ATAC timecourse analyses therefore use Naive, D3, D5 and D8. The
  independent replicate-cohesion check (§4.3) flagged the same two libraries as outliers without
  prior knowledge of the exclusion. The 48 h inhibitor libraries (GSE228193) are a separate
  experiment and were retained.
- No RNA-seq, CUT&RUN or ChIP-seq libraries were excluded.
- Guo *et al.* GSE183616 (Myc-high/Myc-low/naive ATAC) was treated as low-confidence and excluded
  from interpretation (§8.3).

### 4.3 Automated design-integrity audit

A set of reusable checks (`scripts/qc_checks.R`) was applied to every differential contrast in
the project (`extended_analysis/scripts/qc_audit/qc_audit.R`; output `results/extended_analysis/qc_audit/qc_audit_report.{csv,txt}`). Each check
returns PASS / WARN / FAIL:

| Check | Criterion |
|---|---|
| design integrity | Mean within-group vs between-group Pearson correlation of VST values, and each sample's nearest neighbour. FAIL if the within–between gap ≤ 0 or > 25% of samples have a nearest neighbour in another group; WARN if any misfit or gap < 0.005. |
| replicate cohesion | WARN for samples whose mean within-group correlation is < −2 SD of the group. |
| positive control | Expected direction and significance of a known gene (e.g. *Arid1a* down in KO). |
| p-value distribution | WARN if the upper-half density departs strongly from uniform. |
| effect symmetry | WARN if the median log2 fold change is far from zero. |
| dispersion | WARN if dispersion exceeds a reference contrast by > 1.8-fold. |
| library sanity | WARN if library-size CV > 0.75 or size-factor range > 5-fold. |

McDonald contrasts were grouped by the full biological group (timepoint × subset × genotype) plus
per-stratum genotype contrasts. No McDonald contrast failed any check. WARN-level findings were one
outlier replicate in D8 TE (`D8_KO_TE_Exp2_Rep2`) and a small number of nearest-neighbour misfits
in the full-group checks against healthy gaps (ATAC +0.059, RNA +0.066); the RNA misfits were D8
Het samples, which are expected to be intermediate between WT and KO. All FAILs were in the Guo
GSE183618/GSE183616 data (§8). The audit also re-created the pooled-WT GSE183618 design as a
positive control for the design-integrity check, which it correctly failed.

---

## 5. Differential expression and accessibility

Unless stated otherwise, all contrasts used DESeq2, apeglm log2 fold-change shrinkage, and the
thresholds **|log2FC| > 1 (2-fold) and Benjamini–Hochberg FDR < 0.05**. Wald statistics for
ranking were taken from the unshrunken `results()` object (apeglm does not return `stat`).
Each contrast was fitted on only the samples involved (per-stratum subsetting) rather than in one
global model.

### 5.1 RNA-seq (`scripts/core/01_rnaseq_analysis.R`)

Salmon gene counts were rounded to integers and imported with `DESeqDataSetFromMatrix`. Genes with
≥ 10 counts in ≥ 3 samples were retained, both globally and within each contrast. Design `~ genotype`
within each stratum. Contrasts:

- D3: KO vs WT (bulk).
- D8, per subset (TE, EEC, MP): KO vs WT and Het vs WT.
- D8 pseudobulk: counts summed over TE + EEC + MP for each genotype × replicate with all three
  subsets present; KO vs WT and Het vs WT.

No D5 genotype contrast was possible (D5 contains WT only). Gene set enrichment used fgsea
(`nPermSimple = 10000`, `eps = 0`) on genes ranked by the Wald statistic (one entry per gene
symbol, maximum |stat|), against MSigDB Hallmark and C7 ImmuneSigDB collections (msigdbr, mouse).
GO Biological Process over-representation of up- and down-regulated genes used clusterProfiler
`enrichGO` with all tested genes as universe (BH-adjusted p < 0.05).

### 5.2 ATAC-seq (`scripts/core/02_atacseq_analysis.R`)

The consensus featureCounts matrix (§3.1) was filtered to peaks with ≥ 10 counts in ≥ 3 samples
globally and ≥ 10 counts in ≥ 2 samples within each contrast. When a contrast spanned both D8
experiments, `experiment` was included as a covariate (`~ experiment + genotype`); otherwise the
design was `~ genotype` (or `~ treatment`). If apeglm failed to converge, `normal` shrinkage was
used. Contrasts:

- D3 and D5: KO vs WT.
- D8, per subset (TE, EEC, MP): KO vs WT, Het vs WT, Tbet-KO vs WT (batch-adjusted).
- D8 pseudobulk: libraries summed over TE + EEC + MP per genotype × experiment × replicate;
  KO vs WT, Het vs WT, Tbet-KO vs WT (batch-adjusted).
- BAF inhibitors (GSE228193): IL-12 vs untreated; IL-12 + ACBI1 vs IL-12; IL-12 + BRM014 vs IL-12.

Differential peaks were termed *gained* or *lost* by the sign of the shrunken LFC. Consensus peaks
were annotated with ChIPseeker (`annotatePeak`, `tssRegion = c(-3000, 3000)`) against the GENCODE
TxDb, and GO enrichment of genes nearest to gained/lost peaks used `enrichGO` as above.

Because ARID1A loss removes tens of thousands of peaks, median-of-ratios normalisation shifts
unchanged peaks towards positive LFC (in the D8 pseudobulk KO vs WT, the peak-level median LFC is
+0.07 while the mean is −0.21). Effects are therefore summarised as counts of gained and lost
peaks, not as mean LFC.

### 5.3 CUT&RUN and ChIP-seq

Most CUT&RUN antibodies have one replicate per genotype, so count-based differential binding was
not attempted. D5 WT vs KO binding was compared by peak overlap (≥ 1 bp; replicate peak sets
pooled with `GenomicRanges::reduce`) for ARID1A, BATF, ETS1, T-bet and H3K27ac
(`scripts/core/04_cutandrun_analysis.R`). T-bet ChIP-seq peak sets were compared across conditions by
overlap, annotated with ChIPseeker, and scanned with HOMER `findMotifsGenome.pl -size 200 -mask`
against the default GC-matched genomic background (`scripts/core/05_chipseq_analysis.R`).
Signal heatmaps at ATAC peak sets used deepTools 3.5.4 `computeMatrix reference-point`
(peak centre, ± 2 kb, 50-bp bins, lifted mm39 blacklist; `scripts/core/06_deeptools_heatmaps.sh`).

---

## 6. Reproduction of the published panels

Published panels were regenerated from the reprocessed mm39 data rather than from the original
processed files. Where the published analysis used HOMER-specific tools, the equivalent
reprocessed inputs were used (see §10 for the resulting differences). The GEO-deposited bigWigs
(`data/geo_processed/`) were used only as a validation reference.

### Paper panel scripts (scripts/paper/)

Each panel is produced by one script in `scripts/paper/`. Figures are written to `figures/paper/`
and the numbers behind each panel to `results/paper/`. `PANEL_MAP.md` gives the panel-by-panel
correspondence, the status of each panel, and every known difference from the publication.
Thresholds follow the figure legends unless stated: 2-fold change and Benjamini–Hochberg FDR < 0.05
for ATAC-seq and RNA-seq.

*OCR clusters (Fig. 1, 2, 5J).* The published cluster definitions (original mm10 HOMER peaks:
Conserved, Naive, Early Activation, Activation and Late Activation) were lifted to mm39 with UCSC
liftOver (`mm10ToMm39.over.chain.gz`, `-minMatch=0.95`); 3 of 44,221 regions failed to map
(`data/metadata/paper_ocr_clusters/`). Region-based panels (signal, annotation, motifs, overlap)
use these regions directly. For DA-based panels, each region was mapped to the overlapping nf-core
consensus peak, which succeeded for 99.1–99.9% of regions per cluster. As a concordance check, WT
consensus peaks were clustered independently (`figS1_denovo_clusters.R`): peaks DA in any pairwise
WT timepoint contrast were grouped by k-means on row-z-scored VST means (Naive, D3, D5, D8), and
the result was cross-tabulated against the published clusters.

*Signal heatmaps and profiles (Fig. 1A/B, 2D–F, 4C, 5C–E/H/J).* deepTools 3.5.4
`computeMatrix reference-point` (±1 kb around region centres, 10–20-bp bins) was run on nf-core
merged-replicate bigWigs (ATAC: 1e6/mapped-fragment scaling) and CPM-normalised CUT&RUN and ChIP
bigWigs. Conditions spanning several libraries (for example D8 subsets) were averaged with
`bigwigAverage`. D8 genotype comparisons use the Exp2 batch,
the only one containing WT, Het and KO. The 48 h ATAC track (merged replicates) appears in Fig.
1A/1B only and was excluded from all statistical tests (§4.2).

*Genomic annotation (Fig. 1C, 5C, 5F).* HOMER 5.1 `annotatePeaks.pl` with the GENCODE vM35 GTF,
collapsed to promoter, intergenic, intron, exon and other.

*Observed/expected overlap (Fig. 1D).* HOMER `mergePeaks -matrix`, as in the original analysis.
ARID1A peaks were called with MACS2 without a control, because a matched IgG exists only for D5.

*Motif enrichment (Fig. 1F, S1C, 5A).* HOMER 5.1 `findMotifsGenome.pl` on the mm39 FASTA,
`-size 200`, known motifs, GC-matched random background. Motif families were summarised as the
minimum p-value per family, capped at the published display limit.

*RNA-seq panels (Fig. 2G–I, 3E–I).* These use the DESeq2 results of §5.1. DEG counts are reported
at the 2-fold threshold of the legends and at the |log2FC| ≥ 0.585 threshold of the original
Methods. Hallmark GSEA used fgsea (10,000 permutations) on the Wald statistic. For Fig. 3I, the
MP-vs-TE and activated-vs-naive signatures were rebuilt from GSE10239 (MSigDB C7). The legend's
"GSE10739" is an unrelated series and appears to be a typo. The day-4.5 effector contrast was fixed
in advance, not chosen post hoc, and all GSE10239 contrasts are reported.

*Differential panels (Fig. 4A–E, 5G, S5A).* PCA on the 5,000 most variable VST peaks. Lost and gained
OCRs, UpSet intersections (ComplexHeatmap) and Arid1a-KO vs Tbx21-KO overlaps use the §5.2 DA
tables; S5A uses padj < 0.01 as in its legend. For Fig. 4E, TE and MP signature genes were defined
as WT TE vs WT MP DEGs.

*ARID1A-dependent OCRs (Fig. 5C).* OCRs overlapping D5 WT ARID1A MACS2 peaks were split by whether
they were lost in D5 KO.

*BAF inhibitors (Fig. 5D–F).* ACBI1- and BRM014-dependent OCRs were defined among T-bet-bound OCRs
(overlapping IL-12 T-bet ChIP peaks) as lost vs DMSO + IL-12 (2-fold, padj < 0.05). Because the
inhibitor-treated libraries have lower FRiP, DESeq2 was run with total-counted-read size factors
(primary; closest to the original HOMER total-tag normalisation). A median-of-ratios analysis is
reported as a sensitivity check.

*TF PageRank (Fig. 6A, S6A).* This reimplements Taiji (Zhang *et al.* 2019) in R. Nodes are TFs
and genes. TF→gene edges come from JASPAR2020 CORE vertebrate motif matches (motifmatchr,
p < 5e-5) in ATAC consensus peaks. Peaks were linked to genes following Taiji's no-Hi-C rule
(GREAT basal-plus-extension): a promoter link (TSS −5 kb/+1 kb, weight 1), plus distal links to the
nearest upstream and downstream TSS within 50 kb, weighted exp(−d/10 kb). Edge weight was
√(TF expression) × Σ(peak accessibility × link weight) for each genotype; node weight was
exp(z-scored expression). Personalised PageRank (igraph, damping 0.85) was run on the reversed
network for D8 MP WT and KO, and TFs were ranked by the log2 WT/KO PageRank ratio. Replicate
stability was assessed by rerunning with each RNA replicate pair.

*Signal tracks (Fig. S1A, S2D, S4A/B).* bigWig signal in 400 bins over each locus (GENCODE vM35
gene models), drawn with ggplot2.

---

## 7. Extensions beyond the original publication

### 7.1 ARID1A dose-response taxonomy (`extended_analysis/scripts/het_dose_response/`)

The allelic series (WT = 2, Het = 1, KO = 0 copies of *Arid1a*) was used to classify each D8 OCR
by the shape of its dose response. To avoid confounding genotype with batch, only D8 Exp2
libraries were used (WT n = 3, Het n = 3, KO n = 2 per subset; 24 libraries). Each subset was
modelled separately (`~ genotype`; peaks with ≥ 10 counts in ≥ 2 samples):

1. A likelihood-ratio test (`reduced = ~ 1`, `independentFiltering = FALSE`) identified peaks with
   any genotype effect (FDR < 0.05).
2. Wald tests gave Het vs WT, KO vs WT and KO vs Het effects. Unshrunken LFCs were used, because
   the classification compares magnitudes between contrasts and independent shrinkage would
   distort their ratio.
3. Peaks were assigned, in order, to the first matching class (α = 0.05, |LFC| ≥ 0.5):
   - *nonmonotonic*: Het and KO effects of opposite sign, both |LFC| ≥ 0.5, at least one significant;
   - *linear*: Het, KO and KO-vs-Het all significant, same sign, 0.5·|LFC_KO| ≤ |LFC_Het| < |LFC_KO|;
   - *haploinsufficient*: Het and KO significant vs WT, same sign, both |LFC| ≥ 0.5, KO vs Het n.s.;
   - *buffered*: KO significant with |LFC| ≥ 0.5, Het n.s. with |LFC| < 0.5;
   - *other responsive*: LRT-significant, no class matched;
   - *insensitive*: LRT not significant.

   Direction (gained/lost) followed the sign of the KO vs WT effect.

Class features were compared against insensitive peaks among peaks losing accessibility
(`extended_analysis/scripts/het_dose_response/02_feature_enrichment.R`): genomic annotation (ChIPseeker, ± 3 kb promoter),
peak width, baseline WT accessibility, and overlap with D5 WT CUT&RUN MACS2 peaks (ARID1A
replicates pooled; H3K27ac, T-bet, BATF, ETS1). Fisher's exact tests were BH-adjusted across all
84 subset × class × feature tests. D5 CUT&RUN was used as a proxy for D8 binding because D8
factor CUT&RUN was not generated.

Motif enrichment used HOMER 5.1 `findMotifsGenome.pl -size 200 -mask` with explicit backgrounds,
in two schemes: (i) each class's lost peaks vs insensitive-lost peaks of the same subset (TE, EEC,
MP × buffered, linear, haploinsufficient), and (ii) haploinsufficient-lost vs buffered-lost peaks
per subset, which controls for ARID1A occupancy and peak anatomy. De novo motifs were summarised
(`extended_analysis/scripts/het_dose_response/04_denovo_motif_summary.R`) by assigning each to a TF family from its HOMER
best-match annotation; matches with score < 0.6 or to non-vertebrate factors were labelled
"other" and excluded from family summaries. Family-level enrichment is reported as the maximum
−log10 p across that family's de novo motifs.

### 7.2 chromVAR (`extended_analysis/scripts/chromvar/build_chromvar.R`)

chromVAR 1.28.0 was run on the D8 Exp2 WT/Het/KO × TE/EEC/MP libraries using the consensus count
matrix restricted to canonical chromosomes. Peaks were filtered with
`filterPeaks(non_overlapping = TRUE)`, GC bias was computed from the reference FASTA, and
JASPAR2020 CORE vertebrate PWMs were matched with motifmatchr. Deviation z-scores
(`computeDeviations`, `deviationScores`; seed 42) were averaged per motif and genotype. Motifs were
assigned to TF families by name. Dose-buffering was summarised as the fraction of the WT–KO
difference retained at one copy, `(Het − KO) / (WT − KO)`, computed only for motifs with a
monotone decline (WT > Het > KO); values > 0.5 indicate buffering and ≈ 0.5 an additive response.

### 7.3 TOBIAS footprinting (`extended_analysis/scripts/footprinting/run_tobias_footprinting.sh`)

D8 Exp2 filtered BAMs were merged per genotype across subsets (WT 9, Het 9, KO 6 libraries).
TOBIAS 0.16.1 was run as `ATACorrect` (Tn5 bias correction against GRCm39) → `ScoreBigwig
--score footprint` → `BINDetect` with JASPAR2020 CORE vertebrate motifs, over the nf-core
merged-library consensus peaks restricted to chr1–19, X and Y, with conditions WT, Het and KO.
Per-motif het retention was computed from BINDetect mean footprint scores with the same formula
and monotonicity requirement as chromVAR (§7.2), summarised by TF family
(`extended_analysis/scripts/figures/tobias_dose_response.R`), and compared with chromVAR retention by Spearman
correlation across motifs monotone in both methods.

### 7.4 Temporal clustering of WT accessibility (`scripts/core/03_atac_temporal_clustering.R`)

WT ATAC libraries from Naive, D3, D5 and D8 (48 h excluded) were filtered (≥ 10 counts in ≥ 2
samples), variance-stabilised (`varianceStabilizingTransformation`, `blind = TRUE`) and averaged
per timepoint, pooling the three D8 subsets. The 1,500 lowest-variance peaks were set aside as
constitutive; all remaining peaks with non-zero MAD were z-scored across timepoints. Three methods
were applied:

- **k-means**: k chosen over 3–10 by maximum mean silhouette width (within-cluster sum of squares
  and gap statistic, `cluster::clusGap` with B = 50, were also computed); final fit with
  `nstart = 50`, seed 42.
- **Mfuzz** fuzzy c-means with the same number of clusters and fuzzifier *m* estimated by
  `mestimate`; core and strict membership sets were exported.
- **DEGreport `degPatterns`** on the 5,000 highest-MAD peaks, using a replicate-balanced matrix
  (two columns per timepoint; D8 as per-experiment means), `minc = 50`, `reduce = TRUE`,
  `cutoff = 0.7`, `scale = TRUE`.

Clusters were named by the timepoint of maximum centre accessibility, and agreement between
methods was reported.

**WT vs KO comparison.** `temporal_clustering/01_degpatterns_wt_ko.R` ran `degPatterns` (`minc = 100`,
`reduce = TRUE`, `cutoff = 0.7`) separately on WT and KO trajectories (Naive → D3 → D5 → D8 TE;
the WT Naive libraries serve as the shared baseline) on the 5,000 peaks ranked highest by MAD of
timepoint × genotype means from the joint VST.

**Trajectory classes and ARID1A loss** (`extended_analysis/scripts/atac_trajectories/01`–`03`). Dynamic peaks were
called separately in the WT and KO timecourses with timecourse-patterns (github.com/bchick/
timecourse-patterns, commit 087f02b; run by `atac_trajectories/01_timecourse_patterns.sh`). Both arms share the
Naive WT libraries as day 0 (no Naive KO libraries exist), followed by D3, D5 and D8 TE (Exp1 +
Exp2, present in both genotypes). All consensus peaks on primary chromosomes were supplied as raw
counts; peaks with mean count < 10 were removed (68,751 retained). Counts were VST-transformed
(blind), and a peak was dynamic in an arm when a DESeq2 likelihood-ratio test of `~ time` against
`~ 1` gave padj < 0.01 and the range of its per-timepoint means was ≥ 0.5 on the VST scale. That
gave 49,929 dynamic peaks in WT and 34,509 in KO. Sites that differ between the timecourses were
defined as the set difference (`atac_trajectories/02_dynamic_site_sets.R`). 18,200 were dynamic in WT but not KO
("lost"), 2,780 in KO but not WT ("gained") and 31,729 in both. Lost sites were split by why
they failed in KO: KO range < 0.5 ("flat in KO", 10,306) or range ≥ 0.5 with LRT padj ≥ 0.01
("sub-threshold in KO", 7,894). In each arm, and again for the lost sites alone (WT timecourse,
no further selection), `degPatterns` (`reduce = FALSE`; `minc = 50`, scaled to the subsample,
e.g. 15 in the WT arm) clustered a stratified 6,000-peak subsample. The remaining peaks were assigned to the cluster centroid of highest
correlation (r ≥ 0.6, otherwise unassigned). Clusters were named by the workflow's shape rule,
with its second stage splitting "Increasing" into transient, sustained and late
(`superclusters.split.cross = -0.5` for the 0/3/5/8 grid). KO trajectories were compared with WT
classes by projection (`atac_trajectories/03_ko_trajectory_projection.R`). WT and KO libraries were
variance-stabilized together, per-timepoint means were z-scored within genotype, and each peak's
KO profile was assigned to the WT class centroid of highest correlation. KO range < 0.5 was
called "static in KO" and best r < 0.6 "unassigned". As a control, re-assigning WT profiles the
same way recovered the workflow's labels for 77.5% of all dynamic peaks and 70.0% of lost sites.

**Motifs of sites that lose their dynamics** (`atac_trajectories/04`–`05`). Known-motif
enrichment followed the timecourse-patterns T cell example. It used MEME-suite SEA 5.5.4 with a
fixed seed (42) and JASPAR2024 CORE vertebrates non-redundant motifs (879 motifs; sha256
`dd494278…a65b18`). Every region was trimmed to 200 bp around its centre, so differences in peak
width could not pass for enrichment. Two kinds of comparison were run. (i) Each WT trajectory
class of the lost sites was tested against static peaks, those dynamic in neither timecourse
(16,042). (ii) Within each WT class, lost sites were tested against sites that stay dynamic in KO,
in both directions. For this comparison both sets took their class labels from the same WT
clustering. Motifs with SEA q < 1e-5 are reported, keeping the most significant matrix per TF
name. Peaks were called promoters from the ChIPseeker annotation of the consensus peak set
(`consensus_peaks_annotated.csv`, from `core/02_atacseq_analysis.R`), because the motifs separating kept from lost
opening sites include GC-rich CpG-promoter motifs.

**Expression of the TFs behind lost-site motifs** (`atac_trajectories/06_tf_expression_vs_motifs.R`). This asks
whether lost accessibility could follow from lost TF expression. Each JASPAR motif name was split
on "::" and matched case-insensitively to RNA-seq gene symbols, so heterodimers contribute both
genes. RNA-seq KO vs WT contrasts exist at D3 and at D8 TE (sorted, as for ATAC). There are no D5
KO or Naive RNA libraries, so D3 was matched to the early classes (Transient, Transient
Increasing) and D8 TE to the late ones (Sustained and Late Increasing); both timepoints are
reported for every class. A TF was counted as expressed at a timepoint when its mean normalized
count in WT was ≥ 10. For each class and timepoint, the KO vs WT log2 fold changes of TFs whose
motifs were enriched (q < 1e-5) in lost or in kept sites were compared with those of all other
expressed JASPAR TFs (two-sided Wilcoxon). TFs marking lost sites were not down-regulated at the
matched timepoint. At D3 their median log2FC was +0.41 against +0.11 for other TFs (Transient
Increasing, p = 0.013). At D8 TE it was +0.002 against +0.010 (Late Increasing, p = 0.33).

**Joint motif models** (`atac_trajectories/07_glmnet_motif_models.R`). SEA and HOMER test one motif at a time,
so correlated motifs are called together. To ask which motifs predict a site set once all are
fitted jointly, the same 200 bp windows were scanned with motifmatchr (p < 5e-5) for the 879
JASPAR2024 motifs SEA used. That gave a binary window × motif matrix; motifs matched in fewer
than 20 windows of a model were dropped. Elastic-net models (glmnet, alpha = 0.5) were fitted
with three safeguards:
(i) GC fraction, CpG observed/expected and a promoter flag (ChIPseeker, from `core/02_atacseq_analysis.R`) were entered
unpenalized, and a covariates-only model served as the baseline.
(ii) Performance was measured out of sample. Each model was fitted on a stratified 80% split, with
lambda chosen by 10-fold `cv.glmnet` on that split only. It was then scored by AUC on the held-out
20% (one-vs-rest for the multinomial model), with 1,000 stratified bootstrap resamples for 95% CIs
on the AUCs and on ΔAUC.
(iii) Stability selection. Each model was refitted at the training lambda.1se on 100 stratified
half-subsamples, and a motif was called stable at selection frequency ≥ 0.8, keeping the most
frequently selected matrix per TF name.
Model A (multinomial) separated the five WT trajectory classes of lost sites (18,038 windows).
Model B (binomial) separated lost from kept sites within each WT class. HOMER known-motif enrichment
(`findMotifsGenome.pl -size given -nomotif`, its own library and GC-matched background) was run on
the same windows: each lost class against static peaks, and lost against kept in both directions.
In model A, motifs added +0.08 to +0.13 held-out AUC over the covariates in every class (all 95% CIs
excluded zero). Examples are Late Increasing 0.593 → 0.713 and Decreasing 0.696 → 0.823.
In model B, motifs added +0.111 [0.094, 0.129] for Decreasing and +0.038 [0.018, 0.058] for
Transient. They added ≤ 0.02 for Transient Increasing and nothing for Sustained and Late
Increasing, where lambda.1se kept no motif. A lambda.min sensitivity fit gave +0.013 and +0.014
for those two (the Late CI includes zero). Covariates-only AUCs for lost vs kept were 0.54–0.60.

`temporal_clustering/02_kmeans_mfuzz_wt_ko.R` clustered the 30,000 highest-MAD peaks with k-means
(k = 8, chosen from an elbow plot over k = 2–15; `nstart = 25`, `iter.max = 500`) and Mfuzz
(c = 8, `mestimate`) on WT profiles, then projected KO profiles onto the WT centres
(nearest centre for k-means; fuzzy membership for Mfuzz) to quantify disruption of WT patterns.

### 7.5 Multi-omic integration (`extended_analysis/scripts/integration/01_multiomic_integration.R`)

D8 pseudobulk KO vs WT ATAC and RNA results were joined by version-stripped Ensembl gene ID
(ATAC peaks assigned to their ChIPseeker nearest gene; gene-level ATAC LFC = mean over its peaks),
and concordance was assessed by Pearson correlation. Consensus OCRs overlapping D5 WT ARID1A
CUT&RUN peaks were classed as ARID1A-bound, and lost:gained ratios in bound vs unbound OCRs were
compared by Fisher's exact test. Per-subset KO vs WT and Het vs WT chromVAR deviation
differences were computed from §7.2.

### 7.6 Motif grammar (`extended_analysis/scripts/motif_grammar/`)

We tested whether the spatial arrangement of motifs within an enhancer (spacing, orientation,
helical phasing) predicts dose sensitivity beyond motif composition.

*Contrasts.* Primary: buffered vs haploinsufficient peaks (§7.1), restricted to peaks with a
consistent class across subsets and lost in every subset where classified. Control:
cBAF-dependent (D8 pseudobulk KO vs WT, FDR < 0.05 and LFC < −1) vs cBAF-independent
(FDR > 0.5 and |LFC| < 0.25). Both were limited to distal (non-promoter) peaks on canonical
chromosomes, represented as fixed 300-bp windows centred on the peak, and matched 1:1 on
baseMean × GC-content decile strata (seed 42). After matching, the primary contrast had 2,882
windows per class and the control 6,960.

*Motifs.* 115 JASPAR2020 CORE vertebrate PWMs from seven families (ETS, RUNX, T-box, AP-1/bZIP,
NF-κB, TCF/LEF, KLF/SP) were scanned with motifmatchr (p < 1e-4), keeping position and strand.
To remove PWM redundancy, same-family instances within 10 bp were collapsed to the
highest-scoring site.

*Features and model.* Composition features were per-family site counts and maximum PWM scores.
Grammar features were per-family-pair fractions and distances, so that they do not scale with
motif number: fraction of pairs within 50 bp, minimum centre distance, same-strand and convergent
orientation fractions, and helical phasing (10.5-bp periodicity within 100 bp). Elastic-net
logistic regression (glmnet 4.1.10, α = 0.5, inner 5-fold CV, `lambda.min`) was evaluated by
chromosome-held-out 5-fold cross-validation (chromosomes split into five blocks). The test
statistic was ΔAUC = AUC(composition + grammar) − AUC(composition) on pooled out-of-fold predictions.

*Uncertainty.* ΔAUC was given 95% intervals by (i) per-fold ΔAUC, (ii) a paired peak-level
bootstrap (2,000 resamples, both AUCs recomputed on each resample) and (iii) a chromosome-cluster
bootstrap. These intervals resample the evaluation of the fitted models; a label-permutation null
with refitting was not run.

*Unbiased spacing.* MEME suite 5.5.4 SpaMo (`-margin 150 -range 150 -bin 5`) was run with ETS1,
RUNX1, TBX21 and an AP-1 motif as primaries on each of the four class sequence sets (16 runs).

*Composition.* Per family, two measures were compared between classes: mean site count per
window (log2 ratio with 500-resample bootstrap CI) and share, meaning each family's fraction of a
window's sites. Both used Wilcoxon tests with BH correction within each contrast.

*Result.* The result was null. Grammar did not improve prediction for either contrast: primary
ΔAUC = +0.0012 (95% CI, peak bootstrap −0.0053 to +0.0074; chromosome bootstrap −0.0028 to
+0.0045) and control ΔAUC = +0.0001 (−0.0021 to +0.0022; −0.0027 to +0.0027). The sign of the
per-fold ΔAUC varied between folds, and SpaMo found no significant secondary-motif spacing in any
of the 16 runs. Composition-only models reached AUC 0.62 (primary) and 0.68 (control). The only
positively enriched family in buffered windows was ETS, as a share of sites (+7.7 percentage
points, FDR = 4.3 × 10⁻⁴²). In the control, every family was enriched by count, which indicates
residual motif-density imbalance; by share, ETS enrichment was not significant and RUNX and T-box
remained enriched. **Share, not count, is therefore the reported composition metric, and ΔAUC is
always reported with its confidence interval rather than as a bare point estimate.** Full results
are in `docs/reports/motif_grammar.md`.

---

## 8. Cross-study comparison with Guo *et al.* 2022

### 8.1 Data and processing

Guo *et al.* (2022, *Nature* 607:135; GSE183619) profiled *Arid1a*-deficient and c-Myc-deficient
CD8⁺ T cells after *Listeria*-OVA infection. All bulk sub-series were reprocessed on mm39/GENCODE
vM35 with the same pipeline releases as the main data:

- RNA-seq (38 paired-end samples): GSE183615 (*Arid1a* WT/KO 5/5; c-Myc WT/KO 10/10) and GSE199184
  (vehicle vs ARID1A inhibitor BRD-K98645985, 4/4). nf-core/rnaseq 3.22.2, `--aligner star_salmon`,
  strandedness auto-detected, reusing the main STAR index (`scripts/upstream/run_meta_guo_rnaseq.sh`).
- ATAC-seq (60 runs, 30 libraries): GSE183618 (WT, *Arid1a* KO, Myc KO), GSE198894 (vehicle WT,
  vehicle KO, inhibitor) and GSE183616 (Myc-high, Myc-low, naive). nf-core/atacseq 2.1.2 with the
  main-project parameters, `--read_length 50`, reusing the main Bowtie2 index
  (`scripts/upstream/run_meta_guo_atacseq.sh`). All Guo libraries share one consensus peak set.

Baxter *et al.* (2023) was not included: its bulk perturbations target PBAF (*Arid2*, *Pbrm1*)
rather than cBAF, and its remaining data are single-cell multiome.

### 8.2 Differential analysis (`extended_analysis/scripts/guo2022/01_rnaseq_de.R`, `guo2022/02_atacseq_da.R`)

The main-project DESeq2 procedures (§5) were reused, with sample metadata parsed from library
names and without a batch covariate. Each contrast used only samples from one GSE and one
experimental arm: RNA KO vs WT (GSE183615), Myc KO vs Myc WT (GSE183615) and inhibitor vs vehicle
(GSE199184); ATAC KO vs WT and Myc KO vs Myc WT (GSE183618), vehicle KO vs vehicle WT and inhibitor
vs vehicle WT (GSE198894), and Myc-high vs Myc-low and Myc-high vs naive (GSE183616). Thresholds
were |log2FC| > 1 and FDR < 0.05, which are also Guo's published thresholds. Unshrunken LFCs are
reported alongside apeglm estimates because Guo did not apply shrinkage.

### 8.3 The GSE183618 two-arm "WT" split

GSE183618 deposits six WT ATAC libraries under one label, but they are controls for two separate
experiments: three for the *Arid1a* arm (GSM5563555–7) and three for the c-Myc arm
(GSM5563561–3). The assignment was established in three ways:

1. SRR → GSM mapping from the GEO file list.
2. A reads-per-gigabyte batch signature that groups each WT trio with its own KOs (about twofold
   separation between arms).
3. Hierarchical clustering (Jaccard similarity of base-pair overlap, average linkage) of Guo's own
   deposited per-sample peak BEDs. At k = 2 this splits the libraries exactly along the arm boundary.

Note that nf-core replicate numbers follow samplesheet order, not GEO replicate numbering. The c-Myc
arm WT libraries were relabelled `MycWT` before any analysis (guarded by assertions in `guo2022/02_atacseq_da.R`),
and **the two WT groups are never pooled**. Pooling them produces a spurious accessibility
collapse. With the same KO libraries, comparing against the wrong-arm WT alone gives 1,239 gained
and 6,886 lost peaks, whereas the correct within-arm comparison gives 50 gained and 65 lost.

In the correctly paired *Arid1a* arm, KO and WT samples do not separate. The design-integrity check
fails (within-group r 0.9700 vs between-group r 0.9682; three of six samples have a nearest
neighbour in the other group), whereas the c-Myc arm passes with a gap about ten times larger. The
GSE183618 *Arid1a* KO vs WT result is therefore treated as low-confidence. GSE183616 also fails
design integrity (negative gaps; within-condition replicate r ≈ 0.91) and produces no differential
peaks. The GSE198894 contrasts pass all checks and are the primary Guo comparison for ARID1A loss.

### 8.4 Attempted reproduction of Guo Extended Data Fig. 6a

Guo report 11,975 lost and 211 gained peaks for *Arid1a* KO vs WT (n = 3 vs 3; DiffBind 2.16,
`summits = 250`, DESeq2, FDR < 0.05, fold change > 2, mm10). We tested three quantification routes
on the correctly paired arm:

| Route | Regions tested | Gained | Lost |
|---|---:|---:|---:|
| nf-core consensus + featureCounts + DESeq2 (`guo2022/02_atacseq_da.R`) | 179,872 | 50 | 65 |
| Summit-recentred 501-bp windows from the six arm libraries, kept if in ≥ 2 replicates of either condition; featureCounts `-O --fracOverlap 0.2 -p`; same DESeq2 contrast (`guo2022/03_summit_requant.R`) | 87,765 | 34 | 37 |
| DiffBind 2.16.0 in the Bioconductor 3.11 container, `dba.count(summits = 250)`, `DBA_DESEQ2`, all other defaults (full-library-size normalisation, no shrinkage), Guo's thresholds (`guo2022/04_diffbind216.R`) | 92,563 | 47 | 49 |

The DiffBind p-value distribution was close to uniform: 5,450 regions had raw p < 0.05, against
4,628 expected by chance, and only 3,042 regions had |FC| > 2 at any p-value. FRiP was 0.25–0.28
in all six libraries. Thresholds, WT pairing, statistical power, data quality, normalisation,
shrinkage, peak universe and software version were each ruled out as explanations. Genome build
(mm39 vs mm10) was not tested but is not a plausible cause of a 33:1 imbalance arising from a
uniform null. We therefore conclude that Extended Data Fig. 6a does not reproduce from the
deposited data, with the caveat that Guo's DiffBind run itself could not be audited.

---

## 9. Software versions

The authoritative version record is `data/metadata/software_versions.tsv` (workflow, pipeline and
host tool versions), together with `renv.lock` for R packages. The pipeline-internal tool versions
cited in §3 come from each pipeline's `pipeline_info/` output. Principal versions:

- Workflow: Nextflow 25.10.0; nf-core/atacseq 2.1.2, rnaseq 3.22.2, cutandrun 3.2.2, chipseq 2.1.0.
- Host tools: Bowtie2 2.5.2, STAR 2.7.11a, samtools 1.13, bedtools 2.30.0, MACS2 2.2.9.1, MACS3
  3.0.2, deepTools 3.5.4, TOBIAS 0.16.1, HOMER 5.1, MEME suite 5.5.4.
- R 4.4.2 with DESeq2, apeglm, ChIPseeker, clusterProfiler, fgsea, msigdbr, chromVAR 1.28.0,
  motifmatchr 1.28.0, JASPAR2020, TFBSTools, DEGreport, Mfuzz, cluster, glmnet 4.1.10,
  GenomicRanges, ComplexHeatmap and ggplot2 (exact versions in `renv.lock`).
- DiffBind 2.16.0 (Bioconductor 3.11 container `bioconductor-diffbind:2.16.0--r40h5f743cb_2`),
  used only for §8.4.

---

## 10. Deviations from the original publication

| Aspect | Original publication | This reanalysis |
|---|---|---|
| Genome | mm10 | GRCm39 (mm39), GENCODE vM35 |
| Alignment (ATAC, ChIP, CUT&RUN) | STAR, default parameters | Bowtie2 via nf-core pipelines, with the pipelines' duplicate, MAPQ and mitochondrial filtering |
| RNA-seq quantification | STAR 2.5.3a + HOMER `analyzeRepeats.pl -strand both -count exons -condenseGenes -noadj` | STAR 2.7.11b + Salmon 1.10.3 gene counts (nf-core/rnaseq) |
| RNA-seq DE | HOMER `getDiffExpression.pl` (DESeq2), \|log2FC\| ≥ 0.585, FDR < 0.05 | DESeq2 in R, apeglm shrinkage, \|log2FC\| > 1, FDR < 0.05, per-stratum models |
| ATAC-seq peaks | HOMER `findPeaks -style dnase` (4-fold over local) | MACS2 narrow peaks; nf-core merged-replicate consensus quantified with featureCounts |
| Differential accessibility | HOMER `getDifferentialPeaksReplicates.pl` (DESeq2), FC ≥ 2, FDR < 0.05 | DESeq2 on the consensus count matrix, apeglm, with an experiment covariate for D8 contrasts spanning Exp1 and Exp2 |
| ChIP-seq peaks / differential | HOMER `findPeaks -style factor -i input`; `getDiffExpression.pl` FC ≥ 1.5, Poisson p < 1e-4 | MACS3 against input; peak-overlap comparisons only (n = 1 per condition) |
| CUT&RUN peaks | SEACR stringent, `norm`, vs IgG | SEACR 1.3 stringent, `non`, vs IgG, plus MACS2 (q < 0.01) vs IgG; CPM tracks; only IgG-matched (D5) libraries peak-called |
| Motif enrichment | HOMER `findMotifsGenome.pl`, ± 100 bp, GC-matched random background | HOMER 5.1 `-size 200 -mask`; explicit class backgrounds in the dose-response analysis; JASPAR2020 for chromVAR, TOBIAS and grammar analyses |
| TF network ranking | Taiji PageRank | Re-implemented in R (JASPAR2020 motifs, distance-based peak–gene links, personalised PageRank; see §6) |
| OCR clustering | DESeq2-based k-means including the 48 h/D3 column | Paper panels use the published cluster definitions (original mm10 regions, lifted to mm39); an independent re-clustering on Naive/D3/D5/D8 (48 h ATAC excluded for QC) is reported as a concordance check |
| Blacklist | Not stated | Lifted ENCODE mm10 v2 blacklist, applied in deepTools only |
| deepTools | 3.5.1 | 3.5.4 (host), 3.5.1 inside nf-core containers |
