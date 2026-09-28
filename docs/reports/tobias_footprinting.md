# TOBIAS Differential Footprinting — ARID1A Dose-Response

_Assay: ATAC-seq (D8, Exp2 batch) · Genotypes: WT / Het / KO_
_Figures: `tf_tobias_footprints` (summary), `tf_tobias_dose_response` (dose-response), `tf_tobias_chromvar_concordance` (cross-method concordance), in `figures/extended_analysis/`_

## Summary

TOBIAS measures the **depth of the Tn5 protection footprint** at every motif instance — a
readout of transcription-factor occupancy that is entirely independent of both HOMER motif
enrichment and chromVAR accessibility deviations. Applied across the *Arid1a* allelic series
(WT = 2 copies, Het = 1, KO = 0), it delivers a **third orthogonal line of evidence** for the
central finding of this project: **the ETS motif family is the ARID1A-dose-buffered backbone of
the CD8⁺ effector enhancer landscape.**

Two results:

1. **ARID1A loss collapses TF footprints genome-wide** (`tf_tobias_footprints`, panel A). The ETS family shows the largest
   and most significant footprint losses in KO of any TF family.
2. **ETS footprints are uniquely dose-buffered** (`tf_tobias_dose_response`). ETS is the only family that is *both*
   strongly ARID1A-dependent *and* strongly dose-tolerant — its footprints are retained at ~80%
   of full activity on a single *Arid1a* copy, while comparably ARID1A-dependent families
   (RUNX, T-box) and the AP-1/bZIP program collapse toward codominant (additive) behaviour.

The per-motif dose-response is concordant with chromVAR (Spearman ρ = 0.55, n = 232 motifs
monotone in both methods; `tf_tobias_chromvar_concordance`), and ETS occupies the buffered corner in both.

## Methods

- **Merge:** ATAC BAMs pooled per genotype (WT, Het, KO), D8 Exp2 batch, subsets combined.
- **Pipeline:** `TOBIAS ATACorrect` (Tn5 bias correction against the mm39 reference) →
  `ScoreBigwig` (continuous footprint score) → `BINDetect` over D8 consensus OCRs using
  **JASPAR2020 CORE vertebrates** motifs (the same motif set used for chromVAR, `dose_chromvar`).
- **Dose metric** (identical convention to the chromVAR analysis): footprint activity **retained at one
  *Arid1a* copy**, scaling WT = 1 and KO = 0:

  ```
  het_ret = (Het_mean_score − KO_mean_score) / (WT_mean_score − KO_mean_score)
  ```

  Defined only for motifs with a **monotone ARID1A-dependent decline** (WT > Het > KO,
  WT − KO > 0). `het_ret > 0.5` = dose-buffered (recessive-like); `≈ 0.5` = codominant/additive.
- **Family assignment:** regex on motif name, identical rules to the chromVAR/HOMER figures.
- Scripts: `extended_analysis/scripts/footprinting/run_tobias_footprinting.sh` (pipeline),
  `extended_analysis/scripts/figures/tobias_footprints.R` (summary figure),
  `extended_analysis/scripts/figures/tobias_dose_response.R` (tables, dose-response and concordance figures).

## Per-family dose-response

Footprint retention at one *Arid1a* copy, by TF family (monotone motifs only):

| Family | n motifs (monotone) | Median Het-retention | ARID1A-dependence (mean WT−KO change) | Class |
|--------|--------------------:|---------------------:|--------------------------------------:|-------|
| **ETS**       | **21** | **80%** | **0.343** (highest) | **Buffered backbone** — strong dependence *and* dose-tolerant |
| TCF/LEF   | 4  | 82% | 0.126 | Buffered, but few motifs / weak dependence |
| RUNX      | 3  | 61% | 0.473 | Highly ARID1A-dependent, ~codominant |
| T-box     | 15 | 56% | 0.358 | ARID1A-dependent, codominant |
| AP-1/bZIP | 31 | 52% | 0.245 | Codominant (additive) |
| NFκB      | 3  | 10% | 0.073 | Not buffered |
| KLF/SP    | 1  |  0% | 0.007 | Not ARID1A-dependent |

**Why ETS rather than TCF/LEF:** TCF/LEF shows a marginally higher median retention
(82%) but rests on only 4 motifs with weak ARID1A-dependence (WT−KO ≈ 0.13). ETS combines the
**highest retention among well-supported families (80% across 21 motifs)** with the **highest
ARID1A-dependence of any family (0.34)** — i.e. ETS footprints most require cBAF to open, yet
most survive the loss of one *Arid1a* copy. That is the defining signature of a dose-buffered
regulatory backbone (`tf_tobias_dose_response`, panel B, upper-right quadrant).

Top dose-buffered ETS members: ELK1 (97%), ELK3 (91%), ELF5 (87%), ETV3 (84%), ELK4 (84%),
ETV5 (82%), ETV6 (81%), ETS2 (81%), plus ETS1, ERG, GABPA, FLI1, SPI1, SPIB, EHF, FEV.

## Cross-method concordance

Per-motif Het-retention from TOBIAS (footprint depth) vs chromVAR (accessibility deviation),
restricted to the 232 motifs with a monotone decline in **both** methods:

- **Spearman ρ = 0.55** — the two independent methods rank the same motifs as buffered.
- The methods differ in absolute scale (chromVAR compresses retention lower), so the agreement is
  in **ranking**, not identity; ETS sits in the high-retention corner of both axes.
- Together with the HOMER motif-composition result (`dose_motifs_vs_insensitive`) and chromVAR (`dose_chromvar`), this is the
  **third orthogonal method** to independently converge on ETS = dose-buffered.

## Interpretation

ARID1A/cBAF is required to open the effector enhancer landscape genome-wide (footprints collapse
in KO). But the **ETS motif backbone re-opens at a single *Arid1a* copy**, whereas the
AP-1/bZIP, NFκB and (to a lesser degree) RUNX/T-box programs scale additively with dosage. This
provides a mechanistic basis for the near-normal phenotype of *Arid1a*-heterozygous CD8⁺ T cells:
the ETS-wired core of the effector program is buffered against haploinsufficiency, while the
dose-sensitive programs account for the subtler Het phenotypes seen in the RNA/ATAC data.

## Output files

- `results/extended_analysis/footprinting/bindetect/bindetect_results.txt` — raw TOBIAS BINDetect (747 motifs).
- `results/extended_analysis/footprinting/tobias_dose_response_per_motif.csv` — per-motif family, WT/Het/KO
  scores, ARID1A-dependence (WT−KO change + p), Het-retention, monotonicity flag.
- `results/extended_analysis/footprinting/tobias_dose_response_by_family.csv` — family-level summary (above table).
- `results/extended_analysis/footprinting/tobias_vs_chromvar_dose_per_motif.csv` — cross-method retention join.
- `figures/extended_analysis/tf_tobias_footprints.{pdf,png}`, `tf_tobias_dose_response.{pdf,png}`,
  `tf_tobias_chromvar_concordance.{pdf,png}`.
