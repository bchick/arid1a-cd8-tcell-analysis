# Motif grammar of ARID1A dose-sensitivity: does enhancer *syntax* carry information beyond motif *presence*?

**McDonald, Chick et al. 2023 reanalysis — extended analysis**
**Scripts:** `extended_analysis/scripts/motif_grammar/01_region_sets.R` → `motif_grammar/02_grammar_features.R` → `motif_grammar/03_spamo_spacing.R` → `motif_grammar/04_noise_floor.R`
**Run:** `make extensions` (or the four scripts in order, from the repository root)
**Outputs:** `results/extended_analysis/motif_grammar/`, `figures/extended_analysis/motif_grammar/`

---

## TL;DR

Across two independent methods — an elastic-net predictive test and unbiased SpaMo spacing discovery — **motif *arrangement* (grammar) adds nothing to a model that already knows motif *composition***. The grammar gain in held-out AUC is **+0.0012** for the primary contrast and **+0.0001** for the control, and both 95% bootstrap CIs comfortably span zero. SpaMo finds **zero** significantly enriched spacings in all 16 runs.

This is a **substantive negative result, not a failed analysis**: what separates a dose-buffered enhancer from a dose-sensitive one is *which transcription-factor families bind it*, not *how their binding sites are geometrically arranged*. The single strongest predictor is **ETS site share** — buffered enhancers are ETS-dominated — which is consistent with the ETS dose-buffering program established by the footprinting analysis (`docs/reports/tobias_footprinting.md`).

The positive control behaves as a control must (it recovers ETS enrichment), with one caveat quantified below.

---

## 1. The question

The ARID1A dose-response taxonomy (`docs/reports/het_dose_response.md`) splits cBAF-responsive enhancers by how they respond to *gene dosage* — WT (two Arid1a copies) vs Het (one) vs KO (zero):

- **buffered** — one copy suffices; accessibility is maintained at Het and lost only in KO.
- **haploinsufficient** — one copy is not enough; accessibility already erodes at Het.

Composition analysis (which TFs bind) is well trodden. This analysis asks the next question: is dose-sensitivity **encoded in the spatial grammar** of the enhancer — the spacing, orientation, and helical phasing between motifs — over and above which motifs are present? Composite elements (e.g. ETS:RUNX with a fixed 3-bp gap) and helical-phasing constraints are real phenomena, so the hypothesis is not idle.

**Design.** Two contrasts, each reduced to a balanced, confounder-matched two-class problem in fixed 300-bp windows:

| Contrast | Positive class | Negative class | Role |
|---|---|---|---|
| **Primary** | buffered | haploinsufficient | Test question: sequence-encoded dose sensitivity |
| **Control** | cBAF-dependent (lost, LFC < −1) | cBAF-independent (unchanged) | Positive control — the paper's Fig 5C-style split; the pipeline *must* recover the known ETS enrichment |

---

## 2. Why matching was mandatory

Signal strength drives how many motif instances are detectable, and promoters are motif-dense. Before matching, the classes differed sharply on both:

- **buffered vs haploinsufficient:** median baseMean 79.4 vs 30.4 (2.6×); promoter fraction 31.5% vs 20.2%
- **cBAF-dependent vs -independent:** median baseMean 46.1 vs 21.5 (2.1×); promoter fraction 16.1% vs 27.8%

Left uncorrected, any "grammar" signal would be a proxy for read depth and promoter content. Stage 1 therefore restricts to **distal (non-promoter) windows on canonical chromosomes**, drops peaks with inconsistent class assignment, uses **lost-only** peaks for the buffered class (the mixed set pools 15,106 lost with 7,667 gained; the gained peaks likely reflect library-size renormalisation drift after the global loss of accessibility, a different phenomenon), and then **matches classes on baseMean × GC deciles**.

After matching, the confounders are gone:

| Contrast | Class | n | median baseMean | median GC |
|---|---|--:|--:|--:|
| Primary | buffered | 2,882 | 28.9 | 0.487 |
| Primary | haploinsufficient | 2,882 | 28.4 | 0.487 |
| Control | cBAF-dependent | 6,960 | 20.0 | 0.477 |
| Control | cBAF-independent | 6,960 | 19.8 | 0.477 |

Every result below is on these balanced, confounder-matched sets.

---

## 3. Methods in brief

- **Motif scan.** 115 of 746 JASPAR2020 CORE vertebrate PWMs, restricted to seven T-cell–relevant families (ETS, RUNX, T-box, AP-1/bZIP, NFκB, TCF/LEF, KLF/SP), scanned with strand and position on every window. 133,048 primary / 287,843 control instances.
- **PWM-redundancy merge.** JASPAR ships ~24 ETS matrices; one GGAA site matches many, so raw counts measure PWM redundancy, not biology. Instances of the same family within 10 bp are collapsed to the best-scoring representative (**~80% of instances collapse**), giving 25,749 / 56,754 merged sites.
- **Composition features.** Per family: merged-site count and max PWM score. This is the baseline the grammar must beat — the baseline is *not* "no motifs."
- **Grammar features.** All expressed as **fractions or distances**, never counts, so they are ~orthogonal to composition: close-pair fraction (within 50 bp), strand configuration (same/convergent), and helical phasing (10.5-bp periodicity within 100 bp).
- **Model.** Elastic-net logistic regression (α = 0.5), **chromosome-held-out 5-fold CV** — random CV would leak, because near-identical sequences (duplicated peaks, shared repeats) would straddle train/test. Metric: **ΔAUC = AUC(composition + grammar) − AUC(composition)**, the out-of-sample gain from adding syntax.
- **Noise floor (Stage 4).** The ΔAUC point estimate is given an error bar three ways: per-fold ΔAUC spread; a **paired peak-level bootstrap** (2,000 resamples, both AUCs recomputed on each resample so shared model error cancels); and a **chromosome-cluster bootstrap** (resampling whole chromosomes). All operate on the stored out-of-fold predictions — no refitting.
- **SpaMo (Stage 3).** Hypothesis-free spacing discovery: given a primary motif (ETS1, RUNX1, TBX21, JUN), which secondary motif is enriched at which exact bp offset and orientation, against a background of the same sequences. Run per class (4 primaries × 4 class sets = 16 runs).

---

## 4. Results

### 4.1 Headline — grammar adds nothing (both contrasts)

![ΔAUC with bootstrap confidence intervals against zero](../../figures/extended_analysis/motif_grammar/F2_dauc_noise_floor.png)

| Contrast | AUC composition | AUC comp + grammar | **ΔAUC** | 95% CI (peak) | 95% CI (chr-cluster) | Spans 0? |
|---|--:|--:|--:|:--:|:--:|:--:|
| **Primary** | 0.6157 | 0.6169 | **+0.00124** | [−0.0053, +0.0074] | [−0.0028, +0.0045] | **Yes** |
| **Control** | 0.6782 | 0.6783 | **+0.00007** | [−0.0021, +0.0022] | [−0.0027, +0.0027] | **Yes** |

Both ΔAUCs are an order of magnitude smaller than their CI half-widths. Adding the entire grammar feature block moves held-out discrimination by ~0.1 AUC points for the primary and by essentially nothing for the control.

### 4.2 The noise floor makes the null quantitative

![Per-fold AUC, composition vs composition + grammar](../../figures/extended_analysis/motif_grammar/F3_per_fold_auc.png)

Per held-out chromosome block, the grammar model wins on some folds and loses on others (primary ΔAUC range **−0.0080 to +0.0067**; control **−0.0057 to +0.0044**). The sign of the gain is **not stable across folds** — the signature of no real effect. The bootstrap distributions are centered essentially on zero with the observed ΔAUC sitting near the mode:

![Bootstrap distribution of the grammar gain](../../figures/extended_analysis/motif_grammar/F5_bootstrap_distribution.png)

This is why a point estimate alone would have been misleading. +0.0012 *sounds* like a gain; against a bootstrap SD several times its size, it is indistinguishable from zero.

### 4.3 SpaMo agrees — zero enriched spacings

Independent of the elastic net, SpaMo tested every primary motif in every class set for a secondary motif enriched at any spacing/orientation:

**16 / 16 runs executed cleanly; 0 significant secondary spacings** (`spamo_run_manifest.csv`, `spamo_all_results.csv`).

Two methods with entirely different assumptions — a regularized predictive model and a per-offset enrichment scan — **converge on the same negative answer**. That convergence is the strongest evidence in this report.

### 4.4 What *does* separate the classes — composition, and specifically ETS

If grammar is silent, composition should not be. It is not:

![Motif family composition, count-based log2 ratio](../../figures/extended_analysis/motif_grammar/F1_composition_enrichment.png)

**Primary (buffered vs haploinsufficient) — a clean, interpretable composition signal:**

- **ETS is the only positively enriched family** (log2 ratio **+0.327**, FDR = 2.4×10⁻³⁴). Buffered enhancers carry more ETS sites.
- RUNX, T-box, and AP-1/bZIP are all *depleted* in buffered (log2 ratios −0.28 to −0.31, all FDR < 10⁻⁷).

Density-normalization sharpens rather than dissolves this, unlike the control (next):

![Density-normalized composition — each family's share of a window's sites](../../figures/extended_analysis/motif_grammar/F4_composition_share.png)

- Buffered windows devote a **+7.7-point higher share** of their motif content to ETS (FDR = 4.3×10⁻⁴²), with a corresponding drop in every other active family. Buffered enhancers are **ETS-dominated**; haploinsufficient ones spread their sites across RUNX/T-box/AP-1.

**Main point:** the sequence feature that predicts whether one ARID1A copy suffices is *how ETS-centric the enhancer is* — congruent with the footprinting result that ETS occupancy is the dose-buffered backbone (82% footprint retention at one copy; see the footprinting report). It is a composition story, not a grammar story.

---

## 5. The positive control passed — with a caveat worth stating

**Pass:** ETS is significantly enriched in cBAF-dependent (lost) enhancers (log2 ratio **+0.304**, FDR = 1.8×10⁻⁵⁶), recovering the paper's known ETS dependence. The pipeline is not broken.

**Caveat:** in the control, ETS is *not* the top family by count — **RUNX (+0.906) and T-box (+0.546) score higher**, and *all seven* families are enriched (every FDR < 0.05). That pattern is the fingerprint of a residual global motif-density difference that baseMean × GC matching did not fully remove.

The **share analysis resolves it**. Once each family is expressed as a *fraction* of a window's sites, the control's ETS enrichment **collapses to non-significant** (share diff −0.001, FDR = 0.074), while **RUNX (+0.047, FDR = 1.1×10⁻⁹⁴) and T-box (+0.015, FDR = 1.3×10⁻²¹) survive**. So the honest reading of the control is: cBAF-dependent enhancers are *motif-dense*, and among families the ones whose *relative* representation genuinely rises are **RUNX and T-box**, not ETS. ETS rises in absolute count simply because everything does.

This does not undermine the primary result — it strengthens the case for reporting **share, not count**, which is exactly where the primary's ETS signal is strongest and most specific. It is, however, a caveat worth keeping: the count-based control overstates ETS.

---

## 6. Interpretation and limits

**What this shows.** ARID1A dose-sensitivity of an enhancer is predicted by its TF-family *composition* — above all its ETS-centricity — and **not** by the spatial grammar among those motifs, at the resolution these features and SpaMo probe. Two independent methods agree.

**Limits.**
- **A null within a modest ceiling.** Even the full model tops out at AUC ~0.62 (primary) / ~0.68 (control). Sequence composition explains *some* of the taxonomy, not most of it; the rest is presumably trans (cofactor availability, pioneer dynamics) or chromatin context not captured by 300-bp motif content.
- **Grammar is only as good as its featurization.** We tested close-pair fraction, strand configuration, and helical phasing. A more exotic grammar (long-range order, specific ordered triples, exact composite spacings) could carry signal these features smooth over — though SpaMo's per-bp scan, which found nothing, argues against strong pairwise spacing rules specifically.
- **Control caveat (Section 5):** count-based enrichment overstates ETS in the control; prefer the share metric.
- **Noise floor bootstraps the evaluation, not the model fit.** A full permutation null that re-fits under label shuffling (~1 h/contrast) was deferred; for a null claim, the evaluation CI is the directly relevant statement and it already spans zero.

**How to report it.** Report this as a deliberate, well-powered test of the grammar hypothesis that returns a clean negative, and pivot to the positive composition finding: *buffered enhancers are ETS-dominated*, tying the sequence analysis back to the ETS footprinting program. Quote **share**, not count, and quote the **ΔAUC with its CI**, not the bare point estimate.

---

## Appendix — file manifest

**Figures** (`figures/extended_analysis/motif_grammar/`, PDF + PNG @ 300 dpi):
`F1_composition_enrichment`, `F2_dauc_noise_floor`, `F3_per_fold_auc`, `F4_composition_share`, `F5_bootstrap_distribution`

**Tables** (`results/extended_analysis/motif_grammar/`):
`grammar_model_summary_all.csv` (AUCs, ΔAUC, feature-selection counts), `grammar_dauc_noise_floor.csv` (CIs), `grammar_dauc_per_fold.csv`, `grammar_composition_enrichment.csv` (count + share, per family), `spamo_run_manifest.csv` (16 runs), `spamo_all_results.csv` (empty — the null), `{primary,control}_*_coefficients.csv`, `balance_{primary,control}.csv` (matching QC).
