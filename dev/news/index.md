# Changelog

## easyRasch2 (development version)

### Changes

- **`output = "dataframe"` now returns unrounded values** across the
  package; rounding is a presentation concern and now happens only when
  rendering the kable (whose displayed precision is unchanged). Affected
  functions:
  [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md),
  [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md),
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md),
  [`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreBoot.md),
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  (the `$pairs` table),
  [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md),
  [`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md),
  [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md),
  [`RMdimCFA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFA.md)
  (both `$fit` and `$loadings`),
  [`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md),
  [`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md),
  [`RMdifTree()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md)
  (including the stability summary attribute),
  [`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md)
  (the `wle_correlation` element),
  [`RMdimMartinLofResiduals()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLofResiduals.md),
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md)
  (dataframe and file/CSV output), and
  [`RMpersonFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonFit.md).
  [`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md)
  already followed this convention.
- Consequently, `Flagged` labels and sorting are now always computed on
  the exact (unrounded) values. Items or pairs sitting exactly on a
  rounded boundary could in principle change flag relative to earlier
  releases; displayed tables are otherwise identical.
- Deliberately unchanged: the infit *simulation workers* still round the
  per-iteration statistics to 3 decimals before the cutoff computation,
  so simulation-based expected ranges reproduce earlier releases
  exactly.
- Minor layout changes for
  [`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md)
  regarding number of bins if not set by user. Earlier it was sum
  score + 1, now it is sum score divided by 2.
  [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md)
  gets slightly larger points and thicker error bars.

### Improvements

- [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md)
  gains two grouping methods for the observed means, alongside the
  existing quantile grouping (`method = "cut"` is now a legacy alias for
  `method = "quantile"`): `method = "width"` for equal-width intervals
  on the total-score scale, and `method = "manual"` with the new
  `score_breaks` argument – the total scores at which a new group starts
  (the concept of RASCHplot’s `lower.groups`, without the leading zero).
  The documentation of the grouping rules and of the `error_band` (the
  model-implied interval for the observed mean per total score) was
  substantially expanded.
- [`RMplotTile()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotTile.md)
  gains a `text_size` argument for the cell labels (default 4, matching
  [`RMitemCatProb()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatProb.md)),
  and cells with **zero responses are now visually distinct**: they are
  taken off the colour scale and shown unfilled with a light grey
  outline (`zero_fill = "white"`, set to `NULL` for the previous
  behaviour), so the absence of data reads as absence rather than as the
  darkest end of the fill scale, where `n = 0` was nearly
  indistinguishable from `n = 1`. The caption notes the convention when
  empty cells are present.
- [`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md)’s
  `output = "dataframe"` tables now include the `se`, `lower`, and
  `upper` (95% Wald CI) columns for each pair, matching
  [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md);
  the formatted kable output is unchanged.

### Bug fixes

- [`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md)
  is now fully reproducible with the same `seed`: the RMU estimate
  previously varied slightly between identical calls because `mirt`’s
  Metropolis-Hastings plausible-value sampler leaves the R random-number
  stream in a nondeterministic state (the draws themselves are
  reproducible), which perturbed the RMU split-half column assignments
  downstream. The function now re-seeds before the RMU iterations.
- Respondents with **no responses at all** (all-`NA` rows) no longer
  crash the functions that fit a CML model on data retaining missing
  values.
  [`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
  errors on all-NA rows with an opaque “invalid argument type” (and
  [`psychotools::raschmodel()`](https://rdrr.io/pkg/psychotools/man/raschmodel.html)
  segfaults), so such rows are now dropped up front — with the standard
  “`N respondent(s) with no responses dropped.`” message — in
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md),
  [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md),
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md),
  [`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreBoot.md),
  [`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md),
  and the observed-data overlays of
  [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md)
  and
  [`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md).
  [`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md)
  now drops such rows too (jointly with `dif_var`) instead of failing
  with [`eRm::LRtest()`](https://rdrr.io/pkg/eRm/man/LRtest.html)’s
  error. Functions that already used complete cases or called
  `.drop_empty_respondents()` are unaffected.
- [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)’s
  `sample_n` now counts the respondents actually used (all-NA rows
  excluded, incomplete responses still retained) and `sample_n_total`
  the raw input rows, matching
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  and the other `*Cutoff()` objects; previously the two were documented
  as always equal. Captions in the affected functions keep reporting the
  raw total in the `n = X of Y respondents (policy)` form.

## easyRasch2 1.0.0

CRAN release: 2026-07-05

A large consolidation release. The headline change is that the whole
package now runs on a **single estimation engine** — conditional maximum
likelihood (CML) item parameters via `psychotools` and Warm’s
weighted-likelihood (WLE) person locations — replacing the previous mix
of `eRm` (MLE) and `mirt` (MML/EAP). On top of that: new CFA loading
diagnostics, three new functions for person+item parameters and person
fit, bootstrap p-values with multiplicity correction, and a round of
naming/output consistency fixes. After this release `eRm` is used only
by
[`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md)
(Andersen’s LR test) and `mirt` only for the RMU plausible values
([`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)’s
observed-fit overlay moved from an `eRm` fit to
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html),
matching
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md);
values agree to ~1e-6). Accordingly, **`eRm` has moved from Imports to
Suggests**:
[`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md)
errors with an install hint when it is absent.

Many items below shift reported numbers slightly: the underlying
statistics are unchanged, only the estimation engine (WLE locations are
finite at extreme scores, so extreme-score cases are now retained rather
than dropped).

### Unified CML/WLE estimation engine

- **Observed item statistics** —
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md),
  [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md),
  [`RMitemCatProb()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatProb.md),
  [`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemHierarchy.md),
  [`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md)
  — moved to CML (`psychotools`) + WLE. The conditional `iarm`
  statistics are engine-invariant; only minor things move on the common
  grand-mean-zero scale (`Relative_location` via the WLE person mean;
  threshold SEs ~few %; polytomous ICC curve position; the targeting
  person distribution).
  [`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md)
  deliberately stays on
  [`eRm::LRtest()`](https://rdrr.io/pkg/eRm/man/LRtest.html).
- **Parameter / score functions** — the shared `.rasch_fit_cml()` behind
  [`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md),
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md),
  [`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md)
  now fits by CML via `psychotools`. Locations, person estimates and
  score SEs are unchanged (match `eRm` to ~5e-5); only
  [`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md)’s
  threshold SEs shift slightly (psychotools `vcov` vs the eRm delta
  method, ~few %).
- **[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md)**
  — the WLE method now reports the information-based SE
  `1/sqrt(I(theta))` (matching
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md),
  `catR`, `TAM`); point estimates unchanged, SEs differ (larger) at the
  score extremes, and the WLE solver/centring is now shared with
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md).
- **[`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md)**
  — empirical reliability is replaced by a native **marginal
  reliability** (Green 1984; the “Empirical” row is now “Marginal”),
  While
  [`mirt::marginal_rxx()`](https://philchalmers.github.io/mirt/reference/marginal_rxx.html)
  uses `N(0,1)`, this integrates over the estimated latent variance, so
  it is correct on the Rasch logit scale; and the **PSI** is now the
  WLE-based separation reliability (`1 - mean(SEM^2)/Var(theta)`,
  excluding extreme scorers), replacing
  [`eRm::SepRel()`](https://rdrr.io/pkg/eRm/man/SepRel.html) + `mirt`.
  The bootstrap is now fully native (no `mirt` per resample). PSI can
  differ noticeably for scales with many extreme scorers; alpha and RMU
  are unchanged.
- **Simulation-based tools** now both generate and refit on the CML/WLE
  engine:
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)/[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md),
  [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md),
  [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)/[`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md),
  [`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreBoot.md),
  [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md)/[`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCACutoff.md),
  [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md),
  [`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md),
  [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md),
  [`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md).
  The eRm-based internal `extract_item_thresholds()` helper is removed.
  Conditional / scale-invariant statistics (infit/outfit, item-restscore
  classification, the Martin-Löf Monte-Carlo p-value, pooled MI infit
  MSQ/SE) are unchanged; residual-based and simulated quantities shift
  slightly — Q3 values/cut-offs/flags (the statistic is unchanged;
  observed-vs-old-MML correlation is typically \>0.95), PCA eigenvalues
  / variance partition / cut-off, and the gamma cut-offs.
- [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)/[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  keep the previous engine available via `estimator = "MML"` (default
  `"CML"`); the estimator is stored in the cutoff object and reused by
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)/[`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md)
  (a mismatch is overridden with a warning).
  [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)’s
  pooled-SE caveat (Müller 2020) is now documented.

### New functions

- **[`RMpersonFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonFit.md)**
  — per-respondent conditional infit/outfit MSQ and the standardized
  log-likelihood `lz`, computed per response pattern (handles partial
  missingness, no biased person estimate in the residual). Significance
  is by Monte-Carlo resampling under the fitted model, not the
  unreliable asymptotic null (Sinharay 2016; Müller 2020). Output as
  table, data.frame, or a named list of person-fit maps; `zstd = TRUE`
  adds the Wilson-Hilferty transform (comparability only, not for
  inference); `flag = "underfit"` restricts flagging to the
  validity-relevant direction.
- **[`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md)**
  — item difficulty (dichotomous) or Andrich-threshold (polytomous)
  parameters in long/wide format, with optional SEs and Wald CIs.
- **[`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md)**
  — per-respondent theta and SEM computed on each response pattern
  (handles partial missingness), via Warm’s WLE (default) or EAP under a
  normal prior (SD estimated by marginal ML unless fixed).
- Both
  [`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md)
  and
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md)
  support `output = "file"` to write the result table to a CSV at
  `filename` (the data.frame is also returned invisibly).

### CFA dimensionality: loadings + three-function split

CFA now also reports **per-item standardized factor loadings** (observed
vs a simulated expected range), restructured to match the other
simulation tools:

- [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md)
  is now **simulate-only** — returns the simulation object (null
  distributions + cutoffs for fit indices *and* loadings); it no longer
  computes the observed fit, and `output = "kable"` errors with a
  pointer to
  [`RMdimCFA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFA.md).
  **(Breaking.)**
- [`RMdimCFA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFA.md)
  *(new)* takes `data` + a required cutoff object and returns a list of
  two tables, `$fit` (CFI/RMSEA/SRMR) and `$loadings` (observed loading
  vs the two-sided expected range), as kables (default) or data.frames.
- `RMdimCFAPlot(cutoff_res, data)` now returns a **list of two ggplots**
  (`$loadings`, in the
  [`RMitemInfitCutoffPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)
  style; and `$fit`, the faceted fit-index distributions) and requires
  `data`. **(Breaking: previously a single ggplot.)**

### Bootstrap p-values and multiplicity correction

Six simulation-based diagnostics gain optional bootstrap p-values
computed against their simulated null distributions (`p_value = TRUE`;
the default `FALSE` leaves existing output unchanged). Shared mechanics:
each statistic is studentised by its bootstrap mean/SD; marginal
Monte-Carlo p-values have floor `1/(B+1)`; `correction` offers
Westfall-Young studentised-max step-down FWER (default),
Benjamini-Hochberg or Benjamini-Yekutieli FDR, or `"none"`; `Flagged`
then reflects `padj < alpha` while the simulated effect-size bands stay
in the tables. The full `*Cutoff()` object is required, and \>= 1000
cutoff iterations are recommended (warning below that). (Ferreira 2024;
Westfall & Young 1993.) Per function:

- [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md):
  two-sided per item. The asymptotic BH-adjusted p-value and star
  columns from `iarm` are dropped in this mode (one p-value family per
  table); `p_gamma` / `padj_gamma` replace them.
- [`RMdimCFA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFA.md):
  one-sided in the unfavourable direction for CFI / RMSEA / SRMR,
  two-sided for the per-item loadings, corrected as two separate
  families.
- [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md):
  a single one-sided test of the first-contrast eigenvalue, so no
  multiplicity correction is involved.
- [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md):
  two-sided per item.
- [`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md):
  one-sided per pair for excess positive LD (matching
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)),
  computed once per pair in the canonical direction (rest score = total
  − Item2, the simulated direction) and repeated in the direction-2
  table; correction runs over the full pair family before any `n_pairs`
  display filter; the `iarm` BH columns are dropped as in
  [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md).
- [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md):
  one-sided per pair, folded into the new per-pair `$pairs` table
  returned alongside `$matrix` (observed Q3 vs simulated band,
  directional flag, sorted by departure from the per-pair median,
  optional `n_pairs` cap).

[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md),
[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md),
and the newly extended functions share a “Multiple comparisons” help
section.

### Output and argument consistency

- All table and plot captions now report the estimation sample size in a
  common form: `n = X of Y respondents (<policy>)`, where `of Y` appears
  only when respondents were excluded and the policy note
  (`complete cases` / `incomplete responses retained` /
  `missing values imputed`) only when the input contained missing values
  — complete data reads simply `n = X respondents`. Terminology is
  “respondents” throughout (previously a mix of “persons” / “complete
  cases”). Sample size is newly reported by
  [`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md),
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md),
  [`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemHierarchy.md),
  [`RMitemCatProb()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatProb.md),
  [`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md),
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md),
  [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md),
  [`RMpersonFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonFit.md)’s
  kable, and the descriptive plots
  ([`RMplotBar()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotBar.md),
  [`RMplotStackedbar()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotStackedbar.md),
  [`RMplotTile()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotTile.md));
  all other captions are reworded to the common form. The bootstrap-null
  plots
  ([`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md),
  [`RMdifGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaPlot.md),
  [`RMlocdepGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaPlot.md),
  [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md),
  [`RMdimCFAPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFAPlot.md))
  report the simulation sample the same way, suffixed `per dataset`;
  their `*Cutoff()` objects store `sample_n_total` / `sample_has_na`,
  with a graceful fallback for cutoff objects from older versions.
  Respondents with no responses (all-NA rows) are dropped with a
  one-time message rather than triggering a CML fitting error, and
  [`RMplotTile()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotTile.md)
  likewise messages when rows with `NA` in `group` are dropped
  (previously silent; item-level `NA`s are retained in the descriptive
  plots).
- [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)
  no longer prints an
  [`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html)
  result table on every iteration in sequential (`parallel = FALSE`)
  runs; the silencing sink is now also restored safely on iteration
  failure (hardened in
  [`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md)
  too).
- [`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md)
  captions now report the Andersen LR p-value as a plain number (exact
  to three decimals, `p < 0.001` below that) instead of the
  [`format.pval()`](https://rdrr.io/r/base/format.pval.html) scientific
  notation (`<1e-04`), and the test statistic is shown as χ² instead of
  `chi^2`.
- Naming aligned across functions (no deprecation aliases):
  [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md)
  plot output is now `output = "ggplot"` (`"loadings"` kept as a
  backward-compatible alias);
  [`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md)
  `output = "figure"` → `"patchwork"` (matching
  [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md));
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
  `p.adj` → `p_adj`;
  [`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)
  `output` → `statistic`.
- [`RMitemInfitCutoffPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)
  renamed to
  [`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md),
  matching the `<base>Plot` form of the other bootstrap-cutoff plot
  functions
  ([`RMdifGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaPlot.md),
  [`RMlocdepGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaPlot.md),
  [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md),
  [`RMdimCFAPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFAPlot.md)).
  Because the old name shipped in the 0.8.0 CRAN release, it is kept as
  a **deprecated alias** that warns and forwards (unlike the 0.8.0
  renames, which dropped old names outright). Separately,
  [`RMdimCFAPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFAPlot.md)’s
  first argument `cutoff_res` was renamed to `simfit` to match the other
  plot functions (no alias).
- [`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md)
  now defaults to `output = "kable"` (was `"ggplot"`), matching every
  other function that offers both a `kable` and a `ggplot` output
  ([`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md),
  [`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md),
  [`RMpersonFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonFit.md),
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md),
  [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md)).
  Pass `output = "ggplot"` for the previous default.
- The `Flagged` column now labels misfit **direction** (`"overfit"` /
  `"underfit"` / `""`) instead of logical TRUE/FALSE in
  [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
  /
  [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)
  (column type changes from logical to character).
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
  gains a `Flagged` column and drops the redundant significance-star and
  `Location` columns (trimmed headers + explanatory caption). Note the
  value direction differs: a *high* infit is underfit, a *high*
  restscore correlation is overfit.

### Other changes and fixes

- [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)
  and
  [`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md)
  now accept `mids` objects whose items were imputed as ordered factors
  (as required by mice’s `polr` method): the completed data are coerced
  back to numeric responses internally, rather than erroring on the
  factor columns.
- [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
  and
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  gain an experimental `dgp` argument — `"resample"` (default; resample
  WLE locations and simulate, a marginal null) vs `"conditional"`
  (simulate each pattern from the exact Rasch conditional given the
  observed total score, a matched conditional null) — stored in the
  result. See `dev/q3_dgp_comparison.R` / `dev/infit_dgp_comparison.R`.
- Q3 plot and table outputs now share a `$matrix` / `$pairs` structure.
  [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md)
  returns a list of two plots: `$pairs` (the per-pair
  simulated-vs-observed dot-interval, previously the single return
  value) and `$matrix` (a lower-triangle diverging-RdYlBu Q3 tile
  heatmap with above-cutoff pairs outlined, adapted from
  `RASCHplot::ggQ3star()`; needs `data`, else `NULL`).
- [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md)
  is reimplemented on the CML engine, replacing the
  [`iarm::ICCplot()`](https://rdrr.io/pkg/iarm/man/ICCplot.html)
  wrapper. The conditional item curve and its variance now come from
  `psychotools` CML thresholds + the exact conditional distribution
  given the total score. It draws **confidence intervals** on the
  observed class-interval means (`ci`, default on) and an optional model
  band (`error_band`); in DIF mode it shows per-group CIs and annotates
  each panel with the **partial-gamma DIF magnitude**
  ([`iarm::partgam_DIF`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html),
  as in
  [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md)),
  with the full table attached as `attr(., "dif_gamma")`. New arguments
  `ci`, `error_band`, `conf_level`, `min_n` (per-cell floor, 8),
  `items`; `method` / `class_intervals` / `dif_var` / `output` are
  unchanged. The conditional-ICC approach follows Buchardt, Christensen
  & Jensen (2023) and their `RASCHplot` package.
- [`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md):
  corrected the `$direction2` caption / `@return` wording to “total -
  Item2” (with a note that direction 2 lists pairs in reverse order);
  computations unchanged (labelling fix).
- [`RMdifTree()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md):
  fixed a “variable lengths differ” error when a single covariate was
  passed by index or expression (e.g. `covariates = phq9[, 10]`). The
  derived non-syntactic column name is now matched as a literal column
  rather than re-evaluated as an R expression against the original
  (pre-NA-drop) data.
- `RMdifTree(stability = TRUE)` no longer silences the console in
  front-ends (e.g. RStudio) that redirect the message stream. The
  internal stderr capture used to muffle resample-fit noise reset the
  message sink to the default, clobbering the front-end’s sink; it now
  only redirects when no foreign message sink is active.
- `RMitemRestscore(p_adj = "none")` now returns the unadjusted p-values
  instead of `NA` with a “NAs introduced by coercion” warning. With no
  adjustment
  [`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html)
  omits the adjusted-p column, so the p-value is now selected by name
  rather than by a fixed column position (which had picked up the
  significance-stars column). The table header reads “p-value” in this
  case.

## easyRasch2 0.8.0

CRAN release: 2026-06-08

### Initial CRAN release!

### Visual consistency across all plot output

- Every `RM*` function that returns a `ggplot` now applies three shared
  internal theme helpers:
  - `er2_axis_margins()` — extra breathing room around the x and y axis
    titles.
  - `er2_plot_caption()` — left-aligned 9 pt plot caption. When the
    optional `ggtext` package is installed (new in Suggests), the
    caption renders via
    [`ggtext::element_markdown()`](https://wilkelab.org/ggtext/reference/element_markdown.html)
    so the APA-conventional italic “*Note.*” prefix can sit alongside
    roman body text. Without `ggtext`, it falls back to a plain
    `element_text(face = "italic")` and a plain “Note.” prefix — no
    markdown asterisks ever leak through to the rendered text.
  - `er2_caption()` — builds the caption string with the “Note.” prefix
    and wraps long captions at 90 characters via
    [`strwrap()`](https://rdrr.io/r/base/strwrap.html) so they no longer
    run off the right edge of the plot. The body text is wrapped first
    then prefixed, so the prefix is never broken by a line break.
- All three helpers are internal (`@noRd`) and not exported.

### New function

- **[`RMitemCatProb()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatProb.md)**
  plots model-implied category-response probability curves per item,
  similar to
  [`eRm::plotICC()`](https://rdrr.io/pkg/eRm/man/plotICC.html) or `mirt`
  trace plots but with a `ggplot2` / viridis output. Each item gets its
  own facet panel; one curve per response category is coloured from low
  to high using a continuous viridis palette. Polytomous items are fit
  with [`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html); dichotomous
  items fall back to [`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html)
  (recovering the standard two-category logistic ICC). Optional
  descriptive `item_labels` and `category_labels` make the output
  report-ready.
  - Optional **`label_curves = "path"`** mode (requires the optional
    `geomtextpath` package, in Suggests) writes each category’s label
    *along* its own curve in the classic IRT trace-plot style.
    Single-item only — pass the full multi-item dataset and use the
    `item` argument to choose which item’s curves to plot. The model is
    still fit on all items (CML threshold estimation needs the full
    dataset); only the rendering is filtered. Each category’s label is
    positioned at that curve’s modal theta (the peak of its bell for
    middle categories; the relevant edge for the monotone extreme
    categories), clamped a small margin in from the plot edges to
    prevent clipping.

### Function renaming for consistency

This release renames **22 exported functions** under a consistent
*domain-prefix → method → variant-suffix* scheme so that autocompletion
on a prefix (`RMdif`, `RMlocdep`, `RMitem`, `RMdim`, `RMplot`) surfaces
every related function. No semantic changes; only names.

A complete mapping is on the new `?easyRasch2-renaming` help page. The
short version:

| Domain | Prefix | Examples |
|----|----|----|
| DIF | `RMdif` | [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md), [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md), [`RMdifGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaPlot.md), [`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md), [`RMdifTree()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md) |
| Local dependence | `RMlocdep` | [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md), [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md), [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md), [`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md), [`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md), [`RMlocdepGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaPlot.md) |
| Item statistics | `RMitem` | [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md), [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md), [`RMitemInfitCutoffPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md), [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md), [`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md), [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md), [`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreBoot.md), [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md), [`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemHierarchy.md) |
| Dimensionality | `RMdim` | [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md), [`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCACutoff.md), [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md), [`RMdimCFAPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFAPlot.md), [`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md), [`RMdimMartinLofResiduals()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLofResiduals.md) |
| Descriptive plot | `RMplot` | [`RMplotTile()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotTile.md), [`RMplotBar()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotBar.md), [`RMplotStackedbar()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotStackedbar.md) |

Unchanged:
[`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md),
[`RMUreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMUreliability.md),
[`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md),
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md).

**Breaking change.** No deprecation aliases are shipped — existing
scripts will need a search-and-replace. See `?easyRasch2-renaming` for
the full table.

## easyRasch2 0.7.1

- New **`RMlocdepQ3plot()`** for visualising per-pair simulated Q3
  distributions against observed Q3 values, mirroring the design of
  `RMpgLDplot()` and `RMpgDIFplot()`. Supports `items` and `n_pairs`
  filters; with `data` supplied, ranks `n_pairs` by
  `|observed Q3 - median(simulated Q3 per pair)|`.
- `RMlocdepQ3cutoff()` now retains per-pair iteration-level data: the
  returned list gains `pair_results`, `pair_cutoffs`, `item_names`,
  `cutoff_method`, and `hdci_width`. Existing scalar outputs
  (`suggested_cutoff` etc.) are unchanged. New arguments `cutoff_method`
  (`"hdci"` / `"quantile"`) and `hdci_width` (default `0.99`) control
  how per-pair credible intervals are computed.
- [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  can now be called with the full list returned by `RMlocdepQ3cutoff()`
  (it auto-extracts `$suggested_cutoff`), matching the partial-gamma
  family convention. Numeric scalar still works.
- `RMpartgamLD()` and `RMpgLDplot()` gain an `n_pairs` argument to keep
  only the top-N pairs by `|gamma|` (LD) or by
  `|observed - median(sim)|` deviation (plot).
- `RMpartgamLD()`, `RMpartgamDIF()`, and `RMitemrestscore()` now use
  consistent column headers: `Adj. p-value (BH)` and a new
  `p-value sign.` star-string column. `RMitemrestscore()`’s former
  `Absolute_difference` column is now a *signed* `Difference` (observed
  − expected), so over- and underfit are visually distinguishable at a
  glance.

## easyRasch2 0.7.0

- New **`RMbarplot()`** and **`RMstackedbarplot()`** join `RMtileplot()`
  for easy visualization of item response data distributions.
- New **`RMciccPlot()`** for conditional item characteristic curves
  plot, also includes DIF analysis for categorical DIF variables.
- New
  **[`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemHierarchy.md)**,
  which outputs a plot illustrating the item hierarchy with item
  thresholds and confidence intervals.
- Modified `RMcfaCutoff()` to only use .scaled metrics for RMSEA and
  CFI, for stability. See documentation for more details.
- Bug fix for `RMmartinLof()` with dichotomous items.
- Bug fix for
  [`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemHierarchy.md)
  making option `item_labels` work as intended.
- Tests written for all functions, preparing for CRAN submission in the
  nearish future.
- Drafted a package intro/vignette

## easyRasch2 0.6.0

- New **`RMcfaCutoff()`** for testing unidimensionality against
  simulation-based cutoff values for model fit metrics.
  - **`RMcfaPlot()`** shows figure with distribution of simulation
    results and observed model fit values.

## easyRasch2 0.5.8

- New
  **[`RMdifTree()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md)**
  for testing DIF with continuous variables, such as age in years, as
  well as combinations of DIF variables (DIF interactions).
  - Also implements use of `stablelearner`, as recommended by Henninger
    et al.
    2025. to assess stability of tree-based results.
  - Borrows from
    - raschtreeMH (<https://github.com/mirka-henninger/raschtreeMH>)
    - effecttree (<https://github.com/mirka-henninger/effecttree>)
  - See these two papers (and
    [`?RMdifTree`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md))
    for more details:
    - Henninger, M., Debelak, R., & Strobl, C. (2023). A new stopping
      criterion for Rasch trees based on the Mantel-Haenszel effect size
      measure for DIF. Educational Psychological Measurement, 83,
      181-212. <doi:10.1177/00131644221077135>
    - Henninger, M., Radek, J., Debelak, R., & Strobl, C. (2025).
      Partial credit trees meet the partial gamma coefficient for
      quantifying DIF and DSF in polytomous items. Behaviormetrika, 52,
      221-257.

## easyRasch2 0.5.7

- New
  **[`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md)**
  for testing DIF of categorical variables using Andersens’s Likelihood
  Ratio test as implemented in package `eRm`
- Visualization of response data using a tile plot with
  **`RMtileplot()`**.
  - Optional `group` faceting, which is useful for DIF-analyses.
- Added Martin-Löf test of dimensionality **`RMmartinLof()`** based on
  Christensen & Kreiner (2007, doi: 10.1177/0146621605286204), for both
  dichotomous and polytomous (PCM) items.
  - Also added post-processing function **`RMmartinLofresiduals()`**.

## easyRasch2 0.5.6

- Added **`RMresidualPCA()`** for evaluating patterns in the
  standardized residuals from the RM/PCM. Outputs either a table with
  eigenvalues and explained variance or a figure with standardized
  loadings on the first residual contrast and item locations.
- Added **`RMpcaCutoff()`** to determine a simulation-based critical
  value for the largest eigenvalue.

## easyRasch2 0.5.5

- Added **`RMbootRestscore()`** for use with large sample sizes. See
  <https://pgmj.github.io/rasch_itemfit/> for more details.
- Added
  **[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md)**
  that produces a transformation table (or figure) from ordinal sum
  scores to WLE (Weighted Likelihood Estimation) interval scores.
- Modified `RMlocdepQ3cutoff()` to improve speed.

## easyRasch2 0.5.4

- **[`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md)**
  for Wright map style plots.
- New functions for handling missing data with conditional item infit
  MSQ by multiple imputation using package `mice`.
  - **`RMinfitcutoff_mi()`** uses the `mids` object containing multiple
    imputated datasets (output by
    [`mice::mice()`](https://amices.org/mice/reference/mice.html)) to
    run simulations on each datasets and combines the results. Splits
    iterations across imputations (e.g. 250 total / 5 imputations = 50
    each).
  - **`RMiteminfit_mi()`** calculates and pools conditional infit from
    the imputated datasets and optionally uses the `RMinfitcutoff_mi()`
    output for cutoff values.
  - Works with existing **`RMinfitcutoffPlot()`**

## easyRasch2 0.5.3

- New functions **`RMpartgamLD()`** and **`RMpgLDcutoff()`** for
  evaluating local dependency of item pairs in both directions.
  - Plot function **`RMpgLDplot()`**, similar to **`RMpgDIFplot()`**
  - NOTE: The simulation-based cutoffs have not yet been evaluated in a
    systematic way. Default choices may not be sensible, please try
    different numbers of iterations and hdci_width.

## easyRasch2 0.5.2

- New functions **`RMpartgamDIF()`** and **`RMpgDIFcutoff()`** for
  evaluating DIF of categorical external variables.
  - Also with plot function **`RMpgDIFplot()`**, similar to
    `RMinfitcutoffPlot()`
  - NOTE: The simulation-based cutoffs have not yet been evaluated in a
    systematic way. Default choices may not be sensible, please try
    different numbers of iterations and hdci_width.

## easyRasch2 0.5.1

- New function **`RMinfitcutoffPlot()`** to illustrate distribution of
  simulated conditional item infit MSQ values together with the observed
  value.

## easyRasch2 0.5.0

- **`RMinfitcutoff()` gains `cutoff_method` and `hdci_width`
  parameters**: By default (`cutoff_method = "hdci"`,
  `hdci_width = 0.999`), per-item cutoff intervals are now computed
  using the Highest Density Interval via
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html)
  (99.9% HDCI). Set `cutoff_method = "quantile"` to restore the previous
  behaviour (2.5th/97.5th percentiles). The `ggdist` package is only
  required when `cutoff_method = "hdci"` (added to Suggests). The
  returned list now also includes `cutoff_method` and `hdci_width`
  fields.

- **`RMiteminfit()` caption updated**: When `cutoff` is the return value
  of `RMinfitcutoff()`, the kable caption now states the cutoff method,
  e.g.
  `"Cutoff values based on 250 simulation iterations (99.9% HDCI)."` or
  `"Cutoff values based on 250 simulation iterations (2.5th/97.5th percentile)."`.

## easyRasch2 0.4.0

- **`RMiteminfit()` gains optional `cutoff` parameter**: Accepts the
  return value of `RMinfitcutoff()` (or its `$item_cutoffs` data.frame
  directly). When provided, per-item cutoff boundaries (`Infit_low`,
  `Infit_high`) and a logical `Flagged` column are added to both
  `"dataframe"` and `"kable"` output. The kable caption includes the
  number of simulation iterations when available.

- **New `RMinfitcutoff()`**: Simulation-based (parametric bootstrap)
  cutoff determination for \[RMiteminfit()\]. Supports both dichotomous
  and polytomous data. Optional parallel processing via `mirai` (falls
  back to sequential if not installed). Returns per-item 2.5th and
  97.5th percentiles of simulated infit and outfit MSQ distributions
  (`$item_cutoffs`), together with the full iteration-level results
  (`$results`). Requires the `iarm` package (Suggests).

## easyRasch2 0.3.0

- **New `RMiteminfit()`**: Computes conditional infit MSQ statistics for
  each item via
  [`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html),
  enriched with item locations relative to the sample mean person
  location. Supports both dichotomous (Rasch model via
  [`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html)) and polytomous
  (Partial Credit Model via
  [`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html)) data. Requires
  the `iarm` package (Suggests). Output options: `"kable"` (default,
  plain-text
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html)) or
  `"dataframe"`. Optional `sort = "infit"` sorts by infit MSQ
  descending. Only complete cases are used for conditional fit
  calculation.

- **New `RMitemrestscore()`**: Computes observed and model-expected
  item-restscore correlations via
  [`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html),
  enriched with absolute differences between observed and expected
  values, item average locations, and item locations relative to the
  sample mean person location. Supports both dichotomous (Rasch model
  via [`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html)) and polytomous
  (Partial Credit Model via
  [`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html)) data. Requires
  the `iarm` package (Suggests). Output options: `"kable"` (default,
  plain-text
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html)) or
  `"dataframe"`. Optional `sort = "diff"` sorts by absolute difference
  descending.

## easyRasch2 0.2.0

- **[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)**:
  `cutoff` is now optional (default `NULL`). When omitted, the raw Q3
  residual correlation matrix is returned without any dynamic cut-off
  applied. When provided, the dynamic cut-off (mean Q3 + cutoff) is
  shown in the kable caption as before. Supersedes the requirement to
  always supply a cutoff value (PRs
  [\#2](https://github.com/pgmj/easyRasch2/issues/2) and
  [\#3](https://github.com/pgmj/easyRasch2/issues/3)).
- **New `RMlocdepQ3cutoff()`**: Simulation-based (parametric bootstrap)
  cutoff determination for
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md).
  Supports both dichotomous and polytomous data. Optional parallel
  processing via `mirai` (falls back to sequential if not installed).
  The `$suggested_cutoff` is the 99th percentile of the simulated Q3
  max–mean distribution.
- **New internal helpers** in `R/utils-simulation.R`: `sim_poly_item()`,
  `sim_partial_score()`, and `extract_item_thresholds()`.
- **Dependencies**: `eRm`, `psychotools (>= 0.7-3)`, `parallel`, and
  `utils` moved/added to Imports. `mirai` added to Suggests.

## easyRasch2 0.1.0

- Initial package structure.
- Added
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  for Yen’s Q3 residual correlation analysis of local dependence.
