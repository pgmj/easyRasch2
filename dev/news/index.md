# Changelog

## easyRasch2 (development version)

- [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  now computes Yen’s Q3 by **conditional maximum likelihood with Warm’s
  weighted likelihood** person locations by default
  (`estimator = "CML"`), replacing the previous reliance on marginal ML
  / EAP via `mirt`. CML/WLE is true to the Rasch tradition, gives finite
  estimates at extreme scores, and (via `psychotools`) is markedly
  faster in the parametric bootstrap. The previous engine remains
  available with `estimator = "MML"`. **This changes the Q3 values,
  cut-offs, and flags relative to earlier versions** — the statistic is
  the same but the estimator differs (observed-vs-MML off-diagonal
  correlation is typically above 0.95 on real data).
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  gains a matching `estimator` argument and stores it in its return
  value (`$estimator`);
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  and
  [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md)
  read that stored estimator so the observed Q3 and the simulated
  cut-off are always computed the same way (a mismatching `estimator`
  argument is overridden with a warning). The shared CML/WLE/EAP
  person-estimation engine is reusable internally for migrating other
  functions off `mirt`/`eRm` over time.
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)’s
  parametric-bootstrap *generating* distribution now also uses CML item
  parameters and WLE (not MLE) person locations, so the whole pipeline
  is CML/WLE-coherent; this slightly shifts the simulated cut-offs
  relative to earlier development versions.

- [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  given the full
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  object now returns a third list element, `$plot`: a tile heatmap
  (lower triangle, diverging RdYlBu fill) of the observed Q3 matrix,
  with pairs above the global dynamic cut-off outlined. Adapted from the
  `RASCHplot` package’s `ggQ3star()` and styled to match the package’s
  other ggplots. (The returned list is now `$matrix`, `$pairs`, `$plot`;
  the simple `cutoff = NULL` / numeric-cutoff return values are
  unchanged.)

- New
  [`RMpersonFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonFit.md):
  per-respondent person-fit analysis. Reports conditional infit and
  outfit mean-squares (computed on each response pattern via elementary
  symmetric functions, so partial missingness is handled and no biased
  person estimate enters the residual) plus the standardized
  log-likelihood `lz`. Significance is assessed by Monte-Carlo
  resampling under the fitted model — each statistic against the null
  that matches it (MSQ conditional on the total score, lz simulated at
  the estimated location) — rather than by an unreliable asymptotic null
  (Sinharay, 2016; Müller, 2020). The Wilson-Hilferty (ZSTD) transform
  of MSQ is available via `zstd = TRUE` but, following Müller (2020), is
  provided only as a comparability metric and not used for inference.
  Output is a table, a data.frame, or (`output = "ggplot"`) a named list
  of person-fit maps – one per statistic – each captioned with the
  assessed sample size and the proportion of respondents flagged by that
  statistic (split into underfit vs overfit for the MSQ maps), with
  points coloured by fit direction. `flag = "underfit"` restricts
  flagging to the validity-relevant underfit direction (MSQ \> 1) using
  a more powerful one-sided p-value.

- New
  [`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md):
  returns item difficulty (dichotomous) or Andrich-threshold
  (polytomous) parameters in long or wide format, with optional standard
  errors and Wald confidence intervals. Item parameters are estimated by
  CML (via `eRm`) by default, or by MML (via `mirt`) for sparse data;
  polytomous threshold SEs use the delta method.

- New
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md):
  estimates a person location (theta) and its SEM for every respondent,
  working directly on each response pattern so partial missingness is
  handled correctly (unlike the sum-score lookup in
  [`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md)).
  Supports Warm’s WLE (default) and EAP under a normal prior whose SD is
  estimated from the data by marginal maximum likelihood unless fixed
  via `prior_sd`.

- The `Flagged` column now labels the misfit *direction* instead of a
  logical TRUE/FALSE: `"overfit"`, `"underfit"`, or `""` (no misfit), to
  simplify interpretation. In
  [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
  and
  [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)
  infit below the range is `"overfit"` (more predictable than expected)
  and above is `"underfit"` (noisier);
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
  gains a `Flagged` column where over-discrimination (observed \>
  expected, adj. p \< .05, often local dependence) is `"overfit"` and
  under-discrimination is `"underfit"`. Note the value direction is
  opposite between the two (a *high* infit is underfit, a *high*
  restscore correlation is overfit). This changes the type of the
  [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)/[`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)
  `Flagged` column from logical to character.
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
  also drops the now-redundant `Significance` (“p-value sign.”) star
  column – the numeric `p_adjusted` and the new `Flagged` label convey
  the same information more directly – and the `Location` column (the
  relative location is sufficient), and trims the table headers
  (“Observed”, “Expected”) for a narrower table with an explanatory
  caption.

- [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
  gains an optional `p_value` argument. When `TRUE` (and the full
  [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
  object is supplied), it reports per-item bootstrap p-values from the
  simulated null with multiple-comparison correction: Westfall-Young
  studentised-max step-down for the family-wise error rate
  (`correction = "fwer"`, default), or Benjamini-Hochberg /
  Benjamini-Yekutieli FDR. Items are studentised by the bootstrap SD
  (not the ZSTD transform), and p-values are shown alongside the
  simulated effect-size band. At least 1000 cutoff `iterations` are
  recommended (a warning is issued otherwise); default behaviour with
  `p_value = FALSE` is unchanged. Methods follow Ferreira (2024) and
  Westfall & Young (1993).

- [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  now returns two complementary views when the full
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  object is supplied: a `list` with `$matrix` (the Q3
  residual-correlation matrix, flagged against the single global cutoff)
  and `$pairs` (a per-pair table giving each pair’s observed Q3 with its
  simulated low/high reference band, `Flagged` as
  `"above"`/`"below"`/`""`, sorted by departure from the per-pair median
  and optionally truncated to the most extreme `n_pairs`). This matches
  the pairwise Q3 table in the `easyRasch2jmv` jamovi module. With a
  bare numeric `cutoff` (or none) the single matrix is returned as
  before. The optional `p_value` layer now folds one-sided bootstrap
  p-values for excess local dependence (with the same Westfall-Young
  FWER / FDR `correction` options across the *k(k-1)/2* pairs) into
  `$pairs`, which then flags on the adjusted p-value.
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  and
  [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
  share a “Multiple comparisons” help section explaining when to use
  marginal vs corrected p-values.

- [`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md):
  the WLE `method` now reports the information-based standard error
  `1 / sqrt(I(theta))` (matching
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md),
  `catR` and `TAM`) instead of the previous iarm “expected SEM”. Point
  estimates are unchanged; standard errors differ at the score extremes
  (where they are now larger). The WLE path now shares the same solver
  and grand-mean-zero threshold centering as
  [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md),
  so the two functions are fully consistent.

- [`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md):
  corrected the `$direction2` kable caption and the `@return`
  documentation, which described that table’s rest score as “total -
  Item1”.
  [`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html)
  always removes the second item (Item 2) from the rest score in *both*
  directions; the two tables differ only in the order the items of each
  pair are listed (so each pair is tested with each of its items
  removed). The caption now reads “total - Item2” with a note that
  direction 2 lists pairs in the reverse order. Computations are
  unchanged – this is a labelling fix.

- [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
  now computes the observed Q3 matrix with `mirt::residuals(digits = 4)`
  (was `digits = 2`), matching the precision already used for the
  simulated Q3 values in
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  and
  [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md).
  The mean Q3, the dynamic cut-off, and the `above_cutoff` flags are now
  computed from full-precision values (previously the flags were derived
  from values rounded to 2 decimals), and `output = "dataframe"` returns
  4-decimal precision. The kable display remains rounded to 2 decimals.
  Borderline pairs close to the cut-off may flag differently than before
  – the new behaviour is the more accurate one. Aligns with the same fix
  in the easyRasch2jmv jamovi module (1.1.0).

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
| Item statistics | `RMitem` | [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md), [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md), [`RMitemInfitCutoffPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffPlot.md), [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md), [`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md), [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md), [`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreBoot.md), [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md), [`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemHierarchy.md) |
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
