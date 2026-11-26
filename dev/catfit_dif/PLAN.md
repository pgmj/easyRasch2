# Category fit by DIF group: dev plan

Started 2026-11-24. Status: plan only, no code.

## Purpose

RMitemCatFitPlot() shows how an item's response categories work given the
total score. Split by a DIF variable, the same plots should show how the
categories work in each group: a group that avoids a middle category, uses
a step differently, or differs at one end of the score range. RMdifAnchor()
tests whether threshold locations differ. These plots are meant to show the
pattern behind such a difference, and they also apply where RMdifAnchor()
has no answer (three or more groups, DIF that changes across the score
range).

The first version is descriptive: no test and no bands. A test is a later
step (see Roadmap).

## Background and prior art (checked 2026-11-24)

Full texts in `code/shared/articles/`.

- Kreiner and Christensen (2011, Psychometrika, sec. 5.1): no DIF means the
  item is independent of the exogenous variable given the total score.
  Following Holland and Thayer (1988), this holds exactly under the Rasch
  model, not under 2PL/3PL models. That property is what justifies pooled
  predictions given the total score as the no-DIF reference in every group.
- Kreiner (1987, Scand J Stat 14:97-112): exact conditional tests of zero
  partial association, with Monte Carlo p-values from random tables given
  the sufficient margins, slice by slice. A permutation of group labels
  within total scores is this null. Cite it for the later test (step 3).
- Glas (1998, Statistica Sinica 8:647-667): LM tests for DIF with a shift
  per category (observed vs expected focal-group counts in category h) and
  a theta-weighted term for non-uniform DIF. MML-based, not conditional on
  the total score.
- DIF at single steps (differential step functioning, DSF): Penfield NCME
  module, Rasch/PCM trees (Komboz, Zeileis and Strobl), Schoenmakers et al.
  (2024). All parametric.
- DIF that changes across the score range: Hagquist and Andrich (2017)
  residual ANOVA (group by class interval), crossing SIBTEST (Li and Stout,
  1996), smoothed SIBTEST (Douglas, Stout and DiBello, 1996).
- No PAV/CORP calibration used for DIF was found.

## Step 1: group-split RMitemCatFitPlot (descriptive)

### Design decisions (agreed 2026-11-24)

1. **New argument `dif_var`**, as in RMitemICCPlot() and the DIF functions.
   Categorical, two or more groups.
2. **Types: `residuals`, `adjacent`, `calibration`.** `restscore` stops with
   an error in DIF mode for now (see Open questions).
3. **No bands or reference ranges in DIF mode.** Each group's n is shown, and
   the caption says the plot is descriptive and points to RMdifAnchor() /
   RMdifGamma() for tests.
4. **Predictions from one pooled CML fit**, on all respondents used. Groups
   only split the observed data. Reading: group curves that separate = DIF,
   all curves off the reference together = misfit of the item for everyone.
   Separate fits per group are not used, since they would hide that
   distinction.
5. **Groups by colour + line type + shape**, Okabe-Ito, theme_bw(), items as
   facet rows as now. Check the figure in greyscale and with the colour
   validator.

### What each type shows per group

- `calibration`: one PAV curve per group and event \eqn{X_i \ge k},
  observed against the pooled prediction. Points sized or marked by n as in
  the current plot.
- `adjacent`: share in the upper of two adjacent categories, against the
  total score, per group. The current binomial range around the expected
  value is dropped in DIF mode (decision 3).
- `residuals`: observed minus expected counts per category within score
  groups, computed per DIF group from the pooled conditional probabilities.
  Close to Glas's item-by-score table split by group. The score grouping
  (`grouping`, `class_intervals`, `score_breaks`) is shared by all DIF groups,
  so the x positions match.

### Implementation notes (from reading the code 2026-11-24)

- `.catfit_prepare()` drops incomplete rows and extreme total scores but does
  not return which rows it kept. DIF mode needs the kept-row index (or
  filter `data` and `dif_var` jointly first, as RMitemICCPlot() does, then
  map the extreme-score drop). Add a `rows` element to the returned list
  rather than duplicating the filtering.
- Rows with missing `dif_var` are dropped with a message, matching
  RMitemICCPlot() wording. The caption follows the house rule
  `n = X of Y respondents (policy)`, plus n per group.
- `calibration` always runs simulations now (`S` for the band). DIF mode
  must skip `.catfit_simulate()` entirely, so `iterations`, `seed`,
  `parallel` are unused there (say so in the help).
- `fit` with `dif_var`: the `identical(prep$X, fit$prep$X)` check fails when
  `dif_var` has NAs, since the pooled fit is then on fewer rows. Simplest:
  `fit` is not allowed with `dif_var` in step 1 (error). Revisit later.
- `output = "dataframe"` gains a `group` column. Nothing in attributes.
- Helpers `.catplot_residuals()`, `.catplot_adjacent()`,
  `.catplot_calibration()` take X, r, P. Calling them per group on the
  group's rows, with the pooled P, gives the per-group tables, then bind with
  a `group` column. The `_gg` functions need a group aesthetic.
- Sparse groups: PAV and adjacent shares at single scores get noisy. Decide
  whether to drop points below a minimum n per score (RMitemICCPlot() has
  `min_n`). Keep consistent with that function if so.

### Checks before calling step 1 done

- No DIF, equal and unequal group means (impact): group curves overlap up to
  noise in all three types.
- Uniform DIF on one item, DIF at one step, middle-category avoidance in one
  group: the plots show the pattern where expected.
- Unbalanced DIF (several items, same direction): see how DIF-free items look,
  since the total score is contaminated. Document it in the help.
- Container smoke only (small n, few items). Any study-size run is scripted
  for the host and asked about first.

### Docs and housekeeping for step 1

Tests in tests/testthat/test-item_cat_fit.R (or a new file), help text
(param, details on reading and on the total-score caveat), NEWS.md, check
_pkgdown.yml and README, vignette mention if wanted. Log design changes in
this file as work proceeds. jamovi is a separate decision.

## Roadmap after step 1

2. **Dev check of the reference under impact.** No-DIF data with different
   group means and several n, k: quantify how far group curves wander. This
   tells how noisy the plot is without bands.
3. **Permutation bands, then a test.** Permute group labels within total
   scores (Kreiner 1987 null), no refit needed. Bands per group first. Then a
   statistic (for example group-weighted MCB, or a between-group curve
   distance) with max-T over steps and Westfall-Young over items, as in
   RMitemCatFit().
4. **Comparison study.** Item-level detection against partial gamma (same
   conditioning). Step localisation against RMdifAnchor() at threshold and
   spacing level, two groups. Conditions: balanced and unbalanced DIF,
   single-step DIF, middle-category avoidance, crossing DIF, three groups.
   Check whether a flag on cumulative event k points to threshold k, since
   a shifted threshold changes more than one event given the total score.
5. **rasch_book chapter**, with scripts and saved results.

## Validating the permutation bands (agreed 2026-11-25)

A respondent's predicted probability depends only on the total score and the
pooled item parameters. A permutation of group labels within total scores
keeps every group's count at every score. So each group's PAV curve is
evaluated at the same x values in every permutation, and no refit is needed.
Under the Rasch model without DIF, the item is independent of the group given
the total score (Holland and Thayer 1988, via Kreiner and Christensen 2011),
whatever the group ability distributions. The band is then exact for each
point on its own. The studies are therefore about the code, about reading a
whole curve, and about where the null fails.

Band definition. Pointwise: a point is outside when its two-sided permutation
p, `2 * min(#{perm <= obs} + 1, #{perm >= obs} + 1) / (B + 1)`, is at most
alpha. The drawn band is the range of values that would not be outside.
Simultaneous (per curve, that is item x event x group): standardise each
point by its permutation mean and SD, take the maximum |z| over the curve's
points, and compare with the permutation distribution of that maximum.

- **A. Code checks (deterministic).** Group counts and predicted values per
  score unchanged by every permutation. A DIF variable constant within each
  score gives a zero-width band. Same seed, same band. With two groups,
  deviations from the pooled curve mirror each other (weighted by n).
  The fast aggregated PAV agrees with `.pav()`.
- **B. Pointwise or simultaneous (Rasch model, no DIF).** Pointwise
  exceedance overall, by event and by n per point. Share of curves with any
  point outside the pointwise band. Simultaneous exceedance per curve.
  Conditions: impact 0 / 0.5 / 1 logit, groups 2 equal / 2 at 80:20 / 3
  equal, n 250 / 500 / 1000, k 5 / 10, four categories. B = 999.
  Decides which band the plot draws.
- **C. Specificity: no DIF, Rasch model false.** Slope spread x impact x n.
  Curves outside the band without DIF, the Holland and Thayer caveat.
- **D. Unbalanced DIF.** Several items with DIF in one direction, rates for
  the DIF-free items.
- **E. Known patterns.** Uniform DIF, single-step DIF, middle-category
  avoidance, crossing DIF: exceedance at the right event, and spill-over to
  neighbouring events.

A is done in the container, then a B smoke run (two cells, 200 datasets,
4 cores). The full B grid is a host script. C to E wait until B has decided
the band type.

Files: `R/perm_functions.R` (prototype, uses package internals through
`devtools::load_all()`), `R/check_A.R`, `R/run_B.R`, `R/summarise_B.R`.
Results: `results/permband_B_cells/` (per-cell rds, gitignored by the
`*cells*` pattern), `results/permband_B_summary.rds` (tracked), smoke files
`*smoke*.rds` (gitignored), logs tracked.

### Progress log

**2026-11-25, A done** (`R/check_A.R`, log `results/check_A.log`). All checks
pass. Two things found on the way:
- Predictions are only non-decreasing in the score. At extreme scores some
  events are certain (P = 0 or 1), so neighbouring scores tie there. The
  check confirms ties occur only at P = 0 or 1 and that responses agree with
  them, so keeping tied scores as separate points cannot change the PAV
  values. These points cannot vary under permutation and are left out of all
  rates in B.
- The mirror check applies to raw per-score shares
  (n1 (y1 - ybar) = -n2 (y2 - ybar)), not to PAV curves, since PAV is not
  linear. The A list above was worded too broadly.
- Bug fixed before any results: `ifelse()` dropped the matrix shape in the
  simultaneous z, so the simultaneous p was wrong in the first draft.

Permutations reproduce the hypergeometric mean and variance of a group's
count at one score (20,000 permutations). The aggregated PAV equals
`.pav()` on individual data. Timing: 0.3 s (n 250, k 5, B 999) to 1.8 s
(n 1000, k 10, three groups) per dataset.

**2026-11-25, B smoke** (2 cells x 200 datasets, B = 999, 4 cores, 5 min;
`results/run_B_smoke.log`, `results/summarise_B_smoke.log`). All 200
datasets distinct per cell, no failed fits.

| cell | alpha | pointwise | any point outside pointwise band | simultaneous |
|---|---|---|---|---|
| impact 1, 80:20, n 250, k 5 (11 points/curve) | .05 | .011 | .088 | .036 |
| | .10 | .024 | .183 | .079 |
| impact 0.5, three groups, n 500, k 10 (26 points/curve) | .05 | .020 | .275 | .045 |
| | .10 | .042 | .496 | .093 |

Reading, with the caveat that this is two cells: the pointwise p is valid
but well below nominal (discrete PAV values with ties). Reading a curve
against the pointwise band flags 9% to 28% of curves without DIF, rising
with the number of points. The simultaneous band per curve is at or below
nominal in both cells. The naive SEs (.0024, .0015) treat curves as
independent. Curves within a dataset are correlated (with two groups they
mirror), so the real SEs are larger. Even a simultaneous band per curve
leaves multiplicity across the k x events x groups curves in a plot. That
belongs to the test step (max-T across items).

**2026-11-26, full B grid** (host, 54 cells x 1000 datasets, B = 999, 34 min;
`results/run_B.log`, `results/permband_B_summary.rds`). Every cell has 1000
distinct datasets and no failed fits.

At alpha = .05, range over the 54 cells:
- pointwise: .010 to .024 (valid, far below nominal),
- any point of a curve outside the pointwise band: .085 to .370,
- simultaneous per curve: .034 to .049, never above .05.

At alpha = .10 the simultaneous rate is .074 to .101.

Simultaneous rate by factor (mean over cells, alpha = .05): impact 0 / 0.5 /
1 gives .044 / .044 / .043. Groups two equal / 80:20 / three equal gives .044 /
.041 / .045. n 250 / 500 / 1000 gives .040 / .044 / .046. k 5 / 10 gives .042 / .044.
The lowest cells are n 250 with an 80:20 split (.034 to .036), where the
small group's curves are the most discrete.

Any-point rate by n and k (rows n 250, 500, 1000): k 5 gives .11, .16, .20,
k 10 gives .16, .25, .33. It rises with n because the pointwise p becomes less
discrete, and with k because curves have more points.

Conclusions:
1. Impact has no effect, as the theory says (the null holds whatever the
   group distributions under the Rasch model).
2. A pointwise band would mislead: up to a third of DIF-free curves cross it.
3. The simultaneous band per curve holds its level in every cell, at most 1.6
   points below nominal in the smallest unbalanced cells. Drawn as
   permutation mean +- c x permutation SD per point, with c the 1 - alpha
   quantile of the curve's max |z|, a curve crosses the band exactly when
   its simultaneous p is at most alpha (up to ties).
4. Across all curves in a plot there is still multiplicity. That is for the
   test step.

Proposed: the plot draws the simultaneous band (user to decide). Next: C
(specificity, Rasch model false), to be agreed before running.

**2026-11-26, decisions (user).**
- The plot draws the simultaneous per-curve permutation band. This replaces
  decision 3 of step 1 ("no bands").
- Studies C to E are paused. Focus is the visual diagnostic plot, without a
  formal test. Help text states the slope caveat (Kreiner and Christensen
  2011, Holland and Thayer 1988) and the unbalanced-DIF caveat qualitatively,
  and the caption says the band is per curve and not corrected across curves.
- The band must also be checked for `adjacent` and `residuals` before all
  three types get it (option 2). Container smoke, design below.

**Band smoke for adjacent and residuals (design).** A "curve" is what the
plot draws as one series for one group:
- `adjacent`: item x adjacent pair x group, one point per total score.
  The value is the share in the upper category among the group's
  respondents at that score who are in either category. The number in the
  pair varies across permutations, so a point can be empty in a permutation.
  Empty points get z = 0 there, and mean and SD come from the permutations
  where the point exists. Points empty in the observed data are dropped.
- `residuals`: item x group, one point per score group x category (the whole
  row of panels for that item). The value is the standardised residual
  (O - E) / sqrt(V) with E and V from the pooled conditional probabilities and
  the group's counts per score. E and V are fixed under permutation, so only
  O moves.
Same cells as the B smoke (impact 1, 80:20, n 250, k 5 and impact 0.5, three
groups, n 500, k 10), 200 datasets, B = 999. The calibration band is rerun
through the same generic code as a regression check against B.

**2026-11-26, A2 code checks** (`R/check_A2.R`, `results/check_A2.log`), all
pass. The generic code reproduces `perm_bands()` for calibration. Adjacent
shares per group equal `.catplot_adjacent()` on the group's rows with the
pooled P. One-group residuals equal `.catplot_residuals()`. Group E values sum
to the pooled E. A label constant within scores gives fixed points.

Bug found by A2 and fixed in both band functions: the SD as mean of squares
minus squared mean left rounding noise (about 1e-8) where every permutation
gives the same value, so fixed points were reported as varying. Now a
centred SD. Effect on B: negligible (z of order 1e-7 at such points against
critical values near 3). The calibration rerun below reproduces the B smoke
curve rates exactly. check_A.R still passes.

**2026-11-26, smoke for all three types** (`R/run_smoke_types.R`,
`results/run_smoke_types.log`, same datasets as the B smoke, 200 per cell, all
distinct, no failures). Curve rates at alpha .05 (and .10):

| type | cell | points/curve | any point outside pointwise | simultaneous |
|---|---|---|---|---|
| calibration | 80:20, n 250, k 5 | 7.5 | .088 | .036 (.079) |
| calibration | three groups, n 500, k 10 | 18.0 | .275 | .045 (.093) |
| adjacent | 80:20, n 250, k 5 | 6.6 | .058 | .033 (.078) |
| adjacent | three groups, n 500, k 10 | 15.4 | .183 | .044 (.090) |
| residuals | 80:20, n 250, k 5 | 16.0 | .166 | .039 (.096) |
| residuals | three groups, n 500, k 10 | 17.9 | .285 | .050 (.097) |

Points per curve here count only points that vary under permutation, so the
calibration numbers are lower than in the B table (which left out only
certain events). The calibration curve rates equal the B smoke.

Reading: the simultaneous band holds for all three types. Residuals are
closest to nominal (counts over score groups are less discrete). The
.050 cell has 6000 curves (naive SE .003, real SE larger), so it is
consistent with nominal. Any grouping of the residual plot that is a
function of the total score (all `grouping` options) is unchanged by the
permutation, so the argument holds for every grouping, not only the five
quantile groups tested.

**Terminology (2026-11-26).** Pointwise vs simultaneous is standard. For PAV
calibration plots the established term is consistency band (Dimitriadis,
Gneiting and Jordan 2021): the range under the hypothesis, not a confidence
band for the true curve. Säilynoja et al. (2025, articles/) use "95%
point-wise consistency bands" for PAV calibration plots and "simultaneous
confidence bands" for PIT-ECDF plots. Ours: a simultaneous consistency band
per curve. The max-|z| construction is a standard max-T idea, not a
published name. Help text: "simultaneous 95% consistency band (Dimitriadis et
al., 2021), per curve, from permutations of group labels within total scores
(Kreiner, 1987)". Caption in plain words. Never "confidence band".

**2026-11-26, full grid for adjacent and residuals** started in the container
at the user's request (overnight, 4 cores, `R/run_types_full.R`, log
`results/run_types_full.log`, summary by `R/summarise_types_full.R`). This
departs from the container rule in CLAUDE.md (smoke only), by explicit user
request.

**2026-11-26, full grid for adjacent and residuals done** (container, 1.4 h,
54 cells x 1000 datasets, same datasets as B, no failures, all distinct;
`results/run_types_full.log`, `results/permband_types_summary.rds`).
Simultaneous rate per curve, range over cells (mean):

| type | alpha .05 | alpha .10 | any point outside pointwise, .05 |
|---|---|---|---|
| calibration (B grid) | .034 to .049 | .074 to .101 | .085 to .370 |
| adjacent | .023 to .052 (.040) | .065 to .102 | .033 to .381 |
| residuals | .033 to .056 (.047) | .072 to .111 | .145 to .363 |

Mean simultaneous rate at .05 by group split, adjacent / residuals: two
equal .043 / .050, 80:20 .037 / .043, three equal .041 / .048. Impact has no
effect (adjacent .041 / .040 / .039, residuals .048 / .047 / .046 for impact
0 / 0.5 / 1). Lowest cells: n 250 with an 80:20 split (adjacent .023 at k 10).

Cells above .05 (residuals up to .056, adjacent .052) are chance: with the
dataset as cluster (two-group curves mirror, so they are flagged together),
the SE is 1.1 to 1.4 times the naive SE, the largest z against .05 is 1.92
(residuals) and 1.08 (adjacent), and no cell of 108 exceeds z = 1.96 (about 3
expected by chance at exact level). The residual band sits at nominal for two
equal groups (mean .0499). The adjacent band is the most discrete (fewest
points per curve) and the furthest below nominal at small n.

Status: the simultaneous band is validated for all three types on the full
grid. A full grid for
adjacent and residuals is optional (host, about 1 h). Next: implementation
design for RMitemCatFitPlot(dif_var), shown to the user before code.

## Implementation design: RMitemCatFitPlot(dif_var) (draft 2026-11-26)

For review before any package code changes. It supersedes step 1 decision 3
("no bands") and the implementation notes of 2026-11-24 where they differ.
Decisions marked **[D1]** to **[D6]** need the user's call.

### Scope

All changes are behind `dif_var`. Without it, RMitemCatFitPlot() behaves as
now. The only change outside the DIF path: `.catfit_prepare()` also returns
`rows`, the indices of the input rows it kept. That is a shared helper, so
the RMitemCatFit() tests are run as well as the plot tests.

### Arguments

- `dif_var = NULL`: a vector with one value per row of `data` (factor,
  character or numeric, treated as categorical). Placed after `fit` in the
  signature, so existing positional calls do not change.
- Number of permutations: **[D1]** reuse `iterations` (in DIF mode the number
  of permutations), or a new argument `permutations`. Reusing it keeps the
  signature short, and `iterations` already means "datasets drawn under the
  null" throughout the package. In both cases the default in DIF mode is
  1000 (via `missing()`, as `type` already does with `fit`): the band uses
  an extreme quantile of each curve's maximum, B was 999 in the studies, and
  1000 permutations take 0.3 to 2 s. The existing default of 200 for the
  non-DIF simulations is untouched.
- `seed`: as now (`set.seed(seed)` when given). RMdifAnchor() instead uses
  `.with_local_seed()`, which leaves the caller's stream untouched.
  RMitemCatFitPlot() does not, so DIF mode follows RMitemCatFitPlot() for
  consistency within the function.
- `conf_level`: the coverage of the simultaneous band in DIF mode.
- `parallel`, `n_cores`: unused in DIF mode (permutations are cheap and run
  sequentially). Said in the help.
- `grouping`, `class_intervals`, `score_breaks`: as now for `residuals`. Any
  grouping by total score keeps the band valid (see the progress log).

### Input checks and data

1. `type = "restscore"` with `dif_var`: error ("not available with
   `dif_var`", with the reason: the rest score given a category depends on
   the group's ability distribution).
2. `fit` with `dif_var`: error. The `fit` simulations are for item fit, and
   the pooled fit in DIF mode can use fewer rows.
3. `dif_var` through the shared `.prepare_dif_data()` (used by RMdifLR() and
   RMdifAnchor()): length check, rows with missing `dif_var` and empty
   respondents dropped with a message. Then `.catfit_prepare()` as now
   (complete cases, extreme-only categories, CML fit, extreme scores
   dropped), and `dif_var` subset with `prep$rows`.
4. At least two groups after filtering, else error. **[D2]** At most four
   groups, else an error saying to combine groups. Four keeps the ribbons
   readable. The permutation itself has no limit.
5. Pooled CML fit on all rows kept. Predictions are the same for every group.
6. A group with no respondents at some score simply has no point there.
   No minimum n per point (RMitemICCPlot() has `min_n`, but the band already
   reflects sparse points by being wide). Point size shows n.

### What is computed

Ported from `dev/catfit_dif/R/perm_functions.R` into internal functions in
`R/item_cat_fit_plot.R` (or a new `R/item_cat_fit_dif.R`):

- `.catdif_layout_*()` for each type: the fixed points of each group's curve
  (calibration: item x event x group over scores, adjacent: item x pair x
  group over scores, residuals: item x group over score groups x
  categories), and a value function of the group labels.
- `.catdif_bands()`: the generic band code (`perm_bands_generic()`), returning
  per point the observed value, the permutation mean, `lower` and `upper`.
  The band per point is mean +- c x SD, with c the `conf_level` quantile of
  the curve's max |z| over the permutations, clipped to [0, 1] for shares and
  CEPs. A curve leaves its band at a point exactly when that point's |z| > c.
- `.pav_w()`: the weighted PAV on per-score aggregates. It is checked against
  `.pav()` in the tests.
- Permutation within total scores: `permute_within()`.

### How each type is drawn

Groups use colour, line type and shape together (Okabe-Ito, figure
convention), and the legend shows each group with its n, e.g. "Women (n =
312)". **[D3]** RMitemICCPlot() with `dif_var` uses ggplot's default colours
and one shape for all groups, so the two DIF plots would look different.
Options: follow the figure convention here and align RMitemICCPlot() later,
or follow RMitemICCPlot() now. My suggestion is the convention, noting
RMitemICCPlot() for a later fix.

- `calibration`: per group a PAV step line and points (sized by n), and a
  ribbon in the group's colour (alpha about 0.2). Dotted diagonal as now.
  The band is centred on the permutation mean, which is close to the pooled
  curve, not on the diagonal. So an item that misfits for everyone moves all
  bands off the diagonal together, while DIF moves the groups apart.
- `adjacent`: per group points (sized by the number in the pair) and a ribbon.
  The model's expected share is the same for all groups and stays as one
  line. The binomial range of the non-DIF plot is not drawn.
- `residuals`: bars per group, dodged within each category, filled by group
  instead of by sign. The band as a dodged tile behind each bar, as the
  `fit` ranges are drawn now.
- **[D4]** Mark points outside their band: filled points inside, a
  contrasting outline (or larger point) outside. It helps reading many
  panels, but it moves the plot toward a test display. Suggested: yes, since
  the band already defines it.

### Caption (DIF mode)

Built from the existing type caption, with the band sentence replaced, for
example for calibration: "PAV-adjusted conditional event probability for each
group against the predicted probability given the total score, from one fit
to all respondents. Bands: the range of each group's curve expected without
DIF, simultaneous 95% for each curve, from 1000 permutations of group labels
within total scores. Not corrected across curves, so some curves cross their
band by chance. Test DIF with RMdifAnchor() or RMdifGamma()." Then the n
clause (`n = X of Y respondents (complete cases, extreme scores excluded)`)
and n per group. The `range_note` for fewer than 1000 permutations stays.

### Output = "dataframe"

The non-DIF columns, with the non-DIF band columns replaced, plus `group`,
`perm_mean`, `lower`, `upper`, `outside` (point outside its band) and
`iterations_used`. **[D5]** No per-curve p-value: it is an uncorrected test
and the first version is descriptive. It could be added later with the test
step.

### Help text

- `@param dif_var`, the DIF paragraph in `@param iterations`, notes in
  `parallel` and `fit`.
- New details section "**Groups (`dif_var`)**": pooled predictions and how to
  read separation vs common departure, the band (simultaneous consistency
  band, Dimitriadis et al. 2021, permutation within total scores, Kreiner
  1987), not corrected across curves, the Rasch-model condition (conditional
  independence given the total score holds under the Rasch model, not under
  models with unequal slopes, so with unequal slopes and different group
  means curves can separate without DIF: Kreiner and Christensen 2011),
  and unbalanced DIF in other items shifting the total score. Read item fit
  first ([RMitemCatFit()], [RMitemInfit()]).
- References added: Kreiner (1987), Kreiner and Christensen (2011).
  Holland and Thayer (1988) is a book chapter, so it is cited through
  Kreiner and Christensen rather than directly.
- Example: `dif_var` with `pcmdat2` and a simulated group, `seed`.
- **[D6]** Name the band in the help "simultaneous 95% consistency band" and
  in the caption in plain words only (as above).

### Tests

New file `tests/testthat/test-item_cat_fit_dif.R`, ported from `check_A.R`
and `check_A2.R` where fast:
- `.pav_w()` equals `.pav()`.
- Permutation keeps group counts per score.
- A label constant within scores gives zero-width bands.
- Same seed, same output. Non-DIF output unchanged (snapshot of the
  dataframe before and after).
- Adjacent shares per group equal `.catplot_adjacent()` per group, and
  one-group residuals equal `.catplot_residuals()`.
- Errors: restscore, `fit`, one group, five groups, wrong length; NA in
  `dif_var` gives a message and the right n.
- Each type returns a ggplot and a dataframe with the documented columns.
- `.catfit_prepare()$rows` maps back to the input rows.
Then `test-item_cat_fit.R` (shared helper changed).

### Docs and housekeeping

NEWS.md (one short bullet), README and `_pkgdown.yml` checked (no new
function), vignette: a short mention could follow later. Code review log in
dev/code-review/ after the user's /code-review. jamovi is a separate decision.

**2026-11-26, decisions (user, "go ahead" on the suggestions).** D1: reuse
`iterations`, default 1000 in DIF mode. D2: at most four groups. D3: figure
convention (Okabe-Ito, colour + line type + shape), RMitemICCPlot() noted for a
later fix. D4: points outside their band are marked. D5: no per-curve p in the
output. D6: "simultaneous 95% consistency band" in the help, plain words in
the caption. Implementation started the same day.

## Open questions

- `restscore` in DIF mode: the rest score given a category depends on the
  group's ability distribution. Check whether the current expected values are
  computed per respondent from the total score. If so, a group split may be
  valid. Until then, not offered.
- Maximum number of groups for readability (four?).
- Whether to keep the adjacent-plot binomial ranges per group later, or
  replace all bands with the permutation bands of step 3.
- `fit` together with `dif_var` (see Implementation notes).

## Implementation log

**2026-11-26, implemented** (uncommitted).
- `R/item_cat_fit_dif.R` (new): `.catdif_prepare()`, `.catdif_counts()`,
  `.pav_w()`, `.catdif_permute()`, layouts `.catdif_calibration()`,
  `.catdif_adjacent()`, `.catdif_residuals()`, `.catdif_bands()`, the three
  `_gg` functions, `.catdif_scales()`, `.catdif_styles()` and the driver
  `.catdif_run()`.
- `R/item_cat_fit_plot.R`: `dif_var` after `fit`. DIF checks (restscore,
  `fit`) before anything else. `iterations` defaults to 1000 when missing.
  Data through `.prepare_dif_data()`, then `.catfit_prepare()`, then
  `.catdif_run()`. The non-DIF path is unchanged. Help: `dif_var`,
  `iterations`, `conf_level`, `parallel`, `fit`, details "Groups
  (`dif_var`)", return value, Kreiner (1987) and Kreiner and Christensen
  (2011) references (checked against the PDFs), example.
- `R/item_cat_fit.R`: `.catfit_prepare()` also returns `rows`.
- Band rule in the package: c is the (m + 1)-th largest permuted curve
  maximum, m = floor(alpha (B + 1)) - 1, and a point is outside when its |z|
  exceeds c. A curve then has a point outside exactly when its permutation
  p is at most alpha, the rule validated in the studies.
- Calibration points with the same prediction (certain events at extreme
  scores) are merged, n summed.
- Tests: `tests/testthat/test-item_cat_fit_dif.R` (62 pass),
  `test-item_cat_fit.R` (178 pass, 1 skip as before). Help example runs.
- NEWS: the RMitemCatFitPlot() bullet under "New functions" mentions
  `dif_var` (the function is new in this development version).
- Checked by eye on pcmdat2 with a random split and on simulated data with
  DIF in two items: curves separate where DIF was put, ringed points mark
  points outside their band.

Left for later: vignette mention, RMitemICCPlot() group styling (D3), the
user's /code-review, jamovi.

**2026-11-26, code review (high)**, log in
`dev/code-review/2026-11-26-high-R-uncommitted.md`. Fixed: the `≥`
escape in the calibration label (CD1); the help says the band coverage
follows `conf_level` (CD2); in the residual plot, groups 2-4 get a line
pattern over their colour (grid patterns, ggplot2 >= 3.5.0), and `grid` is
now in Suggests (CD4); data is validated before the DIF split, and the caption
says "complete cases" when respondents with no responses were dropped (CD6).
CD5: `.catdif_bands()` computes the permutation moments and maxima one curve
at a time, so the permutation matrix is the only full-size object. Peak
memory at B = 1000 for 20 items x 4 groups (n = 2000) went from 1197 MB to
332 MB, output identical.

**2026-11-26, after the review.** Vignette: a prose paragraph "Category
functioning by group" after the partial gamma DIF section (both
`easyRasch2.Rmd.orig` and the precompiled `easyRasch2.Rmd`, no new chunk),
citing Kreiner & Christensen (2011), added to `vignettes/references.bib`.
D3 check: RMitemICCPlot() with `dif_var` still uses ggplot's default hue
colours and shape 18 for all groups, with no n in the legend, so it remains
a later fix. Greyscale check of the new plots (3 groups, 200 dpi):
calibration and adjacent stay apart by shape and line type. In the
residual plot, a bar outside its band gets a black outline of linewidth
0.8 on a bar of width 0.25, which at normal size fills the bar, so its
group shows only by position (open, see the session report). jamovi:
section G6 added to `easyRasch2jmv/dev/PLAN-easyRasch2-next-release.md`.

**2026-11-26, D3 and residual marks (user).** RMitemICCPlot() with
`dif_var`: points get one shape per group (16, 17, 15, 18, then 3, 4, 8, 7,
recycled), mapped together with colour so the legend stays one; the
caption says "Points (one shape per group)" instead of "Diamonds" in DIF
mode. Colours stay ggplot's default (only the shape was asked for). NEWS
bullet added. Residual plot: outline of bars outside their band thinner
(0.4), plus a black point with a white rim at the bar tip, without a legend
(the caption describes both). Tests: test-cicc_plot.R 59 pass, 1 skip,
test-item_cat_fit_dif.R 62 pass. Checked in greyscale at 200 dpi.
Later the same day (user): RMitemICCPlot() also uses Okabe-Ito colours in
DIF mode, the same first four as RMitemCatFitPlot(), then #56B4E9, #D55E00,
#F0E442, #999999, recycled (black left out, it is the model curve). Tests
pass (59, 1 skip). NEWS and jmv plan G6 updated.
Also: the legend key for the DIF error bars is vertical (`key_glyph = "linerange"`), matching the panels.
ICC shape order changed (user): circle, square, diamond, triangle (16, 15, 18, 17, then 3, 4, 8, 7), triangle last since it reads as an arrow on a vertical bar. Shape order now differs from RMitemCatFitPlot (21, 24, 22, 23), colours match.
2026-11-26 (user): the adjacent DIF plot no longer joins each group's points with a line (back to the design: points and ribbon, as the non-DIF plot). Groups by colour + shape there, line type only in the calibration plot. Tests 62 pass.

**2026-11-26, residual plot after user testing.** The user saw 50+ device
warnings ("Unsupported pattern extend mode; using repeat", from
`.setPattern()`, probably macOS quartz, not reproducible in the container
where the code passes `extend = "repeat"` correctly), a y-axis from -20 to
30 on real data, and groups dodged too far apart. Changes (user approved):
(1) the CD4 patterns are dropped, groups are plain Okabe-Ito fills and
their left-to-right order follows the legend; `grid` is out of Suggests
again. (2a) `.catdif_bands()` cuts each band to the range of the permuted
values and the observed value, for all three types. Cause of the wide
bands: for a cell the model expects to be almost empty, the residual takes
only a few values (e.g. -0.07 or 14 at E = 0.005), and mu +- c * SD
reaches values it cannot take (example: E = 0.07, band -4.0 to 7.3 where
z >= -0.26). The cut never changes `outside` (checked for 2 and 4 groups,
all three types). (2b) Residual y-axis = range of the bars and at least
+-4 (`coord_cartesian`), longer bands cut at the panel edge. (3) Bars of a
category side by side: bar width min(0.3, 0.7 / G), slot bar + 0.03,
tiles bar + 0.02. Help: `perm_mean` is no longer called the band centre,
and the cut is described. Tests 62 + 178 (1 skip) pass.
