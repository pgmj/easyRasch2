# Item impact intervals: package function design

Design note written 2026-10-12. No package code yet. It turns the results of
four dev studies into a function proposal:

- `conditional_rmsd_partA.qmd`: infit MSQ is an n-stable effect size, and
  the large-sample problem is a testing problem, not a measurement problem.
- `fwb_msq_coverage.qmd`: FWB intervals for the observed MSQ have the right
  spread and close to nominal coverage.
- `infit-impact-delta-design.md` and `infit_impact_validation.qmd`: impact
  of misfit on person scores can be simulated against a same-item
  counterfactual, is n-invariant, and depends on test length and mechanism.
  Per-item impact fails for a coherent group of misfitting items.
- The impact-interval section of `infit_impact_validation.qmd`: an FWB
  interval for the MSQ, mapped through the item's impact curve, gives an
  interval for the item's impact in score units.

## 1. Purpose

For each item, report how much its misfit changes the person scores, as an
interval in units a user can interpret: SEM, SD of the trait, or KL. This
answers "does this misfit matter?" instead of "is the model exactly true?".
At large n the second question always has the answer no, and the first is
the one that decides what to do with the item.

The function does not need a tolerance. A tolerance, when given, adds a
verdict (inside, outside, undecided) as a convenience.

## 2. Proposed functions

Following the house pattern of a computing function plus a plot sibling:

```r
RMitemImpact(
  data,
  items = NULL,            # NULL: items whose MSQ interval excludes 1
  metrics = c("kl", "shift", "group", "cbias", "rmse"),
  iterations = 200,        # FWB draws for the MSQ interval
  simultaneous = TRUE,     # one max-type constant over items
  conf_level = 0.95,
  overfit_mechanism = c("auto", "dependence", "slope"),
  n_sim = 5000,            # persons per calibration dataset
  reps = 20,               # calibration replicates per grid point
  tolerance = NULL,        # named list, e.g. list(kl = 0.03, shift = 0.25)
  group_check = TRUE,      # see section 5
  output = c("kable", "dataframe"),
  parallel = TRUE, n_cores = NULL, seed = NULL, verbose = FALSE
)

RMitemImpactPlot(impact, item = NULL, metric = "kl")   # impact curve + intervals
```

Metric names map to the validation study: `rmse` = A, `shift` = B (mean
|score change| in SEM), `kl` = C, `group` = E, `cbias` = F (worst conditional
bias in SEM). Metric D (true error over reported SEM) is added automatically
on the overfit side, where it is the relevant harm under response
dependence.

Placeholder names. `RMitemImpact` fits the planned `RMdifImpact` (DIF impact
function) if both are built.

## 3. What it computes

1. **Fit.** CML thresholds (`.fit_cml_thresholds()`), WLE pool
   (`.wle_theta_pool()`), conditional infit MSQ for every item.
2. **MSQ intervals.** FWB: `iterations` draws of exponential weights scaled
   to mean 1, weighted CML refit (`psychotools::pcmodel(weights = )`),
   weighted conditional infit. Percentile intervals per item, or simultaneous
   intervals from the max-type constant over items (the infit analogue of
   Q3max*). Validated in `fwb_msq_coverage.qmd` (dichotomous, complete data).
3. **Which items to calibrate.** By default the items whose interval
   excludes 1. Calibrating every item costs one to two minutes per item, so
   this keeps the function practical. `items` overrides it.
4. **Impact curve per item.** Simulation from the fitted model, as in
   `infit_impact_validation.R`:
   - persons resampled from the WLE pool, thresholds from the CML fit;
   - upper side (MSQ above 1): the item responds to a second dimension,
     grid of correlations .3 to 1;
   - lower side (MSQ below 1): a steep slope (grid 1 to 3), or uniform
     response dependence on the nearest item (grid 0 to 2);
   - for each grid point, datasets (a) misfitting and (b) the same item
     fitting, with common random numbers, plus (c) the item removed. CML
     refit, WLE, linked through the other items, impact metrics against
     the true trait;
   - the curve is impact against the item's MSQ in (a), with impact 0 at
     MSQ = 1.
   Only the side the item's interval lies on is calibrated, or both if the
   interval contains 1.
5. **Impact interval.** Map the MSQ interval through the curve: the minimum
   and maximum of impact over the interval (0 at the lower end when the
   interval contains 1).
6. **Keep or remove.** RMSE with the item kept over RMSE with it removed, at
   the item's estimated MSQ, read from the same simulation (dataset c).
7. **Verdict (optional).** With `tolerance`, inside if the whole impact
   interval is below it, outside if the whole interval is above it,
   undecided otherwise.

## 4. Output

One row per item (calibrated items have impact columns, the others show only
the MSQ and its interval):

| Item | Infit MSQ | MSQ interval | Side | KL | Shift (SEM) | Group bias (SD) | Worst cond. bias (SEM) | Keep/remove | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| q4 | 1.38 | 1.32 to 1.44 | under | .032 to .041 | 0.12 to 0.15 | .040 to .046 | 0.13 to 0.20 | keep | outside (KL) |

Kable caption, following the house rules: the estimation sample
(`.n_caption()`), the interval type ("simultaneous 95% FWB intervals,
200 draws"), the mechanisms used per side, and the calibration size
(`n_sim`, `reps`). A footnote states what the impact is relative to ("change
in WLE person estimates relative to the same item fitting the Rasch model").

`output = "dataframe"` returns the table plus the impact curves and the FWB
draws as attributes, so `RMitemImpactPlot()` can draw them without
recomputing.

**Plot.** For one item: the impact curve against MSQ, the observed MSQ with
its interval as a horizontal bar on the x axis, and the mapped impact
interval as a bar on the y axis. A dashed line marks the tolerance when one
was given.

## 5. Group check

The joint arm showed that when several items share a second dimension,
each item's MSQ is smaller than alone while the joint impact is far larger
(phq9: RMSE up to 65 percent higher, KL 1.4 to 1.7 times the sum of the
single-item values). Per-item intervals then understate the harm.

With `group_check = TRUE`:

1. If two or more calibrated items are on the underfit side, compute the
   residual correlations (Q3) among them, relative to the mean Q3.
2. If they form a cluster above the `RMlocdepQ3Cutoff()` critical value, or
   if `RMdimResidualPCA()` puts them on the same contrast, add a **group
   row**: the joint impact of the set, from the joint mechanism (all items in
   the set on one shared second dimension), mapped from the set's mean MSQ.
3. Always print a footnote when a group is found: per-item impact understates
   the joint impact, and a dimensionality analysis is the better tool.

The joint mapping is less certain than the per-item one (it assumes one
shared factor), so the footnote should say that the group row is indicative.

## 6. Run time

From the dev runs, on one core:

- FWB intervals: about 0.02 s per draw at n = 500 and 0.1 s at n = 5,000 for
  20 dichotomous items, so 4 to 20 seconds at `iterations = 200`.
- Calibration: 9 grid points x 20 replicates x 3 CML fits at `n_sim = 5000`,
  about one minute per item side.

With mirai on 10 cores, a scale with three flagged items should take well
under a minute. A message at the start should state how many items will be
calibrated and roughly how long it will take.

## 7. What must be validated before release

1. **End-to-end coverage.** The demonstration mapped intervals through a
   curve calibrated on a different person distribution, and its coverage
   only reflected the MSQ interval. The real test: generate data with known
   misfit, run the whole function, and check whether the impact interval
   covers the true impact (computed from the true generating model). This
   includes the error in the curve, which comes from estimated thresholds
   and the WLE pool, not the truth.
2. **Polytomous FWB coverage.** `fwb_msq_coverage.qmd` was dichotomous only.
3. **Missing data.** Weighted CML handles missing responses, but conditional
   infit with missing data needs the per-pattern moments (`.cond_groups()`),
   and the counterfactual simulation needs a missingness model (resampling
   observed missingness patterns, as the Q3 cutoff does).
4. **Mechanism mismatch.** The upper side assumes a second dimension. If the
   real cause is a flat slope, the curve differs. A sensitivity run with the
   slope mechanism on the upper side would show how much.
5. **Group check accuracy.** How often it finds a real group, and how often
   it raises a false alarm when misfitting items are unrelated.
6. **Simultaneous coverage with more draws.** 91 to 95.5 percent at
   B = 200. Recheck at B = 1,000 before deciding the default.

## 8. Decisions

Settled 2026-10-12:

1. **Scope: infit only.** The gamma difference is not included.
2. **Default metrics: all five** (`kl`, `shift`, `group`, `cbias`, `rmse`).
3. **Overfit mechanism default: "auto"** (changed 2026-10-13 from response
   dependence, after the prototype checks in
   `item-impact-validation-design.md` section 12). Dependence when the
   item's largest residual correlation (Q3 minus the mean Q3) exceeds the
   `RMlocdepQ3Cutoff()` value, with that item as the partner. A steep slope
   otherwise. The signature becomes
   `overfit_mechanism = c("auto", "dependence", "slope")`, and the output
   shows the mechanism and partner per item. Metric D (true error over
   reported SEM) is reported on the overfit side.
4. **Complete cases** in the first version, with a message stating how many
   respondents were dropped (`.n_caption()` convention). Section 7 item 3
   is deferred to a later version.

Still open:

5. **jamovi.** A module analysis later, after the package version settles
   (memory: ask before crossing the package boundary).
6. Doc checklist when built: `_pkgdown.yml`, `NEWS.md`, `README.md`.
