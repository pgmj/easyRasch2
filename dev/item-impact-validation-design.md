# Item impact intervals: validation study design

Design note written 2026-10-12. Nothing has been run. This is the
pre-release study for the function proposed in
`item-impact-interval-design.md` (`RMitemImpact()`, infit only, all five
metrics, response dependence as the default overfit mechanism, complete
cases).

## 1. Why this study is needed

The dev studies each validated one step with the truth standing in for the
other steps:

| Study | Validated | Took as known |
|---|---|---|
| `conditional_rmsd_partA.qmd` | infit MSQ is n-stable | |
| `fwb_msq_coverage.qmd` | FWB interval for the MSQ, coverage near nominal | dichotomous items, B = 200 |
| `infit_impact_validation.qmd` | impact curves, n-invariance, mechanisms, joint failure | true thresholds and person distribution |
| impact-interval demo | mapping an interval through a curve | a curve from a different person distribution |

The function estimates the thresholds and the person distribution from the
data, builds the curve from those estimates, and maps an estimated interval
through it. This study asks whether the resulting impact interval covers the
true impact, and where the error comes from when it does not.

## 2. Questions

1. **End-to-end coverage.** Does the impact interval cover the true impact
   when the assumed mechanism is the true one?
2. **Error decomposition.** How much of any miss comes from the MSQ interval
   and how much from the estimated curve?
3. **Polytomous items.** Does the FWB interval for the MSQ hold for
   polytomous items?
4. **Mechanism mismatch.** What happens to the impact interval when the
   assumed mechanism is wrong, on each side?
5. **Partner choice.** Under response dependence, does the largest residual
   correlation find the true partner item?
6. **Cascade items.** Items flagged only because another item misfits fit
   the model in truth, so their true impact is 0. Are their impact intervals
   small enough to clear them?
7. **Group check.** Does it find a real group of items on a shared second
   dimension, and how often does it raise a false alarm?
8. **Simultaneous coverage at B = 1,000**, to set the default number of
   draws.
9. **Run time** of the prototype, for the function's start message and the
   choice of `n_sim` and `reps`.

## 3. Prototype

The study runs on a prototype of the function in
`dev/item_impact_prototype.R`, not on package code. It implements section 3
of the function note: FWB intervals (per item and simultaneous), item
selection, curve calibration per side with the counterfactual design,
mapping, keep or remove, and the group check. The prototype also returns the
estimated curves and the FWB draws so that the error can be decomposed.

Settings at their proposed defaults: `iterations = 200`, `n_sim = 5000`,
`reps = 20`, grids as in `infit_impact_validation.R` (ρ .3 to 1, λ 0 to 2,
slope 1 to 3).

## 4. Base designs

| Base | Items | Categories | True thresholds | True person distribution |
|---|---|---|---|---|
| d20 | 20 | 2 | locations seq(-2, 2), as in Part A | N(0, 1.5²) |
| phq9 | 9 | 4 | CML fit to the phq9 data, complete cases | the phq9 WLE pool |

These are the two bases of the impact validation, so the true impact curves
from that study can be reused as a check on the new truth runs. The phq9
person distribution is skewed with a floor, which is the realistic case for
a public health screening scale.

The target item is the on-target item from the impact validation (d20 item
nearest the person mean, phq9 item 1), plus the off-target item in one
condition.

## 5. Truth

For every condition, the true impact is computed the way the impact
validation computed its curves, with the true thresholds, the true person
distribution and the true mechanism: N = 20,000, 40 replicates, common
random numbers for the misfitting and fitting datasets, the removed-item
dataset, CML refit, WLE, linking through the other items. This gives the
population MSQ of every item and the true value of every metric. Metric E
is defined only under a second dimension, so it has no truth in the slope
and dependence conditions. Metric D is defined everywhere and reported on
the overfit side.

Items that fit in truth (all items in the null condition, all non-target
items elsewhere) have true impact 0 on every metric.

## 6. Arms

### Arm A: coverage, matched mechanism

n = 500, 2,000, 5,000. 100 datasets per cell.

| Condition | Target | Mechanism | Population MSQ (d20 / phq9) |
|---|---|---|---|
| A1 null | none | | 1 / 1 |
| A2 underfit, moderate | on | second dimension, ρ = .7 | 1.25 / 1.41 |
| A3 underfit, strong | on | second dimension, ρ = .5 | 1.40 / 1.65 |
| A4 underfit, off-target | off | second dimension, ρ = .7 | 1.20 / 1.42 |
| A5 overfit | on | response dependence, λ = 1 | 0.91 / 0.83 |

The target item is always calibrated (`items =` set by hand), so coverage
is measured at every n, also when its interval contains 1. The function's
own selection (items whose interval excludes 1) is recorded alongside.

In A2 and A3 at n = 2,000 and 5,000, every flagged item is calibrated, as
the function would do by default. This is where the cascade items appear
(question 6).

### Arm B: mechanism mismatch and partner choice

n = 2,000. 100 datasets per cell. Each dataset is analysed with both
mechanisms on the side its target lies on, so the two curves are compared on
the same data.

| Condition | Truth | Analysed with |
|---|---|---|
| B1 overfit by slope | target slope a = 2 (MSQ .74 / .71) | dependence (default) and slope |
| B2 underfit by flat slope | target slope a < 1, set for MSQ about 1.25 | second dimension (the only upper mechanism) |
| B3 overfit by dependence | λ = 1 on the target, as A5 | partner by largest Q3 and by true partner |

B2 needs a slope below 1 that is not in the earlier grids. A short pilot
sets it. B3 reuses the A5 data at n = 2,000 and adds the analysis with the
true partner, so it costs one extra calibration per dataset.

Coverage under mismatch is not expected. The question is the direction: a
default that overstates impact when it is wrong is safer than one that
understates it.

### Arm C: group check

n = 2,000. 100 datasets per cell.

| Condition | Truth |
|---|---|
| C1 joint group | 3 items on one shared second dimension, ρ = .6 (mean MSQ 1.20 / 1.16) |
| C2 unrelated misfit | 3 items, each on its own second dimension, ρ = .6 |

The A1 and A2 datasets at n = 2,000 serve as further false-alarm checks
(no group to find). The true joint impact for C1 comes from the joint arm of
the impact validation, extended to this item set if needed.

### Arm D: FWB intervals only (no calibration)

Polytomous coverage and the number of draws. phq9 base, null and A2, at
n = 500, 2,000, 5,000, 200 datasets per cell, B = 1,000. The d20 base reuses
the four conditions of `fwb_msq_coverage.R` at B = 1,000. Intervals at
B = 200 are taken from the first 200 draws of the same runs, so both sizes
are compared on the same data.

## 7. Outcomes

Per item, metric and cell:

- **Coverage** of the impact interval, per item and simultaneous over the
  flagged items. Monte Carlo SE about 2.2 points at 100 datasets.
- **Width** of the impact interval, as a ratio to the true impact where
  that is above 0.
- **Decomposition.** MSQ interval coverage of the population MSQ, and curve
  error: the estimated curve evaluated at the population MSQ, minus the true
  impact. Reported as bias and RMSE relative to the true impact.
- **Keep or remove** agreement with the truth.
- **Verdict** classification (inside, outside, undecided) at two example
  tolerances, KL .03 and shift 0.25 SEM, against the true classification.
- **Cascade items** (A2, A3): share whose upper KL and shift bounds fall
  below the example tolerances.
- **Partner choice** (B3): share of datasets where the largest residual
  correlation picks the true partner.
- **Group check** (C, A1, A2): detection rate, false-alarm rate, and
  coverage of the group row.
- **Run time** per dataset, split into FWB and calibration.

## 8. Criteria, set before the run

- Matched-mechanism coverage (A2 to A5) at least 90 percent per item at a
  nominal 95 percent, for every metric reported by default.
- Curve bias within 20 percent of the true impact where the true KL is
  above .01.
- Null (A1): family-wise rate of any item's simultaneous interval excluding
  1 at most 5 percent plus Monte Carlo error.
- Cascade items: upper KL bound below .03 in at least 95 percent of
  datasets.
- Polytomous FWB coverage between 93 and 97 percent per item.
- Group check: detection at least 80 percent in C1, false alarms at most
  5 percent in C2, A1 and A2.

If the end-to-end coverage misses while the MSQ interval covers, the curve
is the cause. Remedies to try, in this order: more `reps` or a larger
`n_sim` (Monte Carlo error), a person pool shrunk toward the estimated
distribution instead of the raw WLE pool (WLE overdispersion), and
propagating curve uncertainty by recalibrating within the FWB draws. Each
remedy is checked on A2 before rerunning the arm.

## 9. Things that could bias the curve

- **The WLE pool is wider than the true person distribution.** The curve is
  calibrated on persons drawn from that pool, so impact in SEM units may be
  biased. The decomposition shows how much.
- **The thresholds are estimated with the misfitting item in the data.** The
  item's own location and the other items' thresholds absorb part of the
  misfit, so the curve is calibrated around a slightly wrong model.
- **Monte Carlo error in the curve** at 20 replicates, larger near MSQ = 1
  where impact is small.

## 10. Run time and staging

From the earlier runs: about 0.2 s per CML fit at n = 5,000 for d20, so one
item side (9 grid points, 20 replicates, 3 fits) costs about 2 core-minutes,
and FWB with 200 draws costs about 0.1 s per draw at n = 5,000. phq9 fits are
faster. Rough totals on 10 cores:

| Arm | Datasets | Calibrations | Hours |
|---|---|---|---|
| A | 3,000 | about 3,400 | 12 to 16 |
| B | 600 | about 1,200 | 4 |
| C | 400 | about 1,600 | 5 |
| D | 1,200 phq9 + 2,400 d20 at B = 1,000 | none | 3 to 5 |

Stages:

1. Prototype and unit checks against the cached impact-validation cells
   (the prototype with true inputs must reproduce the cached curves).
2. Truth runs for every condition (minutes).
3. Timing pilot, 5 datasets per cell, and the B2 slope pilot. This
   replaces the estimates above.
4. Arm D, then A, then B and C, one overnight run each, cached undated and
   dated per cell.
5. Report in `item_impact_validation.qmd`, Typst PDF.

If the pilot shows A is too slow, the first cut is n = 500 to the target
item only (no cascade calibration there anyway), then 75 datasets per cell.

## 11. Files

- `dev/item_impact_prototype.R`: the function internals.
- `dev/item_impact_validation.R`: truth runs, arms, `IMPACT_VAL_SMOKE = 1`
  to define without running.
- `dev/item_impact_validation/`: cache, one file per cell.
- `dev/item_impact_validation.qmd`: report.

## 12. Stage 1 results (prototype and checks, 2026-10-12)

Files: `item_impact_prototype.R`, `item_impact_prototype_check.R`, results in
`item_impact_validation/prototype_check.rds` (+ dated) and
`prototype_check.log`.

**Reproduction.** `sim_rep()` reproduces `one_rep()` exactly (difference 0)
for all four mechanisms on both bases, and one cached replicate at
n = 20,000. Unit-weight FWB infit equals the direct conditional infit
(difference 2e-15).

**Curves at the proposed defaults** (`n_sim = 5000`, `reps = 20`) with the
true model agree with the cached curves (n = 20,000, 40 replicates). MSQ
within .007, KL within about 3 percent, F (cbias) within about 10 percent
(it is a maximum and the noisiest metric). One curve takes 4 to 11 seconds
on 10 cores, about 1 core-minute, half the estimate in section 10.

**Fixes made during the checks.**
- `keep_remove` was set to 1 at the null grid point. Removing a fitting
  item still costs precision (about .975 in d20), so it is now left as
  simulated.
- Intervals often reach past the grid end. The grid is now extended
  (ρ .15 and 0, slope 4 and 5, λ 3 and 4) when the interval reaches past it
  and the MSQ is still moving away from 1.
- Mapping is now along the curve path, so a curve whose MSQ turns back
  toward 1 (dependence in phq9) is handled. Two flags are reported:
  `beyond_grid` (the upper bound is at least the value shown) and
  `turns_back` (the MSQ does not identify the parameter, so the upper bound
  is set by the grid end).

**End-to-end, one dataset each at n = 2,000** (truth from the cached cells):

| Case | Metric | Interval | Truth |
|---|---|---|---|
| underfit, d20, ρ = .5 | KL | .034 to .046 | .037 |
| | shift (SEM) | .126 to .167 | .140 |
| | group (SD) | .048 to .057 | .050 |
| | cbias (SEM) | .155 to .268 | .171 |
| dependence, phq9, λ = 1 | KL | .014 to .042 | .041 |
| | shift (SEM) | .101 to .190 | .199 (missed) |
| | cbias (SEM) | .078 to .113 | .114 (missed) |
| joint group, d20, ρ = .6 | KL | .104 to .132 | .121 |
| | shift (SEM) | .266 to .319 | .297 |
| | cbias (SEM) | .40 to .56 | .47 |

The group check found the three-item group (Q3 cutoff .095) and its row
covered the true joint impact on every metric, while the per-item rows
showed KL of .02 to .03 each. Run time per dataset on 10 cores: 7 s
(null, FWB only) to 58 s (group with joint calibration).

**Finding 1: the dependence curve is biased low when estimated from the
data.** At λ = 1 the curve from the fitted model gives KL .026 against .041
in the true curve, and its MSQ never gets below .846 (true .83). The fitted
thresholds of the dependent item absorb part of the dependence. The
underfit curve shows no such bias (KL .036 against .037 at ρ = .5). This is
the curve error of section 9, and arm A5 will measure it. A possible remedy
is to calibrate the dependence side on the pair's residual correlation
(Q3), which rises with λ, instead of on the MSQ, which turns back.

**Finding 2: dependence is the wrong default for cascade items.** A
cascade item has no real partner. The largest Q3 then picks an arbitrary
item, dependence on that item barely moves the MSQ, and the impact per unit
of MSQ is steep. On the same datasets:

| Item | Partner Q3 (cutoff) | KL, dependence | KL, slope |
|---|---|---|---|
| d20 cascade, V1 | .052 (.093) | .0006 to .0053 | .0007 to .0104 |
| d20 cascade, V20 | .031 (.093) | .0011 to .0386 | .0001 to .0131 |
| d20 joint case, V14 | .061 (.095) | .0014 to .0322 | .0003 to .0094 |
| d20 joint case, V20 | .034 (.095) | .0032 to .0383 | .0003 to .0131 |
| phq9 true dependence, q5 | .366 (.081) | .0135 to .0422 | .0078 to .0214 |

With dependence, two of four cascade items exceed KL .03 at the upper
bound, so they would not be cleared. With slope, all clear by a wide
margin. For the true dependent item, slope understates (upper .021 against
truth .041). The partner's Q3 separates the two cases cleanly here (.37
against .03 to .06, cutoff about .09).

**Decided 2026-10-13:** `overfit_mechanism = "auto"` as the
default: dependence when the item's largest Q3 exceeds the
`RMlocdepQ3Cutoff()` value, slope otherwise, with the choice shown per item
in the output. Arm B would then test "auto" alongside the two fixed
mechanisms.

## 13. Changes after stage 1 (2026-10-13)

- The prototype default is `overfit_mechanism = "auto"`. B1 is analysed
  three ways (auto, dependence, slope). B3 reuses the A5 data at n = 2,000
  with dependence and the true partner forced.
- B2 flat slope set by a pilot at n = 20,000: a = .55 (d20, MSQ 1.25) and
  a = .60 (phq9, MSQ about 1.27). At the same MSQ a flat slope has less
  impact than a second dimension (d20 KL .016 against about .027), so the
  upper-side curve should overstate here.
- Arm D uses this study's generator for d20 (base truth) instead of the
  Part A datasets, so both bases are treated alike.
- Every dataset is analysed with the default selection plus the target
  item(s) (`items_extra`), so cascade items are recorded in every arm.
- Script: `item_impact_validation.R`. Pilot cells in
  `item_impact_validation/pilot/`.

## 14. Timing pilot (2026-10-13)

5 datasets per cell, all 48 cells, no failures. Log in
`item_impact_validation/pilot.log`, cells in `item_impact_validation/pilot/`.

**Projected full run on 10 cores** (100 datasets per cell, 200 for D):

| Arm | Hours |
|---|---|
| A | 17.8 |
| B | 2.2 |
| C | 6.5 |
| D | 4.5 |
| Total | 31 |

Mean core-seconds per dataset: d20 120, phq9 383. Calibration dominates,
about 85 to 90 core-seconds per item side. phq9 is slow because the cascade
flags most of its 9 items at n = 2,000 and 5,000: A2 calibrates 6.2 items
on average at n = 5,000, A3 8.2 at n = 2,000, and C2 all 9. In a short
scale with one strongly misfitting item, the function by default calibrates
nearly every item. On 10 cores that is still about one minute for a user.

**Sanity check of the target items** (5 datasets, not a result):

- Matched second dimension (A2 to A4, both bases): KL and shift covered in
  every dataset.
- B1 steep slope: "auto" chose slope in every dataset and covered. Fixed
  dependence overstated impact and missed (d20 .026 against truth .019),
  which is the conservative direction.
- B2 flat slope analysed as a second dimension: overstated and missed (d20
  KL .021 to .033 against .016), also conservative.
- A5 dependence, phq9: KL covered, shift covered in 1 of 15 target rows
  over the three n. Forcing the true partner (B3) did not help, so the cause
  is the curve bias from section 12 (fitted thresholds absorb the
  dependence), not the partner choice.
- Group check: C1 found in 5 of 5 datasets on both bases, C2 no false alarm
  in 5 of 5.
