# Design: DIF magnitude by impact on person scores

Status: DESIGN, 2026-10-03. Nothing built or run. Two parts. Part A settles the
inference layer of the partial gamma DIF functions, so that detection has the
same defaults as the other cutoff functions. Part B designs a new function
that judges DIF magnitude by its effect on person scores and group
comparisons, which is the last Rasch analysis in the package still interpreted
with fixed rules.

Background: discussion of McNeish and Wolf's dynamic fit index cutoffs
(`articles/simulationbasedcutoffs/`). Their argument against Hu and Bentler
cutoffs applies to DIF magnitude rules. A fixed logit rule (`RMdifLR()`,
`cutoff = 0.5`), the Bjorner et al. (1998) gamma bands (.21/.31, cited in the
`RMdifGamma()` help page) and the ETS A/B/C categories in `RMdifTree()` all
mean different things in different scales.

---

## Part A. Partial gamma DIF inference

### Decided 2026-10-03

1. **Two-sided p-values stay.** DIF has no privileged direction, unlike local
   dependence where `RMlocdepGamma()` is one-sided upper. The asymmetry
   between the two neighbouring functions becomes a stated decision: one
   sentence in each help page saying why.
2. **Westfall-Young p-values become the default** in `RMdifGamma()`, using
   the `p_value = NULL` convention of `RMitemInfit()` and
   `RMitemRestscore()`. `NULL` resolves to `TRUE` when `cutoff` is the full
   `RMdifGammaCutoff()` object and to `FALSE` otherwise. Explicit `TRUE`
   without the full object is an error. `FALSE` with a `cutoff` flags on the
   interval and gives the one-time family-wise error message. This changes
   the default behaviour of an exported function, so it needs a NEWS bullet
   (results move for anyone passing a full cutoff object without
   `p_value`).

### Decided 2026-10-03 (second round)

- **3 approved:** `RMdifGammaCutoff()` defaults become `iterations = 400`,
  `hdci_width = 0.95`, matching the item-fit cutoff pair.
- **4 approved:** the vectorised partial gamma replaces
  `iarm::partgam_DIF()` for the simulated null and the observed statistic
  (validated in `dif_partgam_fast_validation.qmd`). The user-facing SE and
  CI columns may still come from iarm, as on the LD side.
- **5:** one sentence in the `RMdifGamma()` help page saying that iarm's
  adjusted column is a Bonferroni correction over items whatever method name
  it carries.
- **6:** study item 4 RUN 2026-10-03 (`dif_gamma_null_groups.qmd`, R = 300,
  B = 400, 12 conditions). The current random-assignment null fails with a
  minority group 1 SD apart: family-wise .137 (dich) and .147 (poly) at
  480/120, null SD down to 0.75 of the true SD for the easiest item.
  Groupwise, conditional and within-score-stratum permutation all hold
  (.023 to .067). **DECIDED 2026-10-04:** `dgp = c("conditional",
  "permutation")`, conditional default (simulates from the fitted Rasch
  model, as the other cutoff functions do, and matches
  `RMitemRestscoreCutoff()`). Permutation is an option: the Monte Carlo
  form of Kreiner's (1987) exact conditional test, independent of estimated
  thresholds, keeps observed patterns (so it should be more robust to misfit
  elsewhere, untested), and 25 to 60 times faster per iteration (0.6 to
  1.0 ms against 15 to 60 ms, n = 300 to 1000). Random assignment is
  dropped.

No package change yet. Points 1 to 5 can be implemented together when the
hold on package changes is lifted.

### Original proposals

3. **Cutoff defaults match the item-fit pair.** `RMdifGammaCutoff()` has
   `iterations = 250, hdci_width = 0.99`. `RMitemInfitCutoff()` and
   `RMitemRestscoreCutoff()` have `400` and `0.95`. Proposal: `400` and
   `0.95`. At 250 the 0.99 interval is close to the sample range and the
   `m/(B+1)` floor limits the p-values.
4. **Vectorised partial gamma for DIF.** Same Davis (1967) statistic,
   conditioned on the full score instead of a rest score. If it extends, the
   cutoff gets the Q3-side speed-up and higher `iterations` becomes cheap.
   Validate against `iarm::partgam_DIF()` to machine precision first.
   **Validation required before use (2026-10-03).** The package already has
   `.partgam_ld_gamma()` / `.partgam_one()` in `R/ld_partgam.R`, validated
   exactly against `iarm::partgam_LD()`. The DIF case differs in ways that
   the LD validation does not cover:
   - the table is rectangular (m item categories × g groups), while
     `.partgam_one()` builds one m × m `G` and assumes both variables share
     `m`. It needs separate row and column indicator matrices.
   - the conditioning score is the **full** score including the item, not a
     rest score.
   - the group variable is a factor. `iarm::partgam()` treats it as ordered
     by level, so the sign of gamma depends on level order, and with three or
     more groups the groups are treated as ordinal. The fast version must
     reproduce the same coding (factor levels, integer codes, unused levels).
   - NA: `iarm::partgam_DIF()` computes the score with `na.rm = TRUE` on all
     rows and then drops incomplete rows (items or exogenous variable), so it
     is complete-case in effect. The fast path must filter first.
   Validation matrix: as for LD (category counts, unobserved category, sparse
   strata, n = 30 to 1000, dichotomous), plus both level orders, unequal
   group sizes, a group empty within some strata, NA in `dif_var`, three
   groups (for later), and a constant item. Exact agreement with
   `iarm::partgam_DIF()`, brute-force reference where iarm errors.
5. **`iarm` adjusted column. CONFIRMED 2026-10-03 from the source.**
   `partgam_DIF()` calls `p.adjust(pvalue, method = padj, n = l * k)` inside
   the loop on a single p-value. With one value, BH, Holm and Bonferroni all
   return `min(1, p * n)`, so the column is Bonferroni over k items × l
   exogenous variables whatever method is named, and is labelled
   `padj.BH` by default. It only affects the asymptotic path, which loses its
   default role under decision 2.
6. **Null with unequal groups (new question).** `RMdifGammaCutoff()` assigns
   respondents to groups at random, so simulated groups share one latent
   distribution. Under the Rasch model the null itself does not depend on
   the group difference (the response is independent of group given the total
   score), but the null distribution of gamma depends on how each group is
   spread across score strata. Option: simulate each group from its own WLE
   pool (`dgp` argument as in the restscore cutoff). Needs a small null study
   with a group difference of 0, 0.5 and 1 SD before deciding. Part B needs
   this null too, so settle it once.

---

## Part B. `RMdifImpact()` (placeholder name)

### Decided 2026-10-03

| Question | Decision |
|---|---|
| Anchors | Automatic selection by default if the dev study shows it is reliable, plus a user-supplied `anchors` argument |
| Headline metric | Group-mean bias, presented didactically: the biased estimate and its SE, CI and group SDs in a simple two-group comparison |
| Scope v1 | Uniform DIF, one categorical variable, two groups |
| Default tolerances | None. Report impact only, so that no new rule of thumb is created |
| Packaging | New function, `RMdifLR()` unchanged |

### Theory

For uniform DIF of size d on item i, the ML estimate of a person with raw
score r shifts by about d × I_i(θ) / Σ I(θ), the item's share of test
information at θ. Checked against exact score-equation solutions for one
dichotomous item at d = 0.5 logits, other item parameters known:

| k | Max shift (logits) | As fraction of SEM |
|---|---|---|
| 5 | 0.127 | 0.12 |
| 10 | 0.061 | 0.09 |
| 20 | 0.030 | 0.06 |

The first-order formula is within 0.01 logits mid-scale. Consequences:

1. The same d has about four times the effect at k = 5 as at k = 20, and less
   still on an off-target item. A fixed logit rule ignores both.
2. For individuals the shift is small against the SEM. For a group comparison
   it is not, because averaging removes noise but not bias. The bias in the
   group-mean difference is the average shift over each group's score
   distribution, accumulates over items, and cancels for DIF of opposite
   sign (the DTF argument, Chalmers, Counsell & Flora, 2016).
3. Given group-specific item parameters the shift at each raw score is a
   deterministic calculation. Simulation is needed for something else:
   estimated DIF is never zero, so impact computed from estimates is
   inflated by noise.

### Anchors and the split model

The reference is the **split model**: the DIF items get group-specific
parameters, the remaining items are common. The common items are therefore
the anchor set. Choosing which items to split and choosing anchors are the
same decision seen from two sides.

DIF is identified only relative to anchors (Bechger & Maris, 2015,
*Psychometrika*, verify before citing). If most items shift the same way, it
cannot be told apart from a true group difference. Automatic selection is
reliable only under sparse DIF. Plan:

- **Default:** alignment selection (Gini criterion, Strobl et al., 2021,
  verify), implemented on item location differences from separate group
  fits. psychotools 0.7.6 `anchortest()` accepts only `raschmodel` fits,
  so polytomous data needs our own implementation. Alignment works on any
  parameter vector, so this is small.
- **User-supplied:** `anchors = c(...)` overrides selection.
- **Always reported:** the headline impact under the selected anchors next
  to all-item centring (what `RMdifLR()` does now through eRm's per-group
  centring). Disagreement between the two is the signal that DIF is not
  sparse.

### Remedies compared

| Remedy | Bias vs split | Cost |
|---|---|---|
| Keep (pooled model) | the DIF shift | none |
| Remove the item(s) | none, given correct anchors | lost information, larger SEM |
| Split the item(s) | reference | estimation error in group-specific parameters, worst for the smaller group |

The help page says that splitting is a substantive choice (accepted in
clinical measurement, contested where fairness rules apply), not a
recommendation.

### Output (draft)

1. **DIF table:** per item, d̂ in logits with CI from the anchored fits, the
   anchor set used, and the item's information share in each group.
2. **Group comparison table (headline):** for keep, remove and split, the
   group means, SDs, mean difference, its SE and 95% CI, and the
   standardised difference. A raw sum score row is worth considering, since
   that is what many applied users compare.
3. **Shift curve (ggplot):** θ shift by raw score for each group, keep vs
   split, with the SEM band for scale.
4. **Noise floor:** the observed keep-vs-split bias against its distribution
   under no DIF, as a percentile or Monte Carlo p.

No tolerance column and no flag. The user judges.

### Inference

- **Analytic layer:** shift by raw score from the split-model parameters
  (WLE score equation), weighted by each group's observed score distribution,
  gives the expected group-mean bias with no simulation.
- **Uncertainty:** parametric bootstrap from the split model, each group from
  its own θ distribution, refit all three models, recompute all metrics,
  percentile CIs.
- **Noise floor:** parametric bootstrap from the pooled model (no DIF), each
  group from its own θ distribution (Part A question 6), split the same items,
  same pipeline.
- **Selection caveat:** if the split items were chosen because they were
  flagged in the same data, the noise floor understates the no-DIF impact
  (the same lesson as the PCA/t-test pipeline). v1 conditions on the item
  set and says so. A version that repeats the selection step in each
  replicate is a later option.

### Estimation engine (checked 2026-10-03)

The split model is fitted as an item-split dataset: item i becomes `i_g1`
(group 1 responses, NA for group 2) and `i_g2`, all other items common.
One test dataset (8 items, 3 categories, n = 1000, true d = 0.8 on item 1,
group difference 0.5):

- `eRm::PCM()` recovers it: `i1_g2 - i1_g1` = 0.75.
- `psychotools::raschmodel()` handles the structural NA (dichotomised data).
- **`psychotools::pcmodel()` (0.7.6) fails on it**: positive log-likelihood,
  nonsensical parameters, Hessian not invertible. The pooled fit on the same
  data is fine.

**Failure isolated (2026-10-03, scratch tests, no package code).** `pcmodel()`
handles NA in general: 10% MCAR and one item missing for half the sample both
give sensible estimates. It fails when **no row is complete**: complementary
NA on two different items (not a split) fails the same way, and adding five
complete rows makes it work. Column order and matrix vs data.frame make no
difference. The dichotomised split data run through `pcmodel()` did not fail.
This is a minimal example for an upstream report. Adding complete rows is not
a workaround, since a split model never has any.

**Cause found (2026-10-03).** In the missing-data branch, `pcmodel()` builds
the per-pattern parameter lists with
`mapply(split, esf_par, mapply(function(x, y) rep(1:x, y), m_i, oj_max_i))`
in three places (`cloglik`, `agrad`, and the `full` ESF). When every NA
pattern has the same number of items and categories, both `mapply()` calls
simplify their results to matrices instead of lists, and the likelihood is
computed on the wrong structure. With a complete-case pattern present, the
lengths differ and the result stays a list, which is why five complete rows
"fix" it. Dichotomous data survive by accident (one-element items).
Adding `SIMPLIFY = FALSE` to both calls in all three places (patched copy in
the session scratchpad, not installed) gives: DIF 0.755, identical to
`eRm::PCM()`, all thresholds within 1.2e-5 of eRm after centring, and
coefficients identical to unpatched `pcmodel()` on complete data and on 10%
MCAR data. `rsmodel()` has the same construction (three places) and probably
the same bug, not tested. Upstream: psychotools 0.7.6, maintainer Achim
Zeileis, no BugReports field in DESCRIPTION.

The fix does not remove the need for our own CML: `pcmodel()` fits free
parameters only, so it can fit a free-threshold split once fixed, but not the
uniform-shift model below. A fixed `pcmodel()` would serve as a second
validation reference and as the engine for a later DSF extension.

**Decision pending: own CML, eRm as reference only (recommended).**

The model to fit is not the free-threshold split but a **uniform-shift
split**: item i keeps one set of thresholds τ_i, and group 2 gets
τ_i + δ_i. This is the v1 scope, d̂ is the DIF parameter itself (with an SE)
rather than a mean of threshold differences, and it adds one parameter per DIF
item instead of m − 1. It is also the sparse-data fallback: a category empty
in one group leaves δ_i estimable because the thresholds are shared. A free
threshold split (and eRm) fails in that case, since the empty category's
threshold has no information in that group.

eRm cannot fit the shift model directly (LPCM with a design matrix could, but
inherits eRm's behaviour with low cell counts). Implementation:

- Two missingness patterns (one per group), so the conditional log-likelihood
  is the sum over groups of a PCM conditional likelihood in that group's item
  parameters, with τ shared and δ added for group 2 on the DIF items.
- Per group: category counts per item and the raw score distribution as
  sufficient statistics, log γ_r from
  `psychotools::elementary_symmetric_functions()` (order 1 gives the gradient
  of the ESF, so an analytic gradient is cheap). Chain rule from (τ, δ) to the
  per-group threshold vectors.
- `optim()` BFGS or `nlminb()` with the analytic gradient, SEs from a
  numerical Hessian of the gradient. Normalisation matching `pcmodel()` so
  results compare directly.
- Null categories: downcode within an item when a category is empty in
  **both** groups (as `pcmodel(nullcats = "downcode")`). Empty in one group
  only needs no special handling under the shift model, which is the point.
  Extreme scores contribute nothing to a CML and are dropped per group.

Rough size: about 150 lines of core code plus input handling, comparable to
an existing cutoff function. Validation (part of the dev study, item 0):

1. **No DIF items:** must equal `pcmodel()` on complete data to optimiser
   tolerance, dichotomous and polytomous, several k and n, with MCAR NA.
2. **Shift model vs free split:** with one DIF item and a large n, d̂ agrees
   with the mean threshold difference from `eRm::PCM()` on the split data,
   and the log-likelihood of the shift model is at or below the free split.
3. **Recovery:** bias, SE calibration and CI coverage of d̂ by group size,
   including small groups and a category empty in one group.
4. **Gradient check** against numerical derivatives.

Person estimates for all three scorings use the package WLE (dropped item for
"remove", group-specific shifted thresholds for "split"). The simulation DGP
for v1 shifts all thresholds of a DIF item by d, the same model as the fit.
Free-threshold splits (DSF) are a later extension.

### Draft signature

```r
RMdifImpact(
  data,
  dif_var,
  items,                 # items to split (required in v1)
  anchors = NULL,        # NULL = automatic alignment selection
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  seed = NULL,
  verbose = FALSE,
  output = c("kable", "dataframe", "ggplot")
)
```

Open: whether `items` can take an `RMdifGamma()` result, and whether `output`
needs separate kable views for the DIF table and the group comparison table.

---

## Development study (ask before running)

| # | Question | Kind |
|---|---|---|
| 0 | Own uniform-shift CML: validation 1 to 4 above | deterministic + small recovery simulation |
| 0b | Vectorised partial gamma for DIF: validation matrix in Part A | deterministic |
| 1 | First-order shift formula for PCM items, and how close the analytic group-mean bias is to simulated | deterministic, cheap |
| 2 | Anchor selection: how often alignment picks a DIF item, by proportion and balance of DIF items, k, group sizes | simulation |
| 3 | Noise floor: impact under no DIF by size of the smaller group, and whether a correction is needed | simulation |
| 4 | Part A question 6: partial gamma null with unequal group distributions | simulation |
| 5 | Power: partial gamma WY vs item-wise Wald on anchored d̂ (feeds a later sensitivity report) | simulation, later |

Items 1 and the engine check (`pcmodel` minimal example, or the custom CML
validated against eRm) come first, since everything else depends on the
split fit.

### Status 2026-10-03

- **0:** `dif_shift_cml.R` (prototype) + `dif_shift_cml_validation.qmd`.
  Deterministic parts run (code check, not rendered): matches `pcmodel()`
  with no DIF items (max parameter difference ≤ 9e-7), analytic gradient
  matches numerical (≤ 7e-7), one DIF item at n = 2000 + 2000 gives d̂ 0.746
  against 0.745 for both free splits (LR vs free split p = .41), threshold-
  specific DIF is detected by the LR (p = 2e-7), and a category empty in the
  focal group still gives d̂ 0.35 (SE 0.20) where eRm returns NA for the
  threshold and patched `pcmodel()` silently drops it. **Recovery run and
  rendered 2026-10-03** (16 cells × 500 reps, 13 s, cached in
  `dif_shift_cml_recovery.rds`): no failures down to 100/50, bias within
  0.012 for the central item and up to 0.050 for the item at +2 logits at
  100/50 (about 2.5 MCSE), SE ratio 0.93 to 1.10, coverage .928 to .968.
  Empty focal category in 12.4% of reps in the hardest cell, handled.
- **0b:** `dif_partgam_fast.R` + `dif_partgam_fast_validation.qmd`. Run and
  rendered 2026-10-03: exact agreement (difference 0) with
  `iarm::partgam_DIF()` in all 19 cases iarm handles, and with the
  brute-force count in all 20, including reversed levels, character and
  numeric coding, unused level, 90/10 groups, a group absent from strata,
  sparse n = 30, NA in items and grouping variable, three groups. iarm errors
  on a constant item; the fast version returns NA for that item only. 0.5 ms
  against 33 ms per call at n = 600 (about 65x).
- **Bug report:** `psychotools-na-bug/reprex.R` and `email.md`, draft for the
  user to send. `rsmodel()` confirmed to have the same bug and fix.

## Docs checklist when built

`_pkgdown.yml` reference entry, `NEWS.md` (one bullet for the new function,
one for the `RMdifGamma()` default change with "results move"), README if it
lists DIF functions, and a vignette paragraph replacing the hedged ETS
discussion in the `RMdifTree` section.
