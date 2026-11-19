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

### Package function (built 2026-11-18, user: "build RMdifImpact() with tests")

Replaces the 2026-10-03 draft signature (bootstrap, iterations, parallel,
seed): every validated interval is a delta-method interval, so the
function is deterministic and needs no seed.

```r
RMdifImpact(data, dif_var, items, reference = NULL, cut = NULL,
            cut_scale = c("theta", "sum"), target = NULL, conf = 0.95,
            table = c("impact", "dif", "items", "local", "rates",
                      "projection"),
            output = c("kable", "dataframe", "ggplot", "list"))
```

- **Data.** Two groups. Complete cases only (the latent mean is fitted to
  the raw-score distribution and the shift CML has no missing-data
  patterns), incomplete rows dropped with a message and "complete cases"
  in the caption. Full-data category checks only: a category empty in one
  group is allowed, since the shift model shares the item's thresholds.
  At least one item must stay common.
- **Metric.** θ quantities (cut, latent means and SDs, range of the local
  bias, θ at the largest |b|, `target` means) are on the split model's
  scale with the reference group's thresholds centred to mean 0, the scale
  of the validation. The offset to the pooled (RMpersonParameters) scale is
  reported in the settings.
- **Engine.** Port of `dev/dif_shift_cml.R` (shift CML, analytic
  gradient), `dif_impact_theory.R` (WLE by score, equating, MML latent
  moments, local bias, population fits) and the step 3/3b rate and
  projection code, in `R/utils-dif-impact.R`.
- **Intervals.** Joint impact, local bias, projected impact: delta method
  on the shifts (step 2, step 3). Per-item impact: ratio SE (3b). Rates:
  fixed fitted rule, delta method on all split-model parameters with
  (μ, log σ) following through their Jacobian plus their own error
  (sample), or with the target moments fixed (projection) (3b).
- **Output.** `list` returns every table plus settings and the b(θ)
  curve. `kable` and `dataframe` return the table named by `table`
  (default the impact table: Δ keep, Δ split, impact with CI, % of Δ,
  SMD, impact / SE of Δ̂, per estimator). `ggplot` draws b(θ) with its
  pointwise CI and ±SEM(θ) reference lines.
- Remove is not reported (equivalent to split for group comparisons).

**Build notes (2026-11-18).**

- Files: `R/dif_impact.R` (function, tables, kable and ggplot),
  `R/utils-dif-impact.R` (engine), `tests/testthat/test-dif_impact.R`
  (117 expectations, about 15 s), `man/RMdifImpact.Rd`; NEWS, README and
  `_pkgdown.yml` entries. `R CMD check --no-tests --ignore-vignettes`: OK.
- **Agreement with the dev pipeline.** On one simulated dataset (8 items,
  3 categories, 500/500, DIF in I3 and I6), against `estimate_ext()`,
  `estimate_3b()` and `estimate_impact()`: impact, per-item impacts,
  local bias, rates (both scales, sample and projection), the rule and the
  projection agree to 1e-8 or better in the estimates and 5e-7 in the SEs
  (script in the session scratchpad, not kept). This was before the
  sum-rule change below, which moves the sum-scale split focal rule only.
- **Sum-score rule changed from the validation.** Dev rule: an equated
  focal score counts as at or above when `eq >= cut`. With DIF that makes
  the items easier for the focal group, every equated score lies just
  below its raw score, so any DIF, however small, moves the focal cut by
  a whole point (pcmdat2 with random groups, DIF -0.04: focal at or above
  .614 keep against .443 split). The simulations had DIF in the harder
  direction only (or opposite pairs), so they did not show it. Package
  rule: `eq >= cut - 0.5`, the same point that defines the true status.
  The intervals are conditional on whatever rule is fixed, so the 3b
  validation carries over. The Study 3 description of the sum-scale
  results in `dif_impact.qmd` reflects the dev rule. Noted there
  2026-11-19 ("The sum-score rule in RMdifImpact()", Signed rates): at the
  true parameters the package rule moves the focal cut in 8 of the 14
  simulated settings against 2 for the dev rule
  (`dif_impact_sum_rule_check.R` -> `sum_rule_slim.rds`). Coverage not
  rerun with the package rule.
- **Step 3c (written 2026-11-19, user asked).** `dif_impact_ext_sim_3c.R`
  reruns the sum-score rates with the package itself (`.dimp_estimate()`,
  so the package rule and the review fixes), all 42 cells, 500 reps,
  delta SEs only, per-replication L'Ecuyer streams (worker-independent),
  seeds 20261140 + cell. Targets: fitted rule's rates under the truth,
  sample and step 3 target population. Also records the package rule at
  the truth and whether each group's fitted rule differs keep vs split.
  Container smoke (3 cells x 8 reps, 4 workers): 18 s, no errors, no
  warnings. About 2.7 s per rep at k = 10.
  **Results (host 2026-11-19, 1.3 h, 0 errors, 0 warnings;
  `dif_impact_ext_sim_3c_analysis.R` -> `ext3c_results_slim.rds`).** Keep
  and split rates covered .916-.976 (median .950) sample, .920-.976 (.948)
  projected, against .922-.972 / .926-.970 for the 3b rule on the same
  cells; SE ratio median .99. Distinct rows below the 2-MCSE band: 19 of
  375 sample (3b sum: 16 of 291), 14 of 375 projected (3b: 4 of 291), the
  projected ones in correlated `prev` rows of a few k = 10 cells (SE ratio
  .92-.96). Keep - split where the group's rules differ: pooled .954 sample
  / .952 projected over 34 308 reps; per row (>= 100 reps) .924-.984, 3 of
  144 below .93 (all projected). Focal rule differs keep vs split in 71-95%
  (one item 0.8), 92-100% (two same), 2-14% (one item 0.5), <= 3.6%
  (opposite); reference never. Written up in `dif_impact.qmd`, "Sum-score
  rates with the package's rule" (third run).
- **Per-item impacts are not additive.** Their sum exceeded the joint
  impact by about 18% on the test data, the same direction as the 3b
  slope finding. Stated in the help page, tested as sign and size only.
- **Code review fixes (2026-11-18, `dev/code-review/dif_impact.md`).**
  The local bias summaries are now signed at the argmax, with the
  gradient at a fixed argmax, so the interval no longer runs below zero for
  an absolute value. The split rule uses the same WLE and equating solvers
  as the keep rule (the dev pipeline interpolated on a 0.02 grid; a rule
  can only differ for a cut within the interpolation error of a WLE). The
  SMD keeps the unweighted average of the group variances, as in every dev
  script, and the docs now say so. Also: convergence warnings for all main
  fits, `N_incomplete` in the settings, one grid for curve and table, one
  pooled CML fit, warm starts in the projection.
- **rasch_book chapter (draft, 2026-11-19).** `rasch_book/chapters/dif-impact.qmd`,
  in `_quarto.yml` under the DIF part. Reads `results/dif-impact/` (theory,
  sim, ext, ext3b, ext3c slim files). Scripts ported to
  `scripts/dif-impact/` (paths rewritten, internals attached from the
  installed easyRasch2 instead of load_all); all five smoke modes ran from
  the book root against the dev version. Study 3 sum-score rates come from
  3c only, the run history stays here. Audit found four overstated claims
  in `dif_impact.qmd`, corrected in the chapter and (user asked, same day)
  in the dev qmd, both computed inline:
  shift share of parameter variance "over 90%" (median 94%, range 37-111%),
  "remove is never less precise" (split up to 0.5% more precise), "remove
  loses 10 to 15% in RMSE" (4-17%), "keep and split have the same RMSE"
  (ratio 0.96-1.06).
- **Second code review (2026-11-19, `dev/code-review/dif_impact.md`, DI1-DI6).**
  Fixed: a split item at its lowest or highest possible scores in a
  group now stops by name (`.dimp_check_shifts()`). Target rates with
  default means/SDs now include their uncertainty, so they equal the
  sample rates and SEs. Theta ranges widen past +-10 when a threshold lies
  beyond +-5 (`.dimp_span()`). Non-PD MML Hessians are caught. Category
  errors on complete cases say so. Ordinary data give output identical to
  before, so 3c and the chapter stand. DI3 (user chose): the target rates
  keep the rule fitted to these data (as a scoring table used in another
  population), while the projected impact refits the pooled model. The
  code is unchanged and Details say so. The projected impact leaves out the
  uncertainty of default target moments, which adds SE 0.001-0.002
  against 0.066-0.070 on pcmdat2.
- Vignette (2026-11-19, re-knitted by the user, numbers checked against
  the output): "Size of the DIF" (q9) and "When the methods disagree"
  (q2 + q9), tree prose updated for the new RMdifTree defaults.
- Help example (2026-11-19): phq9 by gender, q9 (flagged by RMdifAnchor()),
  plus q2 + q9 for the per-item table; no eRm needed. Runs in about 5 s.
- Open: jamovi module, planned in
  `easyRasch2jmv/dev/PLAN-easyRasch2-next-release.md` (B3, and B4 for the
  RMdifTree defaults).

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

**IMPLEMENTED 2026-10-04 in 1.3.1.9001** (points 1 to 6, uncommitted):
`RMdifGamma()` p_value = NULL with the LD-style notices and captions,
`RMdifGammaCutoff()` 400/.95 + `dgp = c("conditional", "permutation")` +
`requested_iterations`/`dgp` in the return value + no iarm dependency,
`RMdifGammaPlot()` on the fast gamma, `.partgam_dif_gamma()` /
`.partgam_dif_one()` internals, mirai exports. Shared `.check_cutoff_sample()`
group-size message reworded (the null keeps observed groups); consumers drop
unused factor levels. Tests: new `test-dif-flagging-defaults.R`, updated
`test-partgam_dif.R`, DIF case added to `test-parallel-reproducibility.R`.
Full suite 1621 expectations, 0 failures; R CMD check (no tests/vignettes,
UTF-8) clean. NEWS merged into the dev section. Not done: vignette precompile
(DIF numbers in the shipped .Rmd are from the old null), jamovi module.

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

### Restarted 2026-11-17

**Purpose (user, 2026-11-17).** This is the step after `RMdifAnchor()` or
`RMdifTree()` has flagged DIF in one or more items. It answers whether the
DIF is large enough to need action, or whether it can be left as of no
practical relevance. It is a user-facing function. A what-if mode (DIF size
mapped to impact for this scale and this population) is part of the same
function or a sibling function, not just a dev tool.

**Decisions 2026-11-17.**

| Question | Decision |
|---|---|
| Estimators | All three: WLE group mean, latent group mean (CML item parameters fixed, MML group moments), raw sum score |
| What-if mode | User-facing |
| Work order | 1 theory + deterministic check, 2 dev simulation engine, 3 oracle vs estimated anchors, 4 detect-then-act pipelines, 5 package function |

**What changed since 2026-10-03.**

- Anchor selection is settled by the anchor studies (`dif_anchor_tests.qmd`,
  now the rasch_book chapter `dif-anchor-tests.qmd`) and `RMdifAnchor()`
  (committed 79d9bc3): the forward consensus anchor is the default, the
  centred anchor fails under unbalanced DIF. Study item 2 below is replaced
  by these. The "Anchors and the split model" section is kept for the
  reasoning, but its default (Gini alignment) is superseded.
- Part A question 6 is settled (conditional DGP, study 4).
- The large-sample fit thread was reframed as a comparison of scale-level
  actions (`misfit-consequence-comparison-note.md`). DIF is its cleanest
  case: the mechanism is known, keep / remove / split are concrete actions,
  and a group comparison is a concrete use.

**Estimand and counterfactual (2026-11-17).** The headline target is the
group difference Δ = μ_focal − μ_ref. It does not depend on how item
parameters are centred under keep and remove, since both groups are scored
with the same parameters. Every estimator has some bias without DIF (WLE
mean through extreme scores and short scales, sum score through
attenuation and the floor and ceiling), and that bias differs between k and
k − 1 items. The DIF impact of a strategy is therefore its estimate under
the DIF minus the same strategy's estimate on the same persons without DIF
(common random numbers in simulation, the d = 0 population in the
deterministic check). For comparison across estimators, the difference is
also given as a standardised mean difference (SMD) with the estimator's own
pooled within-group SD.

**Hypothesis to check (split vs remove).** For a between-group comparison,
splitting an item should recover little over removing it. The split item
cannot inform Δ (Δ is identified through the common items). It lowers each
person's SEM, which shrinks the within-group variance of θ̂ only slightly
when σ² dominates. And the estimation error in δ̂ is shared by every person
in the focal group, so it enters Δ̂ as an error that does not average out
over persons, scaled by the item's information share. Split may therefore
have a larger SE for Δ̂ than remove when the focal group is small, while
doing better for individual scores.

### Step 1: theory and deterministic check (2026-11-17)

Script `dev/dif_impact_theory.R`, results in `dev/dif_impact_theory/`.
Smoke run (6 settings) done in the container. **Full grid (216 settings)
not yet run, host command in the script header** (about 6 min on 10
workers).

**Method.** No simulation. Each CML fit (pooled for keep, pooled on k − 1
items for remove, uniform-shift for split, all with `.shift_negll()` from
`dif_shift_cml.R`) is made to the *expected* sufficient statistics: category
counts and raw-score distribution per group, from quadrature over the true
N(μ_g, 1). That gives the probability limit of every estimate at the given
group proportions. Then, by raw score:

- WLE: vectorised bisection on Warm's score function, the same function as
  `.theta_wle()`.
- Latent mean: MML for (μ, σ) with item parameters fixed. With Rasch-family
  items the marginal likelihood depends on the pattern only through the raw
  score, so it is fitted to the raw-score distribution.
- Sum score: raw score. Under split, the focal group's raw score is
  true-score equated to the reference items, T_ref(T_focal⁻¹(r)).

Impact = the strategy's Δ under DIF minus the same strategy's Δ at d = 0. SE
at given (n_ref, n_focal) = sampling (within-group variance of the
estimator / n_g, or the MML information for the latent mean) plus item
parameter error (delta method on the CML vcov from the expected information,
with `se_shift` the part from δ̂ alone). The two are assumed independent
because CML conditions on the raw scores. Step 2 checks that assumption.

**Checks (all pass).** Vectorised WLE equals `.theta_wle()` to 2e-11.
Population shift CML recovers d = 0.7 and every threshold to 3e-8. MML with
the true parameters recovers μ = 0.5, σ = 1. Equating with identical
parameters is the identity. One simulated dataset with 100,000 per group
(dichotomous, k = 10, d = 0.5, μ_focal = 0.5) gives keep Δ_WLE = 0.434,
the same as the population value. The latent-mean parameter SE is a
numerical derivative through an MML optimum and is accurate to about 10%.
It is about 3% of the total SE, so this does not matter.

**Grid.** (Changed 2026-11-17 after the first full run, user: d = 0.5, 0.8,
1.2 instead of 0.25, 0.5, 1. The results below are from the first run.)
Dichotomous (locations −2 to 2) and polytomous (4 categories,
locations −1.5 to 1.5, thresholds at location ±1) × k = 5, 10, 20 × DIF item
central or off-target (largest location, made harder) × d = 0.25, 0.5, 1 ×
μ_focal = 0, 0.5 × (n_ref, n_focal) = (500, 500), (1000, 200), (2000, 100).
σ = 1 in both groups.

**Smoke results** (central item, d = 0.5, μ_focal = 0.5, n = 1000/200):

| Scale | Keep impact, WLE (logits) | First-order −d × share | Keep bias / SE | 95% CI coverage, keep | SE(Δ̂) split / remove, WLE |
|---|---|---|---|---|---|
| dich k = 5 | −0.12 | −0.126 | −1.1 | .80–.81 | 0.122 / 0.119 |
| poly k = 5 | −0.12 | −0.114 | −1.3 | .74–.75 | 0.095 / 0.095 |
| dich k = 10 | −0.059 | −0.058 | −0.6 | .91 | 0.100 / 0.100 |
| poly k = 10 | −0.059 | −0.057 | −0.7 | .89 | 0.086 / 0.086 |
| dich k = 20 | −0.030 | −0.030 | −0.35 | .94 | 0.089 / 0.089 |
| poly k = 20 | −0.029 | −0.027 | −0.35 | .94 | 0.081 / 0.081 |

What this shows so far:

1. **Keep.** The bias in Δ is −d times the item's average information
   share. The first-order formula is within 0.01 logits, and with equal
   shares it does not depend on how much of d the pooled item location
   absorbs (`c_absorbed`: .50 at 500/500, .05 at 2000/100). The bias does
   not shrink with n, but the SE does, so the same DIF matters more in
   larger studies.
2. **Remove.** Impact exactly 0, as it must be with the correct anchors.
3. **Split.** Latent mean impact 0 to 1e-7 (correct model, consistent
   estimator). WLE within 0.003 logits, equated sum score within 0.02 SMD.
4. **Split vs remove (hypothesis).** Supported: the SE of Δ̂ is the same
   to two or three decimals, and split is slightly worse at k = 5
   (0.122 vs 0.119). Split's smaller sampling error is cancelled by the
   shared error in δ̂ (`se_shift` 0.05 at k = 5, n_focal = 200). The split
   item carries no information about Δ, so an efficient estimator gains
   nothing from it. For group comparisons, split and remove are equivalent.
   Split helps only for individual scores, which step 2 has to show.
5. **Estimator invariance.** Bias in SE units (and so CI coverage) is
   practically the same for WLE, latent mean and sum score. On the SMD scale
   all three have SE ≈ 0.078 here, the SE of an SMD at 1000/200. The latent
   SMD is larger because its SD excludes measurement error, and so is its
   SE. The consequence for a group comparison does not depend on the choice
   of estimator.

**Full grid run on the host 2026-11-17** (216 settings, 7.2 min,
`theory_results.rds`). All fits converged (max scaled gradient 5e-8), checks
reproduce. Findings over the grid:

1. **First-order formula.** Within 0.015 logits for central items up to
   d = 1, but it overstates the bias for off-target items: ratio
   first-order/exact up to 1.18 at d = 0.5 and 1.39 at d = 1. Use the
   exact calculation (cheap) in the function, the formula only for teaching.
2. **The group-mean bias depends on where the groups sit relative to the
   item.** Central items: keep impact nearly constant over n and μ_focal
   (poly k = 5, d = 0.5: −0.115 to −0.121). Off-target item (poly k = 5,
   d = 0.5): −0.074 with μ_focal = 0, −0.092 to −0.100 with μ_focal = 0.5,
   because the focal group is closer to the item and the item's
   information share is larger there. With μ_focal ≠ 0 the share differs
   between groups, so the absorbed fraction c (0.03 to 0.56 over the n
   settings) moves the bias a little. This is the group-level form of the
   local shift discussed below.
3. **Estimator invariance is approximate.** Bias in SE units, latent/WLE
   0.94 to 1.04, sum/WLE 0.87 to 1.11.
4. **Split vs remove for Δ.** For the latent mean the SE is the same to
   0.5% everywhere, as expected for an efficient estimator when the split
   item has no information about Δ. For the WLE mean and the sum score,
   split has up to 13% (WLE) and 7% (sum) LARGER SE at k = 5, 3% at
   k = 10, under 1% at k = 20. Remove is never worse for a group
   comparison.
5. **Coverage of the naive 95% CI under keep** (latent, poly,
   μ_focal = 0.5) is worst with equal groups, since Δ is most precise
   there: k = 5, central, d = 0.5 gives .63 (500/500), .74 (1000/200), .83
   (2000/100). d = 1 gives .08, .22, .47. At k = 20, d = 0.5 gives
   .93 to .94.
6. **Split residual.** Latent exactly 0. WLE at most 0.009 SMD and the
   equated sum score at most 0.030 SMD (k = 5), 0.009 (k = 10).

**Locality (discussion 2026-11-17).** The person-level shift is local, not
a constant: shift(θ) ≈ −d × I_j(θ)/I(θ). The information share of an
off-target item keeps rising beyond its location (poly k = 5, item at 1.5:
share .07 at θ = −1, .35 at 1.5, .49 at 3), because it is the only
informative item there. An item targeted at a high cut can therefore have a
large shift at the cut and a small group-mean bias. Consequences: the
function needs the shift curve and the shift at a user-supplied cut, not
only the group-mean bias. With raw-score cut-offs the classification
consequence is lumpy (a shift either moves the cut past a raw score or
not). Under keep, the reference group shifts by +c·d·share and the focal
group by −(1 − c)·d·share, so with a small focal group nearly all of the
person-level bias falls on the focal group (a fairness point for the help
page). Δ does not depend on c. Exact vs first-order focal shift, poly
k = 5, central, d = 1: −0.18 vs −0.275 at the top raw score.

**Implication for the function's output.** Impact in logits and SMD with a
CI answers "how large". The bias in SE units and the coverage of a nominal
95% CI for the group difference at the user's own n answer "does it matter
for my comparison" without a tolerance: they are consequences, not cutoffs.
They depend on n, which is correct for a consequence (practical
significance of bias grows with precision) and should be said in the help
page.

### Step 2: simulation engine (2026-11-17)

Script `dev/dif_impact_sim.R` (sources the step 1 helpers), results in
`dev/dif_impact_sim/`. Smoke run done in the container (4 cells, 100 reps,
1 min). **Full grid not yet run** (48 cells × 1000 reps: dich/poly × k 5,
10, 20 × central/off-target × d 0.5, 1 × n 1000/200 and 1000/500,
μ_focal = 0.5, cut at θ = 1; about 1.5 h on 10 workers, command in the
script header). 0.5 to 1.8 s per replication.

**Design.** Each replication draws a dataset with the DIF and a twin
without it from the same persons and the same uniform draws, so only the
DIF item's focal-group responses differ. Fits: keep (`.fit_cml_thresholds`),
remove, split (`shift_cml`). Each fit is linked to the true metric through
the DIF-free items (oracle anchor). Output per replication:

- Group level, for every strategy and estimator: Δ̂ in the DIF dataset and
  its twin, the paired impact, SMDs, and the naive SE a user would compute
  (Welch-type for WLE and sum score, MML information for the latent mean).
- Person level, WLE on the true metric, per group: RMSE and bias against
  true θ, mean SEM, the signed mean shift from the twin in SEM units, and
  classification disagreement and error at the cut. The absolute shift
  columns include response noise (the twin's DIF-item responses differ even
  under split), so only the signed mean and the RMSE are clean.
- `summarise_estimated_impact()`: the impact a real-data function can
  compute, Δ̂(keep) − Δ̂(split) or Δ̂(keep) − Δ̂(remove) within one
  dataset, against the true twin impact, and its distribution on the no-DIF
  twin (noise floor).

**Smoke results** (central item, d = 0.5, n = 1000/200, 100 reps):

1. Mean impact equals the step 1 population value within Monte Carlo
   error in all four cells (dich k = 5 keep WLE −0.122 ± 0.003 vs −0.121).
2. Empirical SE / step 1 SE 0.95 to 1.16, moving together across strategies
   within a cell (cell-level Monte Carlo noise at 100 reps, SD of an SD
   about 7%). To be confirmed at 1000 reps before the independence
   assumption of step 1 counts as checked.
3. The naive SE under split ignores the error in δ̂: naive/population 0.91
   at k = 5, coverage .93 instead of .95. Equal to 0.98 or more at k = 10.
4. Coverage of the naive CI under keep as step 1 predicted: .79 to .81 at
   k = 5 (poly .80 against .75 predicted, within 1.2 MCSE), .88 to .92 at
   k = 10.
5. **Person level: split has keep's precision and remove's lack of bias.**
   Focal RMSE against true θ, poly k = 5: split 0.63, keep 0.63, remove 0.71.
   Focal bias: split 0.003, keep −0.098, remove −0.001. Classification error
   at θ = 1 is about equal for split and keep and higher for remove at k = 5.
   So for group comparisons split and remove are equivalent (step 1), and
   for individual scores split is better.
6. **Estimated impact, the quantity the function will report.**
   Δ̂(keep) − Δ̂(split) recovers the true impact for WLE and latent mean
   (dich k = 5 WLE −0.120 vs −0.122, poly k = 10 −0.058 vs −0.059). Its SD
   equals step 1's `se_shift` (0.051 vs 0.050 at dich k = 5, 0.021 vs 0.022
   at dich k = 10). So the estimated impact is approximately
   −δ̂ × information share, with SE ≈ share × SE(δ̂), which gives an
   analytic CI without a bootstrap. Δ̂(keep) − Δ̂(remove) is worse: it
   carries the difference in estimator bias between k and k − 1 items (null
   mean 0.017 logits for WLE at dich k = 5, 0.03 SMD for the sum score). For
   the sum score, keep − split falls 0.01 to 0.02 SMD short of the true
   impact at k = 5 and 10, from the equating of the split item.
7. Noise floor: the estimated impact on the no-DIF twin has mean about 0 for
   the split reference and SD about equal to its SD under DIF. With d = 0.5
   and n_focal = 200, the estimate is beyond the no-DIF 95% band in about
   70% of datasets for dichotomous scales and in nearly all for polytomous.

### Redesign of the simulation, and local bias (agreed 2026-11-17)

**Agreed design.** The simulations validate the estimator (bias, SE, CI
coverage, effect of estimated anchors). General coverage of scale types
comes from the deterministic step 1, which includes k = 20 and dichotomous
items. Simulations are polytomous only (4 categories), cut at θ = 1, focal
mean 0.5, grid `grid_scenarios()` in `dif_impact_theory.R`, shared by both
scripts so every simulated cell has an exact target:

| Scenario | k | DIF items (location) | Position |
|---|---|---|---|
| one | 5 | 1 | central (0) or near cut (0.75) |
| one | 10 | 1 | central (0.17) or near cut (1.17) |
| two_same | 10 | 2, same direction | central (±0.17) or near cut (0.83, 1.17) |
| two_opposite | 10 | +d near cut (1.17), −d near −cut (−0.83) | |

× d = 0.5, 0.8 (changed from 0.4, 0.8 on 2026-11-17, user: 0.4 was flagged
too rarely for the flagged-set results) × n = 250/250, 500/500, 800/200
(42 cells). Sample sizes changed 2026-11-17 (user): 1000/200 became 800/200
(4:1, same total as 500/500), 250/250 added for small studies, and the two
"detected but small" cells (d = 0.3, 2500/2500) were dropped because step 1
gives that message exactly and the simulation cells only confirmed that
large samples give precise estimates. Steps 2
and 3 are merged in `dif_impact_sim.R`: each replication also runs
`RMdifAnchor()` (defaults) on the DIF dataset and its twin, and the
reported estimates are computed for the oracle items and for the flagged
items.

**Local bias (user, 2026-11-17).** A complement to the group-averaged bias:
how far apart the groups can be locally because of the DIF items.
Implemented in `local_bias()`:

- b(θ) = T_ref⁻¹(T_focal(θ)) − θ, the horizontal distance between the
  groups' test characteristic curves: where a focal person at θ lands when
  scored with the reference parameters. It is the between-group distance at
  equal θ and does not depend on c (to first order).
- gap(θ) = T_focal(θ) − T_ref(θ), the expected raw-score difference at equal
  θ. This is the local form of DTF (Chalmers, Counsell & Flora, 2016; verify
  the exact definitions before citing).
- Summaries over the central 95% of the pooled latent distribution: the
  range (b_lo, b_hi, signed), the unsigned maxima of |b|, |b|/SEM and |gap|,
  and b and gap at the cut. The range has to be restricted because an
  off-target item's share keeps rising into empty tails. For items near a
  high cut the maximum lies at the edge of that range (θ ≈ 2.1), so it
  depends on the range definition, and its location is reported.
- The signed value at the argmax is unusable when DIF goes both ways:
  |b_lo| and |b_hi| are similar and the sign flips between samples (smoke:
  mean −0.026 against a target of −0.070, coverage .68). The range and the
  unsigned maximum are what to report.
- `cut_flip()`: classification at the cut, keep vs split, by WLE. It is
  lumpy: at k = 10 one raw score (about 6% of the focal group) flips in every
  scenario including the central item, because the cut sits near a
  raw-score step. At k = 5 nothing flips. It is reported but should not be
  the headline.

Population values (poly k = 10, d = 0.4, n = 1000/200; from the smoke run before
the change to 0.5):

| Scenario | Group-mean bias | b range | b at cut |
|---|---|---|---|
| one, central | −0.047 | −0.051 to −0.026 | −0.051 |
| one, near cut | −0.043 | −0.080 to −0.011 | −0.055 |
| two_same, near cut | −0.089 | −0.145 to −0.025 | −0.111 |
| two_opposite | −0.012 | −0.067 to +0.057 | −0.033 |

**Step 1 rerun on the host 2026-11-17** (grid A with d = 0.5, 0.8, 1.2 and
grid B with d = 0.3, 0.5, 0.8; 246 settings, 8.2 min, all converged, max
scaled gradient 7e-8; `theory_results.rds` replaces the first run).

1. Earlier findings hold. Split/remove SE ratio for Δ: latent 0.995 to
   1.004, sum up to 1.07, WLE up to 1.14. First-order/exact for the keep
   bias: central 0.91 to 1.11, off-target up to 1.48 at d = 1.2.
2. **Max local bias versus group-mean bias depends on where the item
   sits.** Ratio max |b| / |group-mean bias|: central items 1.1 at every d
   and k. Off-target or near the cut 1.4 to 2.6 (grid A off-target, d = 0.5:
   dich k = 5 0.200 vs 0.087, poly k = 10 0.110 vs 0.049). Two opposite
   items 5.5 to 9.5: group-mean bias −0.010 to −0.026, local range about
   ±0.09 at d = 0.5 and ±0.15 at d = 0.8, with CI coverage for Δ still
   .94 to .95. **This is the case where the group mean alone would clear
   a scale that is biased locally.**
3. **For items off target or near a high cut, the maximum is at the edge
   of the occupied range** (θ = 2.08 at 1000/200, 2.27 at 500/500), so its
   value moves with the range: poly k = 10 near cut, d = 0.5, 0.100 vs 0.105.
   It has to be reported with its location, and the range definition is a
   design choice of the function.
4. **The raw-point gap is an item-level quantity.** Its maximum depends
   only on d and the item's categories (0.12 dichotomous and 0.39
   four-category at d = 0.5, 0.29 and 0.91 at d = 1.2), not on k or
   position. It equals max over θ of the item's own expected-score
   difference. It describes the DIF in raw points, not its impact on the
   scale. Scale dependence enters through b (logits, SEM units). Proposal:
   report it as the item's DIF in raw points, not as an impact measure.
5. **Classification at the cut is too coarse to be a headline.** At k = 10,
   about 6% of the focal group flips (raw cut 22 → 21) in almost every cell,
   for d = 0.5 and 0.8, central, near cut and two opposite items alike. At
   k = 5 nothing flips at d = 0.5, and 11% flips at d = 0.8. Two items in
   the same direction near the cut at d = 0.8 move it two raw scores (12%).
   It depends mainly on how close the cut sits to a raw-score step. The
   reference group never flips here.
6. **"Detected but small" is small in size but not for inference at large
   n.** k = 5, d = 0.3, 2500/2500: group-mean bias −0.072 logits, and the 95%
   CI for Δ covers .42. k = 10: −0.035, coverage .79. A bias that is small
   as an effect size still distorts the group comparison when n is large.
   The help page should separate the two readings (size of the bias, and
   its size relative to the study's own SE).

Grid B, latent mean, keep (selected rows):

| k | Scenario | d | n | Group-mean bias | b range | b at cut | Coverage |
|---|---|---|---|---|---|---|---|
| 5 | one, central | 0.5 | 1000/200 | −0.120 | −0.131 to −0.072 | −0.124 | .74 |
| 5 | one, near cut | 0.5 | 1000/200 | −0.118 | −0.160 to −0.037 | −0.144 | .74 |
| 10 | one, central | 0.5 | 1000/200 | −0.059 | −0.065 to −0.032 | −0.064 | .89 |
| 10 | one, near cut | 0.5 | 1000/200 | −0.054 | −0.100 to −0.013 | −0.067 | .90 |
| 10 | two same, near cut | 0.8 | 500/500 | −0.163 | −0.309 to −0.045 | −0.210 | .34 |
| 10 | two opposite | 0.5 | 500/500 | −0.010 | −0.090 to +0.071 | −0.042 | .95 |
| 10 | two opposite | 0.8 | 500/500 | −0.018 | −0.154 to +0.121 | −0.066 | .94 |

**Estimates the function would report** (`estimate_impact()`): split the
given items (all others common), impact = Δ̂(keep) − Δ̂(split) per
estimator, local bias from the split parameters over the estimated latent
mixture, SEs by the delta method on δ̂ only.

**Smoke run (3 cells, 40 reps, 2.2 min, with d = 0.4 before the change).**
Tentative at 40 reps.

1. Oracle items: estimates near the targets, delta-method SE close to the
   empirical SD, coverage .90 to .97. At k = 5 the impact is slightly
   attenuated (−0.087 vs −0.094, about 1.7 MCSE), to watch in the full run.
2. `RMdifAnchor()` cannot run in 7.5% of k = 5, 1000/200 replications: a
   category of the lowest item is unused in the focal group. Those
   replications keep their oracle results, and the flagged set is recorded
   as unavailable.
3. Detection: k = 5 near cut, d = 0.4, 1000/200 is flagged in 35%. Two
   opposite items at 500/500: both flagged in 30%, at least one in 85%. Then
   the estimate leaves out the cancelling item and coverage of the true
   total impact falls to .35. The function speaks only for the items it is
   given.
4. Selection: a false flag on the no-DIF twin gives a non-zero estimate
   with a CI excluding 0 (the item was flagged because δ̂ was large). The
   size is small (0.03 logits at k = 10), so the practical judgement still
   comes out right, but a nominal CI after selection is not valid. This is
   the selection caveat from 2026-10-03.
5. Detected but small (k = 10, d = 0.3, 2500/2500): flagged in every
   replication, impact −0.034 logits (SE 0.0045) against −0.035. This is the
   use case: a precise demonstration that a detected DIF is small.

**Smoke run 2 (2026-11-17, d = 0.5, new sample sizes, 40 reps).** No
replication failed. `RMdifAnchor()` could not run (category unused in one
group) in 7.5% at k = 5, 800/200, 2.5% at k = 10, 250/250, and 5% in an
extra check at k = 5, 250/250, d = 0.8. Detection at d = 0.5: k = 5 near
cut, 800/200, 68%. k = 10 central, 250/250, 85%. Two opposite items,
500/500: at least one 95%, both 65%. False flags on the no-DIF twin 2.5 to
5%.

**Host runs needed.** The simulation (42 cells × 1000 reps, about 4 to
4.5 h on 10 workers). The simulation computes its own population targets
per cell, so step 1 does not have to be rerun for it. The step 1 results
file still holds the earlier grid B (1000/200, 500/500, 2500/2500); a
rerun (about 8 min) would bring it in line. Commands in the script headers. The full
simulation file will be large (estimate rows: 30 × 1000 × 2 datasets × 2
sets × 10 quantities), so it is gitignored and a slim copy (summaries only)
will be written after the run. **Simulation run on the host 2026-11-17,
results below.**

### Step 2 results: full simulation (host, 2026-11-17)

Write-up for the rasch_book chapter: `dif_impact.qmd` (2026-11-17, renders to Typst PDF from the saved results, nothing recomputed). It reads Study 1 from `theory_results.rds` (general grid only) and Study 2 from `sim_results_slim.rds`, and recomputes the population b(θ) curves inline with a check against the simulation targets (agreement to 1e-11).

42 cells × 1000 replications, 2.8 h on 10 workers, no failed replication.
Full results in `dif_impact_sim/sim_results.rds` (41 MB, gitignored).
`dif_impact_sim_analysis.R` reads it and writes
`dif_impact_sim/sim_results_slim.rds` (150 KB, tracked) with every table
below: targets, strategies, split/remove ratios, person level, estimates,
flagged-set results by class, false flags on the no-DIF twin, detection.
Monte Carlo SE of a coverage near .95 is about .007.

**1. With the right items, the estimates work.** Oracle set (the true DIF
items split), DIF data, all 42 cells:

| Quantity | Max abs bias | SE ratio | Coverage |
|---|---|---|---|
| impact, latent | 0.002 | 0.90 to 1.02 | .92 to .96 |
| impact, WLE | 0.002 | 0.90 to 1.02 | .92 to .96 |
| impact, sum (raw points) | 0.011 | 0.85 to 1.03 | .90 to .96 |
| max abs local bias (b_maxabs) | 0.003, two opposite 0.011 | 0.95 to 1.11, two opposite 1.19 to 1.33 | .93 to .98 |
| same in SEM units (b_sem_maxabs) | 0.006, two opposite 0.027 | 0.95 to 1.36 | .93 to .98 |
| local bias at the cut (b_cut) | 0.003 | 0.88 to 1.06 | .92 to .97 |
| raw-point gap at the cut (gap_cut) | 0.014 | 0.86 to 1.11 | .90 to .97 |
| smallest local bias (b_hi, one direction) | 0.001 | **0.56 to 0.90** | **.72 to .92** |

SE ratio = mean delta-method SE / empirical SD. The lowest coverage for the
impact (.90 to .92) is at d = 0.8 with the item near the cut and equal
groups, where the delta method on δ̂ alone understates the SD by 10 to 15%.
The step 1 shortcut (shift error only) is good enough elsewhere.

The maximum of an absolute curve is biased upward under noise, and that
shows when the curve crosses zero (two opposite items): b_maxabs 0.003 to
0.011 too large against a target of 0.086 to 0.154, the raw-point gap 0.05
to 0.08 too large against 0.31 (d = 0.5) and 0.51 (d = 0.8). The SE is then too large as well,
so coverage stays at .96 to .98.

**2. Do not report the least-affected end of the range with a CI.** For DIF
in one direction, b_hi is the smallest bias, found at the edge of the range
far from the item. The delta method holds the range (from the estimated
mixture) fixed and perturbs only δ̂. At that edge b(θ) changes fast relative
to its small value, so the range uncertainty dominates and is missing.
b_lo, the largest bias, is fine (.93 to .97) even when it sits at the other
edge. Proposal: report max |b| with its θ, the value at the cut, and the
curve. Drop the lo/hi pair, except as the signed ends of the curve when
the sign changes.

**3. Split and remove are equivalent for the group comparison.** Ratio
split/remove of the empirical SE:

| Estimator | Native units | SMD | SE / abs(Δ) |
|---|---|---|---|
| latent | 0.99 to 1.00 | 0.99 to 1.00 | 0.99 to 1.00 |
| WLE | 0.99 to 1.01 | 1.01 to 1.05 | 0.99 to 1.00 |
| sum | 1.11 to 1.30 | 1.00 to 1.04 | 1.00 to 1.02 |

The raw-point ratio for the sum score compares scales of different length
(remove has one item less) and is not a fair comparison. The relative SE
(SE / |Δ|, which sets the power) is the same for both. Step 1 showed up to
13% for the WLE in SMD units, from dichotomous k = 5. Polytomous k = 5 gave
4.8% there and 4.6% here, so the two studies agree.

**4. Person level: the keep bias is small against measurement error.**
500/500, focal group:

| Scale | d | Keep bias | RMSE keep / split / remove | Misclassified at θ = 1, keep / split / remove |
|---|---|---|---|---|
| k = 5, central | 0.5 | −0.064 | 0.63 / 0.63 / 0.71 | .149 / .151 / .167 |
| k = 5, central | 0.8 | −0.101 | 0.62 / 0.63 / 0.71 | .148 / .158 / .167 |
| k = 10, two same, central | 0.8 | −0.096 | 0.44 / 0.44 / 0.51 | .109 / .113 / .128 |

The reference group gets the opposite share (k = 5, d = 0.5: +0.058). Keep
has the precision of split and remove loses 10 to 15% in RMSE. The total
misclassification hides the direction: DIF trades false negatives for
false positives in the focal group, and the total barely moves. A signed
measure (the focal group's false-negative rate at the cut) would show
the systematic part. It was not recorded and is a candidate for the
function's person-level output.

**5. Detection with `RMdifAnchor()` defaults** (DIF data, share of
replications):

| Scenario | d = 0.5: 250/250, 500/500, 800/200 | d = 0.8: same order |
|---|---|---|
| k = 5, one, central | .69, .95, .78 | .94, 1.00, .94 |
| k = 5, one, near cut | .57, .90, .67 | .91, 1.00, .91 |
| k = 10, one, central | .80, .99, .86 | .96, 1.00, .94 |
| k = 10, one, near cut | .55, .88, .64 | .88, 1.00, .93 |
| k = 10, two same, central (both) | .66, .98, .79 | .97, 1.00, .94 |
| k = 10, two same, near cut (both) | .33, .82, .53 | .84, 1.00, .91 |
| k = 10, two opposite (both) | .14, .65, .15 | .56, .98, .56 |

At least one false flag: 5 to 11% of DIF replications at k = 5 (artificial
DIF with four anchors), 3 to 6% at k = 10. On the no-DIF twin any flag
in 3.2 to 6.6% (FWER .05, the top value about 2 MCSE above).
`RMdifAnchor()` could not run (category unused in one group) in 0 to 8% of
replications, and in 15.5% for two opposite items, d = 0.8, 800/200.

**6. With the flagged items, what goes wrong is the item set, not the
estimate.** Flagged replications on the DIF data, split by what was
flagged:

| Flagged set | Share of reps | Impact (latent) | b_maxabs |
|---|---|---|---|
| exact (all true items, nothing else) | .13 to .97 | bias up to 15% (winner's curse at the lowest power), coverage .91 to 1.00 | up to 14% too large, coverage .94 to .99 |
| a true item missed, no false flag | up to .58 (two items, d = 0.5) | coverage .00 to .04 | two same: about half the target, coverage .00 to .06. Two opposite: .64 to .89 |
| at least one false flag | .03 to .11 | k = 5: pulled toward 0 (by up to 120% of the impact), coverage .12 to .81 | coverage .30 to .97 |

- **Missed item.** The estimate covers only the items it is given. For two
  items in the same direction it gives about half the impact. For two
  opposite items it reports −0.02 to −0.08 instead of the true −0.01 to
  −0.03 (the near-cancellation). The local maximum is then about right, because
  the two items act in different θ regions.
- **False flag at k = 5.** The falsely flagged item carries apparent DIF in
  the opposite direction (the artificial-DIF mechanism with few anchors).
  Splitting it as well cancels part of the real impact.
- **Exact set at low power.** The estimate is slightly too large in the
  cells where detection is near 55%, the usual post-selection inflation.

**7. False flags on the no-DIF twin.** In 3 to 7% of replications
`RMdifAnchor()` flags an item with no DIF. The impact estimate is then:

| k | Mean abs impact, latent | 95th percentile | CI excludes 0 | True impact of one item at d = 0.5 |
|---|---|---|---|---|
| 5 | 0.042 to 0.068 | 0.065 to 0.141 | 41 to 76% | 0.114 to 0.121 |
| 10 | 0.023 to 0.038 | 0.034 to 0.061 | 70 to 97% | 0.050 to 0.059 |

The CI excludes 0 because the item was flagged for a large δ̂, and the
impact is about share × δ̂. The impact CI is a rescaled DIF test and adds no
evidence that the DIF exists. What it adds is the size. A false flag at
k = 10 gives about half the impact of a real d = 0.5 item, and the upper
tail at k = 5 reaches the size of a real one. So a judgement of practical
relevance has to rest on the size compared with a stated threshold, not on
whether the CI excludes 0.

**Summary for the design.**

- The engine is fit for purpose. With a correct item set, the impact (all
  three estimators), max |b| and the values at the cut are unbiased with
  delta-method CIs near nominal. A bootstrap is not needed for these.
- The weak point is upstream: which items are split. Missed items are the
  main failure (two items at d = 0.5 and small samples), and the function
  cannot see them. Options to discuss: (a) state it in the output and help
  page only, (b) a sensitivity row that splits every item flagged at an
  unadjusted level as well, (c) report per-item impacts, so the user sees
  how much each split item contributes.
- Relevance has to be judged against a stated size. Candidates: the impact
  relative to the observed Δ, to the pooled SD (SMD), or max |b| in SEM
  units, with a user-supplied smallest relevant value.
- Report max |b| with its θ location and the value at the cut. Drop b_hi.
- Person level: report the focal false-negative and false-positive rates at
  a cut under keep and split (signed), not only total misclassification.

**Consequences for the function design (proposed, not decided).**

- Estimate impact as keep − split, not keep − remove. Split is the
  reference model in any case (it defines the anchor through the common
  items).
- The analytic route (δ̂, its SE and the information share from the fitted
  model and observed score distributions) may be enough for the CI. The
  parametric bootstrap would then be needed only for the noise floor and
  for checking the analytic CI. This would make the function fast.
- Report remove next to split for what the user would actually do, with the
  equivalence for group comparisons stated in the help page, and the
  person-level advantage of split shown.

### Decisions after the simulation (2026-11-17)

| Question | Decision |
|---|---|
| Missed items | (c) per-item impacts: each split item's own keep − split contribution, next to the joint impact |
| Relevance scale | Report all three (impact relative to the observed Δ, as SMD, max \|b\| in SEM units). No user threshold, no verdict. Agrees with "no default tolerances" from 2026-10-03 |
| b_hi | Dropped (user, 2026-11-17, after the explanation below). Report max \|b\| with its θ, the value at the cut, and the curve |
| Projection | In v1 if it is easy (user): `target` with n, mean and SD per group |
| Person level | Signed rates at a user-supplied cut, on the θ or the sum-score scale: focal persons pushed across the cut by the DIF, in each direction, keep vs split |

**Why b_hi under-covers (for the open decision).** b_hi is the signed
maximum of b(θ) over the central 95% of the estimated pooled distribution.
With DIF in one direction it is the least-affected point, and it lies at
the edge of that range, far from the DIF item. Write b(θ) ≈ d·f(θ). The
delta method gives SE ≈ f(θ_edge)·SE(δ̂) and leaves out the term from the
estimated edge, d·f′(θ_edge)·SE(θ_edge). The ratio of the missing to the
included part is about d·[f′/f](θ_edge)·SE(θ_edge)/SE(δ̂). It is large at
the far edge, where f is small and changing fast, and it grows with d.
This fits the simulation: coverage .72–.76 at d = 0.8 near the cut and
.86–.92 at d = 0.5. It does not improve with n (250/250 vs 500/500),
because both parts shrink at the same rate, and it is best at 800/200,
where the pooled range is set mainly by the larger, more precisely
estimated reference group. max |b| is unaffected when it is an interior
maximum (b′ = 0 there) and robust at an edge where f is large, which is
why b_lo covers .93–.97.

**Projection to a target population.** Use case (user): the scale is
tested in n ≈ 600–1000 but intended for a survey with 50 000 respondents.
Of the outputs, b(θ) depends only on item parameters and d. The keep bias
of Δ and the rates at a cut depend on the latent distributions, not on n.
What depends on n is the SE of Δ in the target study, so bias/SE, and the
number of persons pushed across a cut. A projection is therefore a
population calculation (Study 1 machinery) with the estimated split-model
thresholds and user-supplied group sizes, means and SDs: a population CML
keep fit to the expected statistics, Δ(keep) − Δ(split), the sampling SE
of Δ at the target n, and the rates at the cut. Its uncertainty comes from
the calibration sample and does not shrink with the target n. It assumes
that d and the item parameters transfer to the target population.

**Definitions fixed for step 3 (2026-11-17).**

- Per-item impact of item j: Δ(split S \ j) − Δ(split S), i.e. keep item j
  common while the other DIF items stay split. SE by the delta method on
  δ̂_j: |∂Δ(split)/∂δ_j| × SE(δ̂_j). The parts add up to about the joint
  impact (first order).
- Rates at a cut, per group and strategy: prev (share classified at or
  above the cut), FN (true status above, classified below) and FP, as
  shares of the group, plus keep − split. Model-based: response model and
  latent N(μ, σ) per group from the split model. θ-scale cut: WLE ≥ cut,
  keep WLEs moved onto the split metric through the common items.
  Sum-score cut: raw ≥ cut, focal raw scores true-score equated under
  split, true status θ ≥ T_ref⁻¹(cut − 0.5).
- **Classification is a step function of the parameters**: a raw score
  either crosses the cut or not, so keep − split is zero or the mass of one
  or two raw scores, and the delta method does not apply. At k = 10 the θ
  cut of 1 lies 0.03 logits above the reference WLE of raw score 21, so
  which raw score is the cut depends on threshold noise. Intervals are
  therefore percentiles of 300 draws of all split-model parameters (CML
  covariance) and of each group's (μ, log σ) (linearised in the item
  parameters plus own sampling error). The keep WLE table follows the
  reference table in the draws and is linearised in the shifts through
  population keep fits. A first version drew only δ and (μ, σ) and missed
  the threshold noise (smoke test, two-opposite cell).
- A parametric bootstrap would be the alternative for the package (one
  dataset, refits are affordable there), but it cannot be validated by
  simulation overnight.

**Step 3 run** (`dev/dif_impact_ext_sim.R`, results in
`dev/dif_impact_ext_sim/`): validates the per-item impacts, the rates
(sample and projection) and the projected impact. 42 cells × 1000
replications, oracle items, DIF data only. Projection target: n = 40 000 /
10 000, reference N(0, 1), focal N(0.25, 1.25²). Smoke test in the
container (3 cells, 4–12 replications): no errors, interpolated WLEs and
equating within 5e-5 of the exact helpers, about 4 s per replication.
Handed to the host 2026-11-17 for an overnight run.

### Step 3 results (host, 2026-11-18)

Run: cells 1–11 on 10 workers, 12–42 on 6 (restart after a stop, finished
cells read back), 4.9 h for the second part, 0 errors. Each cell has a
fixed seed, but mclapply streams depend on the worker count, so a rerun on
one worker count does not reproduce the file bit for bit. Analysis:
`dev/dif_impact_ext_sim_analysis.R`, tables in
`dev/dif_impact_ext_sim/ext_results_slim.rds`.

| Quantity | Coverage (median) | SE ratio (median) |
|---|---|---|
| Joint impact, sample (3 estimators) | .904–.964 (.94) | .86–1.04 (.96) |
| Projected impact, target n | .910–.960 (.94) | .87–1.03 (.97) |
| Projected bias/SE | .911–.960 (.94) | .88–1.03 (.97) |
| Per-item impact, stored SE | .862–.944 (.90) | .77–.98 (.85) |
| Per-item impact, ratio SE (post hoc) | .926–.982 (.95) | .94–1.19 (.99) |
| Rates, sum-score cut (keep, split) | .921–.984 (.94) | |
| Rates, θ cut, k = 5 | .864–1.000 (.95) | |
| Rates, θ cut, k = 10 | .135–1.000 (.95) | |

**Per-item SE.** The stored SE |∂Δ(split)/∂δ_j| × SE(δ̂_j) is about 15%
too small. Regressing the per-item estimate on both δ̂ (cells 28, 34, 40):
the slope on δ̂_j is −0.136, against −0.119 for the joint impact's
gradient in the same item. Making item j common moves Δ by more than the
local gradient of the split model predicts. The estimate is close to
linear through the origin in δ̂_j (intercept about 0.003, R² .93–.97), so
|impact_j / δ̂_j| × SE(δ̂_j) is the slope of the actual per-item
estimator. Post hoc, with the empirical SD of δ̂_j standing in for its SE,
this gives coverage .926–.982. Needs confirming with the CML SE in a
rerun.

**θ-scale rates: knife-edge.** The undercoverage sits where a population
WLE table point lies on the cut. Cell 34 (k10, two same, central, d 0.8):
focal split WLE of raw 20 is 0.999 against a cut of 1. The sample estimate
averages 1.007 (SD 0.035, 60 replications), so raw 20 counts as above the
cut in about half the samples and below it in the population. Cell 14: keep
WLE of raw 21 is 0.998. Minimum coverage of the θ rates by distance from
the cut to the nearest population WLE (sample context): ≤ 0.01, .64 (median
.76); 0.01–0.05, .81 (.95); above 0.05, .93 (.95). The sum-score rates do
not have the problem because the rule there is a fixed integer. The draws
cannot repair it: the estimand is a step function of the parameters, and
when the step lies within about one SE of the estimate the plug-in sits on
the wrong side of it in a large share of samples.

**Proposal: rates conditional on the fitted classification rule.** The
user will screen with the estimated tables, not the population ones. So the
rule (which raw scores are classified at or above the cut, per group and
strategy) is taken from the fit and held fixed, and the rates of that rule
are estimated: prevalence, FN and FP against true θ ≥ cut, in the sample's
latent distributions or a target population. These are smooth in the
parameters, with no step. The θ cut becomes a raw-score cut per group and
strategy, so the θ and sum-score cuts become one mechanism, and the output
can state the rule directly (for example, "splitting moves the focal cut
from raw 21 to 20"). The keep WLE linearisation (population keep fits for
Jw) is no longer needed. Not adopted yet. A rerun would confirm the rates
and the per-item ratio SE (ask first).

**Step 3b (agreed 2026-11-18, option D).** Rates conditional on the fitted
rule adopted for validation. `dev/dif_impact_ext_sim_3b.R`: fixed-rule
rates (sample and projection) with both a delta-method SE (central
differences in all split-model parameters, μ and σ following through their
Jacobian, plus their own error) and the 300 draws, per-item impacts with the
ratio SE, and the fitted rule as a raw cut against the population rule. 42
cells × 500 replications, seeds 20261130 + cell. Targets per replication:
the fitted rule's rates under the true parameters. Container check (48
replications, the two knife-edge cells 14 and 34): delta and draws SEs
agree within 1–2%, keep/split rate coverage .88–1.00 (median .96), while
the fitted θ rule matched the population rule in only 56–60% of samples.
Per-item ratio SE: coverage .94–1.00. Handed to the host 2026-11-18,
expected 1.5 h on 10 workers.

**Write-up, 2026-11-18.** `dif_impact.qmd` gets a "Size on the reference
scales" part in Study 1, next to the group-mean bias figure (user: the
central question for a practical user). Table and figure of the keep bias
as share of Δ, SMD, largest |b| / SEM and |bias| / SE at 500/500, from
the Study 1 results. Implications point 3 now names the agreed scales
(it said "the study's own SE" instead of the share of Δ). Key reading: at
five items and d = 1.2 a central item removes 57% of Δ = 0.5 (4.0 SE at
500/500) but moves a person by at most 0.53 SEM. From 5 to 20 items the
share of Δ falls fourfold and the SEM-scaled local bias only by half.

**Notation (user, 2026-11-18).** No "d" for the DIF size in anything a user
reads: psychology readers take d for Cohen's d. The qmd writes "DIF = 0.5"
in text and $\text{DIF}$, $\widehat{\text{DIF}}$ in formulas. The package
function, its output columns, plots and help page follow the same rule
(for example a column `dif`, axis "DIF size (logits)"). This log keeps d
as internal shorthand.

### Step 3b results (host, 2026-11-18)

Run: 42 cells × 500 replications on 10 workers, 0 errors, all 500
replications in every cell (cells 1–33 in a first session, 34–42 after a
restart, about 2.6 min per cell). Analysis:
`dev/dif_impact_ext_sim_3b_analysis.R`, tables in
`dev/dif_impact_ext_sim/ext3b_results_slim.rds`. Monte Carlo SE of a .95
coverage: .0097. 468 of 3132 rate rows are degenerate and dropped (error SD
zero, or estimate = target = 0 throughout): keep − split for the reference
group on the sum scale is zero by construction, and on the θ scale it is
zero whenever the two fitted rules agree.

| Quantity | Delta coverage (median) | Delta SE ratio (median) |
|---|---|---|
| Rates, sum cut, keep/split, sample | .922–.972 (.950) | .91–1.09 (1.00) |
| Rates, θ cut, keep/split, sample | .920–.978 (.950) | .86–1.16 (1.00) |
| Rates, sum cut, keep/split, projection | .926–.970 (.948) | .90–1.17 (.99) |
| Rates, θ cut, keep/split, projection | .924–.972 (.950) | .84–1.16 (1.00) |
| Keep − split, rules differ, sum (pooled) | .953 sample, .958 projection | |
| Keep − split, rules differ, θ (pooled) | .956 sample, .966 projection | |
| Per-item impact, ratio SE, CML covariance | .912–.982 (.948) | .90–1.15 (.99) |

**Rates.** Holding the fitted rule fixed removes the knife-edge. The
step 3 failures (θ rates, k10, down to .135) are gone: cells 14 and 34 now
cover .932–.968. Rows at or below .931 (two Monte Carlo SE): 57 of 1008 in
the sample context against about 25 expected, 19 of 1008 in projection.
Many rows are duplicates (keep and split share values where the rules
agree), so the excess is smaller than it looks. A slight undercoverage in
the sample context, plausibly from the (μ, σ) part, not large enough to act
on.

**Delta vs draws.** Mean delta SE / mean draws SE: .962–1.030 (median
1.001) over all 2556 non-degenerate rate rows. Draws percentile coverage
runs about .005 below delta (medians .942–.944). The package can use the
delta method alone and drop the 300 draws.

**Keep − split.** Its distribution is a mixture: zero (with SE
zero) when the two rules classify the same raw scores, otherwise the mass
of one or two score groups. Unconditional coverage (.934–1.00) is inflated
by the exact zeros. Over the replications where the rules differ, pooled
coverage is .953–.966. Per cell (≥ 100 such replications) most rows cover
.93–1.00, but 24 of 390 fall below .93, the lowest .840 (cell 31, k10, two
same, central, DIF 0.5, 250/250, 106 replications), with SE ratios above
1 in some of the same rows. That points to a skewed, non-normal error when
the rules differ by one score group. The difference is the least reliable
of the reported intervals. Worth a remark in the help page, no fix needed.

**Per-item impacts.** The ratio SE with the CML covariance confirms the
post hoc check: .912–.982 (median .948). Lowest: cell 36 (k10, two same,
central, DIF 0.8, 800/200), I06, .912 with SE ratio .90. Adopt.

**Rule agreement.** The fitted θ rule equals the population rule in .44–1.00
of samples per cell (median .85 keep, .91 split focal, .94 split
reference). The sum rule is fixed by construction except for the split
focal equating, .83–1.00. So the rule itself is uncertain, and the rates
are conditional on it. The output has to state the rule (the raw-score cut
per group and strategy) next to the rates.

**For `RMdifImpact()`.** Fixed-rule rates at a θ or sum cut, delta-method
SEs, no draws. Rule shown as raw-score cuts per group and strategy.
Per-item impacts with the ratio SE. Keep − split with its SE, and a note
that it is zero whenever both strategies give the same rule.

**Write-up (2026-11-18, user asked).** `dif_impact.qmd` has a new Study 3
(design, the knife-edge and the fixed rule, coverage table for both runs,
signed rates at θ = 1 for 500/500, sum-scale note) and updated Purpose,
Study 2 person paragraph, Implications 6 and 7, Limitations and Run notes.
For the record: the step 3 θ-rate minimum of .135 is in the projection
context. In the sample the step 3 minimum was .636.
Same day (user): long results tables become figures, short tables stay.
Study 2 estimates (oracle) and flagged-set tables are now figures, and
they and the Study 1 scales table (shown by its figure) moved to
"Appendix: tables" with the Study 3 tables.

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
| 1 | First-order shift formula for PCM items, and how close the analytic group-mean bias is to simulated | deterministic, cheap. **2026-11-17: done as step 1, run on the host** |
| 2 | Anchor selection: how often alignment picks a DIF item, by proportion and balance of DIF items, k, group sizes | **superseded 2026-11-17** by the anchor studies and `RMdifAnchor()` |
| 3 | Noise floor: impact under no DIF by size of the smaller group, and whether a correction is needed | simulation. **2026-11-17: answered by the step 2 simulation (results 1 and 7): no correction needed with a correct item set, false flags give 0.02 to 0.14 logits** |
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
