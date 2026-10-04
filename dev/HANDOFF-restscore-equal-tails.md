# Handoff: an equal-tailed family-wise test for item-restscore

First written 2026-10-02 after the first equal-tailed step-down was built,
tested and withdrawn. Updated 2026-10-03 after the split one-sided
step-down was validated on synthetic nulls, the real restscore null, power
and seed stability. **Nothing here is in the package.** What remains is a
decision, not a research question (see "Open decisions").

## The problem

`RMitemRestscore(cutoff = RMitemRestscoreCutoff(...))` tests each item's
observed minus expected gamma (`Difference`) against a parametric-bootstrap
null. The statistic is studentised by the bootstrap mean and SD, and the
test is two-sided on |t|: marginal p = (1 + #{|t*| >= |t|}) / (B + 1), and
Westfall-Young (WY) studentised-max step-down for the family-wise rate.

Gamma is bounded at 1, so the null of `Difference` is **left-skewed**
(mean per-item skewness −0.28 for 20 dichotomous items at n = 150 and
mistargeted by 1.5 logits, −0.05 for 9 polytomous items at n = 1000).
Folding |t| at the mean therefore gives the long (underfit) tail more than
α/2:

- Marginal, per item, under the null (`restscore_cutoff_validation.qmd`,
  conditional DGP, 500 datasets per cell): underfit 2.6 to 3.4 percent,
  overfit 1.8 to 2.7. Per-item skew correlates −0.79 with the
  underfit-minus-overfit rate.
- Family-wise, over the 2000 null datasets: at least one underfit flag in
  3.70 percent of datasets, at least one overfit flag in 1.55 percent. The
  total family-wise rate is nominal (5.0).

The goal was a family-wise test that keeps the nominal rate **and** splits
it evenly between the directions.

## Summary of the answer

**Two one-sided WY step-downs with pooled moments (`split_t_pooled`).**
Studentise each item by the mean and SD of all B + 1 values (observed plus
simulated), run WY on t for overfit and on −t for underfit, each at α/2,
and flag an item if either run rejects.

- **Level and balance:** 5.0 percent family-wise on the real null, 2.40
  underfit / 2.75 overfit (|t|: 3.70 / 1.55). Exact under exchangeability
  by construction.
- **Power:** 2.9 and 2.2 points less underfit detection than |t| (SE 0.4),
  no measurable overfit gain (+0.2 and +0.5, SE 0.4).
- **Other errors:** fewer false underfit flags when items overfit (4.7
  against 6.2 percent with two overfitting items), 3 to 4 points more clean
  decisions in the overfit cells.
- **Seed stability:** the same as |t| at B = 400 (10.3 against 10.0
  percent of seed pairs disagreeing), better at B = 1000 (4.6 against 6.1).

## History: what failed (2026-10-01 to 10-02)

All in `restscore_eti_null.qmd` (real null run, 2000 datasets, plus a cached
chunk of synthetic nulls in `restscore_eti_null/synthetic.rds`). The
synthetic nulls are exchangeable by construction: B + 1 iid rows, the first
taken as observed, columns `-exp(s * Z)` with `Z` normal and correlated
`rho`; `s = 0.5` (strong skew) and `s = 0.15` (milder, still stronger than
restscore's). 4000 replications per condition, B = 400, MCSE about 0.35.

| variant | marginal | WY family-wise rate | WY tail balance |
|---|---|---|---|
| \|t\| (shipped) | fine overall, unequal tails | nominal | lopsided toward the long tail (e.g. 364/6 low/high, 20 correlated items, s = .15) |
| equal-tailed marginal, 2·min(one-sided p) | **works**: 2.2 to 2.6 percent per tail in the real null | n/a | n/a |
| `first`: normal scores, sims ranked among sims, obs against sims | | **liberal**: 8.65 (k = 9), 9.23 (k = 20); real null 8.85 vs 5.0 for \|t\| | balanced |
| `pooled`: obs + sims ranked together, `qnorm((r − .5)/(B + 1))` | | 3.5 to 4.65 | all upper tail at k = 20 (0/193): floating-point artefact |
| `pooled_sym`: as pooled, |z| from the rank counted from the nearer end | | **0.00 at k = 20**: cannot reject | |
| kernel-smoothed pooled CDF, then `qnorm()` | | about nominal | very lopsided (204/1) |
| `two_piece`: (x − median) / (median − q.025) below, / (q.975 − median) above, from the sims | | 4.8 to 5.9 (slightly liberal at mild skew) | lopsided, but **less than \|t\|** at s = .15 (267/110 vs 364/6) |

Why each failed:

1. **`first`, exchangeability.** Simulated values can never fall outside the
   simulated range; the observed value can. An observed value beyond every
   simulated value in its column beats every simulated maximum, probability
   about 2/(B + 1) per item, so roughly 2k/(B + 1) excess family-wise.
2. **`pooled`, floating point.** `qnorm(400.5/401)` and `-qnorm(0.5/401)`
   differ in the last bits. With many tied maxima, `>=` then resolved every
   tie toward the upper tail.
3. **`pooled_sym`, resolution.** Any rank-based (min-p) step-down needs
   B well above 2k/α. With pooled ranks each item has one value at each end
   of its B + 1 values, so about 2k simulated rows tie at the largest |z|,
   and the smallest attainable adjusted p is about 2k/(B + 1), near .10 for
   k = 20 at B = 400. 800 iterations would be the floor for 20 items, before
   any margin.
4. **Kernel**, fixed bandwidth stretches the long tail again.
5. **`two_piece`**, matching the tails at the 2.5 and 97.5 percent points does
   not match them at the much more extreme quantiles where WY decides
   (roughly α/(2k) per tail).

General lesson: equal tails at the family-wise level means matching the two
tails at extreme quantiles, which either needs many more draws (rank-based)
or a continuous statistic calibrated separately per direction.

## What worked (2026-10-02 to 10-03)

Each step was run in the order the first version of this file prescribed:
synthetic null, real null, power, then (added) seed stability. Every run
reproduced the shipped |t| from the package with zero difference before any
rate was read.

### 1. Synthetic nulls (`restscore_eti_null.qmd`, "Follow-up" sections)

Caches `restscore_eti_null/synthetic_split.rds` and `synthetic_pooled.rds`.
16 conditions (the earlier 12, plus a `mixed` arm whose skew rises from 0.02
to 0.3 across items), 8000 replications in two seed blocks.

- **Splitting balances the directions.** `split_t` gave 2.70 / 2.68 percent
  of replications with a lower / upper flag (12 cells with k ≥ 9), against
  4.94 / 0.09 for |t|.
- **Self-influence makes simulated-only moments liberal.** The mean and SD
  come from the simulated values only, so an extreme simulated value
  enlarges its own scale, which the observed value never does to itself.
  Paired over 96,000 replications, simulated-only versions flagged 331
  (|t|), 368 (`split_t`) and 636 (`split_two_piece`) replications that the
  pooled versions did not, against 0, 0 and 27. Size about 0.35 points,
  roughly (t² − 1)/(2B) of scale. **This also affects the shipped |t|.**
- **Pooled moments make the step-down exact** under exchangeability:
  `split_t_pooled` 2.48 / 2.48 lower and 2.50 / 2.36 upper in the two
  blocks, against 2.49 by construction.
- **Within a direction the imbalance moves to the items.** In the mixed arm
  `split_t_pooled` gave the most skewed third of the items 5.9 times the
  underfit flags of the least skewed third. Two-piece scaling cut that to
  2.1.

### 2. Real restscore null (`restscore_eti_null.qmd`, last section)

Cache `restscore_split_null/`, including the observed differences and the
simulated matrices of all 2000 datasets (`<cell>_sims.rds`), so further
variants can be checked offline.

| variant | family-wise | lower / upper |
|---|---|---|
| \|t\| (shipped) | 5.0 | 3.70 / 1.55 |
| `split_t_pooled` | 5.0 | 2.40 / 2.75 |
| `split_two_piece_pooled` | 4.9 | 2.40 / 2.55 |
| `split_t` (simulated-only) | 5.5 | |
| `split_two_piece` (simulated-only) | 6.0 | |

`dich_tgt_300` was high for every variant (6.4 to 8.0), as for |t| in the
validation study: a bootstrap issue in that cell, not the step-down. The
items' null skewness barely varies within the real cells, so the item
imbalance that two-piece scaling corrects is small in practice.

### 3. Power (`restscore_power_conditional.qmd`, last section)

Caches `restscore_power_split/` (100 datasets per cell) and
`restscore_power_split_ext/` (datasets 101 to 300), both with the simulated
matrices. Figures from 300 datasets per cell, 1800 misfitting items per
direction and targeting, `split_t_pooled` against |t|:

| | difference (SE) | discordant items |
|---|---|---|
| underfit, at −2 | −2.9 (0.4) | 1 / 54 |
| underfit, at 0 | −2.2 (0.4) | 0 / 40 |
| overfit, at −2 | +0.2 (0.4) | 25 / 21 |
| overfit, at 0 | +0.5 (0.4) | 30 / 21 |

- False underfit flags with overfit planted: 1.0 and 4.7 percent (one, two
  items), against 1.6 and 6.2 for |t|.
- Clean decisions in the overfit cells: 55/46/18/12 against 51/42/14/9.
  Underfit cells level or 1 to 3 points behind.
- `split_two_piece_pooled` was behind `split_t_pooled` in three of four
  comparisons (up to 1.6 points). **Dropped.**
- Pooling costs 0.4 to 0.8 points, for |t| and the split alike.

### 4. Seed stability (`restscore_stability.qmd`)

Cache `restscore_stability/`. 5 cells × 20 datasets × 5 seeds × B = 400
and 1000. Share of seed pairs disagreeing on at least one item:

| method | B = 400 | B = 1000 |
|---|---|---|
| WY (shipped) | 10.0 | 6.1 |
| WY, pooled | 9.9 | 6.3 |
| `split_pooled` | 10.3 | 4.6 |
| BH | 13.5 | 11.4 |

The split test's more extreme critical value (about the 10th largest of 400
maxima instead of the 20th) did not make it less stable.

## Open decisions (the user's)

1. **Ship `split_t_pooled` in `RMitemRestscore()`?** The trade: equal tails
   under the null and fewer false underfit flags when items overfit, against
   2 to 3 points of underfit detection and no overfit gain. The simulations
   do not settle which error matters more.
2. **Pooled moments in `.bootstrap_pvalues()` for every cutoff function?**
   Self-influence makes the simulated-only studentisation liberal by about
   0.35 points at B = 400 on synthetic nulls. On the real restscore null
   |t| was still 5.0, so the effect is small in practice. Pooling costs 0.4
   to 0.8 points of power. Package-wide, so separate from decision 1.
3. **Default iterations.** Seed dependence of WY falls from about 10 to 6
   percent of analyses from B = 400 to 1000. Ties in with the pending
   default updates (memory: `easyrasch2-pending-default-updates`).

## If decision 1 is yes: implementation notes

- An option of `.bootstrap_pvalues()` used by `RMitemRestscore()` only.
  Other cutoff functions keep their behaviour unless decision 2 says
  otherwise.
- Moments from `rbind(observed, sim)`, with the existing `na.rm = TRUE` and
  zero-variance guard. The validated arithmetic is `split_variants()` in
  `restscore_power_conditional.qmd`.
- `padj_hi <- .wy_stepdown(t_obs, t_sim, B)`,
  `padj_lo <- .wy_stepdown(-t_obs, -t_sim, B)`. Reporting
  `padj = pmin(1, 2 * pmin(padj_hi, padj_lo))` keeps the `padj < alpha`
  convention, and the direction comes from which run rejected.
- Marginal p: 2·min of the one-sided Monte Carlo p-values, validated in the
  first real-null run (2.2 to 2.6 percent per tail).
- BH has no split counterpart. Decide whether `correction = "fdr_bh"`
  keeps |t| or uses the equal-tailed marginal p.
- Update the help page's paragraph on tail imbalance, NEWS.md, and run the
  doc review checklist.

## Related files

- `restscore_cutoff_validation.qmd`: the |t| test's null validation.
- `restscore_eti_null.qmd`: the withdrawn equal-tailed run, the synthetic
  nulls, the split and pooled follow-ups, and the real split null.
- `restscore_power_conditional.qmd`: power of the shipped test against
  infit and iarm, then the split follow-up and its 300-dataset extension.
- `restscore_stability.qmd`: seed stability at B = 400 and 1000, including
  the split and pooled variants.
- `restscore-cutoff-design.md`: design record, with the withdrawal and a
  2026-10-03 summary of these results.
- Caches with simulated matrices, for checking new variants without a
  rerun: `restscore_split_null/`, `restscore_power_split/`,
  `restscore_power_split_ext/`.

The withdrawn first version is what produced the `cond_eq_*` columns in
`restscore_power_conditional/` and the `eq` columns in
`restscore_eti_null/`; those columns are kept as a record, not as results.
