# Joint studentization in `.bootstrap_pvalues()`

Status 2026-10-23: IMPLEMENTED in 1.4.0.9000 (uncommitted) after user approval.
Decisions: changed default, internal `moments = "simulated"` switch kept for
the paper scripts, below-400 notice reworded as a reproducibility notice.
Origin: Johansson2026fwer APM rev1 (R1.2/E4), `rev1/rerun_results.qmd`,
section "R1.2: the studentization test".

## The change

`.bootstrap_pvalues()` (`R/utils-multiplicity.R`) studentizes the observed
statistic and the B simulated statistics by the mean and SD of the simulated
values only. The simulated values are then internally studentized (each one
contributed to its own scale) and the observed value externally. Across
columns, the max of externally studentized values has a heavier tail, so the
Westfall-Young max is liberal at small B.

Proposal: compute the moments over all B + 1 values, observed included. Every
dataset is then scored on the same scale, so under the null the B + 1 rows of
test statistics are exchangeable. The marginal p-value and the first step of
the step-down become exact Monte Carlo tests. Under the complete null the
step-down rejects only if its first step does, so FWER is exact (at most
`floor(.05 (B + 1)) / (B + 1)`). Under partial nulls validity still rests on
subset pivotality, as now.

Code delta (the only lines that change):

```r
pooled <- rbind(observed, sim)
m <- colMeans(pooled, na.rm = TRUE)
s <- apply(pooled, 2L, stats::sd, na.rm = TRUE)
```

## Null check (synthetic, exchangeable)

`dev/studentization_null_check.R`, results in `dev/studentization_null_check.rds`
and `.log`. 10 000 replications per cell, B + 1 iid rows, columns with
different locations and scales. k = 5, 10, 20 and B = 100, 200, 400, 1000, 2000.
Normal and chi-square(4) noise. Monte Carlo SE of a rate is 0.22 points. An
assertion block first checks that the fast path equals the package helper
(and the candidate code) on 40 datasets for both tails.

Familywise error (%), normal noise. Attainable rate 4.95 at B = 100, 5.00 at 2000.

| B | k | two-sided standard | two-sided joint | upper standard | upper joint |
|---|---|---|---|---|---|
| 100 | 5 | 6.44 | 5.04 | 5.89 | 5.19 |
| 100 | 10 | 7.06 | 5.13 | 6.57 | 5.09 |
| 100 | 20 | 7.68 | 4.99 | 6.87 | 4.90 |
| 200 | 20 | 6.11 | 4.91 | 5.87 | 4.97 |
| 400 | 20 | 5.62 | 4.93 | 5.14 | 4.68 |
| 1000 | 20 | 5.28 | 5.12 | 5.22 | 5.01 |
| 2000 | 20 | 5.58 | 5.45 | 5.21 | 5.05 |

Paired difference standard minus joint at B = 100, k = 20: 2.69 points (SE 0.16)
two-sided, 1.97 (SE 0.14) upper. Chi-square noise gives the same picture
(standard 8.07% at B = 100, k = 20, joint 5.04%).

- Joint: 4.67 to 5.53% over all 60 cell-by-tail rates, which is the spread
  0.22-point SEs predict around an exact 5%.
- Standard: liberal at B = 100 (5.7 to 8.1%), growing with k, still 5.1 to
  5.9% at B = 400 with k = 20, indistinguishable from joint from B = 1000.
- The symmetry property behind the proposed unit test was checked once
  (B = 59, k = 6, both tails): rotating each row into the observed slot gives
  minimum padj values that are a permutation of `(1:(B + 1)) / (B + 1)` under
  joint moments and not under the current code.
- Marginal p, two-sided (normal noise): standard 5.10 to 5.20% at B = 100 to
  200, joint 4.93 to 5.05%. A small excess from centering at a mean that
  excludes the observed value.
- Marginal p, upper tail: identical under both, by construction (the
  studentization is monotone within a column, so the p is a rank of raw
  values). Asserted with `identical()`.
- Caveat on the skewed arm: chi-square(4) noise has a short lower tail, so its
  two-sided and upper columns coincide. Two-sided behaviour is tested by the
  normal arm only.

## Power cost (cached, no new simulation)

From `Johansson2026fwer/rev1/fresh_dich_misfit.rds`, which already stored both
variants. 20 dichotomous items, 250 replications per cell. Detection of the
misfitting items (padj < .05):

- k = 1 misfitting item: 99.6 to 100% under both, no difference.
- k = 3, n = 150: 77.3 vs 74.5 at B = 100 (difference 2.8, SE 0.59), 1.5 at
  B = 400 (SE 0.43), 0.13 from B = 1000.
- k = 3, n = 250: 1.2 at B = 100, 0.7 at 400, 0.3 or less from 1000.

The observed misfit enters its own column's moments and shrinks its own t, at
order t^2 / B. Part of the standard variant's small-B power is the same
liberality that inflates its FWER.

## Design proposal

1. **Changed behaviour, no new user-facing argument.** The standard variant
   has no property worth choosing. From B = 1000 the two agree within Monte
   Carlo noise, so the main visible effect is at the small B where the
   current results are wrong. An exported argument in seven functions for an
   option nobody should pick adds API surface without a use.
2. **Optional internal switch** `moments = c("joint", "simulated")` on
   `.bootstrap_pvalues()` only, not exposed. Reason: `Johansson2026fwer`
   scripts call `easyRasch2:::.bootstrap_pvalues()` as the "standard" arm
   (`rev1/fresh_rerun.R`, `eval_nested()`). After the change, a fresh run of
   that script would silently give joint values in both columns. The caches
   are safe, a re-run is not. The alternative is to copy the old function
   into the paper scripts.
3. **Apply the joint moments to the marginal p as well.** It is exact for the
   two-sided marginal (removes a 0.1 to 0.25-point excess) and changes nothing
   for the upper tail. One set of moments keeps `padj >= p` and feeds BH/BY
   the corrected marginal.

### What moves, per call site

| Function | File | Tail | Marginal p | padj (fwer) | padj (BH/BY) |
|---|---|---|---|---|---|
| `RMitemInfit` (infit and outfit) | conditional_infit.R:480 | two-sided | small | yes | small |
| `RMitemRestscore` | item_restscore.R:458 | two-sided | small | yes | small |
| `RMdifGamma` | dif_partgam.R:388 | two-sided | small | yes | small |
| `RMdimCFA` loadings | cfa_cutoff.R:723 | two-sided | small | yes | small |
| `RMlocdepQ3` | local_dependence.R:701 | upper | none | yes | none |
| `RMlocdepGamma` | ld_partgam.R:475 | upper | none | yes | none |
| `RMdimCFA` fit indices | cfa_cutoff.R:692 | upper | none | yes | none |
| `RMdimResidualPCA` | residual_pca.R:323 | upper, `"none"` | none | n/a | n/a |

`RMdimResidualPCA` is provably unaffected (single upper-tail marginal).
`person_fit.R` has its own p-value code and is unaffected.
`RMitemRestscore` labels direction from the sim-only null mean
(item_restscore.R:468). That is separate from the p-value and can stay.

### Text that becomes wrong

The "below 400 iterations the correction is mildly liberal" claim, which this
change removes as a reason:

- `.notify_low_iterations()` and `.iteration_note()` (utils-multiplicity.R)
- header comment of utils-multiplicity.R ("studentised by the *bootstrap* mean and SD")
- `RMitemInfit` roxygen, conditional_infit.R:117-140 and the caption at :572
- `RMitemRestscore` caption text, item_restscore.R:542
- `RMlocdepQ3` docs, local_dependence.R:105
- `RMdifGamma` docs, dif_partgam.R:90
- `RMlocdepGamma` docs, ld_partgam.R:109
- `RMdimCFA` docs and message, cfa_cutoff.R:509 and :595

The 400 threshold would then rest on reproducibility (seed disagreement) and
p-value resolution only. Whether to keep a console notice below 400, reworded,
or drop it and keep only the caption tier, is open.

### Tests (once approved)

In `tests/testthat/test-infit_pvalue.R`:

- Exact symmetry, deterministic and fast: take B + 1 continuous rows, treat
  each row in turn as observed and the rest as simulated, collect the minimum
  padj. Under joint moments the B + 1 values are a permutation of
  `(1:(B + 1)) / (B + 1)`. The current code fails this.
- Upper tail: marginal p identical to the raw-rank p.
- Existing three helper tests should pass unchanged (check).

Shared helper change, so the full suite runs, plus `devtools::document()`.
NEWS (one bullet, 1.4.0.9000): bootstrap p-values in the seven functions now
studentize by moments of observed and simulated values together, which makes
the family-wise correction exact at any number of iterations. **Corrected
p-values change, mostly below 1000 iterations.**

## Open decisions for the user

1. Change the behaviour (recommended) or add a user-facing argument.
2. Keep an internal `moments` switch for the paper scripts, or copy the old
   function into them.
3. Below-400 notice: reword as a reproducibility notice, or drop it.
