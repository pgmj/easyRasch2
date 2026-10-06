# S-X2 against conditional infit, revised: design

Status: RUN 2026-10-04 to 10-05 (7.6 h, 24 cells, no failures). Results and interpretation in `sx2_comparison_v2.qmd`. Written 2026-10-04. Replaces the comparison in
the July 2026 comparison (`sx2_itemfit_sim.qmd`, removed from dev/ on 2026-10-05, recoverable from git), whose review found seven problems (listed below). Builds
on `sx2_null_diagnostic.qmd`, which found that mirt's S-X2 df is one too small
for Rasch models and that a far-tail excess remains for dichotomous items after
correcting it.

## Question

Does S-X2 add anything to conditional infit with simulation-based cutoffs, when
both are held to the same family-wise error rate? Specifically, does it detect
misfit that infit misses, or flag fewer fitting items when several items
misfit (the overfit cascade)?

Item-restscore is left out. Its comparison with infit is covered by
`restscore_infit_power.qmd` and `restscore_power_conditional.qmd`, and its
p-values depend on the open equal-tailed decision.

## What changes from the original study

| Problem in the original | Change |
|---|---|
| 1. Infit per-item cutoffs without multiplicity correction, S-X2 with BH | Every method under family-wise control, plus a comparison at matched null error rate |
| 2. S-X2 p-values liberal (mirt df) | S-X2 with corrected df, and a parametric-bootstrap S-X2 |
| 3. Conclusion 3 wrong family-wise | Follows from the new null results |
| 4. Perfect-recovery claim with no code, one cell reported | Computed for every cell |
| 5. All replications drawn from one master dataset | Fresh data for every replication |
| 6. One clustered `runif` draw of item locations | Evenly spaced locations |
| 7. easyRasch2 1.1.0, infit at 100/200 iterations | Current easyRasch2 defaults, B = 400. The installed version is recorded in every cell file (1.3.1.9002 at the smoke test) |

## Data generation

- **Items.** 20 dichotomous items, locations evenly spaced on [-2, 2].
- **Persons.** $(\theta_1, \theta_2)$ bivariate normal, both SD 1.5, mean 0,
  correlation $\rho$. Fitting items respond to $\theta_1$, misfitting items to
  $\theta_2$, as `eRm::sim.xdim()` with a 0/1 weight matrix in the original.
  Generated directly with `plogis()`, fresh for every replication.
- **Misfitting items.** k = 1: the item nearest 0 logits. k = 3: the items
  nearest 0, -1 and -2 logits, as in the original.
- **Strength.** $\rho$ = 0.15 (strong, the original main study), 0.50 and 0.75
  (subtle, the original supplement).
- **Sample sizes.** 150, 500, 1000, 2000. The original's 250 is dropped to save
  time, since 150 and 500 bracket it.

Cells:

| k | $\rho$ | n | Replications |
|---|---|---|---|
| 0 | none | 4 | 500 |
| 1 | 0.50, 0.75 | 4 | 200 |
| 3 | 0.15, 0.50, 0.75 | 4 | 200 |

k = 1 at $\rho$ = 0.15 is left out because detection saturates there in the
original. That gives 4 null cells and 20 misfit cells. With 500 replications the
Monte Carlo SE of a 5 percent family-wise rate is 1.0 point. With 200, a
detection rate near 50 percent has an SE of 3.5 points per item, and the
paired comparisons below are more precise than that.

## Methods

All on the same dataset in each replication.

| Label | Statistic | Reference | Multiplicity |
|---|---|---|---|
| `infit_wy` | conditional infit MSQ | `RMitemInfitCutoff()`, B = 400, defaults | WY (`correction = "fwer"`), as the package |
| `infit_bh` | same | same bootstrap p | BH on the marginal bootstrap p |
| `sx2_mirt` | S-X2 | chi-square on mirt's df | BH. Reference only, reproduces the problem |
| `sx2_df1` | S-X2 | chi-square on df + 1 | BH |
| `sx2_boot` | S-X2 | parametric bootstrap, B = 400 | WY |
| `rmsea_boot` | RMSEA.S_X2 | same bootstrap draws | WY |
| `rmsea_05` | RMSEA.S_X2 | fixed .05 | none. Descriptive, as in the original |

**The S-X2 bootstrap.** Fit the Rasch model in mirt (latent variance
estimated). Draw B = 400 datasets from the fitted model, with person locations
from $N(0, \hat\sigma^2)$, refit each and compute S-X2 for every item. The
number of score groups changes between draws because of collapsing, so S-X2 is
not on a common scale across draws. The bootstrap statistic is therefore
$z = \Phi^{-1}(1 - p)$, with $p$ from chi-square on df + 1, which is on a common
scale whatever the df. RMSEA is used as it is. Both go through
`.bootstrap_pvalues(tail = "upper", correction = "fwer")`, the same
step-down the package uses for infit, so the two arms differ only in the
statistic. A draw where the refit fails or an item has no variance is dropped
and counted.

This also tests the untested hypothesis in conclusion 5 of the original study,
that a simulation-based RMSEA cutoff would recover the detection the fixed .05
threshold misses.

## Stored per item and replication

Infit MSQ, `p_infit`, `padj_infit`. S-X2, mirt df, its asymptotic p on mirt's
df and on df + 1, RMSEA, and the bootstrap marginal and WY-adjusted p for
`sx2_boot` and `rmsea_boot`. Number of bootstrap draws dropped. Storing
p-values, not flags, lets any threshold be applied afterwards.

Results are cached per cell in `sx2_comparison_v2/`, so an interrupted run
resumes where it stopped.

## Summaries

1. **Null (k = 0).** Per-item rate and family-wise rate for every method, by n.
   `infit_wy`, `sx2_boot` and `rmsea_boot` should hold 5 percent. The
   diagnostic predicts 7 to 9 percent for `sx2_df1`.
2. **Detection.** Per misfitting item, by location, for every misfit cell.
3. **False flags on fitting items.** Per item, and the share of replications
   with any fitting item flagged. This is the cascade.
4. **Perfect recovery.** The share of replications where every misfitting item
   is flagged and no fitting item is, for every cell.
5. **Matched null error rate.** For each asymptotic method and n, the threshold
   on its adjusted p that gives 5 percent family-wise in the null cell at that
   n, applied to the misfit cells. The threshold is taken from the same design,
   so it is an oracle and favours the method slightly. It shows what a method
   could do if its null were right.
6. **Paired comparison.** `infit_wy` against `sx2_boot` on the same datasets:
   discordant pairs for detecting each misfitting item and for false flags,
   with McNemar-style counts, since that is more precise than comparing two
   rates.

## What would change the package

Stated before the run. S-X2 is worth adding as a complement if `sx2_boot` or
`rmsea_boot` holds 5 percent under the null and, in at least one realistic cell
($\rho$ = 0.50 or 0.75), either detects misfit `infit_wy` misses or flags
clearly fewer fitting items, without a comparable loss elsewhere. If it holds
the null but does neither, the original conclusion that S-X2 is a large-sample
complement does not survive. If only the asymptotic `sx2_df1` would be shipped,
its residual excess from the diagnostic has to be stated wherever it is
reported.

## Cost

Measured on one core, one replication, 20 items: infit with WY 9 s at n = 500
and 20 s at n = 2000, S-X2 bootstrap 27 s and 40 s. At about 40 s on average,
the 6000 replications take about 67 core-hours, 6 to 7 hours on 10 cores. That
is an overnight run. If that is too long, the cheapest cuts are B = 200 for the
S-X2 bootstrap only (about 4.5 hours) or dropping n = 1000 (about 5 hours).

## Not in scope

- Item-restscore, see above.
- Polytomous items. The diagnostic found S-X2 close to nominal for them after
  the df correction, so a polytomous arm is a natural follow-up if the
  dichotomous result is positive.
- A CML-based S-X2.
