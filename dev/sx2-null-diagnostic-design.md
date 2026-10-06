# S-X2 null diagnostic: design

Status: RUN 2026-10-04 (7.3 min), results in `sx2_null_diagnostic.qmd`. Written 2026-10-04. Precedes the revised S-X2
comparison (fixes 1, 5, 6 and 7 in the review of the July 2026 comparison (`sx2_itemfit_sim.qmd`, removed from dev/ on 2026-10-05, recoverable from git)), since
its result decides how S-X2 enters that run.

## The problem

Under a true Rasch model (`k = 0` in the July 2026 comparison (`sx2_itemfit_sim.qmd`, removed from dev/ on 2026-10-05, recoverable from git), 20 dichotomous
items, 200 replications per n), `mirt::itemfit(fit_stats = "S_X2")` rejects
too often:

| n | per item, uncorrected | datasets with any BH flag |
|---|---|---|
| 150 | 6.9 | 10.5 |
| 250 | 7.1 | 10.5 |
| 500 | 8.7 | 13.0 |
| 1000 | 7.1 | 13.5 |
| 2000 | 8.8 | 19.0 |

Per item, pooled over n, the rate runs from 4.4 to 11.7 percent. That is wider
than the binomial SE of about 0.8 points on 1000 item-replications implies, but
replications share a master dataset (fix 5), so the SE is understated and the
spread may be partly noise. This conflicts with the nominal Type I error that
Kang and Chen (2008) report. The cache stores flags only, so nothing further can
be read from it.

## Candidate causes

**H1. One degree of freedom too few.** In mirt 1.46.1, `itemfit()` computes, for
each item,

    df = (non-empty cells) - (rows) - (item parameters) - (latent parameters)

where the last term is `sum(pars[[J + 1]]@est)`, the estimated group
parameters. For `itemtype = "Rasch"` the latent variance is estimated (one
parameter), so every item loses one df beyond its own location parameter. As I
read Orlando and Thissen (2000), df is the number of score groups after
collapsing minus the item's own parameters, with no latent term. The variance
plays the role of a slope shared by all items, and charging it in full to
every item counts it J times. **Verify the definition against Orlando and
Thissen (2000) and Kang and Chen (2008) before relying on this.**

Size, without simulation: if S-X2 follows chi-square on df + 1 and is referred
to df, the rejection rate at .05 is 7.5% at df = 10, 7.1% at df = 14 and 6.9%
at df = 17. mirt's df on 20 dichotomous items was 8 to 11 at n = 150 and 15 to
17 at n = 2000 in a test fit. That matches the observed 6.9 to 8.8% closely
enough to be the main candidate.

**H2. Sparse cells.** `mincell = 1` collapses score groups until every expected
count is at least 1. Expected counts of 1 to 5 still make the chi-square
approximation poor in the tails, which could add to H1 and would matter most at
small n and for items far from the sample.

**H3. Plug-in parameters.** Expected proportions are computed at MML estimates,
and S-X2 makes no allowance for their sampling error beyond the df deduction.

H1 predicts a constant per-item excess at every n. It does not obviously explain
why the BH family-wise rate rises with n, which depends on the far tail of the p
distribution. That part is left open until the diagnostic shows where in the
distribution the excess sits.

## Design

Null only. Fresh data for every replication, no master dataset.

- **Item sets.** (a) 20 dichotomous items evenly spaced on [-2, 2]. (b) The 20
  dichotomous locations of the original study (`runif(20, -2, 2)` with seed
  20250727), to connect to the cached results. (c) 9 polytomous items with four
  categories, locations evenly spaced on [-1.5, 1.5], thresholds at location
  -1.2, 0 and +1.2, matching `restscore_asymptotic_null.qmd`.
- **Persons.** theta ~ N(0, 1.5^2).
- **Sample sizes.** 150, 500, 2000.
- **Replications.** 2000 per cell. Per-item SE of a 5% rate is about 0.5
  points per item, and about 0.5 points for a 5% family-wise rate.
- **Discard rule.** As in the restscore study: a replication with an unobserved
  category is discarded and counted.

Four variants per replication, all from at most two fits:

| Variant | Fit | itemfit call | Tests |
|---|---|---|---|
| V0 | Rasch, variance estimated | default (`mincell = 1`) | reproduces the problem |
| V1 | same fit as V0 | same call, p recomputed on df + 1 | H1 |
| V2 | same fit as V0 | `mincell = 5` | H2 |
| V3 | Rasch, variance fixed at 2.25 | default | H1 against H3 |

V3 has no estimated latent parameter, so mirt's df is the Orlando-Thissen df
with no correction needed, and it also removes the variance from the plug-in
error. If V1 and V3 agree, the df is the cause and estimating the variance does
not matter. If V3 is better than V1, plug-in error in the variance contributes.

**Stored per item and replication:** S-X2, mirt df, p, RMSEA, for each variant.
Number of score groups is recoverable as df + 2 (V0, V2) or df + 1 (V3) for
dichotomous items.

## Summaries

1. Per-item rejection at .05, uncorrected, overall and by item location.
2. Mean S-X2 against mean df. Under a correct reference they are equal. H1
   predicts mean S-X2 near df + 1 for V0.
3. Uniformity of p: rejection at .01, .05 and .10, and a histogram of p for V0
   and V1. The far tail (p < .005) separately, since that drives BH over 20
   items.
4. Family-wise rate with BH and with Bonferroni.
5. Discarded replications per cell.

## What each outcome means

- **V1 near nominal at every n and item set.** The df is the cause. The revised
  comparison uses corrected df for S-X2. Whether to report this to the mirt
  maintainer is your call.
- **V1 improves but the far tail stays heavy, rising with n.** H1 plus something
  else. V2 and V3 show whether it is collapsing or estimation. If neither fixes
  it, S-X2 enters the revised comparison with a simulation-based cutoff (the
  untested hypothesis in conclusion 5 of the original study) instead of its
  asymptotic p-value.
- **V1 no better than V0.** H1 is wrong, and the design note above needs
  correcting before going further.

## Cost

One Rasch fit plus `itemfit()` takes about 0.1 s at n = 150 to 2000. V0 to V2
share a fit and V3 adds one, so about 0.4 s per replication. 3 item sets x 3 n
x 2000 replications is 18,000 replications, about 2 hours on one core and 12 to
15 minutes on 10.

## Not in scope

- A CML-based S-X2. Noted in the original study as a way to drop the normality
  assumption. Separate question.
- Power. That belongs to the revised comparison.
