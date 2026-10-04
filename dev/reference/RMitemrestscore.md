# Item Restscore Analysis

Computes observed and model-expected item-restscore correlations using
[`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html),
and enriches the output with the absolute difference between observed
and expected values, item average locations, and item locations relative
to the sample mean person location.

## Usage

``` r
RMitemRestscore(
  data,
  cutoff = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable",
  sort,
  p_adj = "BH"
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed, but at least one complete case (row with no `NA`) must be
  present.

- cutoff:

  Optional. Default `NULL`, which tests each item against the asymptotic
  normal reference distribution from
  [`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html),
  adjusted by `p_adj`. Can be:

  - The return value of
    [`RMitemRestscoreCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md):
    the parametric bootstrap null distribution of the observed minus
    expected gamma. Adds the columns `Diff_low` and `Diff_high` (the
    descriptive interval), and, with `p_value` resolving to `TRUE`,
    `p_restscore` and `padj_restscore`.

  - The `$item_cutoffs` data.frame from
    [`RMitemRestscoreCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
    directly (columns `Item`, `diff_low`, `diff_high`). Items are then
    flagged against the interval.

  With a `cutoff`, the asymptotic `p_adjusted` column is dropped and
  `p_adj` is ignored.

- p_value:

  Logical or `NULL`. Whether to compute bootstrap p-values from the
  simulated null distribution and flag on them.

  `NULL` (the default) means **use them when they are available**:
  `TRUE` when `cutoff` is the full
  [`RMitemRestscoreCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
  object, which carries the per-item simulated values in its `$results`
  element, and `FALSE` otherwise. Explicit `TRUE` without the full
  object is an error rather than a silent downgrade. With `FALSE` and a
  `cutoff`, items are flagged against the interval instead, and a
  one-time message reports the family-wise error rate that implies.

- correction:

  Character. Multiple-comparison correction applied across items when
  `p_value = TRUE`: `"fwer"` (default) for the Westfall-Young
  studentised-max step-down (family-wise error rate), `"fdr_bh"` /
  `"fdr_by"` for Benjamini-Hochberg / Benjamini-Yekutieli false
  discovery rate control, or `"none"` for uncorrected per-item p-values.

- alpha:

  Numeric in (0, 1). Significance level for the `Flagged` column when
  `p_value = TRUE`. Default `0.05`.

- output:

  Character string controlling the return value. Either `"kable"`
  (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"` for the underlying data.frame.

- sort:

  Optional character string. When `sort = "diff"`, rows are sorted by
  the absolute magnitude of `Difference` in descending order, so that
  both over- and underfitting items appear near the top.

- p_adj:

  Character string specifying the p-value adjustment method passed to
  [`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html).
  Default `"BH"` (Benjamini-Hochberg); use `"none"` for unadjusted
  p-values. Run
  [`?stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for the
  list of available methods. Only used when `cutoff = NULL`.

## Value

- If `output = "kable"`: a `knitr_kable` object (plain text table via
  `format = "pipe"`) with columns for item name, observed and expected
  restscore correlations, the signed difference (observed minus
  expected), adjusted p-value, the `Flagged` misfit label, and item
  location relative to the sample mean person location.

- If `output = "dataframe"`: a data.frame with columns `Item`,
  `Observed`, `Expected`, `Difference`, `p_adjusted`, `Flagged`, and
  `Relative_location`. `Flagged` is `"overfit"` (observed above
  expected, adj. p \< .05), `"underfit"` (below, adj. p \< .05), or `""`
  (not flagged).

With a `cutoff`, `p_adjusted` is replaced by `Diff_low` and `Diff_high`,
followed by `p_restscore` (marginal two-sided bootstrap p-value) and
`padj_restscore` (corrected p-value) when `p_value` resolves to `TRUE`.
`Flagged` then reflects `padj_restscore < alpha`, or the interval when
`p_value = FALSE`.

The `Difference` column is signed (observed minus expected): *positive*
values indicate that the item correlates more strongly with the
rest-score than the Rasch model predicts (over-discrimination /
*overfit*, often associated with local dependence), and *negative*
values indicate weaker-than-expected association (under-discrimination /
*underfit*, often associated with multidimensionality or noise).

## Details

Item-restscore correlations using Goodman-Kruskal's gamma (Kreiner,
2011) measure the association between a person's score on a single item
and their total score on the remaining items (the "restscore"). Under a
correctly fitting Rasch model, observed and model-expected correlations
should agree closely.

Item parameters are estimated by conditional maximum likelihood via
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category PCM); the item-restscore statistic
itself comes from
[`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html)
and is conditional on the total score, so it is invariant to the
estimation engine. Per-item average locations are the means of the CML
thresholds, and the person-location reference is the mean of the Warm
WLE estimates.

Relative item location is defined as the item's average location minus
the sample mean person location, providing a measure of item targeting.

The `iarm` package must be installed (it is in Suggests, not Imports).

**The asymptotic p-value is miscalibrated.** Without a `cutoff`, each
item is tested with `(observed - expected) / SE` against the standard
normal, where the SE is that of the observed gamma alone. The expected
gamma is estimated from the same data and correlates with the observed
one, so the SE is too large for the difference. The difference is also
biased upwards in small samples. Under a true Rasch model the test
therefore flags too many items as overfit in small samples and too few
as underfit at every sample size, and the adjustment chosen in `p_adj`
does not correct either. In simulation, 20 dichotomous items at n = 150
gave at least one BH-flagged item in 13 percent of datasets (30 percent
when 1.5 logits off target), almost all of them overfit, while 9
polytomous items at n = 1000 gave about 1 percent. Pass the object from
[`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
as `cutoff` for p-values from a parametric bootstrap null instead.

**Bootstrap p-values.** When `p_value = TRUE`, each item's observed
`Difference` is compared against its simulated null distribution (from
`cutoff$results`), studentised by the bootstrap mean and SD. Because the
bootstrap mean is subtracted, the small-sample upward bias of the
difference is removed, and because the bootstrap SD is used, the
correlation between observed and expected gamma is accounted for. The
marginal p-value is the two-sided Monte-Carlo p-value
`(1 + #{|t*| >= |t|}) / (B + 1)`, and `correction = "fwer"` applies the
Westfall-Young studentised-max step-down across items (Ferreira, 2024).

The two directions do not get equal shares of alpha. Gamma is bounded at
1, so the null of `Difference` is left-skewed and a two-sided test on
`|t|` rejects more often in the long (underfit) tail. In simulation
under a true Rasch model the total rate was nominal, but the per-item
marginal rate was about 3 percent for underfit and 2 percent for overfit
with 20 dichotomous items at n = 150 and 1.5 logits off target,
narrowing toward equal shares with polytomous items and larger samples.
With `correction = "fwer"`, over four such conditions, at least one item
was flagged as underfit in 3.7 percent of datasets and as overfit in 1.6
percent, 5.0 percent in total. An equal-tailed test, with alpha/2 for
each direction, was evaluated and not adopted: it detected 2 to 3
percentage points fewer underfitting items and no more overfitting ones.

**Reproducibility.** The bootstrap p-values depend on the simulated
null, so two analyses of the same data with different seeds can disagree
about items near the decision boundary. In simulation, in conditions
chosen to include such items, two seeds disagreed about the flag of at
least one item in about 10 percent of analyses at 400 iterations and
about 6 percent at 1000 with `correction = "fwer"`, and in 13 and 11
percent with `correction = "fdr_bh"`. Set `seed` in
[`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
for a reproducible analysis, and use 1000 or more iterations for a final
one.

The direction in `Flagged` follows the studentised value, not the sign
of `Difference`. At small sample sizes the null mean of `Difference` is
slightly positive, so an item can be flagged `"underfit"` while its
`Difference` is still just above zero. This is rare.

## Multiple comparisons

The marginal p-value controls the error rate of a *single* comparison:
for one item (or item pair) decided on in advance it is the relevant
value. But scanning all *k* comparisons and flagging whichever fall
below `alpha` tests *k* hypotheses at once, so the chance of at least
one false flag inflates to roughly \\1 - (1 - \alpha)^k\\ (e.g. about
34% for *k* = 8 at `alpha = 0.05`) – even when every marginal p-value is
correctly calibrated. The corrected (adjusted) p-value controls this:
`correction = "fwer"` bounds the probability of *any* false flag
(strict, lower power), while `"fdr_bh"` / `"fdr_by"` bound the expected
*proportion* of false flags among those raised (a more lenient middle
ground). Rule of thumb: use the marginal p-value for a single
pre-specified comparison, and a corrected p-value when screening the
whole table – the usual workflow.

## Flags depend on the other items

Item infit and item-restscore both judge each item against expectations
from a Rasch model fitted to all items, so misfitting items change what
is expected of the others. A flag in one direction can therefore produce
flags in the other. In simulation with seven polytomous items, two
overfitting items led to about 6 percent of the fitting items being
flagged as underfit, against about 1.5 percent with one overfitting
item. The reverse effect was smaller: with underfitting items planted, 1
to 2 percent of the fitting items were flagged as overfit. When items
are flagged in both directions, consider removing the clearest misfit
and testing again before interpreting the rest.

## References

Kreiner, S. (2011). A Note on Item–Restscore Association in Rasch
Models. *Applied Psychological Measurement, 35*(7), 557–561.
[doi:10.1177/0146621611410227](https://doi.org/10.1177/0146621611410227)

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(6).
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

## See also

[`RMitemRestscoreCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md),
[`RMitemRestscorePlot`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscorePlot.md)

## Examples

``` r
# \donttest{
if (requireNamespace("iarm", quietly = TRUE)) {
  # Simulate binary item response data (8 items, 200 persons)
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
  )
  colnames(sim_data) <- paste0("Item", 1:8)

  # Default kable output
  RMitemRestscore(sim_data)

  # Sorted by absolute difference
  RMitemRestscore(sim_data, sort = "diff")

  # Return as data.frame for further processing
  df <- RMitemRestscore(sim_data, output = "dataframe")

  # Bootstrap null distribution, flagging on Westfall-Young corrected
  # p-values (use 1000 or more iterations in a final analysis)
  if (requireNamespace("ggdist", quietly = TRUE)) {
    cutoff_res <- RMitemRestscoreCutoff(sim_data, iterations = 100,
                                        parallel = FALSE, seed = 42)
    RMitemRestscore(sim_data, cutoff = cutoff_res)
  }
}
#> 
#> 
#> 
#> 
#> 
#> 
#> Table: Item-restscore associations. n = 200 respondents. Two-sided parametric bootstrap p-values for the difference from 100 iterations, conditional DGP, multiplicity correction: Westfall-Young step-down (FWER). p-values cannot be smaller than 1/(100+1) = 0.0099. This is below the calibrated floor of 400, where the correction is mildly liberal (Johansson, 2026). Flagged (adj. p < 0.05): overfit = observed above expected (over-discrimination, often local dependence); underfit = below (under-discrimination, often multidimensionality/noise).
#> 
#> |Item  | Observed| Expected| Difference| Diff low| Diff high|      p| p (adj)|Flagged | Rel. location|
#> |:-----|--------:|--------:|----------:|--------:|---------:|------:|-------:|:-------|-------------:|
#> |Item1 |     0.00|     0.01|     -0.006|   -0.154|     0.147| 0.9604|  1.0000|        |         -0.23|
#> |Item2 |     0.00|     0.01|     -0.016|   -0.148|     0.182| 0.7624|  1.0000|        |          0.13|
#> |Item3 |    -0.01|     0.01|     -0.026|   -0.133|     0.169| 0.5446|  0.9802|        |         -0.01|
#> |Item4 |    -0.05|     0.01|     -0.061|   -0.168|     0.187| 0.4950|  0.9802|        |         -0.05|
#> |Item5 |     0.07|     0.01|      0.058|   -0.160|     0.215| 0.4554|  0.9802|        |          0.09|
#> |Item6 |    -0.04|     0.01|     -0.049|   -0.166|     0.166| 0.4653|  0.9802|        |         -0.09|
#> |Item7 |     0.15|     0.01|      0.136|   -0.159|     0.221| 0.1782|  0.6931|        |         -0.05|
#> |Item8 |    -0.01|     0.01|     -0.021|   -0.148|     0.171| 0.7822|  1.0000|        |         -0.03|
# }
```
