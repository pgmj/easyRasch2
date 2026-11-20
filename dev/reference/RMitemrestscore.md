# Item-restscore association

Computes, for each item, the observed and model-expected association
between the item and the rest score (the total score on the other
items), measured by Goodman-Kruskal's gamma with
[`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html).
Reports their signed difference, the item's location relative to the
mean person location, and a misfit flag. With the result of
[`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md),
items are flagged on a bootstrap p-value corrected across items, which
is recommended over the default asymptotic test (see Details).

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

  A data.frame or matrix of item responses, items scored from 0. Missing
  values (`NA`) are allowed, but the statistic uses complete cases only,
  so at least one is needed.

- cutoff:

  Optional. Default `NULL`, which tests each item with the asymptotic
  normal test of
  [`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html),
  adjusted by `p_adj`. Either of:

  - The full result of
    [`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md).
    This enables the bootstrap p-values (see `p_value`) and adds the
    columns `Diff_low` and `Diff_high`, the interval of the simulated
    differences.

  - Its `$item_cutoffs` data.frame alone (columns `Item`, `diff_low`,
    `diff_high`). Items are then flagged against the interval.

  With a `cutoff`, the asymptotic `p_adjusted` column is dropped and
  `p_adj` is ignored.

- p_value:

  Logical or `NULL`. Whether to flag on bootstrap p-values from the
  simulated values in `cutoff`. `NULL` (default) uses them whenever
  `cutoff` is the full
  [`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
  result. `TRUE` without it is an error. With `FALSE` and a `cutoff`,
  items are flagged against the interval, and a one-time message gives
  the family-wise error rate that implies.

- correction:

  Character. Correction of the bootstrap p-values across items: `"fwer"`
  (default, Westfall-Young step-down, controls the family-wise error
  rate), `"fdr_bh"` or `"fdr_by"` (Benjamini-Hochberg or
  Benjamini-Yekutieli, control the false discovery rate), or `"none"`.

- alpha:

  Numeric in (0, 1). An item is flagged when its corrected bootstrap
  p-value is below `alpha`. Default `0.05`. The asymptotic test without
  a `cutoff` always flags at .05.

- output:

  Character. `"kable"` (default) for a
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"`.

- sort:

  Optional. `sort = "diff"` sorts the rows by the absolute `Difference`,
  largest first, so that over- and underfitting items both come first.
  Other values leave the order unchanged.

- p_adj:

  Character. Adjustment of the asymptotic p-values across items, passed
  to
  [`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html).
  Default `"BH"` (Benjamini-Hochberg). `"none"` gives unadjusted
  p-values, and any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) can be
  used. Only used when `cutoff = NULL`.

## Value

- `output = "kable"`: a `knitr_kable` table with the item name, the
  observed and expected gamma, their difference, the adjusted p-value,
  the `Flagged` label and the relative item location.

- `output = "dataframe"`: a data.frame with the columns `Item`,
  `Observed`, `Expected`, `Difference`, `p_adjusted`, `Flagged` and
  `Relative_location`. `Flagged` is `"overfit"` (observed above
  expected, adjusted p \< .05), `"underfit"` (below, adjusted p \< .05)
  or `""`.

With a `cutoff`, `p_adjusted` is replaced by `Diff_low` and `Diff_high`,
followed by `p_restscore` (uncorrected equal-tailed bootstrap p) and
`padj_restscore` (corrected p) when bootstrap p-values are used.
`Flagged` then marks `padj_restscore < alpha`, or values outside the
interval with `p_value = FALSE`.

`Difference` is observed minus expected gamma. Positive values mean the
item is more strongly associated with the rest score than the model
predicts (overfit, often from local dependence). Negative values mean a
weaker association (underfit, often from multidimensionality or noise).

## Details

**Statistic.** The item-restscore gamma (Kreiner, 2011) measures the
association between the responses to one item and the rest score. Under
the Rasch model the observed and expected gamma agree, apart from
sampling error. The statistic is conditional on the total score.

**Data and model.** The item parameters are estimated by conditional
maximum likelihood with
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category partial credit model). The statistic
uses complete cases, and every response category must be used in them,
or the function stops and names the item. The relative item location is
the mean of the item's CML thresholds minus the mean WLE person
location, a description of targeting. The iarm package must be installed
(it is in Suggests).

**The asymptotic p-value is miscalibrated.** Without a `cutoff`, each
item is tested with `(observed - expected) / SE` against the standard
normal, where the SE is that of the observed gamma alone. The expected
gamma is estimated from the same data and correlates with the observed
one, so the SE is too large for the difference. The difference is also
biased upwards in small samples. When the Rasch model holds, the test
therefore flags too many items as overfit in small samples and too few
as underfit at every sample size, and the adjustment in `p_adj` corrects
neither. In simulation, 20 dichotomous items at n = 150 gave at least
one BH-flagged item in 13% of datasets (30% when 1.5 logits off target),
almost all of them overfit, while 9 polytomous items at n = 1000 gave
about 1%. Pass the result of
[`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
as `cutoff` for p-values from a parametric bootstrap instead. The
simulation study is reported at
<https://pgmj.github.io/rasch_book/chapters/restscore-asymptotic-null.html>.

**Bootstrap p-values.** Each item's observed `Difference` is compared
with its simulated values from `cutoff$results`, studentised by the mean
and SD of the observed and simulated values together. Subtracting that
mean removes the small-sample bias of the difference, and dividing by
that SD accounts for the correlation between the observed and expected
gamma.

The test is equal-tailed: overfit and underfit each get `alpha / 2`.
Gamma is bounded at 1, so the simulated `Difference` is skewed to the
left, and a two-sided test on `|t|` would reject more often in the long
(underfit) tail. The uncorrected p is twice the smaller of the one-sided
Monte Carlo p-values, `(1 + #{t* >= t}) / (B + 1)` and
`(1 + #{t* <= t}) / (B + 1)` for B simulated datasets, so it cannot be
smaller than `2 / (B + 1)`. With `correction = "fwer"`, the
Westfall-Young step-down on the studentised maximum (Ferreira, 2024) is
run once in each direction, and `padj_restscore` is twice the smaller of
the two corrected p-values. `"fdr_bh"` and `"fdr_by"` adjust the
equal-tailed uncorrected p-values.

In simulation with the Rasch model true, over four conditions (20
dichotomous or 9 polytomous items, n = 150 to 1000, on or off target),
`correction = "fwer"` flagged at least one item in 5.0% of datasets,
2.4% underfit and 2.75% overfit. A test on `|t|` also held 5.0%, but
split it 3.7% and 1.6%. The equal-tailed test also gave fewer false
underfit flags when other items overfit (4.7% against 6.2% with two
overfitting items), at a cost of 2 to 3 percentage points of power to
detect underfitting items. Before easyRasch2 1.4.0.9004 the test was on
`|t|`.

**Reproducibility.** The bootstrap p-values depend on the simulations,
so two analyses of the same data with different seeds can disagree about
items near the decision boundary. In simulation, in conditions chosen to
include such items, two seeds disagreed about the flag of at least one
item in about 10% of analyses at 400 iterations and 6% at 1000 with
`correction = "fwer"`, and in 13% and 11% with `"fdr_bh"`. Set `seed` in
[`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
for a reproducible analysis, and use 1000 or more iterations for a final
one.

The direction in `Flagged` follows the studentised value, not the sign
of `Difference`. At small sample sizes the simulated `Difference` has a
slightly positive mean, so an item can, rarely, be flagged `"underfit"`
while its `Difference` is just above zero.

## Multiple comparisons

The uncorrected p-value controls the error rate of a single comparison.
For one item (or item pair) chosen in advance it is the relevant value.
Flagging whichever of k comparisons fall below `alpha` tests k
hypotheses at once, and the chance of at least one false flag rises to
about \\1 - (1 - \alpha)^k\\, for example 34% for k = 8 at
`alpha = 0.05`, even when every uncorrected p-value is calibrated. The
corrected p-value controls this. `correction = "fwer"` bounds the
probability of any false flag (strict, less power), and `"fdr_bh"` and
`"fdr_by"` bound the expected share of false flags among those raised.
Use the uncorrected p-value for a single comparison chosen in advance,
and a corrected one when screening the whole table, which is the usual
case.

## Flags depend on the other items

Item infit and item-restscore both judge each item against expectations
from a Rasch model fitted to all items, so misfitting items change what
is expected of the others. A flag in one direction can therefore produce
flags in the other. Overfitting items can push fitting items into the
underfit range, and underfitting items, for example from a second
dimension, can make fitting items look overfit. How large the effect is
depends on how strong the misfit is and how many items share it. In
simulation with seven polytomous items, two overfitting items led to
about 6% of the fitting items being flagged as underfit, and one or two
mildly multidimensional items to 1% to 2% being flagged as overfit. When
items are flagged in both directions, consider removing the clearest
misfit and testing again before interpreting the rest.

## References

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(1), 6.
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

Kreiner, S. (2011). A Note on Item–Restscore Association in Rasch
Models. *Applied Psychological Measurement, 35*(7), 557–561.
[doi:10.1177/0146621611410227](https://doi.org/10.1177/0146621611410227)

## See also

[`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
for the simulations,
[`RMitemRestscorePlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscorePlot.md)
to plot them,
[`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreBoot.md)
for large samples, and
[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
for a second item fit test.

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
#> Table: Item-restscore associations. n = 200 respondents. Equal-tailed parametric bootstrap p-values for the difference from 100 iterations, conditional DGP, multiplicity correction: Westfall-Young step-down (FWER). p-values cannot be smaller than 2/(100+1) = 0.02. This is below the recommended 400, where decisions depend on the seed (Johansson, 2026). Flagged (adj. p < 0.05): overfit = observed above expected (over-discrimination, often local dependence); underfit = below (under-discrimination, often multidimensionality/noise).
#> 
#> |Item  | Observed| Expected| Difference| Diff low| Diff high|      p| p (adj)|Flagged | Rel. location|
#> |:-----|--------:|--------:|----------:|--------:|---------:|------:|-------:|:-------|-------------:|
#> |Item1 |     0.00|     0.01|     -0.006|   -0.154|     0.147| 0.8911|       1|        |         -0.23|
#> |Item2 |     0.00|     0.01|     -0.016|   -0.148|     0.182| 0.7723|       1|        |          0.13|
#> |Item3 |    -0.01|     0.01|     -0.026|   -0.133|     0.169| 0.6139|       1|        |         -0.01|
#> |Item4 |    -0.05|     0.01|     -0.061|   -0.168|     0.187| 0.4752|       1|        |         -0.05|
#> |Item5 |     0.07|     0.01|      0.058|   -0.160|     0.215| 0.4554|       1|        |          0.09|
#> |Item6 |    -0.04|     0.01|     -0.049|   -0.166|     0.166| 0.4554|       1|        |         -0.09|
#> |Item7 |     0.15|     0.01|      0.136|   -0.159|     0.221| 0.2178|       1|        |         -0.05|
#> |Item8 |    -0.01|     0.01|     -0.021|   -0.148|     0.171| 0.7129|       1|        |         -0.03|
# }
```
