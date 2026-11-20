# Yen's Q3 residual correlations

Computes Yen's \\Q_3\\, the correlations between the model residuals of
each pair of items, to assess local dependence. By default the item
parameters are estimated by conditional maximum likelihood and the
person locations by WLE (`estimator = "CML"`). With the result of
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md),
it also tests each pair against simulations, by default with a bootstrap
p-value corrected across pairs.

## Usage

``` r
RMlocdepQ3(
  data,
  cutoff = NULL,
  output = "kable",
  n_pairs = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  estimator = c("CML", "MML")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0. Missing
  values (`NA`) are allowed. \\Q_3\\ is computed from the respondents
  who answered both items of a pair, and respondents with no responses
  are dropped.

- cutoff:

  Optional. Either of:

  - `NULL` (default): the \\Q_3\\ matrix.

  - A single number: the matrix with a global dynamic cutoff, the mean
    off-diagonal \\Q_3\\ plus this number (for example
    `suggested_cutoff` from
    [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)).

  - The full result of
    [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md):
    a list of two tables, the matrix with the global cutoff and a table
    per pair (see Value).

- output:

  Character. `"kable"` (default) for
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) tables,
  or `"dataframe"`.

- n_pairs:

  Integer or `NULL` (default). With the full cutoff object, shows only
  the `n_pairs` pairs that depart most from their expected range in the
  pair table. `NULL` shows all pairs. The correction always covers all
  pairs.

- p_value:

  Logical or `NULL`. Whether to flag pairs on one-sided bootstrap
  p-values (`p_q3`, `padj_q3`) for excess positive dependence. `NULL`
  (default) uses them whenever `cutoff` is the full
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  result, so a numeric cutoff or none keeps the interval. `FALSE` gives
  the behaviour before version 1.2.0. The interval describes where a
  pair's \\Q_3\\ is expected to fall. Used as a decision rule over m
  pairs, a width `w` gives a family-wise error rate of `1 - w^m`, and
  the number of pairs grows with the square of the number of items
  (Johansson, 2026).

- correction:

  Character. Correction of the p-values across all pairs: `"fwer"`
  (default, Westfall-Young step-down), `"fdr_bh"` or `"fdr_by"`
  (Benjamini-Hochberg or Benjamini-Yekutieli), or `"none"`.

- alpha:

  Numeric in (0, 1). A pair is flagged when its corrected p-value is
  below `alpha`. Default `0.05`. Used only with p-values.

- estimator:

  Character. How the residuals are computed. `"CML"` (default) uses CML
  item parameters (`psychotools`) and WLE person locations, which are
  finite at extreme scores. `"MML"` uses marginal maximum likelihood and
  EAP (`mirt`), kept for older analyses. It must match the estimator of
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md).
  With the full cutoff object, its stored estimator is used, with a
  warning if `estimator` differs.

## Value

With `cutoff = NULL` or a number, one table (`knitr_kable` or
data.frame) with the lower triangle of the \\Q_3\\ matrix. A number adds
a flag for pairs above the dynamic cutoff and a caption describing it.

With the full
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
result, a named list:

- `matrix`:

  the \\Q_3\\ matrix with the global dynamic cutoff (mean off-diagonal
  \\Q_3\\ plus `suggested_cutoff`), as above.

- `pairs`:

  one row per item pair with `Item1`, `Item2`, `Observed` (\\Q_3\\),
  `Low` and `High` (the pair's simulated interval) and `Flagged`:
  `"above"` (above the interval, local dependence), `"below"` or `""`.
  Sorted by the departure from the pair's simulated median and cut to
  `n_pairs`. With p-values, also `p_q3` and `padj_q3`, and `Flagged`
  marks `padj_q3 < alpha`, only as `"above"`.

The \\Q_3\\ heatmap is drawn by
[`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md)
(its `matrix` element), so the tables and plots share the `matrix` and
`pairs` structure.

## Details

**Statistic.** \\Q_3\\ (Yen, 1984) is the correlation between the
residuals of two items after the latent variable is accounted for. Under
local independence it is around \\-1/(k-1)\\ for k items. With a
`cutoff`, the dynamic cutoff is the mean of all off-diagonal \\Q_3\\
plus `cutoff` (Christensen et al., 2017).

**Estimation.** \\Q_3\\ is the correlation matrix of the standardized
residuals \\(x - E)/\sqrt{Var}\\. With `"MML"` the residuals come from
`mirt`. The two estimators give very similar \\Q_3\\ (the off-diagonal
values typically correlate above .95). What matters for inference is
that the observed \\Q_3\\ and the simulations in
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
use the same estimator, which the functions enforce.

**Two views.** With the full
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
result, two tables are returned. `matrix` applies one global cutoff, the
99th percentile of the simulated maximum minus mean \\Q_3\\ (Christensen
et al., 2017), an overview of whether there is any local dependence.
`pairs` compares each observed \\Q_3\\ with its own simulated interval,
so that single dependent pairs can be read off and ranked.

**Bootstrap p-values.** Each observed \\Q_3\\ is tested against its
simulated values (`cutoff$pair_results`), one-sided for excess local
dependence. The statistic is the deviation from the mean of the observed
and simulated values together, divided by their SD. The uncorrected p is
`(1 + #{Q3* >= Q3}) / (B + 1)` for B simulated datasets, and
`correction` adjusts across the \\k(k-1)/2\\ pairs (Westfall-Young
step-down, Ferreira, 2024, or FDR). As for item fit, the family-wise
error rate is controlled at any number of iterations, but decisions
depend on the seed. Below 400 iterations a message says so once per
session, and between 400 and 1000 the caption notes it. A false
discovery rate correction needs far more iterations than the family-wise
one, since a bootstrap p-value cannot fall below `1/(B + 1)` and
Benjamini-Hochberg compares the smallest with `alpha/m`. Over the 36
pairs of a nine-item scale that takes 720 iterations, and over the 190
pairs of a twenty-item scale 3799, against 19 for Westfall-Young
whatever the number of pairs.

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

## References

Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
for Yen's \\Q_3\\: Identification of local dependence in the Rasch
model. *Applied Psychological Measurement, 41*(3), 178–194.
[doi:10.1177/0146621616677520](https://doi.org/10.1177/0146621616677520)

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(1), 6.
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

Johansson, M. (2026). Simulation-based cutoffs for conditional item fit
in Rasch models: Iterations, multiplicity correction, and decision
stability. *PsyArXiv*.
[doi:10.31234/osf.io/7pqz4_v2](https://doi.org/10.31234/osf.io/7pqz4_v2)

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement, 8*(2), 125–145.
[doi:10.1177/014662168400800201](https://doi.org/10.1177/014662168400800201)

## See also

[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
for the simulations,
[`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md)
for the plots,
[`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md)
for partial gamma.

## Examples

``` r
# \donttest{
# Simulate binary item response data (10 items, 200 persons)
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Raw Q3 matrix (no cutoff)
RMlocdepQ3(sim_data)
#> 
#> 
#> Table: Raw Q3 residual correlations (lower triangle). Use RMlocdepQ3Cutoff() to derive a cutoff. n = 200 respondents.
#> 
#> |       |Item1 |Item2 |Item3 |Item4 |Item5 |Item6 |Item7 |Item8 |Item9 |Item10 |
#> |:------|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:------|
#> |Item1  |      |      |      |      |      |      |      |      |      |       |
#> |Item2  |-0.02 |      |      |      |      |      |      |      |      |       |
#> |Item3  |-0.17 |-0.2  |      |      |      |      |      |      |      |       |
#> |Item4  |-0.13 |-0.2  |-0.1  |      |      |      |      |      |      |       |
#> |Item5  |-0.19 |-0.06 |-0.02 |-0.12 |      |      |      |      |      |       |
#> |Item6  |-0.1  |-0.1  |-0.17 |0.03  |-0.19 |      |      |      |      |       |
#> |Item7  |-0.17 |-0.13 |-0.07 |0.03  |-0.04 |-0.11 |      |      |      |       |
#> |Item8  |-0.03 |-0.06 |-0.11 |-0.28 |-0.05 |-0.16 |-0.18 |      |      |       |
#> |Item9  |-0.03 |-0.08 |-0.16 |-0.19 |-0.22 |-0.03 |-0.1  |-0.02 |      |       |
#> |Item10 |-0.16 |-0.14 |0     |-0.09 |-0.1  |-0.18 |-0.15 |-0.12 |-0.09 |       |

# Get the underlying data.frame
q3_df <- RMlocdepQ3(sim_data, output = "dataframe")

# Simulation-based cutoff (use 1000 or more iterations in real analyses)
if (requireNamespace("ggdist", quietly = TRUE)) {
  cutoff_res <- RMlocdepQ3Cutoff(sim_data, iterations = 50, parallel = FALSE)

  # A single number gives the matrix only
  RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)

  # The full object gives two tables, $matrix and $pairs
  res <- RMlocdepQ3(sim_data, cutoff = cutoff_res, output = "dataframe")
  res$pairs

  # Top 5 pairs, with bootstrap p-values (use iterations >= 1000 in practice)
  RMlocdepQ3(sim_data, cutoff = cutoff_res, n_pairs = 5, p_value = TRUE,
             output = "dataframe")$pairs
}
#>   Item1  Item2      Observed        Low         High       p_q3   padj_q3
#> 1 Item4  Item8 -0.2765225048 -0.2377678  0.009937381 1.00000000 1.0000000
#> 2 Item4  Item7  0.0285193666 -0.2183735  0.016194476 0.01960784 0.6274510
#> 3 Item5  Item9 -0.2199294805 -0.2001671  0.032761639 1.00000000 1.0000000
#> 4 Item4  Item6  0.0279341896 -0.2592712 -0.003896209 0.01960784 0.6666667
#> 5 Item3 Item10  0.0005837293 -0.2339065  0.003444371 0.05882353 0.9411765
#>   Flagged
#> 1        
#> 2        
#> 3        
#> 4        
#> 5        
# }
```
