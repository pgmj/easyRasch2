# Item Restscore Analysis

Computes observed and model-expected item-restscore correlations using
[`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html),
and enriches the output with the absolute difference between observed
and expected values, item average locations, and item locations relative
to the sample mean person location.

## Usage

``` r
RMitemRestscore(data, output = "kable", sort, p_adj = "BH")
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed, but at least one complete case (row with no `NA`) must be
  present.

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
  Default `"BH"` (Benjamini-Hochberg). See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
  available methods.

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

## References

Kreiner, S. (2011). A Note on Item–Restscore Association in Rasch
Models. *Applied Psychological Measurement, 35*(7), 557–561.
[doi:10.1177/0146621611410227](https://doi.org/10.1177/0146621611410227)

## Examples

``` r
# \donttest{
# Simulate binary item response data (8 items, 200 persons)
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
)
colnames(sim_data) <- paste0("Item", 1:8)

# Default kable output
RMitemRestscore(sim_data)
#> 
#> 
#> 
#> Table: Item-restscore associations (n = 200 complete cases). Flagged (adj. p < .05): overfit = observed above expected (over-discrimination, often local dependence); underfit = below (under-discrimination, often multidimensionality/noise).
#> 
#> |Item  | Observed| Expected| Difference| Adj. p-value (BH)|Flagged | Rel. location|
#> |:-----|--------:|--------:|----------:|-----------------:|:-------|-------------:|
#> |Item1 |     0.00|     0.01|      -0.01|             0.955|        |         -0.23|
#> |Item2 |     0.00|     0.01|      -0.01|             0.955|        |          0.13|
#> |Item3 |    -0.01|     0.01|      -0.02|             0.955|        |         -0.01|
#> |Item4 |    -0.05|     0.01|      -0.06|             0.955|        |         -0.05|
#> |Item5 |     0.07|     0.01|       0.06|             0.955|        |          0.09|
#> |Item6 |    -0.04|     0.01|      -0.05|             0.955|        |         -0.09|
#> |Item7 |     0.15|     0.01|       0.14|             0.955|        |         -0.05|
#> |Item8 |    -0.01|     0.01|      -0.02|             0.955|        |         -0.03|

# Sorted by absolute difference
RMitemRestscore(sim_data, sort = "diff")
#> 
#> 
#> 
#> Table: Item-restscore associations (n = 200 complete cases). Flagged (adj. p < .05): overfit = observed above expected (over-discrimination, often local dependence); underfit = below (under-discrimination, often multidimensionality/noise).
#> 
#> |Item  | Observed| Expected| Difference| Adj. p-value (BH)|Flagged | Rel. location|
#> |:-----|--------:|--------:|----------:|-----------------:|:-------|-------------:|
#> |Item7 |     0.15|     0.01|       0.14|             0.955|        |         -0.05|
#> |Item4 |    -0.05|     0.01|      -0.06|             0.955|        |         -0.05|
#> |Item5 |     0.07|     0.01|       0.06|             0.955|        |          0.09|
#> |Item6 |    -0.04|     0.01|      -0.05|             0.955|        |         -0.09|
#> |Item3 |    -0.01|     0.01|      -0.02|             0.955|        |         -0.01|
#> |Item8 |    -0.01|     0.01|      -0.02|             0.955|        |         -0.03|
#> |Item1 |     0.00|     0.01|      -0.01|             0.955|        |         -0.23|
#> |Item2 |     0.00|     0.01|      -0.01|             0.955|        |          0.13|

# Return as data.frame for further processing
df <- RMitemRestscore(sim_data, output = "dataframe")
#> 
# }
```
