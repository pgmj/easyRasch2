# Item Restscore Analysis

Computes observed and model-expected item-restscore correlations using
[`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html),
and enriches the output with the absolute difference between observed
and expected values, item average locations, and item locations relative
to the sample mean person location.

## Usage

``` r
RMitemrestscore(data, output = "kable", sort, p.adj = "BH")
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
  `Absolute_difference` in descending order before output.

- p.adj:

  Character string specifying the p-value adjustment method passed to
  [`iarm::item_restscore()`](https://rdrr.io/pkg/iarm/man/item_restscore.html).
  Default `"BH"` (Benjamini-Hochberg). See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
  available methods.

## Value

- If `output = "kable"`: a `knitr_kable` object (plain text table via
  `format = "pipe"`) with columns for item name, observed and expected
  restscore correlations, absolute difference, adjusted p-value,
  significance level, item location, and item location relative to the
  sample mean person location.

- If `output = "dataframe"`: a data.frame with columns `Item`,
  `Observed`, `Expected`, `Absolute_difference`, `p_adjusted`,
  `Significance`, `Location`, and `Relative_location`.

## Details

Item-restscore correlations using Goodman-Kruskal's gamma (Kreiner,
2011) measure the association between a person's score on a single item
and their total score on the remaining items (the "restscore"). Under a
correctly fitting Rasch model, observed and model-expected correlations
should agree closely.

For **dichotomous** data (maximum score = 1), a Rasch model is fitted
via [`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html). Item locations
are the negative beta parameters. Person locations are estimated via
[`eRm::person.parameter()`](https://rdrr.io/pkg/eRm/man/person.parameter.html).

For **polytomous** data (maximum score \> 1), a Partial Credit Model is
fitted via [`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html). Item
average locations are the row-means of the threshold parameter table
returned by
[`eRm::thresholds()`](https://rdrr.io/pkg/eRm/man/thresholds.html).
Person locations are estimated via
[`eRm::person.parameter()`](https://rdrr.io/pkg/eRm/man/person.parameter.html).

Relative item location is defined as the item's average location minus
the sample mean person location, providing a measure of item targeting.

The `iarm` package must be installed (it is in Suggests, not Imports).

## References

Kreiner, S. (2011). A Note on Item–Restscore Association in Rasch
Models. *Applied Psychological Measurement, 35*(7), 557–561.
<https://doi.org/10.1177/0146621611410227>

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate binary item response data (8 items, 200 persons)
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
)
colnames(sim_data) <- paste0("Item", 1:8)

# Default kable output
RMitemrestscore(sim_data)

# Sorted by absolute difference
RMitemrestscore(sim_data, sort = "diff")

# Return as data.frame for further processing
df <- RMitemrestscore(sim_data, output = "dataframe")
} # }
```
