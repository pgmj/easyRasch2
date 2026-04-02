# Conditional Item Infit MSQ

Computes conditional infit mean-square (MSQ) statistics for each item
using
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html),
enriched with item locations relative to the sample mean person
location.

## Usage

``` r
RMiteminfit(data, cutoff = NULL, output = "kable", sort)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed, but at least one complete case (row with no `NA`) must be
  present.

- cutoff:

  Optional. Default `NULL` (no cutoff applied, behaviour is identical to
  the current version). Can be:

  - The return value of
    [`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md)
    (a list with `$item_cutoffs`): the data.frame is extracted
    automatically and the number of simulation iterations,
    `cutoff_method`, and `hdci_width` are included in the kable caption.

  - The `$item_cutoffs` data.frame from
    [`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md)
    directly: must have columns `Item`, `infit_low`, and `infit_high`.
    When provided, adds columns `Infit_low`, `Infit_high`, and `Flagged`
    (logical; `TRUE` when `Infit_MSQ` falls outside the credible range)
    to the result.

- output:

  Character string controlling the return value. Either `"kable"`
  (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"` for the underlying data.frame.

- sort:

  Optional character string. When `sort = "infit"`, rows are sorted by
  `Infit_MSQ` in descending order before output.

## Value

- If `output = "kable"`: a `knitr_kable` object (plain text table via
  `format = "pipe"`) with columns "Item", "Infit MSQ", and "Relative
  location", and a caption showing the number of complete cases. When
  `cutoff` is provided, columns "Infit low", "Infit high", and "Flagged"
  are also included, and the caption notes the simulation-based cutoffs.

- If `output = "dataframe"`: a data.frame with columns `Item`,
  `Infit_MSQ`, and `Relative_location`. When `cutoff` is provided,
  columns `Infit_low`, `Infit_high`, and `Flagged` are also included
  (inserted after `Infit_MSQ`, before `Relative_location`).

## Details

Infit MSQ is a weighted fit statistic that emphasises deviations near
the item location. Values close to 1.0 indicate good fit. Values
substantially above 1.0 suggest underfit (unexpected responses), while
values substantially below 1.0 suggest overfit (overly predictable
responses). The definition of "substantially" depends on several factors
such as sample size, and needs to be determined by simulation using
[`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md).
There is no general rule-of-thumb value that is correct.

Conditional infit MSQ statistics are computed via
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html),
which uses the conditional distribution of the sufficient statistics
(Müller, 2020). Only complete cases (rows without any `NA`) are used in
the conditional fit calculation.

For **dichotomous** data (maximum score = 1), a Rasch model is fitted
via [`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html). Item locations
are the negative beta parameters. Person locations are estimated via
[`eRm::person.parameter()`](https://rdrr.io/pkg/eRm/man/person.parameter.html).

For **polytomous** data (maximum score \> 1), a Partial Credit Model is
fitted via [`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html). Item
average locations are taken from the "Location" column of the threshold
parameter table returned by
[`eRm::thresholds()`](https://rdrr.io/pkg/eRm/man/thresholds.html); if
that column is absent, row means of the threshold columns are used
instead. Person locations are estimated via
[`eRm::person.parameter()`](https://rdrr.io/pkg/eRm/man/person.parameter.html).

Relative item location is defined as the item's average location minus
the sample mean person location, providing a measure of item targeting.

The `iarm` package must be installed (it is in Suggests, not Imports).

## References

Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust
them? *Journal of Statistical Distributions and Applications*, 7(5).
[doi:10.1186/s40488-020-00108-7](https://doi.org/10.1186/s40488-020-00108-7)

## See also

[`RMinfitcutoff`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate binary item response data (5 items, 40 persons)
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 40 * 5, replace = TRUE), nrow = 40, ncol = 5)
)
colnames(sim_data) <- paste0("Item", 1:5)

# Default kable output
RMiteminfit(sim_data)

# Sorted by infit MSQ descending
RMiteminfit(sim_data, sort = "infit")

# Return as data.frame for further processing
df <- RMiteminfit(sim_data, output = "dataframe")

# With simulation-based cutoffs from RMinfitcutoff()
cutoff_res <- RMinfitcutoff(sim_data, iterations = 100, parallel = FALSE,
                            seed = 42)
RMiteminfit(sim_data, cutoff = cutoff_res)
RMiteminfit(sim_data, cutoff = cutoff_res, output = "dataframe")
} # }
```
