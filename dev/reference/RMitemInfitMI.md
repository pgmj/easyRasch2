# Conditional Item Infit MSQ for Multiply Imputed Data

Extends
[`RMitemInfit`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
to work with multiply imputed datasets produced by the `mice` package.
Computes conditional infit MSQ on each imputed dataset and pools the
results using Rubin's rules.

## Usage

``` r
RMitemInfitMI(mids_object, cutoff = NULL, output = "kable", sort)
```

## Arguments

- mids_object:

  A `mids` object (multiply imputed dataset) as returned by
  [`mice::mice()`](https://amices.org/mice/reference/mice.html). Each
  completed dataset must contain only the item response columns to be
  analysed (i.e., no ID or grouping variables). Items must be scored
  starting at 0 (non-negative integers).

- cutoff:

  Optional. Default `NULL` (no cutoff applied). Can be:

  - The return value of
    [`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
    or
    [`RMitemInfitCutoffMI`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md)
    (a list with `$item_cutoffs`): the data.frame is extracted
    automatically and metadata is included in the kable caption.

  - The `$item_cutoffs` data.frame directly: must have columns `Item`,
    `infit_low`, and `infit_high`. When provided, adds columns
    `Infit_low`, `Infit_high`, and `Flagged` to the result.

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
  `format = "pipe"`) with columns "Item", "Infit MSQ", "Infit SE",
  "Relative location", and a caption noting the number of imputations
  and complete cases. When `cutoff` is provided, columns "Infit low",
  "Infit high", and "Flagged" are also included.

- If `output = "dataframe"`: a data.frame with columns `Item`,
  `Infit_MSQ`, `Infit_SE`, and `Relative_location`. When `cutoff` is
  provided, columns `Infit_low`, `Infit_high`, and `Flagged` are also
  included (inserted after `Infit_SE`, before `Relative_location`).

## Details

For each of the `m` imputed datasets, the function:

1.  Fits a Rasch model
    ([`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) for dichotomous
    data or [`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html) for
    polytomous data).

2.  Computes conditional infit MSQ and its standard error via
    [`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html).

3.  Computes item and person locations.

The per-imputation estimates are then pooled using Rubin's rules:

- Pooled MSQ:

  The mean of the `m` infit MSQ point estimates.

- Within-imputation variance:

  The mean of the `m` squared standard errors.

- Between-imputation variance:

  The sample variance of the `m` point estimates.

- Total variance:

  Within + (1 + 1/m) \* Between.

- Pooled SE:

  The square root of the total variance.

Relative item location is the mean of per-imputation relative locations
(item location minus sample mean person location).

Imputed datasets that cause model convergence failures are dropped with
a warning. If all imputations fail, the function stops with an error. At
least two successful imputations are required to estimate
between-imputation variance.

The `mice` and `iarm` packages must be installed (they are in Suggests,
not Imports).

## See also

[`RMitemInfit`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md),
[`RMitemInfitCutoffMI`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md)

## Examples

``` r
# \donttest{
if (requireNamespace("mice", quietly = TRUE)) {
  # Create example data with missing values
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
  )
  colnames(sim_data) <- paste0("Item", 1:8)
  sim_data[sample(length(sim_data), 0.10 * length(sim_data))] <- NA

  # Impute (use more imputations, e.g. m = 5+, in real analyses)
  imp <- mice::mice(sim_data, m = 2, method = "polr", seed = 123,
                    printFlag = FALSE)

  # Pooled infit table (no cutoffs)
  RMitemInfitMI(imp)

  # With simulation-based cutoffs
  # (use more iterations, e.g. 250+, in real analyses)
  cutoff_mi <- RMitemInfitCutoffMI(imp, iterations = 50, parallel = FALSE,
                                seed = 42)
  RMitemInfitMI(imp, cutoff = cutoff_mi)

  # As data.frame
  df <- RMitemInfitMI(imp, cutoff = cutoff_mi, output = "dataframe")
}
# }
```
