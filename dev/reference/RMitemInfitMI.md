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
    `Infit_low`, `Infit_high`, and `Flagged` to the result. `Flagged` is
    a character column labelling the misfit direction: `"overfit"`
    (pooled infit below the range), `"underfit"` (above), or `""`
    (within range).

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
  `Flagged` is a character column (`"overfit"` / `"underfit"` / `""`),
  not the previous logical.

## Details

For each of the `m` imputed datasets, the function:

1.  Fits a Rasch model by CML via
    [`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
    (a dichotomous item is a 2-category partial credit model),
    consistent with
    [`RMitemInfit`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
    and the rest of the package.

2.  Computes conditional infit MSQ and its standard error via
    [`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html).

3.  Computes item locations (mean of the grand-mean-centred CML Andrich
    thresholds) and the mean WLE person location.

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

**Caveat on the pooled SE.** The within-imputation variance is the
squared conditional infit SE from
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html).
Müller (2020) showed that this asymptotic SE is an unreliable measure of
uncertainty for the conditional infit statistic; Rubin's pooled SE
inherits that limitation, so the `Infit_SE`/`Infit SE` column should be
read as an approximate indication of imputation-related variability
rather than a trustworthy inferential standard error. For item misfit
decisions, prefer the simulation-based cutoffs from
[`RMitemInfitCutoffMI`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md).

Imputed datasets that cause model convergence failures are dropped with
a warning. If all imputations fail, the function stops with an error. At
least two successful imputations are required to estimate
between-imputation variance.

The `mice` and `iarm` packages must be installed (they are in Suggests,
not Imports).

## References

Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust
them? *Journal of Statistical Distributions and Applications*, 7(5).
[doi:10.1186/s40488-020-00108-7](https://doi.org/10.1186/s40488-020-00108-7)

## See also

[`RMitemInfit`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md),
[`RMitemInfitCutoffMI`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md)

## Examples

``` r
# \donttest{
if (requireNamespace("mice", quietly = TRUE) &&
    requireNamespace("iarm", quietly = TRUE) &&
    requireNamespace("ggdist", quietly = TRUE)) {
  # Create example data with ~10% MCAR missingness
  set.seed(42)
  mat <- matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
  mat[sample(length(mat), round(0.10 * length(mat)))] <- NA
  sim_data <- as.data.frame(mat)
  colnames(sim_data) <- paste0("Item", 1:8)

  # mice's ordinal method (`polr`) requires the items to be ordered
  # factors, so code them as such before imputing. RMitemInfitMI()
  # converts the completed factors back to numeric internally.
  sim_data[] <- lapply(sim_data, function(x) factor(x, ordered = TRUE))

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
