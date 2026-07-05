# Conditional Item Infit MSQ

Computes conditional infit mean-square (MSQ) statistics for each item
using
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html),
enriched with item locations relative to the sample mean person
location.

## Usage

``` r
RMitemInfit(
  data,
  cutoff = NULL,
  p_value = FALSE,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable",
  sort
)
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
    [`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoff.md)
    (a list with `$item_cutoffs`): the data.frame is extracted
    automatically and the number of simulation iterations,
    `cutoff_method`, and `hdci_width` are included in the kable caption.

  - The `$item_cutoffs` data.frame from
    [`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoff.md)
    directly: must have columns `Item`, `infit_low`, and `infit_high`.
    When provided, adds columns `Infit_low`, `Infit_high`, and `Flagged`
    to the result. `Flagged` labels the misfit direction: `"overfit"`
    (infit below the range – more predictable than the model expects),
    `"underfit"` (above the range – noisier than expected), or `""`
    (within range, no misfit).

- p_value:

  Logical. If `TRUE`, bootstrap p-values are computed from the simulated
  null distribution and added to the output. This requires `cutoff` to
  be the **full**
  [`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoff.md)
  object (which carries the per-item simulated values in its `$results`
  element); the summarised `$item_cutoffs` data.frame is not sufficient.
  Default `FALSE`, in which case behaviour is unchanged.

- correction:

  Character. Multiple-comparison correction applied across items when
  `p_value = TRUE`: `"fwer"` (default) for the Westfall-Young
  studentised-max step-down (family-wise error rate), `"fdr_bh"` /
  `"fdr_by"` for Benjamini-Hochberg / Benjamini-Yekutieli false
  discovery rate control, or `"none"` for uncorrected per-item p-values.

- alpha:

  Numeric in (0, 1). Significance level for the `Flagged` column when
  `p_value = TRUE` (an item is flagged when its corrected p-value is
  below `alpha`). Default `0.05`.

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
  (inserted after `Infit_MSQ`, before `Relative_location`). `Flagged` is
  a character column with values `"overfit"`, `"underfit"`, or `""` (not
  the previous logical). When `p_value = TRUE`, columns `p_infit`
  (marginal two-sided p-value) and `padj_infit` (corrected p-value) are
  added and `Flagged` reflects items with `padj_infit < alpha`
  (direction from `Infit_MSQ` relative to 1).

## Details

Infit MSQ is a weighted fit statistic that emphasises deviations near
the item location. Values close to 1.0 indicate good fit. Values
substantially above 1.0 suggest underfit (unexpected responses), while
values substantially below 1.0 suggest overfit (overly predictable
responses). The definition of "substantially" depends on several factors
such as sample size, and needs to be determined by simulation using
[`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoff.md).
There is no general rule-of-thumb value that is correct.

Conditional infit MSQ statistics are computed via
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html),
which uses the conditional distribution of the sufficient statistics
(Müller, 2020). Only complete cases (rows without any `NA`) are used in
the conditional fit calculation.

Item parameters are estimated by conditional maximum likelihood via
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category PCM); the conditional infit/outfit
MSQ comes from
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html) and
is invariant to the estimation engine. Per-item average locations are
the means of the CML thresholds, and the person-location reference is
the mean of the Warm WLE estimates.

Relative item location is defined as the item's average location minus
the sample mean person location, providing a measure of item targeting.

The `iarm` package must be installed (it is in Suggests, not Imports).

**Bootstrap p-values.** When `p_value = TRUE`, each item's observed
infit is compared against its simulated null distribution (from
`cutoff$results`). The per-item statistic is the residual studentised by
the bootstrap mean and SD – deliberately the empirical SD rather than
the Wilson-Hilferty / ZSTD transform, which is uninformative for
conditional MSQ (Müller, 2020). The marginal p-value is the two-sided
Monte-Carlo p-value `(1 + #{|t*| >= |t|}) / (B + 1)`. For
`correction = "fwer"` the family-wise adjustment uses the Westfall-Young
studentised-max step-down, which exploits the bootstrap dependence among
items and is more powerful than Bonferroni/Holm (Ferreira, 2024); its
validity rests on subset pivotality. `"fdr_bh"`/`"fdr_by"` apply
Benjamini-Hochberg / Benjamini- Yekutieli instead. These are
model-conditional, sample-size-sensitive p-values and are reported
alongside the simulated effect-size band, not in place of it. p-values
can be no smaller than `1 / (B + 1)`, and the studentised-max (FWER)
correction is *liberal* when the simulation is small (the bootstrap
mean/SD used for studentisation are then too noisy). At least 1000
`iterations` in
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoff.md)
are recommended – in simulations the family-wise error rate is then
controlled at the nominal level – and a warning is issued when the
simulation is smaller. The marginal p-values are well calibrated even at
a few hundred iterations.

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

## References

Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust
them? *Journal of Statistical Distributions and Applications*, 7(5).
[doi:10.1186/s40488-020-00108-7](https://doi.org/10.1186/s40488-020-00108-7)

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(6).
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.

## See also

[`RMitemInfitCutoff`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoff.md)

## Examples

``` r
# \donttest{
if (requireNamespace("iarm", quietly = TRUE)) {
  # Simulate binary item response data (5 items, 40 persons)
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 40 * 5, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(sim_data) <- paste0("Item", 1:5)

  # Default kable output
  RMitemInfit(sim_data)

  # Sorted by infit MSQ descending
  RMitemInfit(sim_data, sort = "infit")

  # Return as data.frame for further processing
  df <- RMitemInfit(sim_data, output = "dataframe")

  # Simulation-based cutoffs (100 Monte-Carlo iterations)
  if (requireNamespace("ggdist", quietly = TRUE)) {
    cutoff_res <- RMitemInfitCutoff(sim_data, iterations = 100,
                                    parallel = FALSE, seed = 42)
    RMitemInfit(sim_data, cutoff = cutoff_res)
    RMitemInfit(sim_data, cutoff = cutoff_res, output = "dataframe")

    # Bootstrap p-values with family-wise (Westfall-Young) correction
    # (use iterations >= 1000 in real analyses for stable p-values)
    RMitemInfit(sim_data, cutoff = cutoff_res, p_value = TRUE,
                output = "dataframe")
  }
}
#> Warning: Bootstrap p-values are based on only 100 simulation iterations. With few iterations the studentised-max (FWER) correction is liberal and small p-values are imprecise; use iterations >= 1000 in RMitemInfitCutoff() for reliable p-values.
#>    Item Infit_MSQ Infit_low Infit_high p_infit padj_infit Flagged
#> 1 Item1     1.049     0.607      1.577  0.5842     0.9208        
#> 2 Item2     0.929     0.625      1.519  0.6436     0.9208        
#> 3 Item3     0.828     0.674      1.354  0.1980     0.5347        
#> 4 Item4     1.216     0.643      1.434  0.0792     0.2970        
#> 5 Item5     0.931     0.648      1.353  0.6634     0.9208        
#>   Relative_location
#> 1              0.04
#> 2             -0.48
#> 3             -0.16
#> 4              0.33
#> 5             -0.70
# }
```
