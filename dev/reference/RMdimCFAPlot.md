# Plot one-factor CFA fit and loadings against the simulations

Returns two figures comparing the observed one-factor CFA with the
simulations from
[`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md):
the standardized loadings per item, with each item's simulated values
and the observed loading, in the style of
[`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md),
and the simulated CFI, RMSEA and SRMR with the observed values.

## Usage

``` r
RMdimCFAPlot(simfit, data, percentile = NULL)
```

## Arguments

- simfit:

  The result of
  [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md).

- data:

  The item responses the simulations were made from. Required, since the
  observed values are computed from it.

- percentile:

  Numeric in (50, 100) or `NULL`. When given, the cutoffs and flags are
  recomputed at this percentile from the stored simulations, without
  simulating again. Default `NULL`, the percentile of the
  [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md)
  call.

## Value

A named list of two `ggplot` objects:

- `loadings`:

  Standardized loadings per item: the simulated values (dots), black
  segments for their central 66% and the expected range at `percentile`,
  the median (black dot), and the observed loading as a diamond, red
  when it is outside the expected range.

- `fit`:

  The simulated CFI, RMSEA and SRMR in panels, with the observed value
  and the cutoff.

## Details

The outer segment in the loadings plot is the expected range that
[`RMdimCFA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFA.md)
flags against: the central interval of the simulated loadings covering
`percentile`%. The inner segment is the central 66%. The ranges are per
item and not corrected for testing several items.

## See also

[`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md),
[`RMdimCFA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFA.md)

## Examples

``` r
# \donttest{
if (requireNamespace("lavaan", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("eRm", quietly = TRUE)) {
  data("raschdat1", package = "eRm")
  sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 50,
                        parallel = FALSE, seed = 1)
  plots <- RMdimCFAPlot(sim, data = raschdat1[, 1:8])
  plots$loadings
  plots$fit
}

# }
```
