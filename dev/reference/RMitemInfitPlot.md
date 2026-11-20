# Plot the simulated infit and outfit MSQ values

Plots the simulated conditional item fit values from
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
(or
[`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md))
per item, and optionally the observed values in your data.

## Usage

``` r
RMitemInfitPlot(simfit, data, statistic = "infit")

RMitemInfitCutoffPlot(...)
```

## Arguments

- simfit:

  The result of
  [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
  or
  [`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md).

- data:

  Optional. The item responses the simulations were made from, items
  scored from 0. The observed conditional infit and outfit MSQ are then
  computed (on complete cases) and drawn as orange diamonds.

- statistic:

  Character. `"infit"`, `"outfit"`, or `"both"` side by side. With
  `data` the default is `"infit"`, and `"both"` needs patchwork. Without
  `data` the default is `"both"`, drawn as two facets.

- ...:

  Arguments passed on to `RMitemInfitPlot()`.

## Value

A `ggplot` object, or a `patchwork` object with `statistic = "both"` and
`data`.

## Details

Each item's simulated values are drawn as dots, with the median, an
inner 66% interval and an outer interval. The outer interval is the one
stored in `simfit` and tabulated by
[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md):
the HDCI of width `hdci_width`, or the 2.5th to 97.5th percentiles with
`cutoff_method = "quantile"`. The inner interval uses the same method.

**Without `data`** the intervals are drawn with
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html).

**With `data`**, the model is fitted by CML
([`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html))
and the observed conditional MSQ computed with
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html). The
intervals are black segments with a black dot for the median, and the
observed value is the orange diamond.

The intervals are per item and not corrected for testing several items.
With the full cutoff object,
[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
flags items by adjusted p-values instead, so an item just outside its
interval need not be flagged there. The plot shows no p-values.

The ggplot2 and ggdist packages must be installed, and patchwork for
`statistic = "both"` (all are in Suggests).

`RMitemInfitCutoffPlot()` is a deprecated alias for `RMitemInfitPlot()`,
kept for code written for easyRasch2 0.8.0. It warns and calls
`RMitemInfitPlot()`.

## See also

[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md),
[`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)

## Examples

``` r
# \donttest{
if (requireNamespace("iarm", quietly = TRUE) &&
    requireNamespace("ggdist", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)

  # Run simulation
  cutoff_res <- RMitemInfitCutoff(sim_data, iterations = 100,
                                  parallel = FALSE, seed = 42)

  # Simulated values only (infit and outfit side by side)
  RMitemInfitPlot(cutoff_res)

  # With observed fit overlaid (infit only, the default)
  RMitemInfitPlot(cutoff_res, data = sim_data)

  # Both infit and outfit panels side by side
  if (requireNamespace("patchwork", quietly = TRUE)) {
    RMitemInfitPlot(cutoff_res, data = sim_data, statistic = "both")
  }
}

# }
```
