# Plot the simulated partial gamma local dependence values

Plots the simulated partial gamma of each item pair from
[`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md),
and optionally the observed values in your data.

## Usage

``` r
RMlocdepGammaPlot(simfit, data, items = NULL, n_pairs = NULL)
```

## Arguments

- simfit:

  The result of
  [`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md).

- data:

  Optional. The item responses the simulations were made from, items
  scored from 0. The observed partial gamma (complete cases) is then
  drawn as orange diamonds.

- items:

  Optional character vector of item names. Only pairs with both items in
  it are shown. Default `NULL`, all pairs.

- n_pairs:

  Optional positive integer. Only the `n_pairs` pairs with the largest
  absolute partial gamma are shown, largest first: the observed values
  with `data`, the simulated medians without. Applied after `items`.
  Values above the number of pairs show all pairs.

## Value

A `ggplot` object.

## Details

One row per item pair, labelled "Item1 - Item2". The observed value is
the larger of the pair's two rest-score directions, the statistic
[`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md)
tests, computed with the same internal routine as the simulations.

Each pair's simulated values are drawn as dots, with the median, an
inner 66% interval and an outer interval. The outer interval is the one
in `simfit$pair_cutoffs`, which
[`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md)
tabulates: the HDCI of width `hdci_width`, or the 2.5th to 97.5th
percentiles with `cutoff_method = "quantile"`. The inner interval uses
the same method. Without `data` the intervals are drawn with
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html).
With `data` they are black segments with a black dot for the median, and
the observed value is an orange diamond.

The intervals are per pair and not corrected for testing many pairs.
With the full cutoff object,
[`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md)
flags pairs by adjusted p-values instead, so a pair just outside its
interval need not be flagged there. The plot shows no p-values.

The ggplot2 and ggdist packages must be installed, and iarm with `data`
(all are in Suggests).

## See also

[`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md),
[`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md),
[`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md)

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
  cutoff_res <- RMlocdepGammaCutoff(sim_data, iterations = 100,
                                    parallel = FALSE, seed = 42)

  # Simulated distribution only
  RMlocdepGammaPlot(cutoff_res)

  # With observed partial gamma overlaid
  RMlocdepGammaPlot(cutoff_res, data = sim_data)

  # Plot only a subset of items
  RMlocdepGammaPlot(cutoff_res, data = sim_data,
                    items = c("Item1", "Item2", "Item3"))
}

# }
```
