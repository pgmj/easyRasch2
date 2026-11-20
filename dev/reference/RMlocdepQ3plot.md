# Plot the simulated Q3 residual correlations

Plots the simulated \\Q_3\\ of each item pair from
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md),
and optionally the observed \\Q_3\\ in your data and a heatmap of the
observed matrix.

## Usage

``` r
RMlocdepQ3Plot(simfit, data, items = NULL, n_pairs = NULL)
```

## Arguments

- simfit:

  The result of
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md).

- data:

  Optional. The item responses the simulations were made from, items
  scored from 0. The observed \\Q_3\\ is then computed with the
  estimator of `simfit` and drawn as orange diamonds, and the heatmap is
  returned.

- items:

  Optional character vector of item names. Only pairs with both items in
  it are shown. Default `NULL`, all pairs.

- n_pairs:

  Optional positive integer. Only the `n_pairs` pairs with the largest
  departure are shown: `|observed - simulated median|` with `data`,
  `|simulated median|` without. Applied after `items`. Values above the
  number of pairs show all pairs.

## Value

A named list of two `ggplot` objects, with the same structure as the
tables of
[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md):

- `pairs`:

  one row per item pair, labelled "Item1 - Item2".

- `matrix`:

  a heatmap of the lower triangle of the observed \\Q_3\\ matrix, with
  pairs above the global dynamic cutoff outlined. It needs `data`, so it
  is `NULL` (with a message) without it. With `items`, it shows those
  items only. `n_pairs` does not apply to it.

## Details

\\Q_3\\ is symmetric, so each pair is shown once (unlike partial gamma,
which depends on the direction).

Each pair's simulated values are drawn as dots, with the median, an
inner 66% interval and an outer interval. The outer interval is the one
in `simfit$pair_cutoffs`, which
[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
tabulates: the HDCI of width `hdci_width`, or the 2.5th to 97.5th
percentiles with `cutoff_method = "quantile"`. The inner interval uses
the same method. Without `data` the intervals are drawn with
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html).
With `data` they are black segments with a black dot for the median, and
the observed value is an orange diamond.

The intervals are per pair and not corrected for testing many pairs.
With the full cutoff object,
[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)
flags pairs by adjusted p-values instead, so a pair just outside its
interval need not be flagged there. The plot shows no p-values.

The ggplot2, ggdist and scales packages must be installed, and mirt for
`estimator = "MML"`.

## See also

[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md),
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md),
[`RMlocdepGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaPlot.md)

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("ggdist", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)

  # Run the simulation (use 1000 or more iterations in real analyses)
  cutoff_res <- RMlocdepQ3Cutoff(sim_data, iterations = 50,
                                 parallel = FALSE, seed = 42)

  # Simulated distribution only
  RMlocdepQ3Plot(cutoff_res)

  # With observed Q3 overlaid
  RMlocdepQ3Plot(cutoff_res, data = sim_data)

  # Top 10 pairs by departure from null
  RMlocdepQ3Plot(cutoff_res, data = sim_data, n_pairs = 10)
}
#> The Q3 tile heatmap ($matrix) requires `data`; returning $matrix = NULL.
#> $matrix

#> 
#> $pairs

#> 
# }
```
