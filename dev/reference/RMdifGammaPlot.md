# Plot the simulated partial gamma DIF values

Plots the simulated partial gamma of each item from
[`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md),
and optionally the observed values in your data, computed with the same
routine as the simulations, which reproduces
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).

## Usage

``` r
RMdifGammaPlot(simfit, data, dif_var)
```

## Arguments

- simfit:

  The result of
  [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md).

- data:

  Optional. The item responses the simulations were made from, items
  scored from 0. The observed partial gamma is then drawn as orange
  diamonds.

- dif_var:

  Needed with `data`. The vector of groups used in
  [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md),
  of length `nrow(data)`.

## Value

A `ggplot` object.

## Details

Each item's simulated values are drawn as dots, with the median, an
inner 66% interval and an outer interval. The outer interval is the one
in `simfit$item_cutoffs`, which
[`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md)
tabulates: the HDCI of width `hdci_width`, or the 2.5th to 97.5th
percentiles with `cutoff_method = "quantile"`. The inner interval uses
the same method. Without `data` the intervals are drawn with
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html).
With `data` they are black segments with a black dot for the median, and
the observed value is an orange diamond.

The intervals are per item and not corrected for testing several items.
With the full cutoff object,
[`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md)
flags items by adjusted p-values instead, so an item just outside its
interval need not be flagged there. The plot shows no p-values.

The ggplot2 and ggdist packages must be installed (they are in
Suggests).

## See also

[`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md),
[`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md)

## Examples

``` r
# \donttest{
if (requireNamespace("ggdist", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)
  dif_group <- factor(sample(c("A", "B"), 200, replace = TRUE))

  # Run the simulation
  cutoff_res <- RMdifGammaCutoff(sim_data, dif_var = dif_group,
                                 iterations = 100, parallel = FALSE,
                                 seed = 42)

  # Simulated distribution only
  RMdifGammaPlot(cutoff_res)

  # With observed partial gamma overlaid
  RMdifGammaPlot(cutoff_res, data = sim_data, dif_var = dif_group)
}

# }
```
