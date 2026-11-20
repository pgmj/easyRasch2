# Plot the simulated item-restscore differences

Plots each item's simulated observed minus expected item-restscore gamma
from
[`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md),
and optionally the observed differences in your data.

## Usage

``` r
RMitemRestscorePlot(simfit, data)
```

## Arguments

- simfit:

  The result of
  [`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md).

- data:

  Optional. The item responses the simulations were made from, items
  scored from 0. The observed differences are then computed with
  [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
  and drawn as orange diamonds.

## Value

A `ggplot` object.

## Details

The x-axis is the observed minus expected gamma, the statistic
[`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
tests, with a dashed line at zero. Positive values indicate overfit,
negative values underfit. In small samples the simulated differences are
generally not centred on zero, which is one reason the asymptotic test
is miscalibrated (see
[`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)).

Each item's simulated values are drawn as dots, with the median, an
inner 66% interval and an outer interval. The outer interval is the one
stored in `simfit` and tabulated by
[`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md):
the HDCI of width `hdci_width`, or the 2.5th to 97.5th percentiles with
`cutoff_method = "quantile"`. The inner interval uses the same method.
Without `data` the intervals are drawn with
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html).
With `data` they are black segments with a black dot for the median, and
the observed value is the orange diamond.

The intervals are per item and not corrected for testing several items.
With the full cutoff object,
[`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
flags items by adjusted p-values instead, so an item just outside its
interval need not be flagged there. The plot shows no p-values.

The ggplot2 and ggdist packages must be installed (they are in
Suggests).

## See also

[`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md),
[`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)

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

  cutoff_res <- RMitemRestscoreCutoff(sim_data, iterations = 100,
                                      parallel = FALSE, seed = 42)

  # Simulated differences only
  RMitemRestscorePlot(cutoff_res)

  # With the observed differences overlaid
  RMitemRestscorePlot(cutoff_res, data = sim_data)
}
#> 

# }
```
