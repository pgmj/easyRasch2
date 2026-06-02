# Plot observed CFA fit indices against the simulated null distribution

Faceted histograms of the simulated CFI, RMSEA, and SRMR distributions
(one panel per index), with the observed value overlaid as a vertical
line. The line is coloured red when the observed value falls outside the
chosen percentile of the simulation in the unfavourable direction (CFI
from below; RMSEA / SRMR from above), and grey otherwise.

## Usage

``` r
RMdimCFAPlot(cutoff_res, percentile = NULL)
```

## Arguments

- cutoff_res:

  The list returned by
  [`RMdimCFACutoff`](https://pgmj.github.io/easyRasch2/reference/RMdimCFACutoff.md)
  with `output = "list"`, or the kable returned with `output = "kable"`
  (the underlying list is read from `attr(., "result")`).

- percentile:

  Numeric in (50, 100) or `NULL`. When supplied, the cutoff and the
  flagged status are recomputed at this percentile from the simulated
  distribution stored in `cutoff_res` (no re-simulation needed). When
  `NULL` (default), the percentile that was used by the original
  [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/reference/RMdimCFACutoff.md)
  call is reused.

## Value

A `ggplot` object.

## See also

[`RMdimCFACutoff`](https://pgmj.github.io/easyRasch2/reference/RMdimCFACutoff.md)

## Examples

``` r
# \donttest{
if (requireNamespace("lavaan", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  data("raschdat1", package = "eRm")
  res <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 50,
                     parallel = FALSE, seed = 1, output = "list")
  RMdimCFAPlot(res)                 # use the percentile from RMdimCFACutoff()
  RMdimCFAPlot(res, percentile = 95) # override (no re-simulation needed)
}

# }
```
