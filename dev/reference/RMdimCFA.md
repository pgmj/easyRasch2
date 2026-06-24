# Observed one-factor CFA fit and loadings vs a simulated reference

Fits the observed one-factor categorical CFA to `data` and compares its
fit indices and per-item standardized loadings against the simulated
null distribution produced by
[`RMdimCFACutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md).
Returns a list of two tables: model-fit indices and per-item loadings,
each with the observed value, the expected reference from the
simulation, and a flag.

## Usage

``` r
RMdimCFA(data, cutoff, output = c("kable", "dataframe"))
```

## Arguments

- data:

  A data.frame or matrix of item responses (non-negative integers,
  0-based), the same items used for the cutoff simulation.

- cutoff:

  The list returned by
  [`RMdimCFACutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md).
  Required: observed CFA fit indices are not interpretable without the
  simulated reference, so the function errors if it is missing.

- output:

  Character. `"kable"` (default) returns each table as a
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html);
  `"dataframe"` returns plain data.frames.

## Value

A named list with two elements, `fit` and `loadings`:

- `fit`:

  CFI / RMSEA / SRMR with columns `Index`, `Observed`, `Cutoff`,
  `Direction`, `Flagged` (one-sided, in the unfavourable direction).

- `loadings`:

  One row per item with columns `Item`, `Observed`, `Expected_low`,
  `Expected_high`, `Flagged` (`"below"` / `"above"` / `""`).

Each element is a `knitr_kable` (when `output = "kable"`) or a
data.frame (when `output = "dataframe"`).

## See also

[`RMdimCFACutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md),
[`RMdimCFAPlot`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFAPlot.md)

## Examples

``` r
# \donttest{
if (requireNamespace("lavaan", quietly = TRUE)) {
  data("raschdat1", package = "eRm")
  sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 50,
                        parallel = FALSE, seed = 1)
  tabs <- RMdimCFA(raschdat1[, 1:8], cutoff = sim)
  tabs$fit
  tabs$loadings
}
#> 
#> 
#> Table: Standardized factor loadings (one-factor, lavaan WLSMV; n = 100) vs the simulated expected range under RM unidimensionality (50 iterations at n = 100). Expected range is the two-sided central 99% interval of the simulated loadings; Flagged = below / above that range.
#> 
#> |Item | Observed| Expected low| Expected high|Flagged |
#> |:----|--------:|------------:|-------------:|:-------|
#> |I1   |    0.520|        0.272|         0.938|        |
#> |I2   |    0.441|        0.245|         0.850|        |
#> |I3   |    0.641|        0.249|         0.833|        |
#> |I4   |    0.718|        0.075|         0.873|        |
#> |I5   |    0.303|        0.215|         0.857|        |
#> |I6   |    0.182|        0.209|         0.747|below   |
#> |I7   |    0.419|        0.291|         0.815|        |
#> |I8   |    0.582|        0.250|         0.745|        |
# }
```
