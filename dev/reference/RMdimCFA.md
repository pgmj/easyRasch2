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
RMdimCFA(
  data,
  cutoff,
  p_value = FALSE,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = c("kable", "dataframe")
)
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

- p_value:

  Logical. When `TRUE`, adds bootstrap p-values from the simulated null
  distributions: one-sided in the unfavourable direction for the fit
  indices (CFI low; RMSEA / SRMR high), two-sided for the per-item
  loadings. The `Flagged` columns then reflect `padj < alpha` instead of
  the percentile cutoffs. The fit indices and the loadings are corrected
  as two separate families. Default `FALSE`.

- correction:

  Character. Multiplicity correction for the p-values: `"fwer"`
  (default; Westfall-Young studentised-max step-down), `"fdr_bh"`
  (Benjamini-Hochberg), `"fdr_by"` (Benjamini-Yekutieli), or `"none"`.
  Ignored when `p_value = FALSE`.

- alpha:

  Numeric in (0, 1). Significance level used to flag comparisons on the
  corrected p-value. Default `0.05`. Ignored when `p_value = FALSE`.

- output:

  Character. `"kable"` (default) returns each table as a
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html);
  `"dataframe"` returns plain data.frames.

## Value

A named list with two elements, `fit` and `loadings`:

- `fit`:

  CFI / RMSEA / SRMR with columns `Index`, `Observed`, `Cutoff`,
  `Direction`, `Flagged` (one-sided, in the unfavourable direction).
  With `p_value = TRUE`, columns `p` and `padj` are added and `Flagged`
  reflects `padj < alpha`.

- `loadings`:

  One row per item with columns `Item`, `Observed`, `Expected_low`,
  `Expected_high`, `Flagged` (`"below"` / `"above"` / `""`). With
  `p_value = TRUE`, columns `p_loading` and `padj_loading` are added and
  `Flagged` reflects `padj_loading < alpha` (direction from the sign of
  the deviation from the simulated mean).

Each element is a `knitr_kable` (when `output = "kable"`) or a
data.frame (when `output = "dataframe"`).

## Details

**Bootstrap p-values.** The per-comparison statistic is the residual
studentised by the bootstrap mean and SD. Marginal p-values are
Monte-Carlo, `(1 + count) / (B + 1)`, so they can be no smaller than
`1 / (B + 1)`. `correction = "fwer"` uses the Westfall-Young
studentised-max step-down, which exploits the bootstrap dependence among
the statistics (Ferreira, 2024); it is liberal when the simulation is
small, so at least 1000 `iterations` in
[`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md)
are recommended (a warning is issued below that). These p-values are
model-conditional and sample-size-sensitive and are reported alongside
the simulated expected ranges, not in place of them.

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

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(6).
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.

## See also

[`RMdimCFACutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md),
[`RMdimCFAPlot`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFAPlot.md)

## Examples

``` r
# \donttest{
if (requireNamespace("lavaan", quietly = TRUE) &&
    requireNamespace("eRm", quietly = TRUE)) {
  data("raschdat1", package = "eRm")
  sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 50,
                        parallel = FALSE, seed = 1)
  tabs <- RMdimCFA(raschdat1[, 1:8], cutoff = sim)
  tabs$fit
  tabs$loadings

  # Bootstrap p-values with family-wise (Westfall-Young) correction
  # (use iterations >= 1000 in real analyses for stable p-values)
  RMdimCFA(raschdat1[, 1:8], cutoff = sim, p_value = TRUE)
}
#> Warning: Bootstrap p-values are based on only 50 simulation iterations. With few iterations the studentised-max (FWER) correction is liberal and small p-values are imprecise; use iterations >= 1000 in RMdimCFACutoff() for reliable p-values.
#> $fit
#> 
#> 
#> Table: Rasch Model posterior-predictive CFA fit-index check. Observed CFA fit (one-factor, lavaan WLSMV, ordered = TRUE) vs simulated null under RM unidimensionality (50 iterations at n = 100 simulees). n = 100 respondents. Cutoffs shown at the 99th percentile for reference. One-sided bootstrap p-values in the unfavourable direction (CFI low; RMSEA/SRMR high); multiplicity correction: Westfall-Young step-down (FWER); flagged at padj < 0.05. p-values cannot be smaller than 1/(50+1) = 0.0196.
#> 
#> |Index | Observed| Cutoff|Direction  |      p| p (adj)|Flagged |
#> |:-----|--------:|------:|:----------|------:|-------:|:-------|
#> |CFI   |   0.9894| 0.8380|< 1st pct  | 0.5098|  0.5098|        |
#> |RMSEA |   0.0162| 0.0916|> 99th pct | 0.5098|  0.5098|        |
#> |SRMR  |   0.1158| 0.1417|> 99th pct | 0.2941|  0.3137|        |
#> 
#> $loadings
#> 
#> 
#> Table: Standardized factor loadings (one-factor, lavaan WLSMV) vs the simulated expected range under RM unidimensionality (50 iterations at n = 100 simulees). n = 100 respondents. Expected range shown for reference (two-sided central 99% interval). p/p (adj): two-sided bootstrap p-values; multiplicity correction: Westfall-Young step-down (FWER); flagged at padj < 0.05 (below / above = direction of the deviation).
#> 
#> |Item | Observed| Expected low| Expected high|      p| p (adj)|Flagged |
#> |:----|--------:|------------:|-------------:|------:|-------:|:-------|
#> |I1   |    0.520|        0.272|         0.938| 0.7255|  0.9804|        |
#> |I2   |    0.441|        0.245|         0.850| 0.4314|  0.8235|        |
#> |I3   |    0.641|        0.249|         0.833| 0.6275|  0.9412|        |
#> |I4   |    0.718|        0.075|         0.873| 0.2745|  0.7843|        |
#> |I5   |    0.303|        0.215|         0.857| 0.0588|  0.2157|        |
#> |I6   |    0.182|        0.209|         0.747| 0.0196|  0.0392|below   |
#> |I7   |    0.419|        0.291|         0.815| 0.2745|  0.7451|        |
#> |I8   |    0.582|        0.250|         0.745| 0.8627|  0.9804|        |
#> 
# }
```
