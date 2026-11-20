# One-factor CFA fit and loadings against a simulated reference

Fits a one-factor CFA for ordinal items to `data` and compares its fit
indices and the standardized loading of each item with the simulations
from
[`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md).
Returns two tables, fit indices and loadings, each with the observed
value, the simulated reference and a flag.

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

  A data.frame or matrix of item responses, items scored from 0, the
  same items as in the simulation. Only complete cases are used.

- cutoff:

  The result of
  [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md).
  Required, since the observed fit indices cannot be judged without it.

- p_value:

  Logical. Add bootstrap p-values from the simulations and flag on them
  instead of on the percentile cutoffs. The p-values are one-sided in
  the unfavourable direction for the fit indices (CFI low, RMSEA and
  SRMR high) and two-sided for the loadings. The fit indices and the
  loadings are corrected as two separate families. Default `FALSE`.

- correction:

  Character. Correction of the p-values: `"fwer"` (default,
  Westfall-Young step-down), `"fdr_bh"` (Benjamini-Hochberg), `"fdr_by"`
  (Benjamini-Yekutieli) or `"none"`. Used only with `p_value = TRUE`.

- alpha:

  Numeric in (0, 1). A comparison is flagged when its corrected p-value
  is below `alpha`. Default `0.05`. Used only with `p_value = TRUE`.

- output:

  Character. `"kable"` (default) returns each table as a
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html),
  `"dataframe"` as a data.frame.

## Value

A named list with `fit` and `loadings`, each a `knitr_kable` or a
data.frame by `output`:

- `fit`:

  CFI, RMSEA and SRMR with the columns `Index`, `Observed`, `Cutoff`,
  `Direction` and `Flagged` (one-sided, in the unfavourable direction).
  With `p_value = TRUE`, also `p` and `padj`, and `Flagged` marks
  `padj < alpha`. An index lavaan could not compute for `data` has
  `Observed = NA` and `Flagged = "not tested"`.

- `loadings`:

  One row per item with the columns `Item`, `Observed`, `Expected_low`,
  `Expected_high` and `Flagged` (`"below"`, `"above"` or `""`). With
  `p_value = TRUE`, also `p_loading` and `padj_loading`, and `Flagged`
  marks `padj_loading < alpha`, with the direction from the deviation
  from the simulated mean.

## Details

**Bootstrap p-values.** The statistic for each comparison is the
deviation from the mean of the observed and simulated values together,
divided by their SD. The uncorrected p is `(1 + count) / (B + 1)` for B
simulated datasets, so it can be no smaller than `1 / (B + 1)`.
`correction = "fwer"` uses the Westfall-Young step-down on the
studentised maximum, which uses the dependence between the statistics
(Westfall & Young, 1993; Ferreira, 2024). It controls the family-wise
error rate at any number of iterations, but small p-values are imprecise
and decisions depend on the seed when there are few simulations. Use at
least 1000 `iterations` in
[`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md).
A warning says so below that. The p-values depend on the model and the
sample size, and are reported with the simulated ranges, not instead of
them.

## Multiple comparisons

The uncorrected p-value controls the error rate of a single comparison.
For one item (or item pair) chosen in advance it is the relevant value.
Flagging whichever of k comparisons fall below `alpha` tests k
hypotheses at once, and the chance of at least one false flag rises to
about \\1 - (1 - \alpha)^k\\, for example 34% for k = 8 at
`alpha = 0.05`, even when every uncorrected p-value is calibrated. The
corrected p-value controls this. `correction = "fwer"` bounds the
probability of any false flag (strict, less power), and `"fdr_bh"` and
`"fdr_by"` bound the expected share of false flags among those raised.
Use the uncorrected p-value for a single comparison chosen in advance,
and a corrected one when screening the whole table, which is the usual
case.

## References

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(1), 6.
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.

## See also

[`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md)
for the simulations,
[`RMdimCFAPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFAPlot.md)
for the figures.

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
#> Warning: Bootstrap p-values are based on only 50 simulation iterations. With few iterations small p-values are imprecise and decisions depend on the seed; use iterations >= 1000 in RMdimCFACutoff() for reliable p-values.
#> $fit
#> 
#> 
#> Table: Rasch Model posterior-predictive CFA fit-index check. Observed CFA fit (one-factor, lavaan WLSMV, ordered = TRUE) vs simulated null under RM unidimensionality (50 iterations at n = 100 simulees). n = 100 respondents. Cutoffs shown at the 99th percentile for reference. One-sided bootstrap p-values in the unfavourable direction (CFI low; RMSEA/SRMR high); multiplicity correction: Westfall-Young step-down (FWER); flagged at padj < 0.05. p-values cannot be smaller than 1/(50+1) = 0.02.
#> 
#> |Index | Observed| Cutoff|Direction  |      p| p (adj)|Flagged |
#> |:-----|--------:|------:|:----------|------:|-------:|:-------|
#> |CFI   |   0.9894| 0.8380|< 1st pct  | 0.5098|  0.5294|        |
#> |RMSEA |   0.0162| 0.0916|> 99th pct | 0.5098|  0.5294|        |
#> |SRMR  |   0.1158| 0.1417|> 99th pct | 0.2941|  0.3333|        |
#> 
#> $loadings
#> 
#> 
#> Table: Standardized factor loadings (one-factor, lavaan WLSMV) vs the simulated expected range under RM unidimensionality (50 iterations at n = 100 simulees). n = 100 respondents. Expected range shown for reference (two-sided central 99% interval). p/p (adj): two-sided bootstrap p-values; multiplicity correction: Westfall-Young step-down (FWER); flagged at padj < 0.05 (below / above = direction of the deviation).
#> 
#> |Item | Observed| Expected low| Expected high|      p| p (adj)|Flagged |
#> |:----|--------:|------------:|-------------:|------:|-------:|:-------|
#> |I1   |    0.520|        0.272|         0.938| 0.7255|  0.9804|        |
#> |I2   |    0.441|        0.245|         0.850| 0.4510|  0.8627|        |
#> |I3   |    0.641|        0.249|         0.833| 0.6471|  0.9608|        |
#> |I4   |    0.718|        0.075|         0.873| 0.2941|  0.7843|        |
#> |I5   |    0.303|        0.215|         0.857| 0.0588|  0.3137|        |
#> |I6   |    0.182|        0.209|         0.747| 0.0196|  0.0588|        |
#> |I7   |    0.419|        0.291|         0.815| 0.2745|  0.7647|        |
#> |I8   |    0.582|        0.250|         0.745| 0.8627|  0.9804|        |
#> 
# }
```
