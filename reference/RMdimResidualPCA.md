# PCA of Standardized Rasch Residuals

Fits a Rasch model by CML via `psychotools` (a dichotomous item is a
2-category PCM), extracts the standardized residuals \\(x -
E)/\sqrt{Var}\\ at WLE person locations, and runs an unrotated
principal-component analysis on those residuals via
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html). The function
reports the top `n_components` eigenvalues and their proportions of
unexplained variance, and optionally compares the first-contrast
eigenvalue against a simulation-based bound from
[`RMdimResidualPCACutoff`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCACutoff.md).

## Usage

``` r
RMdimResidualPCA(
  data,
  cutoff = NULL,
  p_value = FALSE,
  n_components = 5L,
  output = "kable"
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Rows with any `NA` are dropped
  before PCA, since [`prcomp()`](https://rdrr.io/r/stats/prcomp.html)
  does not accept missing values.

- cutoff:

  Optional. The list returned by
  [`RMdimResidualPCACutoff`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCACutoff.md)
  (its `suggested_cutoff` is used), or a single numeric value to use as
  the cutoff directly. When provided, the result includes a `Flagged`
  column (logical: is the eigenvalue above the simulated bound?) and the
  kable caption notes the cutoff.

- p_value:

  Logical. When `TRUE`, adds a one-sided bootstrap p-value for the
  first-contrast eigenvalue: the proportion of simulated first-contrast
  eigenvalues at least as large as the observed one,
  `(1 + #\{lambda* >= lambda\}) / (B + 1)`. Requires the **full**
  [`RMdimResidualPCACutoff`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCACutoff.md)
  object as `cutoff` (it carries the simulated eigenvalues in
  `$results`); a bare numeric cutoff is not sufficient. The simulated
  null is the distribution of the *largest* eigenvalue, so the p-value
  applies to PC1 only (`NA` for the other components). This is a single
  test — no multiplicity correction is involved. Default `FALSE`.

- n_components:

  Integer. Number of eigenvalues to report. Capped at the number of
  items. Default `5`.

- output:

  Character. `"kable"` (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table,
  `"dataframe"` for the underlying data.frame, or `"ggplot"` for a
  ggplot of PC1 loadings against item locations. `"loadings"` is
  accepted as a backward-compatible alias for `"ggplot"`.

## Value

- If `output = "kable"`: a `knitr_kable` object with columns Component,
  Eigenvalue, Proportion of variance (and `Flagged` when `cutoff` is
  provided; `p` when `p_value = TRUE`). The caption gives the variance
  partition (% of total observed variance explained by measures vs.
  unexplained), the model fitted, sample size, and cutoff/p-value
  metadata if applicable.

- If `output = "dataframe"`: a data.frame with columns `Component`,
  `Eigenvalue`, `Proportion_of_variance` (and `Flagged` when `cutoff` is
  provided; `p` when `p_value = TRUE`, non-`NA` for PC1 only). The
  variance partition is attached as the `"variance_partition"` attribute
  — a list with elements `total`, `explained`, `unexplained`,
  `pct_explained`, `pct_unexplained`, `n_persons`. Access via
  `attr(result, "variance_partition")`.

- If `output = "ggplot"`: a ggplot showing each item's PC1 loading on
  the x-axis and Rasch item location on the y-axis, with dashed
  reference lines at zero, and the variance partition in the figure
  caption. Item names are labelled via
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
  when `ggrepel` is installed; otherwise plain `geom_text()`.

## Details

Rule-of-thumb thresholds for the first-contrast eigenvalue (e.g., the
"\> 2" heuristic occasionally cited from Winsteps documentation) are not
reliable indicators of multidimensionality; the first-contrast
eigenvalue under a correctly fitting unidimensional model varies
systematically with sample size, test length, and item-parameter spread.
Empirical (simulated) bounds tailored to the data structure should be
used instead — see
[`RMdimResidualPCACutoff`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCACutoff.md),
and Chou & Wang (2010) for the underlying simulation argument.

The PCA is performed on the standardized residuals \\(x -
E)/\sqrt{Var}\\ from the shared CML/WLE engine (CML item parameters via
`psychotools`, WLE person locations). The reported eigenvalues are
unrotated; rotation is appropriate for *interpreting* a multidimensional
solution but obscures the dominant first contrast that dimensionality
assessment is concerned with.

Item locations on the loadings plot are the per-item mean of the CML
Andrich thresholds.

The variance partition follows Linacre's convention: per-item observed
variance is compared to per-item *expected* variance under the fitted
model, summed across items. Expected scores are computed from the CML
item parameters and WLE person locations. WLE is finite at extreme
scores, so all persons are retained (the previous MLE partition dropped
extreme-score cases).

**Bootstrap p-value.** When `p_value = TRUE`, the observed
first-contrast eigenvalue is compared against the simulated null
distribution of largest eigenvalues (from `cutoff$results`), giving the
one-sided Monte-Carlo p-value `(1 + #\{lambda* >= lambda\}) / (B + 1)`.
Because the maximum eigenvalue is a single family-wise statistic, no
multiplicity correction applies. The p-value is model-conditional and
sample-size-sensitive; it is reported alongside the simulated cutoff,
not in place of it, and can be no smaller than `1 / (B + 1)`.

## References

Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item
response models with principal component analysis on standardized
residuals. *Educational and Psychological Measurement, 70*(5), 717-731.
[doi:10.1177/0013164410379322](https://doi.org/10.1177/0013164410379322)

## See also

[`RMdimResidualPCACutoff`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCACutoff.md)

## Examples

``` r
# \donttest{
set.seed(1)
dat <- as.data.frame(
  matrix(sample(0:1, 200 * 12, replace = TRUE), nrow = 200, ncol = 12)
)
colnames(dat) <- paste0("I", 1:12)

# Default kable output
RMdimResidualPCA(dat)
#> 
#> 
#> Table: Rasch model (12 items), n = 200 respondents. Total observed variance: 8.3% explained by measures, 91.7% unexplained.
#> 
#> |Component | Eigenvalue| Proportion of variance|
#> |:---------|----------:|----------------------:|
#> |PC1       |      1.431|                  0.121|
#> |PC2       |      1.407|                  0.119|
#> |PC3       |      1.254|                  0.106|
#> |PC4       |      1.132|                  0.096|
#> |PC5       |      1.105|                  0.093|

# PC1 loadings vs item location plot
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("ggrepel", quietly = TRUE)) {
  RMdimResidualPCA(dat, output = "ggplot")
}


# Simulation-based cutoff (use 250+ iterations in real analyses)
bound <- RMdimResidualPCACutoff(dat, iterations = 50, parallel = FALSE, seed = 1)
RMdimResidualPCA(dat, cutoff = bound)
#> 
#> 
#> Table: Rasch model (12 items), n = 200 respondents. Total observed variance: 8.3% explained by measures, 91.7% unexplained. First-contrast cutoff = 1.623 based on 50 simulation iterations (99th percentile).
#> 
#> |Component | Eigenvalue| Proportion of variance|Flagged |
#> |:---------|----------:|----------------------:|:-------|
#> |PC1       |      1.431|                  0.121|FALSE   |
#> |PC2       |      1.407|                  0.119|FALSE   |
#> |PC3       |      1.254|                  0.106|FALSE   |
#> |PC4       |      1.132|                  0.096|FALSE   |
#> |PC5       |      1.105|                  0.093|FALSE   |

# With the one-sided bootstrap p-value for the first contrast
RMdimResidualPCA(dat, cutoff = bound, p_value = TRUE)
#> 
#> 
#> Table: Rasch model (12 items), n = 200 respondents. Total observed variance: 8.3% explained by measures, 91.7% unexplained. First-contrast cutoff = 1.623 based on 50 simulation iterations (99th percentile). One-sided bootstrap p-value for the first contrast (single test, no multiplicity correction); it cannot be smaller than 1/(50+1) = 0.0196.
#> 
#> |Component | Eigenvalue| Proportion of variance|Flagged |      p|
#> |:---------|----------:|----------------------:|:-------|------:|
#> |PC1       |      1.431|                  0.121|FALSE   | 0.7451|
#> |PC2       |      1.407|                  0.119|FALSE   |     NA|
#> |PC3       |      1.254|                  0.106|FALSE   |     NA|
#> |PC4       |      1.132|                  0.096|FALSE   |     NA|
#> |PC5       |      1.105|                  0.093|FALSE   |     NA|
# }
```
