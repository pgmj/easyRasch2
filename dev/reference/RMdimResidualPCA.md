# PCA of standardized Rasch residuals

Fits a Rasch or partial credit model by CML (`psychotools`), computes
the standardized residuals \\(x - E)/\sqrt{Var}\\ at the WLE person
locations, and runs an unrotated principal component analysis on them
with [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html). Reports
the largest eigenvalues and their shares of the unexplained variance.
With the result of
[`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCACutoff.md),
it compares the largest eigenvalue (the first contrast) with a simulated
bound.

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

  A data.frame or matrix of item responses, items scored from 0. Rows
  with any `NA` are dropped before the PCA, since
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) does not accept
  missing values.

- cutoff:

  Optional. The result of
  [`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCACutoff.md)
  (its `suggested_cutoff`, the 99th percentile, is used), or a single
  number to use as the bound. Adds a logical `Flagged` column (is the
  eigenvalue above the bound?), and the caption reports the bound.

- p_value:

  Logical. Add a one-sided bootstrap p-value for the first-contrast
  eigenvalue, the share of simulated largest eigenvalues at least as
  large as the observed one, `(1 + #\{lambda* >= lambda\}) / (B + 1)`
  for B simulated datasets. Needs the full
  [`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCACutoff.md)
  result as `cutoff`. The simulations are of the largest eigenvalue, so
  the p-value is given for the first component only (`NA` for the
  others). It is a single test, so no correction for multiple testing is
  needed. Default `FALSE`.

- n_components:

  Integer. Number of eigenvalues to report, at most the number of items.
  Default `5`.

- output:

  Character. `"kable"` (default) for a
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table,
  `"dataframe"`, `"list"` for the table and the variance partition, or
  `"ggplot"` for a plot of the first-component loadings against the item
  locations. `"loadings"` is an old name for `"ggplot"`, still accepted.

## Value

- `output = "kable"`: a `knitr_kable` table with the columns Component,
  Eigenvalue and Proportion of variance, plus `Flagged` with a `cutoff`
  and `p` with `p_value = TRUE`. The caption gives the variance
  partition (the percentages of the observed variance explained by the
  model and unexplained), the model, the sample size, and the bound or
  p-value.

- `output = "dataframe"`: a data.frame with the columns `Component`,
  `Eigenvalue` and `Proportion_of_variance`, plus `Flagged` and `p` as
  above. The variance partition is also attached as
  `attr(result, "variance_partition")`, kept for compatibility. Use
  `output = "list"` to get it as a table.

- `output = "list"`: a list with `components`, the data.frame above, and
  `variance`, a one-row data.frame with `total`, `explained`,
  `unexplained`, `pct_explained`, `pct_unexplained` and `n_persons`, and
  `kable`, the table from `output = "kable"`.

- `output = "ggplot"`: each item's first-component loading (x-axis)
  against its location (y-axis), with dashed lines at zero and the
  variance partition in the caption. Item names are placed with
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
  when ggrepel is installed.

## Details

Fixed thresholds for the first-contrast eigenvalue, such as "above 2",
do not indicate multidimensionality reliably. Under a unidimensional
model that fits, the eigenvalue varies with the sample size, the number
of items and the spread of the item parameters, so the bound should be
simulated for the data at hand (Chou & Wang, 2010).

The residuals come from the package's CML and WLE engine. The
eigenvalues are unrotated. Rotation helps interpret a multidimensional
solution, but it hides the dominant first contrast that the assessment
of dimensionality is about. The item locations in the plot are the means
of each item's CML thresholds.

The variance partition follows Linacre's convention: each item's
observed variance is compared with its expected variance under the
model, summed over items. The expected scores use the CML item
parameters and the WLE person locations. WLE estimates are finite at
extreme scores, so all respondents are kept.

**Bootstrap p-value.** With `p_value = TRUE`, the observed
first-contrast eigenvalue is compared with the simulated largest
eigenvalues in `cutoff$results`. The largest eigenvalue is already a
single statistic over all components, so no correction for multiple
testing applies. The p-value depends on the model and the sample size,
is reported with the simulated bound, not instead of it, and can be no
smaller than `1 / (B + 1)`.

## References

Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item
response models with principal component analysis on standardized
residuals. *Educational and Psychological Measurement, 70*(5), 717-731.
[doi:10.1177/0013164410379322](https://doi.org/10.1177/0013164410379322)

## See also

[`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCACutoff.md)
for the simulated bound,
[`RMdimCFA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFA.md)
and
[`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md)
for other checks of unidimensionality.

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
#> Table: Rasch model (12 items), n = 200 respondents. Total observed variance: 8.3% explained by measures, 91.7% unexplained. First-contrast cutoff = 1.623 based on 50 simulation iterations (99th percentile). One-sided bootstrap p-value for the first contrast (single test, no multiplicity correction); it cannot be smaller than 1/(50+1) = 0.02.
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
