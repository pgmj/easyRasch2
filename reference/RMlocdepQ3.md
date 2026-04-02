# Q3 Residual Correlations for Local Dependence Assessment

Computes Yen's Q3 residual correlations between item pairs using a Rasch
model fitted via Marginal Maximum Likelihood (MML) in `mirt`. High
correlations (above the dynamic cut-off) indicate potential local
dependence between items. See
[`RMlocdepQ3cutoff`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3cutoff.md)
for how to determine the appropriate dynamic cut-off for your data.

## Usage

``` r
RMlocdepQ3(data, cutoff = NULL, output = "kable")
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed.

- cutoff:

  Optional single numeric value added to the mean off-diagonal Q3
  correlation to produce the dynamic cut-off threshold. When `NULL`
  (default), the raw Q3 residual correlation matrix is returned without
  any dynamic cut-off applied. Use
  [`RMlocdepQ3cutoff`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3cutoff.md)
  to derive a simulation-based cutoff value.

- output:

  Character string controlling the return value. Either `"kable"`
  (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"` for the underlying numeric data.frame.

## Value

- If `output = "kable"`: a `knitr_kable` object showing the lower
  triangle of the Q3 correlation matrix. When `cutoff` is provided, a
  footnote describing the dynamic cut-off is included and an extra
  `above_cutoff` column marks rows containing at least one value above
  the threshold.

- If `output = "dataframe"`: a data.frame (rounded to 2 decimal places)
  with the lower triangle of the Q3 correlation matrix; the upper
  triangle and diagonal are set to `NA`. When `cutoff` is provided, an
  additional logical column `above_cutoff` indicates whether the row
  contains any value exceeding the dynamic cut-off.

## Details

The Q3 statistic (Yen, 1984) is the correlation between residuals of
pairs of items after accounting for the latent trait. Under local
independence, Q3 values are expected to be around \\-1/(k-1)\\ where
\\k\\ is the number of items. When `cutoff` is supplied, the dynamic
cut-off is the mean of all off-diagonal Q3 values plus `cutoff`,
following the approach of Christensen et al. (2017). Use
[`RMlocdepQ3cutoff`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3cutoff.md)
to obtain a simulation-based cutoff recommendation.

`mirt` is used for model fitting here because Q3 requires model-based
expected responses, which are most readily available from MML
estimation.

## References

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement*, 8(2), 125–145.

Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
for Yen's Q3: Identification of local dependence in the Rasch model.
*Applied Psychological Measurement*, 41(3), 178–194.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate binary item response data (10 items, 200 persons)
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Raw Q3 matrix (no cutoff)
RMlocdepQ3(sim_data)

# Get the underlying data.frame
q3_df <- RMlocdepQ3(sim_data, output = "dataframe")

# With a simulation-based cutoff
cutoff_res <- RMlocdepQ3cutoff(sim_data, iterations = 100)
RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
} # }
```
