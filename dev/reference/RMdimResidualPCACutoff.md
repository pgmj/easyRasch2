# Simulation-based bound for the first-contrast eigenvalue

Simulates datasets from a Rasch or partial credit model that fits and
resembles your data, and records the largest eigenvalue (the first
contrast) of the PCA of standardized residuals in each. Returns several
upper percentiles of these eigenvalues. The 99th is the suggested bound
for
[`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md),
as in
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md).
Fixed thresholds do not work, since the eigenvalue depends on the sample
size, the number of items and the spread of the item parameters (Chou &
Wang, 2010).

## Usage

``` r
RMdimResidualPCACutoff(
  data,
  iterations = 250,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0. Only
  complete cases are used.

- iterations:

  Integer. Number of simulated datasets. Default `250`.

- parallel:

  Logical. Run the simulations in parallel via mirai when it is
  installed. Default `TRUE`.

- n_cores:

  Integer or `NULL`. Number of parallel workers. When `NULL`,
  `getOption("mc.cores")` is used. If that is unset too, a warning is
  given and the simulations run sequentially.

- verbose:

  Logical. Show a progress bar. Default `FALSE`.

- seed:

  Integer or `NULL`. Random seed for reproducibility. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md)
  for what this guarantees and how it interacts with `parallel`.

## Value

A list with:

- `results`:

  data.frame with `iteration` and `eigenvalue`.

- `p95`, `p99`, `p995`, `p999`:

  Percentiles of the simulated eigenvalues.

- `max`:

  The largest simulated eigenvalue.

- `suggested_cutoff`:

  The 99th percentile (`p99`). Pass the whole list to
  [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md)
  as `cutoff`, or use this value directly.

- `suggested_cutoff_percentile`:

  The percentile behind `suggested_cutoff`, always `99`.

- `actual_iterations`:

  Number of simulated datasets that could be used.

- `requested_iterations`:

  The `iterations` requested.

- `sample_n`:

  Number of complete cases used.

- `item_names`:

  Item names.

## Details

The generating model has the CML item parameters (`psychotools`) and the
WLE person locations of `data`. In each iteration, person locations are
drawn with replacement from the WLE estimates, responses are simulated
from the model, the standardized residuals are computed with the same
engine as for the observed data, and the largest eigenvalue of their PCA
is recorded.

A simulated dataset that cannot be refitted, for example because a
response category is empty, is left out, and
[`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md)
reports how many in its caption. `actual_iterations` gives the number
that remained. The iarm package is not needed.

## References

Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item
response models with principal component analysis on standardized
residuals. *Educational and Psychological Measurement, 70*(5), 717-731.
[doi:10.1177/0013164410379322](https://doi.org/10.1177/0013164410379322)

## See also

[`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md)

## Examples

``` r
# \donttest{
set.seed(1)
dat <- as.data.frame(
  matrix(sample(0:1, 200 * 12, replace = TRUE), nrow = 200, ncol = 12)
)
colnames(dat) <- paste0("I", 1:12)

# Few iterations for a fast example. Use 250 or more in real analyses.
bound <- RMdimResidualPCACutoff(dat, iterations = 50, parallel = FALSE, seed = 1)
bound$suggested_cutoff
#> [1] 1.623358

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
# }
```
