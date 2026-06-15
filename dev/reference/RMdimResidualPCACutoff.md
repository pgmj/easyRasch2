# Simulation-based Cutoff for First-Contrast Eigenvalue

Parametric bootstrap producing an empirical upper percentile for the
largest eigenvalue from a PCA of standardized residuals under a
correctly fitting Rasch model. For each iteration the function resamples
theta values from the data-based estimates, simulates response data,
refits the appropriate Rasch model, and extracts the first-contrast
eigenvalue. Several upper-tail percentiles of the resulting distribution
are returned; the 99th percentile is reported as the suggested cutoff
(matching the convention used by
[`RMlocdepQ3Cutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)).

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

  A data.frame or matrix of item responses (0-based, non-negative
  integers).

- iterations:

  Integer. Number of simulation iterations. Default `250`.

- parallel:

  Logical. Use parallel processing via `mirai`. Default `TRUE`.

- n_cores:

  Integer or `NULL`. Number of parallel workers. When `NULL`,
  `getOption("mc.cores")` is checked first; if neither is set,
  `parallel = TRUE` falls back to sequential with a warning.

- verbose:

  Logical. Show a progress bar. Default `FALSE`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

## Value

A list with components:

- `results`:

  data.frame: `iteration`, `eigenvalue`.

- `p95`, `p99`, `p995`, `p999`:

  Empirical percentiles of `eigenvalue`.

- `max`:

  The largest simulated eigenvalue.

- `suggested_cutoff`:

  The 99th percentile (`p99`) — pass this list back into
  [`RMdimResidualPCA`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md)
  via `cutoff = `, or use `suggested_cutoff` directly.

- `suggested_cutoff_percentile`:

  The percentile used for `suggested_cutoff`, currently always `99`.

- `actual_iterations`:

  Number of successful iterations.

- `sample_n`:

  Number of complete cases used.

- `item_names`:

  Character vector of item names.

## Details

Rule-of-thumb cutoffs for the first-contrast eigenvalue depend strongly
on sample size, test length, and item-parameter spread; simulation-based
cutoffs tailored to the data are more defensible (Chou & Wang, 2010).

Per iteration: theta values are sampled with replacement from the
WLE/MLE estimates derived from the fitted full-sample model; response
data are simulated under the model
([`psychotools::rrm`](https://rdrr.io/pkg/psychotools/man/rrm.html) for
dichotomous, partial-credit simulator for polytomous); the model is
refitted with [`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) /
[`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html); standardized
residuals are extracted via `eRm::itemfit()$st.res`;
[`prcomp()`](https://rdrr.io/r/stats/prcomp.html) is run; the
first-contrast eigenvalue is recorded.

Iterations that fail (e.g., due to a degenerate simulated dataset where
some category isn't represented) are silently dropped. The `iarm`
package is **not** required.

## References

Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item
response models with principal component analysis on standardized
residuals. *Educational and Psychological Measurement, 70*(5), 717-731.
[doi:10.1177/0013164410379322](https://doi.org/10.1177/0013164410379322)

## See also

[`RMdimResidualPCA`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md)

## Examples

``` r
# \donttest{
set.seed(1)
dat <- as.data.frame(
  matrix(sample(0:1, 200 * 12, replace = TRUE), nrow = 200, ncol = 12)
)
colnames(dat) <- paste0("I", 1:12)

# Few iterations for a fast example; use 250+ in real analyses
bound <- RMdimResidualPCACutoff(dat, iterations = 50, parallel = FALSE, seed = 1)
bound$suggested_cutoff
#> [1] 1.656817

RMdimResidualPCA(dat, cutoff = bound)
#> 
#> 
#> Table: Rasch model (200 complete cases, 12 items). Total observed variance: 9.7% explained by measures, 90.3% unexplained
#> (basis for PCA; n = 200 non-extreme cases). First-contrast cutoff = 1.657 based on 50 simulation iterations (99th percentile).
#> 
#> |Component | Eigenvalue| Proportion of variance|Flagged |
#> |:---------|----------:|----------------------:|:-------|
#> |PC1       |      1.464|                  0.121|FALSE   |
#> |PC2       |      1.425|                  0.118|FALSE   |
#> |PC3       |      1.281|                  0.106|FALSE   |
#> |PC4       |      1.155|                  0.096|FALSE   |
#> |PC5       |      1.125|                  0.093|FALSE   |
# }
```
