# Partial Gamma DIF Analysis

Computes partial gamma coefficients for Differential Item Functioning
(DIF) using
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).
Each item is tested for association with a single categorical DIF
variable, controlling for the total score.

## Usage

``` r
RMpartgamDIF(data, dif_var, cutoff = NULL, output = "kable")
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed, but at least one complete case must exist after combining
  `data` and `dif_var`.

- dif_var:

  A vector (factor or character) of the same length as `nrow(data)`,
  representing the grouping variable for DIF analysis.

- cutoff:

  Optional. Default `NULL` (no cutoff applied). Can be:

  - The return value of
    [`RMpgDIFcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgDIFcutoff.md)
    (a list with `$item_cutoffs`): the data.frame is extracted
    automatically and simulation metadata is included in the kable
    caption.

  - The `$item_cutoffs` data.frame from
    [`RMpgDIFcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgDIFcutoff.md)
    directly: must have columns `Item`, `gamma_low`, `gamma_high`. When
    provided, adds columns `Gamma_low`, `Gamma_high`, and `Flagged`
    (logical; `TRUE` when the observed partial gamma falls outside the
    credible range) to the result.

- output:

  Character string controlling the return value. Either `"kable"`
  (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"` for the underlying data.frame.

## Value

- If `output = "kable"`: a `knitr_kable` object with columns "Item",
  "Partial gamma", "SE", "Lower CI", "Upper CI", and "Adjusted p-value
  (BH)". When `cutoff` is provided, additional columns "Gamma low",
  "Gamma high", and "Flagged" are included.

- If `output = "dataframe"`: a data.frame with columns `Item`, `gamma`,
  `se`, `lower`, `upper`, `padj_bh`. When `cutoff` is provided, columns
  `gamma_low`, `gamma_high`, and `flagged` are also included.

## Details

Partial gamma (Bjorner et al., 1998) measures the association between
item response and an exogenous grouping variable, controlling for the
total score. Values near 0 indicate no DIF. Recommended interpretive
thresholds (Bjorner et al., 1998):

- **No or negligible DIF**: gamma within \\\[-0.21, 0.21\]\\, *or* gamma
  not significantly different from 0.

- **Slight to moderate DIF**: gamma within \\\[-0.31, 0.31\]\\ (and
  outside \\\[-0.21, 0.21\]\\), *or* not significantly outside
  \\\[-0.21, 0.21\]\\.

- **Moderate to large DIF**: gamma outside \\\[-0.31, 0.31\]\\, **and**
  significantly outside \\\[-0.21, 0.21\]\\.

The `iarm` package must be installed (it is in Suggests, not Imports).

## References

Bjorner, J., Kreiner, S., Ware, J., Damsgaard, M. and Bech, P.
Differential item functioning in the Danish translation of the SF-36.
*Journal of Clinical Epidemiology*, 51(11), 1998, pp. 1189-1202.

## See also

[`RMpgDIFcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgDIFcutoff.md)

## Examples

``` r
# \donttest{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)
dif_group <- factor(sample(c("A", "B"), 200, replace = TRUE))

# Default kable output
RMpartgamDIF(sim_data, dif_group)
#> 
#> 
#> Table: Partial gamma DIF analysis (n = 200 complete cases). Positive gamma indicates higher scores in higher DIF group levels.
#> 
#> |Item   | Partial gamma|    SE| Lower CI| Upper CI| Adjusted p-value (BH)|
#> |:------|-------------:|-----:|--------:|--------:|---------------------:|
#> |Item1  |         0.029| 0.161|   -0.287|    0.345|                     1|
#> |Item2  |        -0.168| 0.156|   -0.474|    0.138|                     1|
#> |Item3  |         0.111| 0.168|   -0.219|    0.441|                     1|
#> |Item4  |        -0.030| 0.155|   -0.333|    0.273|                     1|
#> |Item5  |         0.076| 0.159|   -0.235|    0.387|                     1|
#> |Item6  |         0.042| 0.156|   -0.264|    0.349|                     1|
#> |Item7  |        -0.146| 0.160|   -0.460|    0.169|                     1|
#> |Item8  |         0.019| 0.162|   -0.298|    0.336|                     1|
#> |Item9  |        -0.141| 0.161|   -0.457|    0.175|                     1|
#> |Item10 |         0.178| 0.156|   -0.127|    0.483|                     1|

# Return as data.frame
RMpartgamDIF(sim_data, dif_group, output = "dataframe")
#>      Item  gamma    se  lower upper padj_bh
#> 1   Item1  0.029 0.161 -0.287 0.345       1
#> 2   Item2 -0.168 0.156 -0.474 0.138       1
#> 3   Item3  0.111 0.168 -0.219 0.441       1
#> 4   Item4 -0.030 0.155 -0.333 0.273       1
#> 5   Item5  0.076 0.159 -0.235 0.387       1
#> 6   Item6  0.042 0.156 -0.264 0.349       1
#> 7   Item7 -0.146 0.160 -0.460 0.169       1
#> 8   Item8  0.019 0.162 -0.298 0.336       1
#> 9   Item9 -0.141 0.161 -0.457 0.175       1
#> 10 Item10  0.178 0.156 -0.127 0.483       1
# }
if (FALSE) { # \dontrun{
# Simulation-based cutoffs (slow): 100+ Monte-Carlo iterations
cutoff_res <- RMpgDIFcutoff(sim_data, dif_var = dif_group,
                                  iterations = 100, parallel = FALSE,
                                  seed = 42)
RMpartgamDIF(sim_data, dif_group, cutoff = cutoff_res)
} # }
```
