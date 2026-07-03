# Partial Gamma DIF Analysis

Computes partial gamma coefficients for Differential Item Functioning
(DIF) using
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).
Each item is tested for association with a single categorical DIF
variable, controlling for the total score.

## Usage

``` r
RMdifGamma(
  data,
  dif_var,
  cutoff = NULL,
  p_value = FALSE,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable"
)
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
    [`RMdifGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)
    (a list with `$item_cutoffs`): the data.frame is extracted
    automatically and simulation metadata is included in the kable
    caption.

  - The `$item_cutoffs` data.frame from
    [`RMdifGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)
    directly: must have columns `Item`, `gamma_low`, `gamma_high`. When
    provided, adds columns `Gamma_low`, `Gamma_high`, and `Flagged`
    (logical; `TRUE` when the observed partial gamma falls outside the
    credible range) to the result.

- p_value:

  Logical. When `TRUE`, adds two-sided bootstrap p-values (`p_gamma`,
  `padj_gamma`) comparing each item's observed partial gamma against its
  simulated null distribution, and `flagged` reflects
  `padj_gamma < alpha` instead of the credible range. The asymptotic
  BH-adjusted p-value and star columns from
  [`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html)
  are **dropped** in this mode (two p-value families in one table would
  invite double-reading); the simulated `gamma_low` / `gamma_high` band
  is kept as the effect-size reference. Requires the **full**
  [`RMdifGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)
  object as `cutoff` (it carries the simulated distributions in
  `$results`). Default `FALSE`.

- correction:

  Character. Multiplicity correction for the bootstrap p-values:
  `"fwer"` (default; Westfall-Young studentised-max step-down),
  `"fdr_bh"`, `"fdr_by"`, or `"none"`. Ignored when `p_value = FALSE`.

- alpha:

  Numeric in (0, 1). Significance level used to flag items on the
  corrected p-value. Default `0.05`. Ignored when `p_value = FALSE`.

- output:

  Character string controlling the return value. Either `"kable"`
  (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"` for the underlying data.frame.

## Value

- If `output = "kable"`: a `knitr_kable` object with columns "Item",
  "Partial gamma", "SE", "Lower CI", "Upper CI", "Adj. p-value (BH)",
  and "p-value sign." (a star-string indicator from
  [`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html)).
  When `cutoff` is provided, additional columns "Gamma low", "Gamma
  high", and "Flagged" are included. With `p_value = TRUE`, the
  asymptotic p-value columns are replaced by bootstrap "p" and "p
  (adj)".

- If `output = "dataframe"`: a data.frame with columns `Item`, `gamma`,
  `se`, `lower`, `upper`, `padj_bh`, `Significance`. When `cutoff` is
  provided, columns `gamma_low`, `gamma_high`, and `flagged` are also
  included. With `p_value = TRUE`, `padj_bh` and `Significance` are
  replaced by `p_gamma` and `padj_gamma`.

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

**Bootstrap p-values.** When `p_value = TRUE`, each item's observed
partial gamma is compared against its simulated null distribution (from
`cutoff$results`, where the DIF variable is random by construction). The
per-item statistic is the residual studentised by the bootstrap mean and
SD; the marginal p-value is the two-sided Monte-Carlo p-value
`(1 + #\{|t*| >= |t|\}) / (B + 1)`, so it can be no smaller than
`1 / (B + 1)`. `correction = "fwer"` uses the Westfall-Young
studentised-max step-down, which exploits the bootstrap dependence among
items (Ferreira, 2024); it is liberal when the simulation is small, so
at least 1000 `iterations` in
[`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)
are recommended (a warning is issued below that). Unlike the asymptotic
p-values from
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html),
these are calibrated against the *simulated Rasch null* rather than the
asymptotic SE; they are model-conditional and sample-size-sensitive, and
are reported alongside the simulated effect-size band, not in place of
it.

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

Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
(1998). Differential item functioning in the Danish translation of the
SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189–1202.
[doi:10.1016/S0895-4356(98)00111-5](https://doi.org/10.1016/S0895-4356%2898%2900111-5)

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(6).
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.

## See also

[`RMdifGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)

## Examples

``` r
# \donttest{
if (requireNamespace("iarm", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)
  dif_group <- factor(sample(c("A", "B"), 200, replace = TRUE))

  # Default kable output
  RMdifGamma(sim_data, dif_group)

  # Return as data.frame
  RMdifGamma(sim_data, dif_group, output = "dataframe")

  # Simulation-based cutoffs (100 Monte-Carlo iterations)
  if (requireNamespace("ggdist", quietly = TRUE)) {
    cutoff_res <- RMdifGammaCutoff(sim_data, dif_var = dif_group,
                                   iterations = 100, parallel = FALSE,
                                   seed = 42)
    RMdifGamma(sim_data, dif_group, cutoff = cutoff_res)

    # Bootstrap p-values with family-wise (Westfall-Young) correction
    # (use iterations >= 1000 in real analyses for stable p-values)
    RMdifGamma(sim_data, dif_group, cutoff = cutoff_res, p_value = TRUE,
               output = "dataframe")
  }
}
#> Warning: Bootstrap p-values are based on only 100 simulation iterations. With few iterations the studentised-max (FWER) correction is liberal and small p-values are imprecise; use iterations >= 1000 in RMdifGammaCutoff() for reliable p-values.
#>      Item  gamma    se  lower upper gamma_low gamma_high p_gamma padj_gamma
#> 1   Item1  0.029 0.161 -0.287 0.345    -0.392      0.419  1.0000     1.0000
#> 2   Item2 -0.168 0.156 -0.474 0.138    -0.334      0.297  0.2376     0.9010
#> 3   Item3  0.111 0.168 -0.219 0.441    -0.404      0.486  0.4554     0.9802
#> 4   Item4 -0.030 0.155 -0.333 0.273    -0.338      0.395  0.9010     1.0000
#> 5   Item5  0.076 0.159 -0.235 0.387    -0.368      0.379  0.5743     1.0000
#> 6   Item6  0.042 0.156 -0.264 0.349    -0.347      0.354  0.9010     1.0000
#> 7   Item7 -0.146 0.160 -0.460 0.169    -0.391      0.273  0.4158     0.9802
#> 8   Item8  0.019 0.162 -0.298 0.336    -0.363      0.325  0.9208     1.0000
#> 9   Item9 -0.141 0.161 -0.457 0.175    -0.372      0.419  0.5347     0.9802
#> 10 Item10  0.178 0.156 -0.127 0.483    -0.386      0.354  0.2178     0.9010
#>    flagged
#> 1    FALSE
#> 2    FALSE
#> 3    FALSE
#> 4    FALSE
#> 5    FALSE
#> 6    FALSE
#> 7    FALSE
#> 8    FALSE
#> 9    FALSE
#> 10   FALSE
# }
```
