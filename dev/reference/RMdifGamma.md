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
#>      Item       gamma        se      lower     upper  gamma_low gamma_high
#> 1   Item1  0.02929936 0.1612158 -0.2866779 0.3452766 -0.3918869  0.4190476
#> 2   Item2 -0.16753927 0.1561349 -0.4735580 0.1384794 -0.3341646  0.2969502
#> 3   Item3  0.11111111 0.1681747 -0.2185053 0.4407275 -0.4043210  0.4864595
#> 4   Item4 -0.03045685 0.1545831 -0.3334342 0.2725205 -0.3377265  0.3954457
#> 5   Item5  0.07616708 0.1587018 -0.2348828 0.3872169 -0.3684211  0.3787375
#> 6   Item6  0.04239401 0.1562001 -0.2637526 0.3485406 -0.3467337  0.3541153
#> 7   Item7 -0.14560440 0.1603616 -0.4599073 0.1686985 -0.3906899  0.2727273
#> 8   Item8  0.01897019 0.1618599 -0.2982693 0.3362097 -0.3626374  0.3248998
#> 9   Item9 -0.14088398 0.1614016 -0.4572253 0.1754573 -0.3724247  0.4188101
#> 10 Item10  0.17777778 0.1556159 -0.1272237 0.4827792 -0.3863216  0.3540313
#>      p_gamma padj_gamma flagged
#> 1  1.0000000  1.0000000   FALSE
#> 2  0.2376238  0.9009901   FALSE
#> 3  0.4554455  0.9801980   FALSE
#> 4  0.9009901  1.0000000   FALSE
#> 5  0.5742574  1.0000000   FALSE
#> 6  0.9009901  1.0000000   FALSE
#> 7  0.4158416  0.9801980   FALSE
#> 8  0.9207921  1.0000000   FALSE
#> 9  0.5346535  0.9801980   FALSE
#> 10 0.2178218  0.9009901   FALSE
# }
```
