# Partial gamma DIF

Computes the partial gamma coefficient of each item, the association
between the item and a grouping variable given the total score, to
assess differential item functioning (DIF). With the result of
[`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md),
it tests each item against simulations without DIF, by default with a
bootstrap p-value corrected across items.

## Usage

``` r
RMdifGamma(
  data,
  dif_var,
  cutoff = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable"
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0. Only
  complete cases (in `data` and `dif_var`) are used, so at least one is
  needed.

- dif_var:

  A vector (factor or character) of length `nrow(data)` with the groups,
  two or more.

- cutoff:

  Optional. Default `NULL`, the asymptotic test of
  [`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).
  Either of:

  - The full result of
    [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md).
    This enables the bootstrap p-values (see `p_value`), and the caption
    reports the simulation settings.

  - Its `$item_cutoffs` data.frame alone, with columns `Item`,
    `gamma_low` and `gamma_high`. Items are flagged against the
    interval.

  Either adds the columns `Gamma_low`, `Gamma_high` and a logical
  `Flagged` (`TRUE` when the partial gamma is outside the simulated
  interval, or with p-values when the item is flagged).

- p_value:

  Logical or `NULL`. Whether to flag on two-sided bootstrap p-values
  (`p_gamma`, `padj_gamma`) from the simulations. The asymptotic p-value
  and star columns are then dropped, so that the table has one family of
  p-values, and the simulated interval is kept as a reference. `TRUE`
  needs the full
  [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)
  result. `NULL` (default) uses the p-values whenever `cutoff` is that
  result, so a call without a cutoff keeps the asymptotic table. `FALSE`
  flags against the interval, a poor decision rule, since over k items a
  width `w` gives a family-wise error rate of `1 - w^k` (Johansson,
  2026).

- correction:

  Character. Correction of the bootstrap p-values across items: `"fwer"`
  (default, Westfall-Young step-down), `"fdr_bh"`, `"fdr_by"` or
  `"none"`. Used only with p-values.

- alpha:

  Numeric in (0, 1). An item is flagged when its corrected p-value is
  below `alpha`. Default `0.05`. Used only with p-values.

- output:

  Character. `"kable"` (default) for a
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"`.

## Value

- `output = "kable"`: a `knitr_kable` table with the columns "Item",
  "Partial gamma", "SE", "Lower CI", "Upper CI", "Adj. p-value (BH)" and
  "p-value sign." (stars for the adjusted p-value, as iarm shows them),
  plus "Gamma low", "Gamma high" and "Flagged" with a `cutoff`. With
  p-values, the asymptotic columns are replaced by "p" and "p (adj)".

- `output = "dataframe"`: a data.frame with the columns `Item`, `gamma`,
  `se`, `lower`, `upper`, `padj_bh` and `Significance`, plus
  `gamma_low`, `gamma_high` and `flagged` with a `cutoff`. With
  p-values, `padj_bh` and `Significance` are replaced by `p_gamma` and
  `padj_gamma`.

## Details

**Statistic.** The partial gamma (Bjorner et al., 1998) measures the
association between the item response and an external grouping variable
given the total score. Values near 0 indicate no DIF. Bjorner et al.
(1998) also proposed classes of DIF size, which
[`RMdifTree()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md)
uses. They are conventions, not values calibrated to the sample and
items at hand, and this function does not use them. The iarm package
must be installed (it is in Suggests).

**Asymptotic p-values.** The adjusted p-value is a Benjamini-Hochberg
adjustment, made by easyRasch2, of the p-values from
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html)
over the items. iarm's own adjusted column is a Bonferroni correction
whatever method is named, since iarm adjusts one p-value at a time, so
it is not used.

**Bootstrap p-values.** Each item's observed partial gamma is compared
with its simulated values (`cutoff$results`, simulated without DIF). The
statistic is the deviation from the mean of the observed and simulated
values together, divided by their SD. The uncorrected p is two-sided,
`(1 + #\{|t*| >= |t|\}) / (B + 1)` for B simulated datasets, so it can
be no smaller than `1 / (B + 1)`. It is two-sided since DIF can favour
either group, unlike
[`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md),
which tests one-sided for excess positive local dependence.
`correction = "fwer"` uses the Westfall-Young step-down, which uses the
dependence between the items (Westfall & Young, 1993, Ferreira, 2024).
It controls the family-wise error rate at any number of iterations in
[`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md),
but below 400 decisions depend on the seed and a message says so. Unlike
the asymptotic p-values, these are calibrated against simulations from
the Rasch model. They depend on the model and the sample size, and are
reported with the simulated interval, not instead of it. The observed
partial gamma they test is computed by the same internal routine as the
simulations, which reproduces
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html)
exactly.

**Real and artificial DIF.** Partial gamma conditions on the total
score, which includes every item. Real DIF in one item then induces
artificial DIF in the other items, in the opposite direction, since at a
given total score the item scores must add up to that score (Hagquist &
Andrich, 2017). An item can therefore be flagged only because another
item has DIF. The simulations in
[`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)
do not guard against this, since they have no DIF in any item. Resolve
DIF items one at a time, starting with the largest, and repeat the
analysis before deciding on the next (Hagquist & Andrich, 2017).
[`RMdifAnchor()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifAnchor.md)
avoids the problem by linking the groups through anchor items.

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

Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
(1998). Differential item functioning in the Danish translation of the
SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189–1202.
[doi:10.1016/S0895-4356(98)00111-5](https://doi.org/10.1016/S0895-4356%2898%2900111-5)

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(1), 6.
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

Hagquist, C., & Andrich, D. (2017). Recent advances in analysis of
differential item functioning in health research using the Rasch model.
*Health and Quality of Life Outcomes, 15*, 181.
[doi:10.1186/s12955-017-0755-0](https://doi.org/10.1186/s12955-017-0755-0)

Johansson, M. (2026). Simulation-based cutoffs for conditional item fit
in Rasch models: Iterations, multiplicity correction, and decision
stability. *PsyArXiv*.
[doi:10.31234/osf.io/7pqz4_v2](https://doi.org/10.31234/osf.io/7pqz4_v2)

Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.

## See also

[`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)
for the simulations,
[`RMdifGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaPlot.md)
for the plot,
[`RMdifAnchor()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifAnchor.md)
for anchored tests,
[`RMdifImpact()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifImpact.md)
for the size of the impact.

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

  # Simulation-based cutoffs (100 iterations to keep the example quick)
  if (requireNamespace("ggdist", quietly = TRUE)) {
    cutoff_res <- RMdifGammaCutoff(sim_data, dif_var = dif_group,
                                   iterations = 100, parallel = FALSE,
                                   seed = 42)
    RMdifGamma(sim_data, dif_group, cutoff = cutoff_res)

    # The full cutoff object switches on bootstrap p-values with the
    # family-wise (Westfall-Young) correction. Use 1000 or more
    # iterations in real analyses.
    RMdifGamma(sim_data, dif_group, cutoff = cutoff_res,
               output = "dataframe")
  }
}
#> Registered S3 methods overwritten by 'vcdExtra':
#>   method      from
#>   pairs.table vcd 
#>   print.Kappa vcd 
#> Bootstrap p-values are based on 100 iterations, below the recommended 400.
#> ℹ The family-wise error rate is controlled at any number of iterations. With few iterations the p-values are coarse, and which comparisons get flagged depends on the seed.
#> ℹ See Johansson (2026), doi:10.31234/osf.io/7pqz4_v2.
#> ℹ Raise `iterations` in RMdifGammaCutoff().
#> This message is displayed once per session.
#>      Item       gamma        se      lower     upper  gamma_low gamma_high
#> 1   Item1  0.02929936 0.1612158 -0.2866779 0.3452766 -0.2791519  0.3475936
#> 2   Item2 -0.16753927 0.1561349 -0.4735580 0.1384794 -0.2672522  0.3369713
#> 3   Item3  0.11111111 0.1681747 -0.2185053 0.4407275 -0.4172185  0.2781790
#> 4   Item4 -0.03045685 0.1545831 -0.3334342 0.2725205 -0.3276904  0.2719986
#> 5   Item5  0.07616708 0.1587018 -0.2348828 0.3872169 -0.3374844  0.1909814
#> 6   Item6  0.04239401 0.1562001 -0.2637526 0.3485406 -0.2731989  0.3440197
#> 7   Item7 -0.14560440 0.1603616 -0.4599073 0.1686985 -0.3110048  0.3402597
#> 8   Item8  0.01897019 0.1618599 -0.2982693 0.3362097 -0.2474490  0.3723958
#> 9   Item9 -0.14088398 0.1614016 -0.4572253 0.1754573 -0.3283396  0.2865475
#> 10 Item10  0.17777778 0.1556159 -0.1272237 0.4827792 -0.2857143  0.3036849
#>      p_gamma padj_gamma flagged
#> 1  1.0000000   1.000000   FALSE
#> 2  0.2871287   0.950495   FALSE
#> 3  0.5247525   1.000000   FALSE
#> 4  0.9009901   1.000000   FALSE
#> 5  0.5247525   1.000000   FALSE
#> 6  0.8415842   1.000000   FALSE
#> 7  0.2970297   0.970297   FALSE
#> 8  0.9207921   1.000000   FALSE
#> 9  0.3960396   0.970297   FALSE
#> 10 0.2970297   0.950495   FALSE
# }
```
