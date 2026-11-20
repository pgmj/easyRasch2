# Simulation-based partial gamma DIF cutoffs

Simulates the partial gamma DIF coefficient (Bjorner et al., 1998,
Henninger et al., 2025) of each item when there is no DIF. Passed to
[`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md),
the result gives bootstrap p-values corrected across items, and an
interval per item. Every simulated dataset keeps each respondent's
observed group membership and total score, so the null reflects the
observed group sizes and any difference between the groups' latent
distributions.

## Usage

``` r
RMdifGammaCutoff(
  data,
  dif_var,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95,
  dgp = c("conditional", "permutation")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0. Only
  complete cases (in `data` and `dif_var`) are used.

- dif_var:

  A vector (factor, character or integer) of length `nrow(data)` with
  the groups. Group membership is held fixed in every simulated dataset,
  so the group sizes and each group's distribution of total scores are
  those observed.

- iterations:

  Integer. Number of simulated datasets. Default `400`. Use 1000 to 2000
  for a final analysis.

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

- cutoff_method:

  Character. How the interval per item is computed: `"hdci"` (default),
  the highest-density continuous interval from
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html)
  (needs ggdist), or `"quantile"`, the 2.5th to 97.5th percentiles.

- hdci_width:

  Numeric. Width of the interval with `cutoff_method = "hdci"`. Default
  `0.95`, as in
  [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md).
  Ignored with `"quantile"`.

- dgp:

  Character. How the datasets without DIF are simulated. `"conditional"`
  (default) draws each respondent's response pattern from the model
  given their observed total score, with the CML item thresholds fixed,
  as the conditional option of
  [`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreCutoff.md)
  does. `"permutation"` keeps the observed responses and permutes group
  membership among respondents with the same total score, the Monte
  Carlo form of an exact conditional test of item-group independence
  given the score (Kreiner, 1987). It needs no item parameters, keeps
  any misfit elsewhere in the observed responses, and is 25 to 60 times
  faster per iteration. Both held the nominal error rate in simulation,
  including when the groups differed by 1 SD.

## Value

A list with:

- `results`:

  data.frame with `iteration`, `Item` and `gamma`, one row per item and
  simulated dataset.

- `item_cutoffs`:

  data.frame with the interval per item: `Item`, `gamma_low`,
  `gamma_high`, computed by `cutoff_method`.

- `actual_iterations`:

  Number of simulated datasets that could be used.

- `requested_iterations`:

  Number of iterations asked for.

- `sample_n`:

  Number of complete cases used.

- `sample_n_total`:

  Number of respondents in the raw input data, before removing rows with
  `NA` in `data` or `dif_var`.

- `sample_has_na`:

  Logical. Whether `data` or `dif_var` contained any missing values.

- `sample_summary`:

  Summary statistics of the WLE person locations, or `NULL` when the
  model could not be fitted with `dgp = "permutation"`.

- `item_names`:

  Character vector of item names from data.

- `dif_group_sizes`:

  Integer vector of group sizes, held fixed in every simulated dataset.

- `cutoff_method`:

  The method used to compute cutoffs (`"hdci"` or `"quantile"`).

- `hdci_width`:

  The HDCI width used (only meaningful when `cutoff_method = "hdci"`).

- `dgp`:

  The data-generating process used.

## Details

Partial gamma conditions on the total score, and under the Rasch model
the response to an item is independent of group membership given the
total score. Both ways of simulating keep each respondent's group and
total score, so the simulations have the observed distribution of groups
over the total scores. Assigning simulated respondents to groups at
random, as an earlier version did, gives both groups the same latent
distribution. When a minority group differed from the majority by 1 SD,
that was too narrow and flagged at least one item in about 14% of
datasets without DIF, against a nominal 5%.

In each iteration one dataset is simulated (see `dgp`) and the partial
gamma of every item computed. The computation reproduces
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html)
exactly but is vectorised, so iarm is not needed here. Datasets that
cannot be used are left out, and `actual_iterations` gives the number
used.

The conditional data-generating process uses CML item thresholds from
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category PCM). The group order follows the
levels of `dif_var` (alphabetical for character vectors), which sets the
sign of gamma as in
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).

The simulations have no DIF in any item, so the cutoffs do not guard
against artificial DIF, which real DIF in one item induces in the others
(Hagquist & Andrich, 2017). See "Real and artificial DIF" in
[`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md).

Parallel processing needs the optional mirai package.

## References

Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
(1998). Differential item functioning in the Danish translation of the
SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189–1202.
[doi:10.1016/S0895-4356(98)00111-5](https://doi.org/10.1016/S0895-4356%2898%2900111-5)

Hagquist, C., & Andrich, D. (2017). Recent advances in analysis of
differential item functioning in health research using the Rasch model.
*Health and Quality of Life Outcomes, 15*, 181.
[doi:10.1186/s12955-017-0755-0](https://doi.org/10.1186/s12955-017-0755-0)

Henninger, M., Radek, J., Sengewald, M.-A., & Strobl, C. (2025). Partial
credit trees meet the partial gamma coefficient for quantifying DIF and
DSF in polytomous items. *Behaviormetrika, 52*, 221–257.
[doi:10.1007/s41237-024-00252-3](https://doi.org/10.1007/s41237-024-00252-3)

Kreiner, S. (1987). Analysis of multidimensional contingency tables by
exact conditional tests: Techniques and strategies. *Scandinavian
Journal of Statistics, 14*(2), 97–112.

## See also

[`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md)
for the test,
[`RMdifGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaPlot.md)
for the plot,
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).

## Examples

``` r
# \donttest{
if (requireNamespace("ggdist", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)
  dif_sex <- sample(c("male", "female"), 200, replace = TRUE)

  # 100 iterations, sequentially, for a quick example
  cutoff_res <- RMdifGammaCutoff(sim_data, dif_var = dif_sex,
                                 iterations = 100, parallel = FALSE,
                                 seed = 42)
  cutoff_res$item_cutoffs

  # The permutation null needs no item parameters and is much faster
  perm_res <- RMdifGammaCutoff(sim_data, dif_var = dif_sex,
                               iterations = 100, parallel = FALSE,
                               seed = 42, dgp = "permutation")
  perm_res$item_cutoffs
}
#>      Item  gamma_low gamma_high
#> 1   Item1 -0.2858958  0.3232323
#> 2   Item2 -0.2365457  0.2855346
#> 3   Item3 -0.3102493  0.3402597
#> 4   Item4 -0.3560976  0.2978986
#> 5   Item5 -0.3096695  0.2663317
#> 6   Item6 -0.2935323  0.2971014
#> 7   Item7 -0.2694938  0.2739362
#> 8   Item8 -0.3387534  0.2286220
#> 9   Item9 -0.3145161  0.2383489
#> 10 Item10 -0.2345310  0.3630881
# }
```
