# Conditional item infit MSQ

Computes the conditional infit mean-square (MSQ) of each item with
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html), and
each item's location relative to the mean person location. Given the
result of
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md),
it also flags misfitting items, by default on a bootstrap p-value
corrected across items.

## Usage

``` r
RMitemInfit(
  data,
  cutoff = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable",
  sort,
  statistic = "infit"
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0. Missing
  values (`NA`) are allowed, but the infit statistic uses complete cases
  only (see Details), so at least one is needed.

- cutoff:

  Optional. Default `NULL`, no flagging. Either of:

  - The full result of
    [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md).
    This enables the p-values (see `p_value`), and the caption reports
    the number of iterations and the interval method. A result of
    [`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md)
    gives an error, since its simulations come from imputed data. Use it
    with
    [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md).

  - Its `$item_cutoffs` data.frame alone, with columns `Item`,
    `infit_low` and `infit_high` (or `outfit_low` and `outfit_high` with
    `statistic = "outfit"`). Items are then flagged against the
    interval.

  Either adds the columns `Infit_low`, `Infit_high` and `Flagged`.
  `Flagged` gives the direction of misfit: `"overfit"` (below 1, or
  below the interval), `"underfit"` (above) or `""` (no misfit).

- p_value:

  Logical or `NULL`. Whether to flag on bootstrap p-values computed from
  the simulated values in `cutoff`. `NULL` (default) uses them whenever
  `cutoff` is the full
  [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
  result, and flags against the interval otherwise. `TRUE` without the
  full result is an error.

  With `FALSE`, items are flagged against the interval and a one-time
  message gives the family-wise error rate that implies. The interval
  describes where a fitting item's statistic is expected to fall. Used
  as a decision rule over `k` items, a width `w` gives a family-wise
  error rate of `1 - w^k`, 37% for the default 95% interval over nine
  items. The corrected p-value keeps the error rate at `alpha` and needs
  far fewer iterations to do so (Johansson, 2026). Versions up to 1.1.1
  flagged against the interval by default, so older scripts that pass
  the full cutoff object may now flag different items.

- correction:

  Character. Correction of the p-values across items when they are used:
  `"fwer"` (default, Westfall-Young step-down, controls the family-wise
  error rate), `"fdr_bh"` or `"fdr_by"` (Benjamini-Hochberg or
  Benjamini-Yekutieli, control the false discovery rate), or `"none"`.

- alpha:

  Numeric in (0, 1). An item is flagged when its corrected p-value is
  below `alpha`. Default `0.05`. Used only with p-values.

- output:

  Character. `"kable"` (default) for a
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"`.

- sort:

  Optional. `sort = "infit"` (or `"outfit"`, matching `statistic`) sorts
  the rows by MSQ, largest first. Other values leave the order
  unchanged.

- statistic:

  Character. Which conditional fit statistic to report and test:
  `"infit"` (default) or `"outfit"`. Both come from the same
  [`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html)
  call and the same simulations, and the column names, interval,
  p-values and `Flagged` follow the choice. The correction is applied
  across the k items for the chosen statistic. Testing both statistics
  and reporting whichever flags is a family of 2k tests, larger than
  either call corrects for.

## Value

- `output = "kable"`: a `knitr_kable` table with the columns "Item",
  "Infit MSQ" and "Relative location", and a caption with the sample
  size. With `cutoff`, also "Infit low", "Infit high" and "Flagged", and
  with p-values "p" and "p (adj)".

- `output = "dataframe"`: a data.frame with the columns `Item`,
  `Infit_MSQ` and `Relative_location`. With `cutoff`, also `Infit_low`,
  `Infit_high` and `Flagged` (after `Infit_MSQ`). With p-values, also
  `p_infit` (uncorrected two-sided p) and `padj_infit` (corrected p).
  `Flagged` then marks items with `padj_infit < alpha`, with the
  direction from `Infit_MSQ` relative to 1.

With `statistic = "outfit"` the columns are named after it:
`Outfit_MSQ`, `Outfit_low`, `Outfit_high`, `p_outfit`, `padj_outfit`.

## Details

**Statistic.** Infit MSQ is a weighted fit statistic that emphasises
unexpected responses near the item's location. Values near 1 indicate
fit. Values above 1 indicate underfit (responses less predictable than
the model expects), and values below 1 overfit (responses more
predictable). How far from 1 is too far depends on the sample size and
the items, so it is judged against simulations from
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md).
No fixed value is right in general. The conditional statistic from
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html) uses
the distribution of the responses given the total score (Müller, 2020).

**Data and model.** The item parameters are estimated by conditional
maximum likelihood with
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category partial credit model). The
conditional infit uses complete cases only, and every response category
must be used in them, or the function stops and names the item. The
relative item location is the mean of the item's CML thresholds minus
the mean WLE person location, a description of targeting. The `iarm`
package must be installed (it is in Suggests).

**Bootstrap p-values.** Each item's observed statistic is compared with
its simulated values from `cutoff$results`. The test statistic is the
deviation from the mean of the observed and simulated values together,
divided by their SD. The ZSTD (Wilson-Hilferty) transform is not used,
since it does not work for the conditional MSQ (Müller, 2020). The
uncorrected p is two-sided, `(1 + #{|t*| >= |t|}) / (B + 1)` for B
simulated datasets, so it can be no smaller than `1 / (B + 1)`. With
`"fwer"` the correction is the Westfall-Young step-down on the
studentised maximum, which uses the dependence between the items'
statistics and has more power than Bonferroni or Holm (Westfall & Young,
1993; Ferreira, 2024). Its validity rests on subset pivotality.
`"fdr_bh"` and `"fdr_by"` apply Benjamini-Hochberg or
Benjamini-Yekutieli. The p-values depend on the model and the sample
size, and are reported with the simulated interval, not instead of it.

**Iterations.** The observed value is part of the mean and SD it is
studentised by, so the observed and simulated statistics are
exchangeable under the model, and the family-wise error rate is at the
nominal level at any number of iterations in
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md).
More iterations make decisions reproducible. Two analysts using
different seeds disagree about at least one item roughly 10% of the time
at 400 iterations and 4% at 2000 (Johansson, 2026). A message says so
below 400 iterations, and the caption below 1000. Use 1000 to 2000 for a
final analysis. The uncorrected p-values are well calibrated even at a
few hundred iterations.

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

## Flags depend on the other items

Item infit and item-restscore both judge each item against expectations
from a Rasch model fitted to all items, so misfitting items change what
is expected of the others. A flag in one direction can therefore produce
flags in the other. Overfitting items can push fitting items into the
underfit range, and underfitting items, for example from a second
dimension, can make fitting items look overfit. How large the effect is
depends on how strong the misfit is and how many items share it. In
simulation with seven polytomous items, two overfitting items led to
about 6% of the fitting items being flagged as underfit, and one or two
mildly multidimensional items to 1% to 2% being flagged as overfit. When
items are flagged in both directions, consider removing the clearest
misfit and testing again before interpreting the rest.

## References

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(1), 6.
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

Johansson, M. (2026). Simulation-based cutoffs for conditional item fit
in Rasch models: Iterations, multiplicity correction, and decision
stability. *PsyArXiv*.
[doi:10.31234/osf.io/7pqz4_v2](https://doi.org/10.31234/osf.io/7pqz4_v2)

Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust
them? *Journal of Statistical Distributions and Applications*, 7(5).
[doi:10.1186/s40488-020-00108-7](https://doi.org/10.1186/s40488-020-00108-7)

Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.

## See also

[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md)
for the simulations,
[`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)
to plot them with the observed values,
[`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md)
for multiply imputed data, and
[`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
for a second item fit test.

## Examples

``` r
# \donttest{
if (requireNamespace("iarm", quietly = TRUE)) {
  # Simulate binary item response data (5 items, 40 persons)
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 40 * 5, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(sim_data) <- paste0("Item", 1:5)

  # Default kable output
  RMitemInfit(sim_data)

  # Sorted by infit MSQ descending
  RMitemInfit(sim_data, sort = "infit")

  # Return as data.frame for further processing
  df <- RMitemInfit(sim_data, output = "dataframe")

  # Simulation-based cutoffs (100 Monte-Carlo iterations)
  if (requireNamespace("ggdist", quietly = TRUE)) {
    cutoff_res <- RMitemInfitCutoff(sim_data, iterations = 100,
                                    parallel = FALSE, seed = 42)
    RMitemInfit(sim_data, cutoff = cutoff_res)
    RMitemInfit(sim_data, cutoff = cutoff_res, output = "dataframe")

    # Bootstrap p-values with family-wise (Westfall-Young) correction
    # (use iterations >= 1000 in real analyses for stable p-values)
    RMitemInfit(sim_data, cutoff = cutoff_res, p_value = TRUE,
                output = "dataframe")

    # The same for conditional outfit MSQ
    RMitemInfit(sim_data, cutoff = cutoff_res, p_value = TRUE,
                statistic = "outfit", output = "dataframe")
  }
}
#>    Item Outfit_MSQ Outfit_low Outfit_high  p_outfit padj_outfit Flagged
#> 1 Item1  1.1294533      0.558       1.366 0.3861386   0.7821782        
#> 2 Item2  0.9381845      0.661       1.399 0.7128713   0.7821782        
#> 3 Item3  0.8125048      0.694       1.381 0.3465347   0.6534653        
#> 4 Item4  1.2755158      0.583       1.379 0.1683168   0.4851485        
#> 5 Item5  0.9035259      0.517       1.457 0.5742574   0.7821782        
#>   Relative_location
#> 1        0.03566442
#> 2       -0.47593330
#> 3       -0.16351968
#> 4        0.32749828
#> 5       -0.69983314
# }
```
