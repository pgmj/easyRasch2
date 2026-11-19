# Category Fit Test Given the Total Score

Tests, for each item, whether its cumulative response events \\X_i \ge
k\\ depend on the total score as the partial credit model predicts. The
prediction is \\P(X_i \ge k \mid r)\\, which under the model depends on
the item parameters alone, so no person locations are estimated. The
test detects categories that respondents do not distinguish, categories
used in the wrong order, and steps that separate respondents more or
less than the scoring assumes, including when the threshold estimates
are ordered.
[`RMitemCatFitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFitPlot.md)
shows where an item departs.

## Usage

``` r
RMitemCatFit(
  data,
  iterations = 1000L,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  item_labels = NULL,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  output = c("kable", "dataframe", "list")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0.
  Respondents with missing responses are left out (complete cases),
  since the total score must cover the same items for everyone.

- iterations:

  Integer. Number of datasets simulated from the model. Default `1000L`.
  The p-values are counts out of the number of simulated datasets that
  could be refitted plus one (`iterations_used + 1`).

- correction:

  Character. Correction of the item p-values across items: `"fwer"`
  (default, Westfall-Young step-down, controls the family-wise error
  rate), `"fdr_bh"`, `"fdr_by"` (false discovery rate) or `"none"`.

- item_labels:

  Optional character vector of item labels, one per column of `data`.
  Default: column names.

- parallel, n_cores:

  Logical / integer. Run the simulations in parallel via mirai when
  available. Default sequential. When `n_cores` is `NULL`,
  `getOption("mc.cores")` is used, and two workers when that is unset
  too.

- seed:

  Optional integer for reproducible simulations. See
  [easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md).

- output:

  Character. `"kable"` (default), `"dataframe"`, or `"list"`. The list
  also holds the simulations, so that
  [`RMitemCatFitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFitPlot.md)
  can draw every plot from it (`fit` argument) without simulating again.

## Value

- `output = "kable"`: a
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table,
  one row per item, with the event p-values, the item p and the
  corrected item p.

- `output = "dataframe"`: a data.frame with columns `item`, `mcb_1`,
  ..., `mcb_m` (observed MCB per event, `NA` for events an item does not
  have), `p_1`, ..., `p_m` (uncorrected event p), `p` (item p) and
  `padj` (item p corrected across items).

- `output = "list"`: a list of class `RMitemCatFit` with `table` (the
  data.frame above), `correction`, `iterations` (requested),
  `iterations_used` (refitted), `n_used`, `n_total`, `qualifiers`,
  `item_labels`, and `prep` and `sims`, the data and simulations that
  [`RMitemCatFitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFitPlot.md)
  uses.

## Details

**Statistic.** For each item i and event \\X_i \ge k\\, the observed
events are regressed on the predicted probability \\p = P(X_i \ge k \mid
r)\\ by pool-adjacent-violators (PAV). The statistic is the
miscalibration component of the Brier score in the CORP decomposition
(Dimitriadis et al., 2021), \$\$\mathrm{MCB} = \frac{1}{n}\sum_v (p_v -
y_v)^2 - \frac{1}{n}\sum_v (\hat c_v - y_v)^2,\$\$ where \\\hat c_v\\ is
the PAV-recalibrated probability. MCB is zero when the PAV curve lies on
the diagonal of the calibration plot in
[`RMitemCatFitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFitPlot.md)
and grows with its distance from it. Respondents with the lowest or
highest possible total score are left out, since their responses are
fixed by the score.

**Reference distribution.** Response patterns are simulated from the
model given each respondent's total score, at the item parameters
estimated by conditional maximum likelihood, and the parameters are
re-estimated in every simulated dataset. The fit of the prediction to
the same data is therefore part of the reference distribution. Datasets
in which some category of some item is empty cannot be refitted and are
dropped.

**p-values.** For each event, p is the share of the observed MCB and the
simulated ones that lie at or above the observed value. To combine the
events of an item, each event's MCB is studentised by the mean and SD of
the observed and simulated values together, and the item statistic is
the largest of them. The item p is the share of the observed and
simulated datasets whose item statistic is at or above the observed one.
With `"fwer"` the item p-values are corrected across items by the
Westfall-Young step-down on the item statistics (Ge et al., 2003).

In a simulation study (5 or 10 items with five categories, n = 250 to
1000, 199 simulated datasets) the item p rejected at a rate between .039
and .051 at .05 when the model held, and the family-wise rate over items
was .030 to .074 per condition (.041 over all of them). A check of the
condition with the highest rate (n = 1000, 5 items) with 999 simulated
datasets gave .052 with ordered and .045 with disordered thresholds
(Monte Carlo SE .007). The test detected interchangeable or swapped
categories and a step that discriminated twice as strongly as the
others, also when the thresholds were ordered. It rarely detected "don't
know" use in up to 10% of responses or an extreme response style shared
by all items.

## References

Dimitriadis, T., Gneiting, T., & Jordan, A. I. (2021). Stable
reliability diagrams for probabilistic classifiers. *Proceedings of the
National Academy of Sciences, 118*(8), e2016191118.
[doi:10.1073/pnas.2016191118](https://doi.org/10.1073/pnas.2016191118)

Ge, Y., Dudoit, S., & Speed, T. P. (2003). Resampling-based multiple
testing for microarray data analysis. *Test, 12*(1), 1–77.
[doi:10.1007/BF02595811](https://doi.org/10.1007/BF02595811)

## See also

[`RMitemCatFitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFitPlot.md)
for the plots behind the test.

## Examples

``` r
# \donttest{
if (requireNamespace("eRm", quietly = TRUE)) {
  data(pcmdat2, package = "eRm")
  fit <- RMitemCatFit(pcmdat2, iterations = 200, seed = 1,
                      output = "list")
  fit$table
  RMitemCatFitPlot(pcmdat2, type = "calibration", fit = fit)
}
#> Bootstrap p-values are based on 200 iterations, below the recommended 400.
#> ℹ The family-wise error rate is controlled at any number of iterations. With few iterations the p-values are coarse, and which comparisons get flagged depends on the seed.
#> ℹ See Johansson (2026), doi:10.31234/osf.io/7pqz4_v2.
#> ℹ Raise `iterations` in RMitemCatFit().
#> This message is displayed once per session.

# }
```
