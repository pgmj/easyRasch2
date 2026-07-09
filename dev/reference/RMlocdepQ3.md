# \\Q_3\\ Residual Correlations for Local Dependence Assessment

Computes Yen's \\Q_3\\ residual correlations between item pairs. By
default the Rasch model is fitted by conditional maximum likelihood with
WLE person locations (`estimator = "CML"`); marginal ML via `mirt` is
available with `estimator = "MML"`. High correlations (above the dynamic
cut-off) indicate potential local dependence between items. See
[`RMlocdepQ3Cutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
for how to determine the appropriate dynamic cut-off for your data.

## Usage

``` r
RMlocdepQ3(
  data,
  cutoff = NULL,
  output = "kable",
  n_pairs = NULL,
  p_value = FALSE,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  estimator = c("CML", "MML")
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed.

- cutoff:

  Optional. `NULL` (default) returns the raw \\Q_3\\matrix. A single
  numeric value returns the \\Q_3\\ matrix with a global dynamic cut-off
  (the value added to the mean off-diagonal \\Q_3\\). The **full list**
  returned by
  [`RMlocdepQ3Cutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  returns a *list of two tables* (see Value): the cut-off matrix plus a
  per-pair table.

- output:

  Character string controlling the return value. Either `"kable"`
  (default) for formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table(s),
  or `"dataframe"` for the underlying data.frame(s).

- n_pairs:

  Integer or `NULL` (default). When the full cutoff object is supplied,
  limits the per-pair table to the `n_pairs` pairs with the largest
  departure from their expected range. `NULL` shows all pairs.

- p_value:

  Logical. If `TRUE` (requires the full
  [`RMlocdepQ3Cutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  object), the per-pair table also reports one-sided bootstrap p-values
  (`p_q3`, `padj_q3`) and flags `above` pairs (only) on
  `padj_q3 < alpha` instead of the expected range. Default `FALSE`.

- correction:

  Character. Multiple-comparison correction across item pairs when
  `p_value = TRUE`: `"fwer"` (default) for the Westfall-Young
  studentised-max step-down, `"fdr_bh"` / `"fdr_by"` for
  Benjamini-Hochberg / Benjamini-Yekutieli, or `"none"`.

- alpha:

  Numeric in (0, 1). Significance level for `Flagged` when
  `p_value = TRUE`. Default `0.05`.

- estimator:

  Character. Estimation engine for the \\Q_3\\ residual correlations.
  `"CML"` (default) fits item parameters by conditional maximum
  likelihood (`psychotools`) and person locations by Warm's weighted
  likelihood (WLE) – true to the Rasch tradition and finite at extreme
  scores. `"MML"` uses the marginal-ML / EAP engine (`mirt`), retained
  for backward compatibility. The estimator must match the one used by
  [`RMlocdepQ3Cutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md);
  when the full cutoff object is supplied, its stored estimator takes
  precedence (a mismatching `estimator` argument is overridden with a
  warning).

## Value

With `cutoff = NULL` or a bare numeric cut-off, a single object (`kable`
or data.frame) holding the lower triangle of the \\Q_3\\ matrix; a
numeric cut-off adds a `Flagged` row flag and a caption describing the
dynamic cut-off.

With the **full
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
object**, a named list of two:

- `$matrix`:

  the \\Q_3\\ lower-triangle matrix with the *global* dynamic cut-off
  (mean off-diagonal \\Q_3\\ + suggested cut-off), as above.

- `$pairs`:

  one row per item pair: `Item1`, `Item2`, `Observed` (\\Q_3\\),
  `Low`/`High` (the per-pair expected range, i.e. the simulated bounds),
  and `Flagged` – `"above"` (\\Q_3\\ above the upper bound, indicating
  local dependence), `"below"` (below the lower bound), or `""`. Sorted
  by absolute departure from the per-pair simulated median and truncated
  to `n_pairs`. With `p_value = TRUE`, columns `p_q3` and `padj_q3` are
  added and `Flagged` reflects `padj_q3 < alpha` and **only flags
  `"above"`**.

The \\Q_3\\ tile heatmap that earlier versions returned as `$plot` is
now produced by
[`RMlocdepQ3Plot`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md)
(as its `$matrix` element), so the table and plot outputs share the same
`$matrix`/`$pairs` structure.

## Details

The \\Q_3\\ statistic (Yen, 1984) is the correlation between residuals
of pairs of items after accounting for the latent trait. Under local
independence, \\Q_3\\ values are expected to be around \\-1/(k-1)\\
where \\k\\ is the number of items. When `cutoff` is supplied, the
dynamic cut-off is the mean of all off-diagonal \\Q_3\\ values plus
`cutoff`, following the approach of Christensen et al. (2017). Use
[`RMlocdepQ3Cutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
to obtain a simulation-based cutoff recommendation.

\\Q_3\\ is the column-wise correlation matrix of the model standardized
residuals \\(x - E)/\sqrt{Var}\\. By default (`estimator = "CML"`) item
parameters are estimated by conditional maximum likelihood
(`psychotools`) and person locations by Warm's weighted likelihood;
`estimator = "MML"` instead uses `mirt`'s marginal-ML model and its
built-in \\Q_3\\ residuals. The two estimators give very similar \\Q_3\\
values (off-diagonal correlation typically \> 0.95); what matters for
inference is that the observed \\Q_3\\ and the simulated cut-off in
[`RMlocdepQ3Cutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
use the *same* estimator, which the functions enforce.

**Two views of local dependence.** Given the full
[`RMlocdepQ3Cutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
object, two complementary tables are returned. The `$matrix` applies a
single *global* cut-off (the Christensen et al. approach: the 99th
percentile of the simulated max-minus-mean \\Q_3\\) – a family-wise "is
there any local dependence" overview. The `$pairs` table is the
per-comparison view: each observed \\Q_3\\ against its own simulated
expected range (the `Low`/`High` bounds), so individual dependent pairs
can be read off and ranked.

**Bootstrap p-values.** When `p_value = TRUE`, the `$pairs` table also
tests each observed \\Q_3\\ against its simulated null (from
`cutoff$pair_results`) with a one-sided (upper-tail) test for excess
local dependence. The pair statistic is studentised by the bootstrap
mean and SD; the marginal p-value is `(1 + #{Q3* >= Q3}) / (B + 1)`, and
`correction` applies the family-wise (Westfall-Young step-down) or FDR
adjustment across the \\k(k-1)/2\\ pairs. As for item fit, the
family-wise correction is liberal when the simulation is small, so \>=
1000 `iterations` in
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
are recommended (a warning is issued otherwise).

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

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement, 8*(2), 125–145.
[doi:10.1177/014662168400800201](https://doi.org/10.1177/014662168400800201)

Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
for Yen's \\Q_3\\: Identification of local dependence in the Rasch
model. *Applied Psychological Measurement, 41*(3), 178–194.
[doi:10.1177/0146621616677520](https://doi.org/10.1177/0146621616677520)

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(6).
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

## Examples

``` r
# \donttest{
# Simulate binary item response data (10 items, 200 persons)
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Raw Q3 matrix (no cutoff)
RMlocdepQ3(sim_data)
#> 
#> 
#> Table: Raw Q3 residual correlations (lower triangle). Use RMlocdepQ3Cutoff() to derive a cutoff. n = 200 respondents.
#> 
#> |       |Item1 |Item2 |Item3 |Item4 |Item5 |Item6 |Item7 |Item8 |Item9 |Item10 |
#> |:------|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:-----|:------|
#> |Item1  |      |      |      |      |      |      |      |      |      |       |
#> |Item2  |-0.02 |      |      |      |      |      |      |      |      |       |
#> |Item3  |-0.17 |-0.2  |      |      |      |      |      |      |      |       |
#> |Item4  |-0.13 |-0.2  |-0.1  |      |      |      |      |      |      |       |
#> |Item5  |-0.19 |-0.06 |-0.02 |-0.12 |      |      |      |      |      |       |
#> |Item6  |-0.1  |-0.1  |-0.17 |0.03  |-0.19 |      |      |      |      |       |
#> |Item7  |-0.17 |-0.13 |-0.07 |0.03  |-0.04 |-0.11 |      |      |      |       |
#> |Item8  |-0.03 |-0.06 |-0.11 |-0.28 |-0.05 |-0.16 |-0.18 |      |      |       |
#> |Item9  |-0.03 |-0.08 |-0.16 |-0.19 |-0.22 |-0.03 |-0.1  |-0.02 |      |       |
#> |Item10 |-0.16 |-0.14 |0     |-0.09 |-0.1  |-0.18 |-0.15 |-0.12 |-0.09 |       |

# Get the underlying data.frame
q3_df <- RMlocdepQ3(sim_data, output = "dataframe")

# Simulation-based cutoff (use 500+ iterations in real analyses)
if (requireNamespace("ggdist", quietly = TRUE)) {
  cutoff_res <- RMlocdepQ3Cutoff(sim_data, iterations = 50, parallel = FALSE)

  # Bare numeric cutoff -> just the matrix
  RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)

  # Full object -> list of two tables: $matrix and $pairs
  res <- RMlocdepQ3(sim_data, cutoff = cutoff_res, output = "dataframe")
  res$pairs

  # Top 5 pairs, with bootstrap p-values (use iterations >= 1000 in practice)
  RMlocdepQ3(sim_data, cutoff = cutoff_res, n_pairs = 5, p_value = TRUE,
             output = "dataframe")$pairs
}
#> Warning: Bootstrap p-values are based on only 50 simulation iterations. With few iterations the studentised-max (FWER) correction is liberal and small p-values are imprecise; use iterations >= 1000 in RMlocdepQ3Cutoff() for reliable p-values.
#>   Item1  Item2      Observed        Low       High       p_q3   padj_q3 Flagged
#> 1 Item4  Item8 -0.2765225048 -0.2700737 0.01430429 1.00000000 1.0000000        
#> 2 Item4  Item7  0.0285193666 -0.2640942 0.01619448 0.01960784 0.4509804        
#> 3 Item5  Item9 -0.2199294805 -0.2160826 0.03315679 1.00000000 1.0000000        
#> 4 Item4  Item6  0.0279341896 -0.2592712 0.02734526 0.01960784 0.5098039        
#> 5 Item3 Item10  0.0005837293 -0.2430223 0.02814042 0.05882353 0.8823529        
# }
```
