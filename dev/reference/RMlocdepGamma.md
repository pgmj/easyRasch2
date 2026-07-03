# Partial Gamma Local Dependence Analysis

Computes partial gamma coefficients for Local Dependence (LD) assessment
using
[`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html).
Each pair of items is tested for residual association, controlling for
the rest score (total score minus one of the items in the pair).

## Usage

``` r
RMlocdepGamma(
  data,
  cutoff = NULL,
  p_value = FALSE,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable",
  n_pairs = NULL
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed, but at least one complete case must exist.

- cutoff:

  Optional. Default `NULL` (no cutoff applied). Can be:

  - The return value of
    [`RMlocdepGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md)
    (a list with `$pair_cutoffs`): the data.frame is extracted
    automatically and simulation metadata is included in the kable
    caption.

  - The `$pair_cutoffs` data.frame from
    [`RMlocdepGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md)
    directly: must have columns `Item1`, `Item2`, `gamma_low`,
    `gamma_high`. When provided, adds columns `Gamma_low`, `Gamma_high`,
    and `Flagged` (logical; `TRUE` when the observed partial gamma falls
    outside the credible range) to the result.

- p_value:

  Logical. When `TRUE`, adds one-sided bootstrap p-values for *excess
  positive* local dependence (`p_gamma`, `padj_gamma`), matching the
  `p_value` semantics of
  [`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md),
  and `flagged` reflects `padj_gamma < alpha` (positive deviations only)
  instead of the credible range. One test per item pair: the p-value is
  computed in the canonical direction (direction 1, rest score = total -
  Item2, the direction that was simulated) and repeated in the
  direction-2 table for the same pair. The asymptotic BH-adjusted
  p-value and star columns from
  [`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html)
  are **dropped** in this mode; the simulated `gamma_low` / `gamma_high`
  band is kept as the effect-size reference. Requires the **full**
  [`RMlocdepGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md)
  object as `cutoff` (it carries the simulated distributions in
  `$results`). Default `FALSE`.

- correction:

  Character. Multiplicity correction for the bootstrap p-values, applied
  over the family of all item pairs (before any `n_pairs` display
  filter): `"fwer"` (default; Westfall-Young studentised-max step-down),
  `"fdr_bh"`, `"fdr_by"`, or `"none"`. Ignored when `p_value = FALSE`.

- alpha:

  Numeric in (0, 1). Significance level used to flag pairs on the
  corrected p-value. Default `0.05`. Ignored when `p_value = FALSE`.

- output:

  Character string controlling the return value. Either `"kable"`
  (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"` for the underlying data.frame.

- n_pairs:

  Optional positive integer. When supplied, only the `n_pairs` item
  pairs with the largest absolute partial-gamma values (i.e., strongest
  residual dependence in either direction) are retained per rest-score
  direction, sorted by `|gamma|` descending. When `NULL` (default), all
  pairs are returned in `iarm`'s native ordering. Values larger than the
  total number of pairs are silently capped at that total.

## Value

- If `output = "kable"`: an object of class `"RMlocdepGamma"`.
  Internally a list with two `knitr_kable` elements, `$direction1` and
  `$direction2`. In both, the rest score is the total score minus Item 2
  (the second column); the two elements list each item pair in the two
  possible orders, so together they cover both rest-score directions for
  every pair. Each has columns "Item 1", "Item 2", "Partial gamma",
  "Adj. p-value (BH)", and "p-value sign." (a star-string indicator from
  [`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html)).
  When `cutoff` is provided, additional columns "Gamma low", "Gamma
  high", and "Flagged" are included.

  The object has custom [`print()`](https://rdrr.io/r/base/print.html)
  and
  [`knitr::knit_print()`](https://rdrr.io/pkg/knitr/man/knit_print.html)
  methods: in the R console it prints the two tables stacked vertically;
  in a Quarto / R Markdown chunk it renders as two distinct pipe tables.
  Access the individual tables explicitly as `result$direction1` and
  `result$direction2` if needed.

- If `output = "dataframe"`: a named list of two data.frames
  (`$direction1`, `$direction2`) with columns `Item1`, `Item2`, `gamma`,
  `padj_bh`, `Significance`. When `cutoff` is provided, columns
  `gamma_low`, `gamma_high`, and `flagged` are also included. With
  `p_value = TRUE`, `padj_bh` and `Significance` are replaced by
  `p_gamma` and `padj_gamma` (identical for a pair in both directions).

## Details

Partial gamma (Christensen, Kreiner & Mesbah, 2013) measures the
residual association between pairs of items after controlling for the
rest score (total score minus one item). Because it matters which item
is subtracted, calculations are done for each pair in both directions,
yielding two data.frames.

Values near 0 indicate no local dependence. Large positive values
suggest positive LD (items share variance beyond the latent trait),
while large negative values suggest negative LD.

The `iarm` package must be installed (it is in Suggests, not Imports).

**Bootstrap p-values.** When `p_value = TRUE`, each pair's observed
partial gamma (canonical direction) is compared against its simulated
null distribution (from `cutoff$results`, simulated under local
independence). The per-pair statistic is the residual studentised by the
bootstrap mean and SD; the marginal p-value is the one-sided Monte-Carlo
p-value `(1 + #\{t* >= t\}) / (B + 1)` for excess *positive* LD
(redundancy, the diagnostic target — matching
[`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)),
so it can be no smaller than `1 / (B + 1)`. The band still shows both
bounds for reference. `correction = "fwer"` uses the Westfall-Young
studentised-max step-down over the family of all pairs, which exploits
the bootstrap dependence among them (Ferreira, 2024); it is liberal when
the simulation is small, so at least 1000 `iterations` in
[`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md)
are recommended (a warning is issued below that). Unlike the asymptotic
p-values from
[`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html),
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

Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013). *Rasch
Models in Health*, pp. 133–135. ISTE & Wiley.
[doi:10.1002/9781118574454](https://doi.org/10.1002/9781118574454)

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously. *Journal of Statistical Theory and
Practice, 19*(6).
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.

## See also

[`RMlocdepGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md),
[`RMlocdepGammaPlot`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaPlot.md)

## Examples

``` r
# \donttest{
if (requireNamespace("iarm", quietly = TRUE)) {
  set.seed(42)
  sim_data <- as.data.frame(
    matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
  )
  colnames(sim_data) <- paste0("Item", 1:10)

  # Default kable output
  RMlocdepGamma(sim_data)

  # Return as data.frame list
  RMlocdepGamma(sim_data, output = "dataframe")

  # Simulation-based cutoffs (slow): 100+ Monte-Carlo iterations
  if (requireNamespace("ggdist", quietly = TRUE)) {
    cutoff_res <- RMlocdepGammaCutoff(sim_data, iterations = 100,
                                      parallel = FALSE, seed = 42)
    RMlocdepGamma(sim_data, cutoff = cutoff_res)

    # Bootstrap p-values with family-wise (Westfall-Young) correction
    # (use iterations >= 1000 in real analyses for stable p-values)
    RMlocdepGamma(sim_data, cutoff = cutoff_res, p_value = TRUE,
                  output = "dataframe")
  }
}
#> Warning: Bootstrap p-values are based on only 100 simulation iterations. With few iterations the studentised-max (FWER) correction is liberal and small p-values are imprecise; use iterations >= 1000 in RMlocdepGammaCutoff() for reliable p-values.
#> $direction1
#>    Item1  Item2  gamma gamma_low gamma_high p_gamma padj_gamma flagged
#> 1  Item1  Item2  0.320    -0.425      0.320  0.0198     0.6238   FALSE
#> 2  Item1  Item3 -0.143    -0.473      0.440  0.8119     1.0000   FALSE
#> 3  Item1  Item4 -0.069    -0.446      0.326  0.7030     1.0000   FALSE
#> 4  Item1  Item5 -0.153    -0.419      0.446  0.8020     1.0000   FALSE
#> 5  Item1  Item6 -0.031    -0.389      0.337  0.6238     1.0000   FALSE
#> 6  Item1  Item7 -0.163    -0.324      0.453  0.9109     1.0000   FALSE
#> 7  Item1  Item8  0.191    -0.426      0.475  0.1881     1.0000   FALSE
#> 8  Item1  Item9  0.124    -0.320      0.400  0.1980     1.0000   FALSE
#> 9  Item1 Item10 -0.039    -0.321      0.317  0.6139     1.0000   FALSE
#> 10 Item2  Item3 -0.206    -0.342      0.369  0.9109     1.0000   FALSE
#> 11 Item2  Item4 -0.202    -0.282      0.451  0.9109     1.0000   FALSE
#> 12 Item2  Item5  0.058    -0.393      0.418  0.4257     1.0000   FALSE
#> 13 Item2  Item6  0.025    -0.382      0.336  0.5149     1.0000   FALSE
#> 14 Item2  Item7 -0.230    -0.561      0.398  0.8614     1.0000   FALSE
#> 15 Item2  Item8  0.150    -0.431      0.390  0.1980     1.0000   FALSE
#> 16 Item2  Item9  0.154    -0.411      0.375  0.1683     1.0000   FALSE
#> 17 Item2 Item10  0.024    -0.353      0.361  0.3663     1.0000   FALSE
#> 18 Item3  Item4  0.021    -0.370      0.317  0.4257     1.0000   FALSE
#> 19 Item3  Item5  0.190    -0.455      0.434  0.1386     0.9901   FALSE
#> 20 Item3  Item6 -0.143    -0.417      0.446  0.7327     1.0000   FALSE
#> 21 Item3  Item7  0.086    -0.414      0.442  0.2970     1.0000   FALSE
#> 22 Item3  Item8 -0.042    -0.435      0.418  0.5941     1.0000   FALSE
#> 23 Item3  Item9 -0.074    -0.379      0.347  0.6634     1.0000   FALSE
#> 24 Item3 Item10  0.221    -0.351      0.458  0.1089     0.9901   FALSE
#> 25 Item4  Item5 -0.024    -0.360      0.417  0.5248     1.0000   FALSE
#> 26 Item4  Item6  0.224    -0.334      0.456  0.0990     0.9901   FALSE
#> 27 Item4  Item7  0.333    -0.397      0.344  0.0198     0.5743   FALSE
#> 28 Item4  Item8 -0.406    -0.369      0.348  1.0000     1.0000   FALSE
#> 29 Item4  Item9 -0.153    -0.397      0.334  0.8218     1.0000   FALSE
#> 30 Item4 Item10 -0.084    -0.364      0.346  0.6931     1.0000   FALSE
#> 31 Item5  Item6 -0.136    -0.387      0.332  0.7921     1.0000   FALSE
#> 32 Item5  Item7  0.157    -0.462      0.450  0.1584     1.0000   FALSE
#> 33 Item5  Item8  0.130    -0.371      0.438  0.1485     1.0000   FALSE
#> 34 Item5  Item9 -0.348    -0.370      0.407  0.9901     1.0000   FALSE
#> 35 Item5 Item10  0.014    -0.405      0.320  0.3762     1.0000   FALSE
#> 36 Item6  Item7  0.004    -0.450      0.375  0.4356     1.0000   FALSE
#> 37 Item6  Item8 -0.133    -0.250      0.367  0.8812     1.0000   FALSE
#> 38 Item6  Item9  0.062    -0.382      0.390  0.3465     1.0000   FALSE
#> 39 Item6 Item10 -0.243    -0.307      0.484  0.9604     1.0000   FALSE
#> 40 Item7  Item8 -0.198    -0.443      0.366  0.8713     1.0000   FALSE
#> 41 Item7  Item9 -0.005    -0.365      0.364  0.5446     1.0000   FALSE
#> 42 Item7 Item10  0.021    -0.377      0.298  0.4653     1.0000   FALSE
#> 43 Item8  Item9  0.301    -0.410      0.316  0.0297     0.7426   FALSE
#> 44 Item8 Item10  0.045    -0.417      0.400  0.3861     1.0000   FALSE
#> 45 Item9 Item10  0.036    -0.358      0.349  0.3762     1.0000   FALSE
#> 
#> $direction2
#>     Item1 Item2  gamma gamma_low gamma_high p_gamma padj_gamma flagged
#> 1   Item2 Item1  0.300    -0.425      0.320  0.0198     0.6238   FALSE
#> 2   Item3 Item1 -0.158    -0.473      0.440  0.8119     1.0000   FALSE
#> 3   Item3 Item2 -0.230    -0.342      0.369  0.9109     1.0000   FALSE
#> 4   Item4 Item1 -0.141    -0.446      0.326  0.7030     1.0000   FALSE
#> 5   Item4 Item2 -0.263    -0.282      0.451  0.9109     1.0000   FALSE
#> 6   Item4 Item3 -0.042    -0.370      0.317  0.4257     1.0000   FALSE
#> 7   Item5 Item1 -0.152    -0.419      0.446  0.8020     1.0000   FALSE
#> 8   Item5 Item2  0.039    -0.393      0.418  0.4257     1.0000   FALSE
#> 9   Item5 Item3  0.201    -0.455      0.434  0.1386     0.9901   FALSE
#> 10  Item5 Item4  0.021    -0.360      0.417  0.5248     1.0000   FALSE
#> 11  Item6 Item1 -0.102    -0.389      0.337  0.6238     1.0000   FALSE
#> 12  Item6 Item2 -0.063    -0.382      0.336  0.5149     1.0000   FALSE
#> 13  Item6 Item3 -0.208    -0.417      0.446  0.7327     1.0000   FALSE
#> 14  Item6 Item4  0.206    -0.334      0.456  0.0990     0.9901   FALSE
#> 15  Item6 Item5 -0.201    -0.387      0.332  0.7921     1.0000   FALSE
#> 16  Item7 Item1 -0.116    -0.324      0.453  0.9109     1.0000   FALSE
#> 17  Item7 Item2 -0.195    -0.561      0.398  0.8614     1.0000   FALSE
#> 18  Item7 Item3  0.211    -0.414      0.442  0.2970     1.0000   FALSE
#> 19  Item7 Item4  0.428    -0.397      0.344  0.0198     0.5743   FALSE
#> 20  Item7 Item5  0.218    -0.462      0.450  0.1584     1.0000   FALSE
#> 21  Item7 Item6  0.093    -0.450      0.375  0.4356     1.0000   FALSE
#> 22  Item8 Item1  0.237    -0.426      0.475  0.1881     1.0000   FALSE
#> 23  Item8 Item2  0.220    -0.431      0.390  0.1980     1.0000   FALSE
#> 24  Item8 Item3  0.024    -0.435      0.418  0.5941     1.0000   FALSE
#> 25  Item8 Item4 -0.319    -0.369      0.348  1.0000     1.0000   FALSE
#> 26  Item8 Item5  0.152    -0.371      0.438  0.1485     1.0000   FALSE
#> 27  Item8 Item6 -0.022    -0.250      0.367  0.8812     1.0000   FALSE
#> 28  Item8 Item7 -0.222    -0.443      0.366  0.8713     1.0000   FALSE
#> 29  Item9 Item1  0.193    -0.320      0.400  0.1980     1.0000   FALSE
#> 30  Item9 Item2  0.199    -0.411      0.375  0.1683     1.0000   FALSE
#> 31  Item9 Item3 -0.008    -0.379      0.347  0.6634     1.0000   FALSE
#> 32  Item9 Item4 -0.071    -0.397      0.334  0.8218     1.0000   FALSE
#> 33  Item9 Item5 -0.256    -0.370      0.407  0.9901     1.0000   FALSE
#> 34  Item9 Item6  0.218    -0.382      0.390  0.3465     1.0000   FALSE
#> 35  Item9 Item7 -0.001    -0.365      0.364  0.5446     1.0000   FALSE
#> 36  Item9 Item8  0.357    -0.410      0.316  0.0297     0.7426   FALSE
#> 37 Item10 Item1 -0.071    -0.321      0.317  0.6139     1.0000   FALSE
#> 38 Item10 Item2 -0.023    -0.353      0.361  0.3663     1.0000   FALSE
#> 39 Item10 Item3  0.176    -0.351      0.458  0.1089     0.9901   FALSE
#> 40 Item10 Item4 -0.088    -0.364      0.346  0.6931     1.0000   FALSE
#> 41 Item10 Item5 -0.029    -0.405      0.320  0.3762     1.0000   FALSE
#> 42 Item10 Item6 -0.235    -0.307      0.484  0.9604     1.0000   FALSE
#> 43 Item10 Item7 -0.040    -0.377      0.298  0.4653     1.0000   FALSE
#> 44 Item10 Item8 -0.047    -0.417      0.400  0.3861     1.0000   FALSE
#> 45 Item10 Item9 -0.069    -0.358      0.349  0.3762     1.0000   FALSE
#> 
# }
```
