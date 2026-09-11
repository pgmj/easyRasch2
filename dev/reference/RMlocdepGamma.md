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
  p_value = NULL,
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

  Logical or `NULL`. When `TRUE`, adds one-sided bootstrap p-values for
  *excess positive* local dependence (`p_gamma`, `padj_gamma`), matching
  the `p_value` semantics of
  [`RMlocdepQ3`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md),
  and `flagged` reflects `padj_gamma < alpha` (positive deviations only)
  instead of the credible range. One test per item pair: the p-value is
  computed on `gamma_pair`, the larger of the two conditioning
  directions, and repeated in the direction-2 table for the same pair.
  The asymptotic adjusted p-value and star columns from
  [`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html)
  are **dropped** in this mode, and the simulated `gamma_low` /
  `gamma_high` band is kept as the effect-size reference. `NULL`, the
  default, means `TRUE` when `cutoff` is the **full**
  [`RMlocdepGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md)
  object (it carries the simulated distributions in `$results`) and
  `FALSE` otherwise, so calling with no cutoff keeps the asymptotic
  table unchanged. Pass `FALSE` for the pre-1.2.0 behaviour. The
  interval makes a poor decision rule, since its width sets a
  family-wise error rate of `1 - width^m` over all \\m\\ pairs at once
  (Johansson, 2026).

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
  `se`, `lower`, `upper` (95% Wald CI), `padj_bh`, `Significance`. When
  `cutoff` is provided, columns `gamma_low`, `gamma_high`, and `flagged`
  are also included. With `p_value = TRUE`, `padj_bh` and `Significance`
  are replaced by `p_gamma` and `padj_gamma` (identical for a pair in
  both directions).

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

**Bootstrap p-values.** When `p_value = TRUE`, each pair is tested once,
on the larger of its two rest-score directions (the `gamma_pair`
column), against its simulated null distribution (from `cutoff$results`,
simulated under local independence and built from that same maximum).
The per-pair statistic is the residual studentised by the bootstrap mean
and SD; the marginal p-value is the one-sided Monte-Carlo p-value
`(1 + #\{t* >= t\}) / (B + 1)` for excess *positive* LD (redundancy, the
diagnostic target — matching
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
#> Bootstrap p-values are based on 100 iterations, below the calibrated floor of 400.
#> ℹ Below 400 the Westfall-Young correction is mildly liberal under the null, so the family-wise error rate is above the nominal level.
#> ℹ See Johansson (2026), doi:10.31234/osf.io/7pqz4_v2.
#> ℹ Raise `iterations` in RMlocdepGammaCutoff().
#> This message is displayed once per session.
#> $direction1
#>    Item1  Item2        gamma   gamma_pair        se       lower       upper
#> 1  Item1  Item2  0.319681456  0.319681456 0.1432970  0.03882458  0.60053833
#> 2  Item1  Item3 -0.143187067 -0.143187067 0.1606749 -0.45810418  0.17173004
#> 3  Item1  Item4 -0.069478908 -0.069478908 0.1591124 -0.38133356  0.24237574
#> 4  Item1  Item5 -0.152585120 -0.151515152 0.1603455 -0.46685651  0.16168627
#> 5  Item1  Item6 -0.031210986 -0.031210986 0.1609887 -0.34674305  0.28432107
#> 6  Item1  Item7 -0.163240629 -0.116219668 0.1588380 -0.47455731  0.14807605
#> 7  Item1  Item8  0.190697674  0.236714976 0.1542705 -0.11166703  0.49306238
#> 8  Item1  Item9  0.124413146  0.193026152 0.1601356 -0.18944683  0.43827312
#> 9  Item1 Item10 -0.039190898 -0.039190898 0.1606385 -0.35403655  0.27565475
#> 10 Item2  Item3 -0.205542725 -0.205542725 0.1585150 -0.51622650  0.10514105
#> 11 Item2  Item4 -0.202469136 -0.202469136 0.1546746 -0.50562571  0.10068744
#> 12 Item2  Item5  0.057692308  0.057692308 0.1588853 -0.25371708  0.36910170
#> 13 Item2  Item6  0.024937656  0.024937656 0.1596142 -0.28790040  0.33777571
#> 14 Item2  Item7 -0.230024213 -0.194936709 0.1584601 -0.54060033  0.08055190
#> 15 Item2  Item8  0.150057274  0.219927096 0.1603992 -0.16431942  0.46443397
#> 16 Item2  Item9  0.153846154  0.198589894 0.1591964 -0.15817299  0.46586530
#> 17 Item2 Item10  0.023839398  0.023839398 0.1612274 -0.29216056  0.33983935
#> 18 Item3  Item4  0.020737327  0.020737327 0.1650756 -0.30280487  0.34427953
#> 19 Item3  Item5  0.190424374  0.200878156 0.1569296 -0.11715197  0.49800072
#> 20 Item3  Item6 -0.142857143 -0.142857143 0.1577938 -0.45212740  0.16641312
#> 21 Item3  Item7  0.085653105  0.211469534 0.1656989 -0.23911086  0.41041707
#> 22 Item3  Item8 -0.042222222  0.023752969 0.1687135 -0.37289469  0.28845025
#> 23 Item3  Item9 -0.074398249 -0.008206331 0.1659882 -0.39972914  0.25093264
#> 24 Item3 Item10  0.220689655  0.220689655 0.1571460 -0.08731083  0.52869014
#> 25 Item4  Item5 -0.024390244  0.020656136 0.1543578 -0.32692592  0.27814543
#> 26 Item4  Item6  0.223609535  0.223609535 0.1476403 -0.06576015  0.51297922
#> 27 Item4  Item7  0.332624867  0.428246014 0.1397878  0.05864578  0.60660395
#> 28 Item4  Item8 -0.405895692 -0.318595579 0.1385430 -0.67743505 -0.13435633
#> 29 Item4  Item9 -0.152542373 -0.070631970 0.1561449 -0.45858075  0.15349601
#> 30 Item4 Item10 -0.083839611 -0.083839611 0.1570041 -0.39156194  0.22388271
#> 31 Item5  Item6 -0.135678392 -0.135678392 0.1558843 -0.44120610  0.16984931
#> 32 Item5  Item7  0.157417894  0.218116806 0.1534992 -0.14343499  0.45827078
#> 33 Item5  Item8  0.129807692  0.151960784 0.1592747 -0.18236496  0.44198034
#> 34 Item5  Item9 -0.348182884 -0.255689424 0.1462168 -0.63476249 -0.06160328
#> 35 Item5 Item10  0.013836478  0.013836478 0.1581814 -0.29619341  0.32386636
#> 36 Item6  Item7  0.004484305  0.092682927 0.1586812 -0.30652521  0.31549382
#> 37 Item6  Item8 -0.133409350 -0.021879022 0.1567139 -0.44056304  0.17374434
#> 38 Item6  Item9  0.061728395  0.217503218 0.1577577 -0.24747098  0.37092777
#> 39 Item6 Item10 -0.243309002 -0.234932349 0.1478988 -0.53318525  0.04656724
#> 40 Item7  Item8 -0.197916667 -0.197916667 0.1660299 -0.52332928  0.12749595
#> 41 Item7  Item9 -0.004926108 -0.001236094 0.1665351 -0.33132891  0.32147669
#> 42 Item7 Item10  0.021276596  0.021276596 0.1652290 -0.30256636  0.34511955
#> 43 Item8  Item9  0.301478953  0.357058126 0.1460257  0.01527382  0.58768409
#> 44 Item8 Item10  0.044973545  0.044973545 0.1601356 -0.26888653  0.35883362
#> 45 Item9 Item10  0.035761589  0.035761589 0.1655990 -0.28880658  0.36032976
#>     gamma_low gamma_high    p_gamma padj_gamma flagged
#> 1  -0.2777321  0.3198758 0.01980198  0.6534653   FALSE
#> 2  -0.2571429  0.3843537 0.88118812  1.0000000   FALSE
#> 3  -0.2265372  0.3262195 0.72277228  1.0000000   FALSE
#> 4  -0.3109244  0.3537519 0.83168317  1.0000000   FALSE
#> 5  -0.2680723  0.3297003 0.70297030  1.0000000   FALSE
#> 6  -0.2089041  0.2947559 0.87128713  1.0000000   FALSE
#> 7  -0.3167702  0.3052264 0.16831683  0.9900990   FALSE
#> 8  -0.2711268  0.3800000 0.14851485  0.9900990   FALSE
#> 9  -0.2714777  0.2960000 0.66336634  1.0000000   FALSE
#> 10 -0.2605364  0.3693694 0.91089109  1.0000000   FALSE
#> 11 -0.2141653  0.3376906 0.96039604  1.0000000   FALSE
#> 12 -0.1967865  0.4179104 0.55445545  1.0000000   FALSE
#> 13 -0.3203540  0.2933912 0.62376238  1.0000000   FALSE
#> 14 -0.3580264  0.4016620 0.82178218  1.0000000   FALSE
#> 15 -0.2684564  0.3580705 0.13861386  0.9900990   FALSE
#> 16 -0.3421927  0.2688498 0.13861386  0.9900990   FALSE
#> 17 -0.2885662  0.2979177 0.46534653  1.0000000   FALSE
#> 18 -0.2500000  0.3216561 0.52475248  1.0000000   FALSE
#> 19 -0.2911392  0.2585670 0.12871287  0.9900990   FALSE
#> 20 -0.3357143  0.3543860 0.80198020  1.0000000   FALSE
#> 21 -0.3416537  0.3750000 0.20792079  0.9900990   FALSE
#> 22 -0.3453355  0.3333333 0.52475248  1.0000000   FALSE
#> 23 -0.2324723  0.3308271 0.59405941  1.0000000   FALSE
#> 24 -0.2238806  0.4440559 0.13861386  0.9900990   FALSE
#> 25 -0.2923077  0.3902848 0.44554455  1.0000000   FALSE
#> 26 -0.2087912  0.4170854 0.11881188  0.9900990   FALSE
#> 27 -0.3563636  0.2890365 0.01980198  0.1881188   FALSE
#> 28 -0.2226402  0.3563579 0.99009901  1.0000000   FALSE
#> 29 -0.3491311  0.2613241 0.69306931  1.0000000   FALSE
#> 30 -0.2740741  0.3381555 0.74257426  1.0000000   FALSE
#> 31 -0.2629969  0.3025335 0.86138614  1.0000000   FALSE
#> 32 -0.2527472  0.3256151 0.10891089  0.9900990   FALSE
#> 33 -0.3303965  0.3283358 0.18811881  1.0000000   FALSE
#> 34 -0.2671233  0.3651877 0.96039604  1.0000000   FALSE
#> 35 -0.2531876  0.3196481 0.45544554  1.0000000   FALSE
#> 36 -0.2467190  0.3983051 0.30693069  1.0000000   FALSE
#> 37 -0.2130584  0.3502235 0.67326733  1.0000000   FALSE
#> 38 -0.2437886  0.3966102 0.10891089  0.9900990   FALSE
#> 39 -0.2382609  0.3498350 0.97029703  1.0000000   FALSE
#> 40 -0.3460803  0.2953216 0.89108911  1.0000000   FALSE
#> 41 -0.2693498  0.3362319 0.58415842  1.0000000   FALSE
#> 42 -0.2040201  0.3078261 0.53465347  1.0000000   FALSE
#> 43 -0.2593918  0.3160813 0.01980198  0.4950495   FALSE
#> 44 -0.2667877  0.3750000 0.47524752  1.0000000   FALSE
#> 45 -0.2635379  0.3019197 0.48514851  1.0000000   FALSE
#> 
#> $direction2
#>     Item1 Item2        gamma   gamma_pair        se       lower       upper
#> 1   Item2 Item1  0.300448430  0.319681456 0.1445592  0.01711755  0.58377931
#> 2   Item3 Item1 -0.157775255 -0.143187067 0.1605048 -0.47235893  0.15680842
#> 3   Item3 Item2 -0.229563270 -0.205542725 0.1562271 -0.53576281  0.07663627
#> 4   Item4 Item1 -0.140893471 -0.069478908 0.1563189 -0.44727290  0.16548596
#> 5   Item4 Item2 -0.263397948 -0.202469136 0.1494361 -0.55628727  0.02949137
#> 6   Item4 Item3 -0.042162162  0.020737327 0.1614684 -0.35863443  0.27431010
#> 7   Item5 Item1 -0.151515152 -0.151515152 0.1616437 -0.46833104  0.16530074
#> 8   Item5 Item2  0.038961039  0.057692308 0.1611175 -0.27682347  0.35474555
#> 9   Item5 Item3  0.200878156  0.200878156 0.1537484 -0.10046312  0.50221943
#> 10  Item5 Item4  0.020656136  0.020656136 0.1569318 -0.28692454  0.32823681
#> 11  Item6 Item1 -0.101851852 -0.031210986 0.1579840 -0.41149479  0.20779109
#> 12  Item6 Item2 -0.062713797  0.024937656 0.1573442 -0.37110274  0.24567514
#> 13  Item6 Item3 -0.208287895 -0.142857143 0.1535339 -0.50920875  0.09263296
#> 14  Item6 Item4  0.205816555  0.223609535 0.1502536 -0.08867500  0.50030811
#> 15  Item6 Item5 -0.200929152 -0.135678392 0.1506350 -0.49616836  0.09431005
#> 16  Item7 Item1 -0.116219668 -0.116219668 0.1659714 -0.44151768  0.20907835
#> 17  Item7 Item2 -0.194936709 -0.194936709 0.1632885 -0.51497633  0.12510292
#> 18  Item7 Item3  0.211469534  0.211469534 0.1597576 -0.10164965  0.52458872
#> 19  Item7 Item4  0.428246014  0.428246014 0.1340087  0.16559385  0.69089817
#> 20  Item7 Item5  0.218116806  0.218116806 0.1526861 -0.08114249  0.51737611
#> 21  Item7 Item6  0.092682927  0.092682927 0.1596508 -0.22022694  0.40559279
#> 22  Item8 Item1  0.236714976  0.236714976 0.1523104 -0.06180795  0.53523790
#> 23  Item8 Item2  0.219927096  0.219927096 0.1554129 -0.08467656  0.52453075
#> 24  Item8 Item3  0.023752969  0.023752969 0.1650663 -0.29977106  0.34727699
#> 25  Item8 Item4 -0.318595579 -0.318595579 0.1504600 -0.61349172 -0.02369943
#> 26  Item8 Item5  0.151960784  0.151960784 0.1559234 -0.15364344  0.45756501
#> 27  Item8 Item6 -0.021879022 -0.021879022 0.1610349 -0.33750169  0.29374365
#> 28  Item8 Item7 -0.222222222 -0.197916667 0.1610302 -0.53783558  0.09339114
#> 29  Item9 Item1  0.193026152  0.193026152 0.1589923 -0.11859311  0.50464542
#> 30  Item9 Item2  0.198589894  0.198589894 0.1594713 -0.11396802  0.51114781
#> 31  Item9 Item3 -0.008206331 -0.008206331 0.1676919 -0.33687651  0.32046385
#> 32  Item9 Item4 -0.070631970 -0.070631970 0.1628262 -0.38976550  0.24850156
#> 33  Item9 Item5 -0.255689424 -0.255689424 0.1569567 -0.56331899  0.05194014
#> 34  Item9 Item6  0.217503218  0.217503218 0.1563077 -0.08885421  0.52386065
#> 35  Item9 Item7 -0.001236094 -0.001236094 0.1664800 -0.32753090  0.32505872
#> 36  Item9 Item8  0.357058126  0.357058126 0.1457329  0.07142697  0.64268928
#> 37 Item10 Item1 -0.070904645 -0.039190898 0.1597779 -0.38406353  0.24225423
#> 38 Item10 Item2 -0.022754491  0.023839398 0.1604775 -0.33728468  0.29177570
#> 39 Item10 Item3  0.176079734  0.220689655 0.1560401 -0.12975315  0.48191262
#> 40 Item10 Item4 -0.088270859 -0.083839611 0.1560097 -0.39404435  0.21750263
#> 41 Item10 Item5 -0.028915663  0.013836478 0.1560887 -0.33484397  0.27701264
#> 42 Item10 Item6 -0.234932349 -0.234932349 0.1488034 -0.52658173  0.05671703
#> 43 Item10 Item7 -0.040000000  0.021276596 0.1595171 -0.35264769  0.27264769
#> 44 Item10 Item8 -0.047044632  0.044973545 0.1605876 -0.36179047  0.26770120
#> 45 Item10 Item9 -0.069047619  0.035761589 0.1606941 -0.38400232  0.24590708
#>     gamma_low gamma_high    p_gamma padj_gamma flagged
#> 1  -0.2777321  0.3198758 0.01980198  0.6534653   FALSE
#> 2  -0.2571429  0.3843537 0.88118812  1.0000000   FALSE
#> 3  -0.2605364  0.3693694 0.91089109  1.0000000   FALSE
#> 4  -0.2265372  0.3262195 0.72277228  1.0000000   FALSE
#> 5  -0.2141653  0.3376906 0.96039604  1.0000000   FALSE
#> 6  -0.2500000  0.3216561 0.52475248  1.0000000   FALSE
#> 7  -0.3109244  0.3537519 0.83168317  1.0000000   FALSE
#> 8  -0.1967865  0.4179104 0.55445545  1.0000000   FALSE
#> 9  -0.2911392  0.2585670 0.12871287  0.9900990   FALSE
#> 10 -0.2923077  0.3902848 0.44554455  1.0000000   FALSE
#> 11 -0.2680723  0.3297003 0.70297030  1.0000000   FALSE
#> 12 -0.3203540  0.2933912 0.62376238  1.0000000   FALSE
#> 13 -0.3357143  0.3543860 0.80198020  1.0000000   FALSE
#> 14 -0.2087912  0.4170854 0.11881188  0.9900990   FALSE
#> 15 -0.2629969  0.3025335 0.86138614  1.0000000   FALSE
#> 16 -0.2089041  0.2947559 0.87128713  1.0000000   FALSE
#> 17 -0.3580264  0.4016620 0.82178218  1.0000000   FALSE
#> 18 -0.3416537  0.3750000 0.20792079  0.9900990   FALSE
#> 19 -0.3563636  0.2890365 0.01980198  0.1881188   FALSE
#> 20 -0.2527472  0.3256151 0.10891089  0.9900990   FALSE
#> 21 -0.2467190  0.3983051 0.30693069  1.0000000   FALSE
#> 22 -0.3167702  0.3052264 0.16831683  0.9900990   FALSE
#> 23 -0.2684564  0.3580705 0.13861386  0.9900990   FALSE
#> 24 -0.3453355  0.3333333 0.52475248  1.0000000   FALSE
#> 25 -0.2226402  0.3563579 0.99009901  1.0000000   FALSE
#> 26 -0.3303965  0.3283358 0.18811881  1.0000000   FALSE
#> 27 -0.2130584  0.3502235 0.67326733  1.0000000   FALSE
#> 28 -0.3460803  0.2953216 0.89108911  1.0000000   FALSE
#> 29 -0.2711268  0.3800000 0.14851485  0.9900990   FALSE
#> 30 -0.3421927  0.2688498 0.13861386  0.9900990   FALSE
#> 31 -0.2324723  0.3308271 0.59405941  1.0000000   FALSE
#> 32 -0.3491311  0.2613241 0.69306931  1.0000000   FALSE
#> 33 -0.2671233  0.3651877 0.96039604  1.0000000   FALSE
#> 34 -0.2437886  0.3966102 0.10891089  0.9900990   FALSE
#> 35 -0.2693498  0.3362319 0.58415842  1.0000000   FALSE
#> 36 -0.2593918  0.3160813 0.01980198  0.4950495   FALSE
#> 37 -0.2714777  0.2960000 0.66336634  1.0000000   FALSE
#> 38 -0.2885662  0.2979177 0.46534653  1.0000000   FALSE
#> 39 -0.2238806  0.4440559 0.13861386  0.9900990   FALSE
#> 40 -0.2740741  0.3381555 0.74257426  1.0000000   FALSE
#> 41 -0.2531876  0.3196481 0.45544554  1.0000000   FALSE
#> 42 -0.2382609  0.3498350 0.97029703  1.0000000   FALSE
#> 43 -0.2040201  0.3078261 0.53465347  1.0000000   FALSE
#> 44 -0.2667877  0.3750000 0.47524752  1.0000000   FALSE
#> 45 -0.2635379  0.3019197 0.48514851  1.0000000   FALSE
#> 
# }
```
