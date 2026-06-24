# Simulation-Based Partial Gamma DIF Cutoff Determination

Uses parametric bootstrap simulation to determine appropriate cutoff
values for partial gamma DIF analysis via
[`partgam_DIF`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html). Under a
correctly fitting Rasch model where the DIF variable is unrelated to
item responses (i.e., no true DIF), this function generates the expected
distribution of absolute partial gamma values per item, providing
empirical critical values.

## Usage

``` r
RMdifGammaCutoff(
  data,
  dif_var,
  iterations = 250,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.99
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Only complete cases (rows
  without any `NA`) are used.

- dif_var:

  A vector (factor, character, or integer) defining group membership for
  DIF analysis. Must have the same length as `nrow(data)`. The actual
  group labels are used to determine the number of groups and their
  relative sizes; during simulation, respondents are randomly assigned
  to groups with the same proportions, so there is no true DIF by
  construction.

- iterations:

  Integer. Number of simulation iterations (default 250).

- parallel:

  Logical. Use parallel processing via `mirai` if available (default
  `TRUE`).

- n_cores:

  Integer or `NULL`. Number of parallel workers. When `NULL`,
  `getOption("mc.cores")` is checked first. If neither is set and
  `parallel = TRUE`, a warning is issued and execution falls back to
  sequential (single core) processing.

- verbose:

  Logical. Show a progress bar (default `FALSE`).

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- cutoff_method:

  Character string specifying how cutoff intervals are computed. Either
  `"hdci"` (default) for the Highest Density Interval via
  [`ggdist::hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html),
  or `"quantile"` for the 2.5th/97.5th percentiles via
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html).

- hdci_width:

  Numeric. Width of the HDCI when `cutoff_method = "hdci"`. Default is
  `0.99` (99\\ `cutoff_method = "quantile"`.

## Value

A list with components:

- `results`:

  data.frame with columns `iteration`, `Item`, and `gamma` (one row per
  item per successful iteration).

- `item_cutoffs`:

  data.frame with per-item cutoff summaries: `Item`, `gamma_low`,
  `gamma_high`. Bounds are computed using the method specified by
  `cutoff_method`.

- `actual_iterations`:

  Number of successful iterations.

- `sample_n`:

  Number of complete cases used.

- `sample_summary`:

  Summary statistics of estimated person parameters.

- `item_names`:

  Character vector of item names from data.

- `dif_group_sizes`:

  Named integer vector of group sizes used in the simulation (matches
  proportions in the observed `dif_var`).

- `cutoff_method`:

  The method used to compute cutoffs (`"hdci"` or `"quantile"`).

- `hdci_width`:

  The HDCI width used (only meaningful when `cutoff_method = "hdci"`).

## Details

For each simulation iteration the function:

1.  Resamples person parameters (thetas) with replacement from the WLE
    person locations.

2.  Simulates item response data under a Rasch model (dichotomous via
    [`psychotools::rrm()`](https://rdrr.io/pkg/psychotools/man/rrm.html)
    or polytomous via an internal partial credit simulator).

3.  Creates a random DIF variable by sampling group labels with the same
    proportions as the observed `dif_var`, so there is **no true DIF**
    by construction.

4.  Computes partial gamma DIF statistics via
    [`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).

The distribution of partial gamma values across iterations provides
empirical critical values per item. Values from real data that fall
outside these bounds suggest DIF that exceeds what would be expected by
chance under a correctly fitting Rasch model. Failed iterations (e.g.,
due to convergence issues or degenerate data) are silently discarded.

The generating model uses CML item thresholds via
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category PCM) and WLE person locations,
consistent with the rest of the package; responses are simulated with
[`psychotools::rrm()`](https://rdrr.io/pkg/psychotools/man/rrm.html)
(dichotomous) or an internal partial credit score simulator
(polytomous).

Parallel processing is provided by the `mirai` package (optional).
Install it with `install.packages("mirai")` to enable parallelisation.

The `iarm` package must be installed (it is in Suggests, not Imports).

## References

Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
(1998). Differential item functioning in the Danish translation of the
SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189–1202.
[doi:10.1016/S0895-4356(98)00111-5](https://doi.org/10.1016/S0895-4356%2898%2900111-5)

Henninger, M., Radek, J., Debelak, R., & Strobl, C. (2025). Partial
credit trees meet the partial gamma coefficient for quantifying DIF and
DSF in polytomous items. *Behaviormetrika, 52*, 221–257.
[doi:10.1007/s41237-024-00252-3](https://doi.org/10.1007/s41237-024-00252-3)

## See also

[`partgam_DIF`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html)

## Examples

``` r
# \donttest{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)
dif_sex <- sample(c("male", "female"), 200, replace = TRUE)

# Run 100 iterations sequentially for a quick demo
cutoff_res <- RMdifGammaCutoff(sim_data, dif_var = dif_sex,
                                  iterations = 100, parallel = FALSE,
                                  seed = 42)
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1298 0.1745 0.4571  1.0000     -0.4718 0.2123
#> 2   Item2 random_dif -0.1160 0.1651 0.4821  1.0000     -0.4396 0.2075
#> 3   Item3 random_dif  0.1935 0.1584 0.2216  1.0000     -0.1168 0.5039
#> 4   Item4 random_dif -0.2428 0.1674 0.1471  1.0000     -0.5709 0.0854
#> 5   Item5 random_dif -0.2049 0.1658 0.2165  1.0000     -0.5297 0.1200
#> 6   Item6 random_dif  0.3467 0.1501 0.0208  0.2085      0.0526 0.6408
#> 7   Item7 random_dif  0.1003 0.1707 0.5566  1.0000     -0.2342 0.4349
#> 8   Item8 random_dif -0.0845 0.1640 0.6063  1.0000     -0.4060 0.2370
#> 9   Item9 random_dif  0.0487 0.1809 0.7878  1.0000     -0.3059 0.4033
#> 10 Item10 random_dif  0.0514 0.1721 0.7653  1.0000     -0.2859 0.3887
#>      Item        Var   gamma     se pvalue padj.BH  sig   lower   upper
#> 1   Item1 random_dif -0.2632 0.1589 0.0977  0.9766      -0.5746  0.0482
#> 2   Item2 random_dif  0.1141 0.1670 0.4945  1.0000      -0.2132  0.4413
#> 3   Item3 random_dif -0.6088 0.1114 0.0000  0.0000  *** -0.8272 -0.3905
#> 4   Item4 random_dif  0.2176 0.1655 0.1886  1.0000      -0.1068  0.5419
#> 5   Item5 random_dif  0.1780 0.1620 0.2719  1.0000      -0.1395  0.4956
#> 6   Item6 random_dif -0.1446 0.1702 0.3956  1.0000      -0.4783  0.1891
#> 7   Item7 random_dif  0.4585 0.1360 0.0007  0.0075   **  0.1920  0.7251
#> 8   Item8 random_dif  0.2355 0.1591 0.1387  1.0000      -0.0763  0.5473
#> 9   Item9 random_dif -0.1544 0.1572 0.3262  1.0000      -0.4626  0.1538
#> 10 Item10 random_dif  0.0246 0.1657 0.8821  1.0000      -0.3002  0.3493
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1791 0.1591 0.2603       1     -0.4909 0.1327
#> 2   Item2 random_dif  0.1227 0.1573 0.4355       1     -0.1857 0.4311
#> 3   Item3 random_dif -0.1212 0.1641 0.4601       1     -0.4428 0.2004
#> 4   Item4 random_dif -0.0574 0.1644 0.7269       1     -0.3797 0.2648
#> 5   Item5 random_dif -0.0183 0.1662 0.9123       1     -0.3441 0.3075
#> 6   Item6 random_dif -0.1720 0.1626 0.2904       1     -0.4907 0.1468
#> 7   Item7 random_dif  0.0758 0.1696 0.6549       1     -0.2567 0.4083
#> 8   Item8 random_dif  0.1516 0.1615 0.3478       1     -0.1649 0.4681
#> 9   Item9 random_dif  0.0087 0.1716 0.9596       1     -0.3276 0.3450
#> 10 Item10 random_dif  0.1867 0.1595 0.2420       1     -0.1260 0.4994
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2186 0.1643 0.1833  1.0000     -0.5405 0.1034
#> 2   Item2 random_dif  0.1705 0.1617 0.2919  1.0000     -0.1465 0.4875
#> 3   Item3 random_dif  0.1500 0.1796 0.4037  1.0000     -0.2021 0.5021
#> 4   Item4 random_dif  0.0780 0.1633 0.6329  1.0000     -0.2420 0.3980
#> 5   Item5 random_dif -0.0133 0.1704 0.9376  1.0000     -0.3474 0.3207
#> 6   Item6 random_dif -0.0247 0.1665 0.8819  1.0000     -0.3512 0.3017
#> 7   Item7 random_dif -0.2727 0.1552 0.0788  0.7881     -0.5769 0.0314
#> 8   Item8 random_dif -0.1111 0.1735 0.5220  1.0000     -0.4512 0.2290
#> 9   Item9 random_dif  0.0076 0.1685 0.9639  1.0000     -0.3226 0.3379
#> 10 Item10 random_dif  0.2570 0.1652 0.1197  1.0000     -0.0667 0.5808
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0109 0.1650 0.9476  1.0000     -0.3125  0.3342
#> 2   Item2 random_dif -0.0767 0.1735 0.6586  1.0000     -0.4168  0.2635
#> 3   Item3 random_dif  0.2355 0.1638 0.1506  1.0000     -0.0856  0.5566
#> 4   Item4 random_dif  0.0092 0.1637 0.9553  1.0000     -0.3117  0.3301
#> 5   Item5 random_dif -0.3191 0.1593 0.0452  0.4519     -0.6315 -0.0068
#> 6   Item6 random_dif -0.1014 0.1758 0.5642  1.0000     -0.4459  0.2432
#> 7   Item7 random_dif  0.1948 0.1626 0.2308  1.0000     -0.1238  0.5135
#> 8   Item8 random_dif -0.0141 0.1709 0.9341  1.0000     -0.3491  0.3208
#> 9   Item9 random_dif  0.1516 0.1629 0.3519  1.0000     -0.1676  0.4708
#> 10 Item10 random_dif -0.1119 0.1751 0.5227  1.0000     -0.4550  0.2312
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0397 0.1728 0.8183  1.0000     -0.3784  0.2990
#> 2   Item2 random_dif  0.0870 0.1832 0.6346  1.0000     -0.2719  0.4460
#> 3   Item3 random_dif  0.0864 0.1627 0.5954  1.0000     -0.2325  0.4052
#> 4   Item4 random_dif -0.3638 0.1500 0.0153  0.1531     -0.6578 -0.0697
#> 5   Item5 random_dif  0.0295 0.1770 0.8675  1.0000     -0.3173  0.3763
#> 6   Item6 random_dif -0.1090 0.1623 0.5018  1.0000     -0.4272  0.2092
#> 7   Item7 random_dif  0.1056 0.1562 0.4992  1.0000     -0.2006  0.4118
#> 8   Item8 random_dif  0.1745 0.1549 0.2600  1.0000     -0.1291  0.4781
#> 9   Item9 random_dif -0.0664 0.1676 0.6918  1.0000     -0.3949  0.2620
#> 10 Item10 random_dif  0.0981 0.1715 0.5673  1.0000     -0.2380  0.4343
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0493 0.1706 0.7726  1.0000     -0.3838 0.2851
#> 2   Item2 random_dif  0.2491 0.1606 0.1210  1.0000     -0.0658 0.5639
#> 3   Item3 random_dif -0.1314 0.1678 0.4336  1.0000     -0.4602 0.1975
#> 4   Item4 random_dif  0.2623 0.1518 0.0840  0.8395     -0.0352 0.5597
#> 5   Item5 random_dif  0.0018 0.1666 0.9914  1.0000     -0.3247 0.3283
#> 6   Item6 random_dif  0.0037 0.1727 0.9828  1.0000     -0.3347 0.3421
#> 7   Item7 random_dif -0.1965 0.1579 0.2133  1.0000     -0.5059 0.1129
#> 8   Item8 random_dif -0.0600 0.1709 0.7256  1.0000     -0.3949 0.2749
#> 9   Item9 random_dif -0.0551 0.1629 0.7354  1.0000     -0.3744 0.2643
#> 10 Item10 random_dif -0.0400 0.1722 0.8163  1.0000     -0.3775 0.2975
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0669 0.1776 0.7063  1.0000     -0.4150 0.2811
#> 2   Item2 random_dif -0.0346 0.1690 0.8378  1.0000     -0.3659 0.2967
#> 3   Item3 random_dif  0.1196 0.1734 0.4906  1.0000     -0.2204 0.4595
#> 4   Item4 random_dif  0.0345 0.1716 0.8407  1.0000     -0.3018 0.3708
#> 5   Item5 random_dif -0.0510 0.1734 0.7685  1.0000     -0.3909 0.2888
#> 6   Item6 random_dif -0.2224 0.1628 0.1718  1.0000     -0.5414 0.0966
#> 7   Item7 random_dif -0.2036 0.1698 0.2305  1.0000     -0.5364 0.1292
#> 8   Item8 random_dif  0.3209 0.1598 0.0446  0.4460      0.0077 0.6341
#> 9   Item9 random_dif -0.1561 0.1684 0.3539  1.0000     -0.4862 0.1740
#> 10 Item10 random_dif  0.2918 0.1651 0.0772  0.7715     -0.0318 0.6154
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0444 0.1667 0.7898       1     -0.3713 0.2824
#> 2   Item2 random_dif  0.0709 0.1719 0.6802       1     -0.2661 0.4078
#> 3   Item3 random_dif  0.0982 0.1633 0.5474       1     -0.2218 0.4182
#> 4   Item4 random_dif  0.0860 0.1641 0.6000       1     -0.2356 0.4076
#> 5   Item5 random_dif -0.0731 0.1670 0.6615       1     -0.4004 0.2542
#> 6   Item6 random_dif  0.0315 0.1651 0.8487       1     -0.2920 0.3550
#> 7   Item7 random_dif  0.0750 0.1698 0.6588       1     -0.2578 0.4077
#> 8   Item8 random_dif  0.0107 0.1814 0.9531       1     -0.3449 0.3663
#> 9   Item9 random_dif -0.1022 0.1692 0.5457       1     -0.4338 0.2293
#> 10 Item10 random_dif -0.1669 0.1708 0.3285       1     -0.5018 0.1679
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1528 0.1639 0.3512       1     -0.4740 0.1684
#> 2   Item2 random_dif -0.0875 0.1707 0.6080       1     -0.4220 0.2470
#> 3   Item3 random_dif -0.0132 0.1667 0.9367       1     -0.3399 0.3134
#> 4   Item4 random_dif  0.0729 0.1681 0.6644       1     -0.2565 0.4023
#> 5   Item5 random_dif  0.1148 0.1713 0.5030       1     -0.2211 0.4506
#> 6   Item6 random_dif  0.0775 0.1696 0.6478       1     -0.2550 0.4099
#> 7   Item7 random_dif -0.0148 0.1654 0.9286       1     -0.3390 0.3094
#> 8   Item8 random_dif -0.1387 0.1681 0.4092       1     -0.4682 0.1907
#> 9   Item9 random_dif  0.1410 0.1580 0.3721       1     -0.1687 0.4507
#> 10 Item10 random_dif -0.0088 0.1725 0.9595       1     -0.3469 0.3294
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2023 0.1689 0.2310  1.0000     -0.5333 0.1287
#> 2   Item2 random_dif  0.3076 0.1532 0.0447  0.4473      0.0072 0.6079
#> 3   Item3 random_dif  0.2136 0.1559 0.1708  1.0000     -0.0921 0.5192
#> 4   Item4 random_dif  0.1180 0.1760 0.5025  1.0000     -0.2269 0.4630
#> 5   Item5 random_dif  0.0080 0.1703 0.9627  1.0000     -0.3258 0.3417
#> 6   Item6 random_dif -0.0980 0.1676 0.5589  1.0000     -0.4265 0.2305
#> 7   Item7 random_dif -0.1243 0.1638 0.4479  1.0000     -0.4453 0.1967
#> 8   Item8 random_dif -0.0435 0.1769 0.8059  1.0000     -0.3902 0.3033
#> 9   Item9 random_dif -0.2207 0.1637 0.1776  1.0000     -0.5417 0.1002
#> 10 Item10 random_dif -0.0018 0.1647 0.9913  1.0000     -0.3247 0.3211
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0469 0.1643 0.7752  1.0000     -0.3690 0.2752
#> 2   Item2 random_dif -0.1563 0.1617 0.3337  1.0000     -0.4732 0.1606
#> 3   Item3 random_dif -0.0063 0.1801 0.9721  1.0000     -0.3592 0.3466
#> 4   Item4 random_dif  0.0111 0.1654 0.9464  1.0000     -0.3130 0.3352
#> 5   Item5 random_dif  0.3797 0.1486 0.0106  0.1062      0.0884 0.6710
#> 6   Item6 random_dif -0.0237 0.1631 0.8846  1.0000     -0.3434 0.2961
#> 7   Item7 random_dif  0.0018 0.1642 0.9912  1.0000     -0.3201 0.3237
#> 8   Item8 random_dif -0.1357 0.1706 0.4264  1.0000     -0.4700 0.1987
#> 9   Item9 random_dif  0.0241 0.1768 0.8916  1.0000     -0.3224 0.3706
#> 10 Item10 random_dif -0.0471 0.1611 0.7700  1.0000     -0.3628 0.2686
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1604 0.1692 0.3428  1.0000     -0.1711 0.4920
#> 2   Item2 random_dif  0.0940 0.1730 0.5869  1.0000     -0.2451 0.4330
#> 3   Item3 random_dif  0.2805 0.1654 0.0900  0.8999     -0.0438 0.6047
#> 4   Item4 random_dif -0.1018 0.1663 0.5406  1.0000     -0.4277 0.2242
#> 5   Item5 random_dif -0.1681 0.1687 0.3189  1.0000     -0.4988 0.1625
#> 6   Item6 random_dif -0.2872 0.1531 0.0606  0.6059     -0.5872 0.0128
#> 7   Item7 random_dif -0.0352 0.1695 0.8354  1.0000     -0.3674 0.2969
#> 8   Item8 random_dif  0.0314 0.1669 0.8509  1.0000     -0.2957 0.3584
#> 9   Item9 random_dif  0.1479 0.1590 0.3522  1.0000     -0.1637 0.4595
#> 10 Item10 random_dif -0.0949 0.1700 0.5768  1.0000     -0.4282 0.2384
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0446 0.1696 0.7926       1     -0.3771 0.2879
#> 2   Item2 random_dif -0.0559 0.1645 0.7342       1     -0.3783 0.2666
#> 3   Item3 random_dif  0.1027 0.1594 0.5193       1     -0.2098 0.4152
#> 4   Item4 random_dif  0.1672 0.1595 0.2944       1     -0.1453 0.4798
#> 5   Item5 random_dif -0.1527 0.1613 0.3437       1     -0.4689 0.1634
#> 6   Item6 random_dif  0.2185 0.1713 0.2020       1     -0.1172 0.5542
#> 7   Item7 random_dif -0.1560 0.1644 0.3426       1     -0.4783 0.1662
#> 8   Item8 random_dif -0.2059 0.1604 0.1993       1     -0.5202 0.1085
#> 9   Item9 random_dif  0.1781 0.1718 0.2998       1     -0.1585 0.5147
#> 10 Item10 random_dif -0.0538 0.1687 0.7500       1     -0.3845 0.2769
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0517 0.1666 0.7563  1.0000     -0.3782 0.2748
#> 2   Item2 random_dif -0.0505 0.1734 0.7707  1.0000     -0.3904 0.2893
#> 3   Item3 random_dif -0.2444 0.1624 0.1323  1.0000     -0.5627 0.0738
#> 4   Item4 random_dif -0.1986 0.1663 0.2325  1.0000     -0.5245 0.1274
#> 5   Item5 random_dif  0.1610 0.1621 0.3206  1.0000     -0.1567 0.4786
#> 6   Item6 random_dif  0.0293 0.1686 0.8622  1.0000     -0.3012 0.3598
#> 7   Item7 random_dif  0.0866 0.1646 0.5989  1.0000     -0.2360 0.4092
#> 8   Item8 random_dif  0.2923 0.1568 0.0624  0.6235     -0.0151 0.5997
#> 9   Item9 random_dif -0.1261 0.1705 0.4594  1.0000     -0.4603 0.2080
#> 10 Item10 random_dif  0.0687 0.1712 0.6882  1.0000     -0.2669 0.4043
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1627 0.1678 0.3323  1.0000     -0.1662 0.4916
#> 2   Item2 random_dif -0.1601 0.1603 0.3180  1.0000     -0.4742 0.1541
#> 3   Item3 random_dif -0.0632 0.1701 0.7102  1.0000     -0.3966 0.2702
#> 4   Item4 random_dif  0.1417 0.1753 0.4186  1.0000     -0.2018 0.4853
#> 5   Item5 random_dif  0.1275 0.1637 0.4360  1.0000     -0.1933 0.4482
#> 6   Item6 random_dif  0.1581 0.1647 0.3371  1.0000     -0.1647 0.4809
#> 7   Item7 random_dif -0.1345 0.1793 0.4533  1.0000     -0.4858 0.2169
#> 8   Item8 random_dif -0.1368 0.1778 0.4415  1.0000     -0.4853 0.2116
#> 9   Item9 random_dif -0.2965 0.1551 0.0560  0.5599     -0.6006 0.0076
#> 10 Item10 random_dif  0.1856 0.1683 0.2701  1.0000     -0.1443 0.5155
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1352 0.1648 0.4118  1.0000     -0.4582 0.1877
#> 2   Item2 random_dif  0.0017 0.1744 0.9921  1.0000     -0.3400 0.3435
#> 3   Item3 random_dif  0.2869 0.1586 0.0705  0.7052     -0.0240 0.5978
#> 4   Item4 random_dif  0.0079 0.1642 0.9615  1.0000     -0.3139 0.3298
#> 5   Item5 random_dif  0.0078 0.1657 0.9625  1.0000     -0.3170 0.3326
#> 6   Item6 random_dif -0.0693 0.1781 0.6974  1.0000     -0.4184 0.2799
#> 7   Item7 random_dif  0.0865 0.1684 0.6076  1.0000     -0.2435 0.4165
#> 8   Item8 random_dif -0.0017 0.1763 0.9921  1.0000     -0.3472 0.3437
#> 9   Item9 random_dif  0.0600 0.1688 0.7224  1.0000     -0.2709 0.3908
#> 10 Item10 random_dif -0.2747 0.1631 0.0922  0.9224     -0.5944 0.0451
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3121 0.1532 0.0416  0.4164     -0.6125 -0.0118
#> 2   Item2 random_dif -0.1261 0.1613 0.4342  1.0000     -0.4423  0.1900
#> 3   Item3 random_dif  0.0982 0.1640 0.5490  1.0000     -0.2231  0.4196
#> 4   Item4 random_dif -0.0619 0.1653 0.7079  1.0000     -0.3858  0.2620
#> 5   Item5 random_dif  0.1792 0.1616 0.2673  1.0000     -0.1374  0.4959
#> 6   Item6 random_dif  0.1751 0.1675 0.2958  1.0000     -0.1531  0.5033
#> 7   Item7 random_dif -0.1889 0.1690 0.2638  1.0000     -0.5201  0.1424
#> 8   Item8 random_dif  0.1975 0.1700 0.2454  1.0000     -0.1357  0.5308
#> 9   Item9 random_dif  0.0659 0.1714 0.7007  1.0000     -0.2700  0.4017
#> 10 Item10 random_dif -0.0200 0.1660 0.9043  1.0000     -0.3454  0.3055
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0058 0.1765 0.9737       1     -0.3518 0.3401
#> 2   Item2 random_dif -0.1259 0.1733 0.4677       1     -0.4656 0.2138
#> 3   Item3 random_dif  0.0164 0.1654 0.9210       1     -0.3078 0.3406
#> 4   Item4 random_dif  0.1088 0.1648 0.5089       1     -0.2141 0.4318
#> 5   Item5 random_dif  0.0710 0.1748 0.6846       1     -0.2716 0.4137
#> 6   Item6 random_dif -0.0252 0.1679 0.8808       1     -0.3543 0.3040
#> 7   Item7 random_dif -0.1204 0.1697 0.4780       1     -0.4529 0.2121
#> 8   Item8 random_dif  0.1848 0.1675 0.2699       1     -0.1435 0.5130
#> 9   Item9 random_dif -0.1798 0.1704 0.2914       1     -0.5137 0.1542
#> 10 Item10 random_dif  0.0938 0.1745 0.5909       1     -0.2483 0.4359
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0151 0.1687 0.9286  1.0000     -0.3458 0.3156
#> 2   Item2 random_dif  0.0962 0.1605 0.5491  1.0000     -0.2184 0.4107
#> 3   Item3 random_dif -0.2655 0.1473 0.0714  0.7142     -0.5541 0.0231
#> 4   Item4 random_dif -0.0124 0.1632 0.9395  1.0000     -0.3322 0.3075
#> 5   Item5 random_dif  0.0626 0.1604 0.6963  1.0000     -0.2518 0.3770
#> 6   Item6 random_dif -0.3076 0.1617 0.0572  0.5716     -0.6245 0.0094
#> 7   Item7 random_dif -0.0115 0.1706 0.9465  1.0000     -0.3459 0.3229
#> 8   Item8 random_dif  0.2810 0.1545 0.0688  0.6883     -0.0217 0.5838
#> 9   Item9 random_dif  0.1411 0.1619 0.3836  1.0000     -0.1763 0.4584
#> 10 Item10 random_dif  0.0256 0.1715 0.8812  1.0000     -0.3106 0.3618
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.3344 0.1581 0.0344  0.3442      0.0245 0.6443
#> 2   Item2 random_dif -0.1655 0.1627 0.3092  1.0000     -0.4844 0.1534
#> 3   Item3 random_dif  0.0508 0.1638 0.7562  1.0000     -0.2701 0.3718
#> 4   Item4 random_dif -0.0954 0.1611 0.5538  1.0000     -0.4111 0.2203
#> 5   Item5 random_dif -0.0482 0.1708 0.7776  1.0000     -0.3830 0.2865
#> 6   Item6 random_dif -0.1879 0.1554 0.2267  1.0000     -0.4925 0.1167
#> 7   Item7 random_dif -0.0016 0.1664 0.9922  1.0000     -0.3278 0.3245
#> 8   Item8 random_dif -0.1362 0.1564 0.3838  1.0000     -0.4428 0.1703
#> 9   Item9 random_dif  0.1328 0.1664 0.4248  1.0000     -0.1933 0.4589
#> 10 Item10 random_dif  0.1688 0.1682 0.3154  1.0000     -0.1608 0.4985
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0633 0.1774 0.7213  1.0000     -0.2844  0.4110
#> 2   Item2 random_dif  0.0176 0.1706 0.9178  1.0000     -0.3168  0.3521
#> 3   Item3 random_dif -0.3641 0.1547 0.0186  0.1857     -0.6672 -0.0610
#> 4   Item4 random_dif -0.0423 0.1729 0.8070  1.0000     -0.3812  0.2967
#> 5   Item5 random_dif -0.2032 0.1591 0.2016  1.0000     -0.5150  0.1087
#> 6   Item6 random_dif  0.1458 0.1674 0.3836  1.0000     -0.1823  0.4739
#> 7   Item7 random_dif  0.1847 0.1641 0.2603  1.0000     -0.1369  0.5062
#> 8   Item8 random_dif  0.0550 0.1716 0.7486  1.0000     -0.2813  0.3912
#> 9   Item9 random_dif  0.0070 0.1687 0.9670  1.0000     -0.3237  0.3376
#> 10 Item10 random_dif  0.1641 0.1722 0.3404  1.0000     -0.1733  0.5016
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1693 0.1598 0.2895       1     -0.1440 0.4826
#> 2   Item2 random_dif -0.0968 0.1594 0.5437       1     -0.4091 0.2156
#> 3   Item3 random_dif -0.0850 0.1647 0.6060       1     -0.4078 0.2379
#> 4   Item4 random_dif  0.1517 0.1573 0.3349       1     -0.1566 0.4601
#> 5   Item5 random_dif  0.0245 0.1679 0.8838       1     -0.3046 0.3536
#> 6   Item6 random_dif -0.1828 0.1599 0.2529       1     -0.4961 0.1305
#> 7   Item7 random_dif  0.0000 0.1606 1.0000       1     -0.3148 0.3148
#> 8   Item8 random_dif -0.0046 0.1608 0.9773       1     -0.3198 0.3107
#> 9   Item9 random_dif -0.1748 0.1605 0.2761       1     -0.4892 0.1397
#> 10 Item10 random_dif  0.1975 0.1572 0.2090       1     -0.1106 0.5055
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2508 0.1687 0.1372       1     -0.5815 0.0799
#> 2   Item2 random_dif  0.1124 0.1689 0.5059       1     -0.2187 0.4435
#> 3   Item3 random_dif  0.0247 0.1725 0.8862       1     -0.3135 0.3629
#> 4   Item4 random_dif  0.0163 0.1682 0.9226       1     -0.3133 0.3460
#> 5   Item5 random_dif  0.0330 0.1759 0.8513       1     -0.3118 0.3777
#> 6   Item6 random_dif -0.0548 0.1765 0.7561       1     -0.4009 0.2912
#> 7   Item7 random_dif -0.0396 0.1717 0.8174       1     -0.3761 0.2968
#> 8   Item8 random_dif  0.0311 0.1642 0.8499       1     -0.2908 0.3530
#> 9   Item9 random_dif  0.1711 0.1728 0.3223       1     -0.1677 0.5098
#> 10 Item10 random_dif -0.0667 0.1825 0.7149       1     -0.4243 0.2910
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1638 0.1812 0.3661  1.0000     -0.5190 0.1914
#> 2   Item2 random_dif -0.1082 0.1720 0.5290  1.0000     -0.4453 0.2288
#> 3   Item3 random_dif  0.1572 0.1781 0.3775  1.0000     -0.1919 0.5062
#> 4   Item4 random_dif  0.0184 0.1774 0.9176  1.0000     -0.3293 0.3661
#> 5   Item5 random_dif  0.2751 0.1607 0.0869  0.8694     -0.0399 0.5901
#> 6   Item6 random_dif -0.1603 0.1691 0.3432  1.0000     -0.4916 0.1711
#> 7   Item7 random_dif  0.0061 0.1684 0.9709  1.0000     -0.3239 0.3361
#> 8   Item8 random_dif -0.0762 0.1742 0.6621  1.0000     -0.4177 0.2654
#> 9   Item9 random_dif  0.1916 0.1615 0.2355  1.0000     -0.1250 0.5082
#> 10 Item10 random_dif -0.1622 0.1695 0.3386  1.0000     -0.4945 0.1700
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0512 0.1790 0.7749  1.0000     -0.2997  0.4020
#> 2   Item2 random_dif -0.0051 0.1605 0.9749  1.0000     -0.3196  0.3095
#> 3   Item3 random_dif -0.3301 0.1478 0.0256  0.2557     -0.6199 -0.0403
#> 4   Item4 random_dif -0.1248 0.1622 0.4417  1.0000     -0.4427  0.1931
#> 5   Item5 random_dif -0.0113 0.1622 0.9446  1.0000     -0.3291  0.3065
#> 6   Item6 random_dif  0.0645 0.1777 0.7166  1.0000     -0.2838  0.4128
#> 7   Item7 random_dif  0.1950 0.1611 0.2262  1.0000     -0.1208  0.5109
#> 8   Item8 random_dif  0.0050 0.1618 0.9753  1.0000     -0.3121  0.3221
#> 9   Item9 random_dif  0.1863 0.1690 0.2704  1.0000     -0.1450  0.5176
#> 10 Item10 random_dif  0.0360 0.1676 0.8301  1.0000     -0.2925  0.3645
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0232 0.1792 0.8972       1     -0.3282 0.3745
#> 2   Item2 random_dif -0.0840 0.1726 0.6264       1     -0.4223 0.2543
#> 3   Item3 random_dif  0.0288 0.1718 0.8668       1     -0.3078 0.3654
#> 4   Item4 random_dif  0.0146 0.1727 0.9326       1     -0.3239 0.3531
#> 5   Item5 random_dif  0.0561 0.1658 0.7352       1     -0.2689 0.3811
#> 6   Item6 random_dif -0.0797 0.1586 0.6152       1     -0.3905 0.2311
#> 7   Item7 random_dif  0.2272 0.1658 0.1706       1     -0.0977 0.5521
#> 8   Item8 random_dif  0.1326 0.1614 0.4114       1     -0.1838 0.4489
#> 9   Item9 random_dif -0.1230 0.1676 0.4631       1     -0.4514 0.2055
#> 10 Item10 random_dif -0.1767 0.1604 0.2707       1     -0.4911 0.1377
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1788 0.1573 0.2554  1.0000     -0.4871  0.1294
#> 2   Item2 random_dif  0.1204 0.1685 0.4749  1.0000     -0.2099  0.4507
#> 3   Item3 random_dif  0.4043 0.1429 0.0047  0.0467   *  0.1242  0.6844
#> 4   Item4 random_dif -0.1005 0.1667 0.5466  1.0000     -0.4273  0.2263
#> 5   Item5 random_dif -0.0093 0.1613 0.9538  1.0000     -0.3255  0.3068
#> 6   Item6 random_dif  0.1129 0.1666 0.4979  1.0000     -0.2136  0.4394
#> 7   Item7 random_dif  0.1401 0.1615 0.3857  1.0000     -0.1764  0.4566
#> 8   Item8 random_dif -0.2045 0.1528 0.1807  1.0000     -0.5041  0.0950
#> 9   Item9 random_dif -0.3056 0.1475 0.0383  0.3832     -0.5947 -0.0164
#> 10 Item10 random_dif  0.0556 0.1685 0.7416  1.0000     -0.2747  0.3858
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1144 0.1653 0.4887       1     -0.4383 0.2095
#> 2   Item2 random_dif  0.2610 0.1688 0.1220       1     -0.0698 0.5918
#> 3   Item3 random_dif -0.1993 0.1629 0.2210       1     -0.5186 0.1199
#> 4   Item4 random_dif -0.1786 0.1598 0.2638       1     -0.4918 0.1346
#> 5   Item5 random_dif -0.0137 0.1836 0.9405       1     -0.3736 0.3462
#> 6   Item6 random_dif  0.1804 0.1638 0.2708       1     -0.1406 0.5014
#> 7   Item7 random_dif  0.0909 0.1750 0.6034       1     -0.2520 0.4339
#> 8   Item8 random_dif -0.0270 0.1757 0.8778       1     -0.3714 0.3174
#> 9   Item9 random_dif -0.1471 0.1698 0.3862       1     -0.4798 0.1856
#> 10 Item10 random_dif  0.1696 0.1716 0.3231       1     -0.1668 0.5060
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2189 0.1666 0.1888  1.0000     -0.5453 0.1076
#> 2   Item2 random_dif  0.1915 0.1657 0.2477  1.0000     -0.1332 0.5162
#> 3   Item3 random_dif  0.1155 0.1611 0.4733  1.0000     -0.2002 0.4312
#> 4   Item4 random_dif  0.0386 0.1654 0.8155  1.0000     -0.2855 0.3627
#> 5   Item5 random_dif  0.0175 0.1668 0.9162  1.0000     -0.3093 0.3444
#> 6   Item6 random_dif -0.0782 0.1690 0.6434  1.0000     -0.4095 0.2530
#> 7   Item7 random_dif  0.1325 0.1643 0.4203  1.0000     -0.1897 0.4546
#> 8   Item8 random_dif -0.2712 0.1522 0.0747  0.7467     -0.5695 0.0270
#> 9   Item9 random_dif  0.1639 0.1606 0.3074  1.0000     -0.1508 0.4786
#> 10 Item10 random_dif -0.1280 0.1714 0.4553  1.0000     -0.4640 0.2080
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1326 0.1684 0.4312  1.0000     -0.1975 0.4627
#> 2   Item2 random_dif -0.2870 0.1619 0.0763  0.7626     -0.6043 0.0303
#> 3   Item3 random_dif  0.1965 0.1611 0.2225  1.0000     -0.1192 0.5122
#> 4   Item4 random_dif  0.0207 0.1644 0.8998  1.0000     -0.3014 0.3428
#> 5   Item5 random_dif -0.1640 0.1621 0.3117  1.0000     -0.4818 0.1538
#> 6   Item6 random_dif  0.0310 0.1679 0.8537  1.0000     -0.2982 0.3601
#> 7   Item7 random_dif  0.0684 0.1637 0.6759  1.0000     -0.2524 0.3893
#> 8   Item8 random_dif -0.2743 0.1551 0.0769  0.7693     -0.5783 0.0297
#> 9   Item9 random_dif  0.0676 0.1640 0.6801  1.0000     -0.2538 0.3890
#> 10 Item10 random_dif  0.2305 0.1696 0.1743  1.0000     -0.1020 0.5630
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0146 0.1714 0.9320  1.0000     -0.3506  0.3213
#> 2   Item2 random_dif  0.0855 0.1723 0.6197  1.0000     -0.2522  0.4233
#> 3   Item3 random_dif -0.3248 0.1532 0.0340  0.3399     -0.6251 -0.0245
#> 4   Item4 random_dif -0.0627 0.1704 0.7130  1.0000     -0.3967  0.2713
#> 5   Item5 random_dif -0.1217 0.1679 0.4686  1.0000     -0.4508  0.2074
#> 6   Item6 random_dif -0.1918 0.1677 0.2528  1.0000     -0.5205  0.1369
#> 7   Item7 random_dif  0.2727 0.1542 0.0770  0.7699     -0.0295  0.5750
#> 8   Item8 random_dif -0.0616 0.1694 0.7163  1.0000     -0.3936  0.2705
#> 9   Item9 random_dif  0.1070 0.1699 0.5290  1.0000     -0.2261  0.4400
#> 10 Item10 random_dif  0.3109 0.1623 0.0554  0.5536     -0.0071  0.6290
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.3431 0.1594 0.0313  0.3134      0.0307 0.6554
#> 2   Item2 random_dif -0.0243 0.1765 0.8905  1.0000     -0.3702 0.3216
#> 3   Item3 random_dif -0.1727 0.1726 0.3172  1.0000     -0.5110 0.1657
#> 4   Item4 random_dif  0.1680 0.1740 0.3342  1.0000     -0.1730 0.5089
#> 5   Item5 random_dif -0.2087 0.1738 0.2298  1.0000     -0.5492 0.1319
#> 6   Item6 random_dif -0.1721 0.1784 0.3345  1.0000     -0.5217 0.1774
#> 7   Item7 random_dif  0.0137 0.1653 0.9337  1.0000     -0.3102 0.3377
#> 8   Item8 random_dif -0.0833 0.1780 0.6396  1.0000     -0.4321 0.2655
#> 9   Item9 random_dif  0.1429 0.1752 0.4147  1.0000     -0.2004 0.4862
#> 10 Item10 random_dif -0.0187 0.1733 0.9139  1.0000     -0.3583 0.3209
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1254 0.1668 0.4522  1.0000     -0.4524 0.2016
#> 2   Item2 random_dif -0.2970 0.1610 0.0652  0.6517     -0.6126 0.0187
#> 3   Item3 random_dif  0.1144 0.1683 0.4968  1.0000     -0.2155 0.4443
#> 4   Item4 random_dif  0.0859 0.1710 0.6153  1.0000     -0.2492 0.4210
#> 5   Item5 random_dif  0.2565 0.1673 0.1251  1.0000     -0.0713 0.5844
#> 6   Item6 random_dif  0.0341 0.1706 0.8414  1.0000     -0.3002 0.3685
#> 7   Item7 random_dif  0.0017 0.1714 0.9919  1.0000     -0.3342 0.3377
#> 8   Item8 random_dif  0.0523 0.1781 0.7692  1.0000     -0.2968 0.4014
#> 9   Item9 random_dif -0.0052 0.1775 0.9767  1.0000     -0.3532 0.3428
#> 10 Item10 random_dif -0.0933 0.1703 0.5837  1.0000     -0.4272 0.2405
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0578 0.1682 0.7314  1.0000     -0.2720  0.3875
#> 2   Item2 random_dif  0.3384 0.1617 0.0364  0.3641      0.0214  0.6554
#> 3   Item3 random_dif -0.0671 0.1706 0.6942  1.0000     -0.4015  0.2673
#> 4   Item4 random_dif -0.0934 0.1661 0.5740  1.0000     -0.4190  0.2322
#> 5   Item5 random_dif -0.2526 0.1584 0.1107  1.0000     -0.5631  0.0578
#> 6   Item6 random_dif  0.2241 0.1638 0.1713  1.0000     -0.0969  0.5451
#> 7   Item7 random_dif  0.1895 0.1696 0.2640  1.0000     -0.1430  0.5219
#> 8   Item8 random_dif  0.0303 0.1721 0.8603  1.0000     -0.3071  0.3677
#> 9   Item9 random_dif  0.0060 0.1662 0.9713  1.0000     -0.3197  0.3316
#> 10 Item10 random_dif -0.3901 0.1444 0.0069  0.0692   . -0.6732 -0.1070
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0107 0.1631 0.9476  1.0000     -0.3090 0.3305
#> 2   Item2 random_dif  0.2844 0.1621 0.0793  0.7933     -0.0333 0.6020
#> 3   Item3 random_dif  0.0983 0.1669 0.5558  1.0000     -0.2288 0.4255
#> 4   Item4 random_dif  0.0261 0.1684 0.8767  1.0000     -0.3039 0.3561
#> 5   Item5 random_dif -0.1079 0.1612 0.5031  1.0000     -0.4238 0.2080
#> 6   Item6 random_dif -0.1134 0.1815 0.5320  1.0000     -0.4692 0.2423
#> 7   Item7 random_dif -0.1061 0.1580 0.5018  1.0000     -0.4157 0.2035
#> 8   Item8 random_dif  0.0766 0.1611 0.6341  1.0000     -0.2390 0.3923
#> 9   Item9 random_dif -0.0438 0.1641 0.7896  1.0000     -0.3655 0.2779
#> 10 Item10 random_dif -0.1241 0.1674 0.4585  1.0000     -0.4522 0.2040
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0897 0.1600 0.5749  1.0000     -0.2238  0.4033
#> 2   Item2 random_dif -0.2534 0.1626 0.1192  1.0000     -0.5722  0.0653
#> 3   Item3 random_dif -0.0929 0.1586 0.5578  1.0000     -0.4038  0.2179
#> 4   Item4 random_dif  0.3058 0.1523 0.0447  0.4466      0.0073  0.6043
#> 5   Item5 random_dif -0.0904 0.1601 0.5724  1.0000     -0.4041  0.2234
#> 6   Item6 random_dif  0.3443 0.1463 0.0186  0.1862      0.0575  0.6311
#> 7   Item7 random_dif -0.0948 0.1623 0.5590  1.0000     -0.4129  0.2233
#> 8   Item8 random_dif -0.0883 0.1718 0.6074  1.0000     -0.4251  0.2485
#> 9   Item9 random_dif -0.3101 0.1453 0.0328  0.3280     -0.5949 -0.0254
#> 10 Item10 random_dif  0.1809 0.1656 0.2747  1.0000     -0.1437  0.5055
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1607 0.1640 0.3271       1     -0.4820 0.1607
#> 2   Item2 random_dif -0.2442 0.1558 0.1169       1     -0.5495 0.0611
#> 3   Item3 random_dif  0.0571 0.1600 0.7212       1     -0.2565 0.3707
#> 4   Item4 random_dif  0.1331 0.1655 0.4214       1     -0.1913 0.4575
#> 5   Item5 random_dif  0.0841 0.1747 0.6300       1     -0.2582 0.4265
#> 6   Item6 random_dif -0.1140 0.1612 0.4793       1     -0.4299 0.2019
#> 7   Item7 random_dif  0.0420 0.1691 0.8041       1     -0.2895 0.3734
#> 8   Item8 random_dif  0.0989 0.1772 0.5768       1     -0.2485 0.4463
#> 9   Item9 random_dif  0.1619 0.1674 0.3332       1     -0.1661 0.4900
#> 10 Item10 random_dif -0.0273 0.1689 0.8714       1     -0.3583 0.3036
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0073 0.1696 0.9656  1.0000     -0.3398 0.3252
#> 2   Item2 random_dif  0.0588 0.1633 0.7186  1.0000     -0.2612 0.3788
#> 3   Item3 random_dif  0.2027 0.1669 0.2245  1.0000     -0.1244 0.5298
#> 4   Item4 random_dif -0.1936 0.1550 0.2116  1.0000     -0.4973 0.1101
#> 5   Item5 random_dif -0.1897 0.1660 0.2533  1.0000     -0.5151 0.1357
#> 6   Item6 random_dif  0.3245 0.1509 0.0315  0.3152      0.0287 0.6203
#> 7   Item7 random_dif  0.0305 0.1726 0.8596  1.0000     -0.3078 0.3688
#> 8   Item8 random_dif -0.0203 0.1625 0.9008  1.0000     -0.3388 0.2983
#> 9   Item9 random_dif -0.0312 0.1629 0.8481  1.0000     -0.3505 0.2881
#> 10 Item10 random_dif -0.1712 0.1608 0.2870  1.0000     -0.4864 0.1440
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1176 0.1637 0.4722   1.000     -0.4384 0.2031
#> 2   Item2 random_dif -0.1776 0.1670 0.2877   1.000     -0.5050 0.1498
#> 3   Item3 random_dif -0.1490 0.1615 0.3561   1.000     -0.4655 0.1675
#> 4   Item4 random_dif  0.1544 0.1542 0.3166   1.000     -0.1478 0.4566
#> 5   Item5 random_dif -0.1910 0.1684 0.2566   1.000     -0.5211 0.1390
#> 6   Item6 random_dif  0.2679 0.1511 0.0761   0.761     -0.0281 0.5640
#> 7   Item7 random_dif -0.0381 0.1625 0.8148   1.000     -0.3565 0.2803
#> 8   Item8 random_dif  0.0250 0.1622 0.8773   1.000     -0.2928 0.3429
#> 9   Item9 random_dif  0.1913 0.1570 0.2232   1.000     -0.1165 0.4991
#> 10 Item10 random_dif -0.0196 0.1676 0.9069   1.000     -0.3481 0.3089
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1244 0.1666 0.4555       1     -0.2023 0.4510
#> 2   Item2 random_dif -0.2133 0.1640 0.1934       1     -0.5347 0.1081
#> 3   Item3 random_dif -0.0093 0.1751 0.9578       1     -0.3525 0.3340
#> 4   Item4 random_dif  0.0368 0.1711 0.8298       1     -0.2986 0.3721
#> 5   Item5 random_dif  0.1541 0.1607 0.3376       1     -0.1609 0.4691
#> 6   Item6 random_dif -0.1652 0.1690 0.3282       1     -0.4963 0.1660
#> 7   Item7 random_dif  0.0657 0.1676 0.6949       1     -0.2628 0.3943
#> 8   Item8 random_dif -0.0017 0.1682 0.9918       1     -0.3313 0.3279
#> 9   Item9 random_dif  0.1664 0.1631 0.3076       1     -0.1532 0.4860
#> 10 Item10 random_dif -0.2015 0.1699 0.2357       1     -0.5345 0.1315
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1129 0.1774 0.5246       1     -0.4607 0.2349
#> 2   Item2 random_dif -0.1179 0.1605 0.4623       1     -0.4324 0.1965
#> 3   Item3 random_dif -0.0542 0.1676 0.7463       1     -0.3826 0.2742
#> 4   Item4 random_dif  0.2036 0.1655 0.2186       1     -0.1208 0.5281
#> 5   Item5 random_dif -0.0316 0.1640 0.8473       1     -0.3530 0.2898
#> 6   Item6 random_dif  0.1135 0.1637 0.4882       1     -0.2074 0.4344
#> 7   Item7 random_dif -0.0441 0.1676 0.7923       1     -0.3726 0.2843
#> 8   Item8 random_dif -0.0273 0.1672 0.8702       1     -0.3550 0.3004
#> 9   Item9 random_dif  0.0267 0.1670 0.8728       1     -0.3006 0.3541
#> 10 Item10 random_dif  0.0389 0.1720 0.8210       1     -0.2982 0.3760
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1092 0.1643 0.5064       1     -0.2129 0.4312
#> 2   Item2 random_dif  0.0669 0.1655 0.6861       1     -0.2575 0.3913
#> 3   Item3 random_dif  0.0937 0.1717 0.5854       1     -0.2429 0.4303
#> 4   Item4 random_dif -0.1812 0.1637 0.2684       1     -0.5020 0.1396
#> 5   Item5 random_dif -0.1856 0.1653 0.2616       1     -0.5096 0.1384
#> 6   Item6 random_dif  0.2567 0.1622 0.1136       1     -0.0613 0.5747
#> 7   Item7 random_dif  0.0073 0.1691 0.9654       1     -0.3240 0.3387
#> 8   Item8 random_dif -0.2444 0.1633 0.1346       1     -0.5645 0.0758
#> 9   Item9 random_dif  0.0588 0.1654 0.7220       1     -0.2653 0.3829
#> 10 Item10 random_dif  0.0196 0.1661 0.9060       1     -0.3059 0.3451
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0970 0.1633 0.5525  1.0000     -0.2230  0.4170
#> 2   Item2 random_dif -0.0185 0.1637 0.9101  1.0000     -0.3394  0.3024
#> 3   Item3 random_dif  0.1184 0.1659 0.4753  1.0000     -0.2067  0.4435
#> 4   Item4 random_dif  0.2383 0.1570 0.1291  1.0000     -0.0694  0.5459
#> 5   Item5 random_dif -0.2751 0.1534 0.0729  0.7292     -0.5758  0.0256
#> 6   Item6 random_dif -0.4007 0.1463 0.0062  0.0616   . -0.6874 -0.1139
#> 7   Item7 random_dif  0.0162 0.1598 0.9193  1.0000     -0.2970  0.3294
#> 8   Item8 random_dif  0.0556 0.1705 0.7445  1.0000     -0.2785  0.3897
#> 9   Item9 random_dif  0.2884 0.1556 0.0638  0.6377     -0.0165  0.5934
#> 10 Item10 random_dif -0.0903 0.1660 0.5866  1.0000     -0.4156  0.2351
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1492 0.1778 0.4014       1     -0.1992 0.4976
#> 2   Item2 random_dif  0.1115 0.1630 0.4940       1     -0.2080 0.4309
#> 3   Item3 random_dif -0.1519 0.1674 0.3642       1     -0.4801 0.1762
#> 4   Item4 random_dif -0.1359 0.1579 0.3894       1     -0.4454 0.1736
#> 5   Item5 random_dif  0.1184 0.1635 0.4689       1     -0.2021 0.4389
#> 6   Item6 random_dif  0.0476 0.1708 0.7804       1     -0.2871 0.3823
#> 7   Item7 random_dif -0.1329 0.1684 0.4301       1     -0.4629 0.1972
#> 8   Item8 random_dif -0.0364 0.1655 0.8258       1     -0.3608 0.2880
#> 9   Item9 random_dif  0.0458 0.1665 0.7831       1     -0.2806 0.3722
#> 10 Item10 random_dif -0.0017 0.1734 0.9921       1     -0.3415 0.3381
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0753 0.1821 0.6794       1     -0.2817 0.4322
#> 2   Item2 random_dif -0.2292 0.1695 0.1762       1     -0.5613 0.1029
#> 3   Item3 random_dif -0.0310 0.1654 0.8515       1     -0.3552 0.2933
#> 4   Item4 random_dif -0.1667 0.1667 0.3175       1     -0.4935 0.1601
#> 5   Item5 random_dif  0.2042 0.1675 0.2230       1     -0.1242 0.5325
#> 6   Item6 random_dif  0.1094 0.1672 0.5132       1     -0.2184 0.4371
#> 7   Item7 random_dif  0.1184 0.1738 0.4959       1     -0.2223 0.4590
#> 8   Item8 random_dif -0.1832 0.1730 0.2896       1     -0.5224 0.1559
#> 9   Item9 random_dif  0.0280 0.1738 0.8720       1     -0.3126 0.3686
#> 10 Item10 random_dif  0.0970 0.1903 0.6102       1     -0.2760 0.4701
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1130 0.1725 0.5126       1     -0.4511 0.2252
#> 2   Item2 random_dif -0.1185 0.1666 0.4768       1     -0.4451 0.2080
#> 3   Item3 random_dif -0.0167 0.1673 0.9205       1     -0.3446 0.3112
#> 4   Item4 random_dif -0.0793 0.1696 0.6403       1     -0.4117 0.2532
#> 5   Item5 random_dif  0.1200 0.1688 0.4771       1     -0.2108 0.4508
#> 6   Item6 random_dif  0.2526 0.1617 0.1183       1     -0.0643 0.5695
#> 7   Item7 random_dif -0.0808 0.1782 0.6503       1     -0.4302 0.2685
#> 8   Item8 random_dif -0.0225 0.1652 0.8917       1     -0.3464 0.3014
#> 9   Item9 random_dif  0.0150 0.1751 0.9319       1     -0.3283 0.3582
#> 10 Item10 random_dif  0.0239 0.1730 0.8902       1     -0.3153 0.3631
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1321 0.1632 0.4181  1.0000     -0.4520 0.1877
#> 2   Item2 random_dif  0.3299 0.1562 0.0347  0.3473      0.0237 0.6361
#> 3   Item3 random_dif -0.2178 0.1628 0.1810  1.0000     -0.5370 0.1013
#> 4   Item4 random_dif -0.3046 0.1557 0.0504  0.5042     -0.6098 0.0006
#> 5   Item5 random_dif -0.0063 0.1615 0.9691  1.0000     -0.3227 0.3102
#> 6   Item6 random_dif  0.1155 0.1604 0.4714  1.0000     -0.1988 0.4298
#> 7   Item7 random_dif  0.0781 0.1616 0.6289  1.0000     -0.2386 0.3948
#> 8   Item8 random_dif  0.1328 0.1685 0.4307  1.0000     -0.1975 0.4630
#> 9   Item9 random_dif -0.2032 0.1617 0.2088  1.0000     -0.5201 0.1137
#> 10 Item10 random_dif  0.2105 0.1594 0.1865  1.0000     -0.1018 0.5229
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1801 0.1669 0.2803  1.0000     -0.5072  0.1469
#> 2   Item2 random_dif  0.0079 0.1633 0.9613  1.0000     -0.3122  0.3281
#> 3   Item3 random_dif -0.1256 0.1675 0.4531  1.0000     -0.4539  0.2026
#> 4   Item4 random_dif -0.1165 0.1668 0.4848  1.0000     -0.4433  0.2103
#> 5   Item5 random_dif  0.1009 0.1641 0.5387  1.0000     -0.2207  0.4224
#> 6   Item6 random_dif  0.0307 0.1669 0.8540  1.0000     -0.2965  0.3579
#> 7   Item7 random_dif  0.3228 0.1486 0.0298  0.2978      0.0317  0.6140
#> 8   Item8 random_dif  0.0559 0.1692 0.7411  1.0000     -0.2758  0.3876
#> 9   Item9 random_dif -0.3629 0.1486 0.0146  0.1458     -0.6541 -0.0717
#> 10 Item10 random_dif  0.2420 0.1596 0.1296  1.0000     -0.0709  0.5549
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1569 0.1750 0.3699  1.0000     -0.4998 0.1860
#> 2   Item2 random_dif  0.1231 0.1692 0.4671  1.0000     -0.2086 0.4547
#> 3   Item3 random_dif -0.0833 0.1803 0.6439  1.0000     -0.4367 0.2701
#> 4   Item4 random_dif  0.0365 0.1678 0.8280  1.0000     -0.2924 0.3653
#> 5   Item5 random_dif -0.0604 0.1668 0.7174  1.0000     -0.3872 0.2665
#> 6   Item6 random_dif -0.0778 0.1785 0.6630  1.0000     -0.4276 0.2721
#> 7   Item7 random_dif  0.0405 0.1668 0.8080  1.0000     -0.2863 0.3674
#> 8   Item8 random_dif -0.1210 0.1646 0.4621  1.0000     -0.4436 0.2015
#> 9   Item9 random_dif  0.3724 0.1500 0.0130  0.1302      0.0785 0.6664
#> 10 Item10 random_dif -0.1314 0.1739 0.4498  1.0000     -0.4721 0.2094
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.2127 0.1597 0.1829  1.0000     -0.1003  0.5256
#> 2   Item2 random_dif -0.0676 0.1675 0.6865  1.0000     -0.3958  0.2606
#> 3   Item3 random_dif -0.1214 0.1725 0.4817  1.0000     -0.4595  0.2167
#> 4   Item4 random_dif -0.0510 0.1695 0.7633  1.0000     -0.3832  0.2811
#> 5   Item5 random_dif -0.0530 0.1684 0.7528  1.0000     -0.3830  0.2770
#> 6   Item6 random_dif -0.1694 0.1748 0.3327  1.0000     -0.5120  0.1733
#> 7   Item7 random_dif  0.1875 0.1669 0.2612  1.0000     -0.1396  0.5146
#> 8   Item8 random_dif  0.3344 0.1479 0.0237  0.2375      0.0446  0.6243
#> 9   Item9 random_dif -0.3391 0.1471 0.0211  0.2114     -0.6274 -0.0508
#> 10 Item10 random_dif  0.0366 0.1728 0.8322  1.0000     -0.3020  0.3752
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0430 0.1656 0.7953  1.0000     -0.3675 0.2816
#> 2   Item2 random_dif -0.0974 0.1679 0.5618  1.0000     -0.4266 0.2317
#> 3   Item3 random_dif  0.0578 0.1627 0.7224  1.0000     -0.2611 0.3767
#> 4   Item4 random_dif -0.1908 0.1637 0.2437  1.0000     -0.5117 0.1300
#> 5   Item5 random_dif -0.0276 0.1717 0.8723  1.0000     -0.3642 0.3090
#> 6   Item6 random_dif -0.2432 0.1542 0.1146  1.0000     -0.5454 0.0589
#> 7   Item7 random_dif  0.2676 0.1571 0.0885  0.8845     -0.0403 0.5756
#> 8   Item8 random_dif  0.2277 0.1609 0.1571  1.0000     -0.0877 0.5430
#> 9   Item9 random_dif  0.0300 0.1638 0.8549  1.0000     -0.2911 0.3510
#> 10 Item10 random_dif  0.0273 0.1645 0.8684  1.0000     -0.2952 0.3498
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0412 0.1710 0.8096       1     -0.2939 0.3763
#> 2   Item2 random_dif  0.0444 0.1699 0.7939       1     -0.2885 0.3773
#> 3   Item3 random_dif -0.0099 0.1684 0.9530       1     -0.3400 0.3202
#> 4   Item4 random_dif  0.2389 0.1661 0.1504       1     -0.0867 0.5645
#> 5   Item5 random_dif  0.1891 0.1691 0.2634       1     -0.1423 0.5205
#> 6   Item6 random_dif -0.1640 0.1758 0.3508       1     -0.5086 0.1805
#> 7   Item7 random_dif  0.0438 0.1716 0.7988       1     -0.2927 0.3802
#> 8   Item8 random_dif -0.1548 0.1643 0.3461       1     -0.4769 0.1672
#> 9   Item9 random_dif -0.0265 0.1755 0.8802       1     -0.3705 0.3176
#> 10 Item10 random_dif -0.2582 0.1741 0.1382       1     -0.5995 0.0831
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0033 0.1717 0.9845  1.0000     -0.3332 0.3399
#> 2   Item2 random_dif -0.0081 0.1685 0.9618  1.0000     -0.3383 0.3221
#> 3   Item3 random_dif -0.0744 0.1704 0.6625  1.0000     -0.4084 0.2597
#> 4   Item4 random_dif  0.2067 0.1707 0.2259  1.0000     -0.1278 0.5411
#> 5   Item5 random_dif  0.1974 0.1641 0.2289  1.0000     -0.1242 0.5190
#> 6   Item6 random_dif -0.2995 0.1671 0.0732  0.7316     -0.6271 0.0281
#> 7   Item7 random_dif  0.0605 0.1660 0.7157  1.0000     -0.2650 0.3859
#> 8   Item8 random_dif -0.1012 0.1758 0.5648  1.0000     -0.4458 0.2434
#> 9   Item9 random_dif  0.0032 0.1683 0.9847  1.0000     -0.3266 0.3331
#> 10 Item10 random_dif -0.0051 0.1735 0.9767  1.0000     -0.3451 0.3349
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0360 0.1719 0.8342       1     -0.3009 0.3729
#> 2   Item2 random_dif -0.0924 0.1629 0.5707       1     -0.4117 0.2270
#> 3   Item3 random_dif  0.0420 0.1747 0.8098       1     -0.3003 0.3844
#> 4   Item4 random_dif -0.0137 0.1685 0.9350       1     -0.3440 0.3165
#> 5   Item5 random_dif -0.2128 0.1626 0.1907       1     -0.5316 0.1059
#> 6   Item6 random_dif -0.0205 0.1672 0.9022       1     -0.3482 0.3071
#> 7   Item7 random_dif  0.0249 0.1716 0.8846       1     -0.3114 0.3613
#> 8   Item8 random_dif  0.1256 0.1627 0.4401       1     -0.1933 0.4446
#> 9   Item9 random_dif  0.1319 0.1674 0.4307       1     -0.1961 0.4599
#> 10 Item10 random_dif -0.0141 0.1704 0.9339       1     -0.3482 0.3199
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0000 0.1651 1.0000       1     -0.3236 0.3236
#> 2   Item2 random_dif  0.0604 0.1789 0.7356       1     -0.2902 0.4110
#> 3   Item3 random_dif -0.0361 0.1692 0.8311       1     -0.3677 0.2956
#> 4   Item4 random_dif -0.1132 0.1624 0.4857       1     -0.4316 0.2051
#> 5   Item5 random_dif -0.0524 0.1672 0.7540       1     -0.3802 0.2754
#> 6   Item6 random_dif -0.1784 0.1657 0.2817       1     -0.5033 0.1464
#> 7   Item7 random_dif -0.0902 0.1617 0.5769       1     -0.4072 0.2268
#> 8   Item8 random_dif  0.0437 0.1675 0.7944       1     -0.2847 0.3720
#> 9   Item9 random_dif  0.1749 0.1602 0.2750       1     -0.1391 0.4889
#> 10 Item10 random_dif  0.1967 0.1807 0.2763       1     -0.1575 0.5508
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2455 0.1608 0.1268       1     -0.0696 0.5606
#> 2   Item2 random_dif -0.2238 0.1615 0.1658       1     -0.5402 0.0927
#> 3   Item3 random_dif  0.0903 0.1649 0.5840       1     -0.2328 0.4134
#> 4   Item4 random_dif  0.1525 0.1642 0.3529       1     -0.1692 0.4742
#> 5   Item5 random_dif -0.2500 0.1594 0.1168       1     -0.5624 0.0624
#> 6   Item6 random_dif -0.0348 0.1638 0.8316       1     -0.3560 0.2863
#> 7   Item7 random_dif  0.0838 0.1632 0.6076       1     -0.2360 0.4036
#> 8   Item8 random_dif -0.0916 0.1716 0.5935       1     -0.4280 0.2448
#> 9   Item9 random_dif  0.1699 0.1640 0.3002       1     -0.1515 0.4913
#> 10 Item10 random_dif -0.1517 0.1683 0.3672       1     -0.4815 0.1781
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0094 0.1810 0.9587       1     -0.3454 0.3641
#> 2   Item2 random_dif -0.2122 0.1702 0.2125       1     -0.5458 0.1214
#> 3   Item3 random_dif -0.0968 0.1730 0.5760       1     -0.4359 0.2424
#> 4   Item4 random_dif  0.0820 0.1635 0.6161       1     -0.2384 0.4023
#> 5   Item5 random_dif -0.0699 0.1625 0.6671       1     -0.3885 0.2486
#> 6   Item6 random_dif -0.0915 0.1638 0.5765       1     -0.4126 0.2296
#> 7   Item7 random_dif  0.2043 0.1574 0.1943       1     -0.1042 0.5127
#> 8   Item8 random_dif -0.1973 0.1638 0.2285       1     -0.5184 0.1238
#> 9   Item9 random_dif  0.2267 0.1626 0.1632       1     -0.0919 0.5453
#> 10 Item10 random_dif  0.1268 0.1696 0.4547       1     -0.2056 0.4591
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1882 0.1699 0.2680       1     -0.1448 0.5211
#> 2   Item2 random_dif -0.1387 0.1661 0.4035       1     -0.4642 0.1868
#> 3   Item3 random_dif  0.1138 0.1731 0.5110       1     -0.2255 0.4530
#> 4   Item4 random_dif  0.2426 0.1615 0.1330       1     -0.0739 0.5592
#> 5   Item5 random_dif  0.0181 0.1826 0.9210       1     -0.3397 0.3759
#> 6   Item6 random_dif -0.1242 0.1683 0.4607       1     -0.4541 0.2057
#> 7   Item7 random_dif -0.2036 0.1743 0.2429       1     -0.5452 0.1381
#> 8   Item8 random_dif -0.1144 0.1700 0.5008       1     -0.4476 0.2187
#> 9   Item9 random_dif -0.0158 0.1761 0.9287       1     -0.3610 0.3294
#> 10 Item10 random_dif  0.0369 0.1736 0.8316       1     -0.3033 0.3771
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3313 0.1533 0.0307  0.3070     -0.6318 -0.0308
#> 2   Item2 random_dif  0.0735 0.1590 0.6437  1.0000     -0.2380  0.3851
#> 3   Item3 random_dif  0.0108 0.1641 0.9474  1.0000     -0.3107  0.3324
#> 4   Item4 random_dif  0.3365 0.1539 0.0288  0.2876      0.0349  0.6381
#> 5   Item5 random_dif -0.1060 0.1682 0.5285  1.0000     -0.4358  0.2237
#> 6   Item6 random_dif -0.0693 0.1693 0.6823  1.0000     -0.4011  0.2625
#> 7   Item7 random_dif -0.0808 0.1631 0.6203  1.0000     -0.4006  0.2389
#> 8   Item8 random_dif  0.0252 0.1601 0.8750  1.0000     -0.2886  0.3390
#> 9   Item9 random_dif  0.2989 0.1661 0.0720  0.7205     -0.0268  0.6245
#> 10 Item10 random_dif -0.1563 0.1659 0.3461  1.0000     -0.4814  0.1688
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0228 0.1710 0.8941       1     -0.3125 0.3580
#> 2   Item2 random_dif -0.1973 0.1635 0.2276       1     -0.5177 0.1232
#> 3   Item3 random_dif  0.1536 0.1630 0.3461       1     -0.1659 0.4731
#> 4   Item4 random_dif  0.1115 0.1699 0.5118       1     -0.2216 0.4445
#> 5   Item5 random_dif -0.1149 0.1676 0.4931       1     -0.4433 0.2136
#> 6   Item6 random_dif -0.0903 0.1669 0.5885       1     -0.4174 0.2368
#> 7   Item7 random_dif -0.0673 0.1651 0.6835       1     -0.3909 0.2563
#> 8   Item8 random_dif  0.1217 0.1657 0.4626       1     -0.2030 0.4465
#> 9   Item9 random_dif  0.0695 0.1691 0.6809       1     -0.2618 0.4009
#> 10 Item10 random_dif -0.0175 0.1717 0.9189       1     -0.3540 0.3190
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3424 0.1572 0.0293  0.2935     -0.6504 -0.0344
#> 2   Item2 random_dif  0.0380 0.1827 0.8353  1.0000     -0.3200  0.3960
#> 3   Item3 random_dif  0.0578 0.1658 0.7275  1.0000     -0.2672  0.3828
#> 4   Item4 random_dif -0.2155 0.1719 0.2100  1.0000     -0.5525  0.1215
#> 5   Item5 random_dif  0.0428 0.1701 0.8012  1.0000     -0.2905  0.3762
#> 6   Item6 random_dif  0.1315 0.1634 0.4208  1.0000     -0.1887  0.4517
#> 7   Item7 random_dif -0.0810 0.1764 0.6462  1.0000     -0.4267  0.2647
#> 8   Item8 random_dif  0.0963 0.1745 0.5810  1.0000     -0.2457  0.4382
#> 9   Item9 random_dif -0.0752 0.1699 0.6579  1.0000     -0.4082  0.2577
#> 10 Item10 random_dif  0.3183 0.1574 0.0432  0.4318      0.0098  0.6268
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1279 0.1669 0.4435  1.0000     -0.1992  0.4549
#> 2   Item2 random_dif  0.0349 0.1702 0.8375  1.0000     -0.2986  0.3684
#> 3   Item3 random_dif -0.1083 0.1635 0.5076  1.0000     -0.4287  0.2121
#> 4   Item4 random_dif -0.0959 0.1714 0.5758  1.0000     -0.4318  0.2400
#> 5   Item5 random_dif  0.1537 0.1633 0.3465  1.0000     -0.1663  0.4738
#> 6   Item6 random_dif -0.1266 0.1650 0.4428  1.0000     -0.4500  0.1967
#> 7   Item7 random_dif  0.2197 0.1667 0.1874  1.0000     -0.1070  0.5463
#> 8   Item8 random_dif  0.1018 0.1701 0.5495  1.0000     -0.2316  0.4353
#> 9   Item9 random_dif -0.4521 0.1408 0.0013  0.0133   * -0.7282 -0.1761
#> 10 Item10 random_dif  0.1647 0.1675 0.3255  1.0000     -0.1636  0.4929
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3365 0.1537 0.0285  0.2854     -0.6378 -0.0353
#> 2   Item2 random_dif  0.0601 0.1752 0.7316  1.0000     -0.2833  0.4035
#> 3   Item3 random_dif -0.3534 0.1479 0.0169  0.1689     -0.6433 -0.0635
#> 4   Item4 random_dif  0.2924 0.1593 0.0665  0.6647     -0.0199  0.6046
#> 5   Item5 random_dif -0.0120 0.1732 0.9447  1.0000     -0.3514  0.3274
#> 6   Item6 random_dif  0.2742 0.1612 0.0889  0.8889     -0.0417  0.5901
#> 7   Item7 random_dif  0.1630 0.1623 0.3151  1.0000     -0.1550  0.4810
#> 8   Item8 random_dif -0.0518 0.1717 0.7630  1.0000     -0.3884  0.2848
#> 9   Item9 random_dif -0.0497 0.1764 0.7780  1.0000     -0.3955  0.2960
#> 10 Item10 random_dif  0.0363 0.1757 0.8364  1.0000     -0.3081  0.3806
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.4190 0.1378 0.0024  0.0236   * -0.6892 -0.1489
#> 2   Item2 random_dif -0.1250 0.1625 0.4418  1.0000     -0.4435  0.1935
#> 3   Item3 random_dif  0.1131 0.1709 0.5079  1.0000     -0.2218  0.4480
#> 4   Item4 random_dif -0.1689 0.1598 0.2906  1.0000     -0.4821  0.1443
#> 5   Item5 random_dif -0.0748 0.1645 0.6492  1.0000     -0.3972  0.2476
#> 6   Item6 random_dif  0.2447 0.1639 0.1355  1.0000     -0.0766  0.5660
#> 7   Item7 random_dif  0.1075 0.1631 0.5099  1.0000     -0.2122  0.4272
#> 8   Item8 random_dif -0.0234 0.1575 0.8819  1.0000     -0.3321  0.2853
#> 9   Item9 random_dif  0.1486 0.1600 0.3529  1.0000     -0.1649  0.4621
#> 10 Item10 random_dif  0.2231 0.1538 0.1469  1.0000     -0.0783  0.5245
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2132 0.1594 0.1811  1.0000     -0.5257 0.0993
#> 2   Item2 random_dif  0.3288 0.1559 0.0350  0.3497      0.0232 0.6343
#> 3   Item3 random_dif  0.0164 0.1677 0.9221  1.0000     -0.3122 0.3450
#> 4   Item4 random_dif -0.2021 0.1661 0.2235  1.0000     -0.5276 0.1234
#> 5   Item5 random_dif -0.2049 0.1762 0.2448  1.0000     -0.5503 0.1404
#> 6   Item6 random_dif  0.2445 0.1542 0.1127  1.0000     -0.0576 0.5467
#> 7   Item7 random_dif  0.1978 0.1664 0.2343  1.0000     -0.1282 0.5239
#> 8   Item8 random_dif -0.0036 0.1674 0.9827  1.0000     -0.3317 0.3245
#> 9   Item9 random_dif  0.0698 0.1659 0.6737  1.0000     -0.2552 0.3949
#> 10 Item10 random_dif -0.3140 0.1676 0.0610  0.6098     -0.6426 0.0145
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1357 0.1780 0.4460  1.0000     -0.2132  0.4846
#> 2   Item2 random_dif -0.0272 0.1668 0.8704  1.0000     -0.3541  0.2997
#> 3   Item3 random_dif -0.1084 0.1873 0.5628  1.0000     -0.4755  0.2587
#> 4   Item4 random_dif  0.0377 0.1745 0.8290  1.0000     -0.3044  0.3798
#> 5   Item5 random_dif -0.3787 0.1485 0.0108  0.1077     -0.6698 -0.0876
#> 6   Item6 random_dif -0.1511 0.1781 0.3963  1.0000     -0.5001  0.1980
#> 7   Item7 random_dif  0.2521 0.1615 0.1185  1.0000     -0.0644  0.5686
#> 8   Item8 random_dif -0.1048 0.1775 0.5551  1.0000     -0.4527  0.2432
#> 9   Item9 random_dif  0.1069 0.1695 0.5283  1.0000     -0.2253  0.4391
#> 10 Item10 random_dif  0.1980 0.1620 0.2218  1.0000     -0.1196  0.5155
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1141 0.1706 0.5037  1.0000     -0.2203 0.4484
#> 2   Item2 random_dif  0.0734 0.1727 0.6707  1.0000     -0.2651 0.4119
#> 3   Item3 random_dif -0.1027 0.1665 0.5374  1.0000     -0.4291 0.2236
#> 4   Item4 random_dif  0.0977 0.1684 0.5617  1.0000     -0.2324 0.4279
#> 5   Item5 random_dif -0.0712 0.1662 0.6685  1.0000     -0.3969 0.2546
#> 6   Item6 random_dif  0.1514 0.1670 0.3647  1.0000     -0.1760 0.4788
#> 7   Item7 random_dif  0.1076 0.1611 0.5041  1.0000     -0.2082 0.4234
#> 8   Item8 random_dif -0.2172 0.1639 0.1849  1.0000     -0.5384 0.1039
#> 9   Item9 random_dif  0.1451 0.1776 0.4139  1.0000     -0.2030 0.4933
#> 10 Item10 random_dif -0.3180 0.1715 0.0637  0.6366     -0.6541 0.0181
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0891 0.1656 0.5906  1.0000     -0.4137 0.2355
#> 2   Item2 random_dif  0.0756 0.1718 0.6600  1.0000     -0.2611 0.4122
#> 3   Item3 random_dif -0.0774 0.1681 0.6451  1.0000     -0.4070 0.2521
#> 4   Item4 random_dif  0.0468 0.1709 0.7842  1.0000     -0.2882 0.3818
#> 5   Item5 random_dif  0.3571 0.1570 0.0229  0.2288      0.0495 0.6648
#> 6   Item6 random_dif -0.0558 0.1617 0.7299  1.0000     -0.3727 0.2610
#> 7   Item7 random_dif -0.0333 0.1678 0.8425  1.0000     -0.3622 0.2955
#> 8   Item8 random_dif -0.1363 0.1689 0.4198  1.0000     -0.4674 0.1948
#> 9   Item9 random_dif  0.1738 0.1627 0.2852  1.0000     -0.1450 0.4927
#> 10 Item10 random_dif -0.2417 0.1626 0.1372  1.0000     -0.5604 0.0770
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2463 0.1552 0.1125  1.0000     -0.0579 0.5505
#> 2   Item2 random_dif -0.1835 0.1686 0.2764  1.0000     -0.5140 0.1470
#> 3   Item3 random_dif  0.0094 0.1768 0.9575  1.0000     -0.3371 0.3559
#> 4   Item4 random_dif  0.1424 0.1662 0.3918  1.0000     -0.1834 0.4682
#> 5   Item5 random_dif  0.0915 0.1639 0.5764  1.0000     -0.2297 0.4128
#> 6   Item6 random_dif -0.1012 0.1672 0.5450  1.0000     -0.4289 0.2265
#> 7   Item7 random_dif -0.1339 0.1680 0.4253  1.0000     -0.4632 0.1953
#> 8   Item8 random_dif -0.1128 0.1749 0.5189  1.0000     -0.4557 0.2300
#> 9   Item9 random_dif -0.2680 0.1626 0.0993  0.9933     -0.5866 0.0507
#> 10 Item10 random_dif  0.2860 0.1582 0.0706  0.7064     -0.0241 0.5960
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.4028 0.1408 0.0042  0.0424   * -0.6788 -0.1267
#> 2   Item2 random_dif  0.1429 0.1684 0.3961  1.0000     -0.1871  0.4728
#> 3   Item3 random_dif -0.1085 0.1714 0.5267  1.0000     -0.4445  0.2274
#> 4   Item4 random_dif -0.0083 0.1737 0.9619  1.0000     -0.3487  0.3321
#> 5   Item5 random_dif  0.0146 0.1702 0.9315  1.0000     -0.3189  0.3481
#> 6   Item6 random_dif -0.0915 0.1718 0.5944  1.0000     -0.4283  0.2453
#> 7   Item7 random_dif -0.1647 0.1684 0.3281  1.0000     -0.4948  0.1654
#> 8   Item8 random_dif  0.2852 0.1680 0.0896  0.8964     -0.0441  0.6146
#> 9   Item9 random_dif  0.2421 0.1714 0.1578  1.0000     -0.0938  0.5781
#> 10 Item10 random_dif  0.1693 0.1671 0.3109  1.0000     -0.1582  0.4967
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2812 0.1607 0.0801  0.8008     -0.5962 0.0337
#> 2   Item2 random_dif -0.0098 0.1661 0.9528  1.0000     -0.3354 0.3157
#> 3   Item3 random_dif  0.2424 0.1622 0.1350  1.0000     -0.0754 0.5602
#> 4   Item4 random_dif  0.0305 0.1685 0.8563  1.0000     -0.2997 0.3607
#> 5   Item5 random_dif  0.1641 0.1683 0.3297  1.0000     -0.1658 0.4940
#> 6   Item6 random_dif -0.0333 0.1660 0.8408  1.0000     -0.3586 0.2919
#> 7   Item7 random_dif -0.0206 0.1691 0.9030  1.0000     -0.3521 0.3108
#> 8   Item8 random_dif -0.0052 0.1703 0.9756  1.0000     -0.3390 0.3285
#> 9   Item9 random_dif  0.0172 0.1706 0.9195  1.0000     -0.3171 0.3516
#> 10 Item10 random_dif -0.1239 0.1685 0.4621  1.0000     -0.4542 0.2064
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1384 0.1696 0.4143       1     -0.4707 0.1939
#> 2   Item2 random_dif -0.0878 0.1732 0.6123       1     -0.4273 0.2517
#> 3   Item3 random_dif  0.0069 0.1699 0.9677       1     -0.3261 0.3398
#> 4   Item4 random_dif -0.0184 0.1711 0.9145       1     -0.3538 0.3171
#> 5   Item5 random_dif  0.0045 0.1644 0.9784       1     -0.3177 0.3266
#> 6   Item6 random_dif -0.0643 0.1688 0.7033       1     -0.3952 0.2666
#> 7   Item7 random_dif  0.0874 0.1681 0.6029       1     -0.2420 0.4169
#> 8   Item8 random_dif  0.0334 0.1746 0.8483       1     -0.3088 0.3756
#> 9   Item9 random_dif  0.0304 0.1702 0.8582       1     -0.3032 0.3641
#> 10 Item10 random_dif  0.1341 0.1639 0.4132       1     -0.1871 0.4553
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0350 0.1676 0.8344  1.0000     -0.3635  0.2934
#> 2   Item2 random_dif  0.1063 0.1738 0.5407  1.0000     -0.2343  0.4469
#> 3   Item3 random_dif -0.0227 0.1675 0.8921  1.0000     -0.3510  0.3056
#> 4   Item4 random_dif -0.4271 0.1482 0.0040  0.0396   * -0.7176 -0.1366
#> 5   Item5 random_dif -0.0084 0.1681 0.9601  1.0000     -0.3379  0.3211
#> 6   Item6 random_dif -0.1259 0.1699 0.4587  1.0000     -0.4588  0.2071
#> 7   Item7 random_dif -0.0525 0.1631 0.7477  1.0000     -0.3722  0.2672
#> 8   Item8 random_dif -0.0016 0.1666 0.9925  1.0000     -0.3281  0.3250
#> 9   Item9 random_dif  0.2174 0.1670 0.1931  1.0000     -0.1100  0.5447
#> 10 Item10 random_dif  0.3863 0.1625 0.0174  0.1745      0.0678  0.7048
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.4407 0.1483 0.0030  0.0296   *  0.1501  0.7313
#> 2   Item2 random_dif  0.0671 0.1636 0.6817  1.0000     -0.2536  0.3878
#> 3   Item3 random_dif  0.0433 0.1774 0.8070  1.0000     -0.3043  0.3910
#> 4   Item4 random_dif -0.2682 0.1544 0.0823  0.8228     -0.5707  0.0343
#> 5   Item5 random_dif -0.3787 0.1499 0.0115  0.1149     -0.6724 -0.0850
#> 6   Item6 random_dif  0.2791 0.1563 0.0743  0.7425     -0.0273  0.5855
#> 7   Item7 random_dif -0.0144 0.1727 0.9334  1.0000     -0.3530  0.3241
#> 8   Item8 random_dif  0.0840 0.1689 0.6188  1.0000     -0.2470  0.4150
#> 9   Item9 random_dif -0.1617 0.1593 0.3101  1.0000     -0.4740  0.1506
#> 10 Item10 random_dif -0.0742 0.1784 0.6774  1.0000     -0.4239  0.2754
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1325 0.1549 0.3926       1     -0.1712 0.4361
#> 2   Item2 random_dif  0.0032 0.1575 0.9837       1     -0.3054 0.3119
#> 3   Item3 random_dif -0.1791 0.1664 0.2819       1     -0.5053 0.1471
#> 4   Item4 random_dif  0.2356 0.1612 0.1437       1     -0.0803 0.5515
#> 5   Item5 random_dif -0.0739 0.1640 0.6520       1     -0.3953 0.2475
#> 6   Item6 random_dif -0.0638 0.1800 0.7228       1     -0.4166 0.2889
#> 7   Item7 random_dif -0.0553 0.1732 0.7493       1     -0.3948 0.2841
#> 8   Item8 random_dif -0.2462 0.1534 0.1086       1     -0.5469 0.0546
#> 9   Item9 random_dif  0.2067 0.1578 0.1901       1     -0.1025 0.5159
#> 10 Item10 random_dif  0.0019 0.1689 0.9909       1     -0.3292 0.3330
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0253 0.1652 0.8783  1.0000     -0.2986 0.3491
#> 2   Item2 random_dif  0.0733 0.1807 0.6852  1.0000     -0.2809 0.4275
#> 3   Item3 random_dif -0.0110 0.1747 0.9497  1.0000     -0.3535 0.3314
#> 4   Item4 random_dif  0.1386 0.1604 0.3875  1.0000     -0.1758 0.4530
#> 5   Item5 random_dif  0.2430 0.1591 0.1267  1.0000     -0.0689 0.5548
#> 6   Item6 random_dif  0.0104 0.1739 0.9522  1.0000     -0.3304 0.3512
#> 7   Item7 random_dif -0.1080 0.1695 0.5241  1.0000     -0.4403 0.2243
#> 8   Item8 random_dif -0.2771 0.1641 0.0913  0.9135     -0.5988 0.0446
#> 9   Item9 random_dif -0.0667 0.1719 0.6981  1.0000     -0.4036 0.2702
#> 10 Item10 random_dif -0.0405 0.1787 0.8209  1.0000     -0.3907 0.3098
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3046 0.1515 0.0443  0.4433     -0.6016 -0.0077
#> 2   Item2 random_dif -0.2384 0.1599 0.1360  1.0000     -0.5519  0.0750
#> 3   Item3 random_dif  0.0948 0.1633 0.5613  1.0000     -0.2251  0.4148
#> 4   Item4 random_dif  0.3377 0.1478 0.0223  0.2229      0.0481  0.6274
#> 5   Item5 random_dif  0.3094 0.1458 0.0338  0.3382      0.0237  0.5951
#> 6   Item6 random_dif -0.0473 0.1651 0.7746  1.0000     -0.3709  0.2763
#> 7   Item7 random_dif  0.0680 0.1608 0.6723  1.0000     -0.2472  0.3833
#> 8   Item8 random_dif -0.0526 0.1670 0.7526  1.0000     -0.3799  0.2746
#> 9   Item9 random_dif -0.0388 0.1578 0.8056  1.0000     -0.3481  0.2705
#> 10 Item10 random_dif -0.2215 0.1768 0.2102  1.0000     -0.5680  0.1250
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1164 0.1666 0.4848       1     -0.2101 0.4428
#> 2   Item2 random_dif -0.1300 0.1675 0.4375       1     -0.4583 0.1982
#> 3   Item3 random_dif  0.2066 0.1566 0.1871       1     -0.1004 0.5137
#> 4   Item4 random_dif  0.0403 0.1761 0.8190       1     -0.3048 0.3854
#> 5   Item5 random_dif  0.0030 0.1647 0.9854       1     -0.3198 0.3259
#> 6   Item6 random_dif  0.0955 0.1736 0.5822       1     -0.2447 0.4357
#> 7   Item7 random_dif -0.1555 0.1558 0.3181       1     -0.4608 0.1498
#> 8   Item8 random_dif -0.0427 0.1606 0.7906       1     -0.3574 0.2721
#> 9   Item9 random_dif -0.0558 0.1586 0.7249       1     -0.3667 0.2550
#> 10 Item10 random_dif -0.0401 0.1631 0.8056       1     -0.3597 0.2795
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0875 0.1670 0.6002  1.0000     -0.2398 0.4149
#> 2   Item2 random_dif -0.2895 0.1567 0.0647  0.6473     -0.5967 0.0177
#> 3   Item3 random_dif -0.1339 0.1590 0.3999  1.0000     -0.4455 0.1778
#> 4   Item4 random_dif  0.0332 0.1682 0.8434  1.0000     -0.2964 0.3629
#> 5   Item5 random_dif  0.1051 0.1709 0.5386  1.0000     -0.2299 0.4401
#> 6   Item6 random_dif -0.0681 0.1757 0.6983  1.0000     -0.4124 0.2762
#> 7   Item7 random_dif  0.0101 0.1695 0.9523  1.0000     -0.3220 0.3423
#> 8   Item8 random_dif  0.1761 0.1662 0.2895  1.0000     -0.1497 0.5019
#> 9   Item9 random_dif -0.0175 0.1721 0.9188  1.0000     -0.3549 0.3199
#> 10 Item10 random_dif  0.1429 0.1789 0.4246  1.0000     -0.2078 0.4935
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0215 0.1698 0.8992       1     -0.3114 0.3544
#> 2   Item2 random_dif  0.0435 0.1723 0.8008       1     -0.2942 0.3812
#> 3   Item3 random_dif  0.0083 0.1778 0.9626       1     -0.3401 0.3567
#> 4   Item4 random_dif  0.1024 0.1678 0.5417       1     -0.2265 0.4313
#> 5   Item5 random_dif  0.0591 0.1758 0.7370       1     -0.2856 0.4037
#> 6   Item6 random_dif  0.0531 0.1732 0.7593       1     -0.2864 0.3925
#> 7   Item7 random_dif -0.2073 0.1577 0.1886       1     -0.5163 0.1017
#> 8   Item8 random_dif -0.0609 0.1838 0.7405       1     -0.4210 0.2993
#> 9   Item9 random_dif -0.0056 0.1674 0.9731       1     -0.3338 0.3225
#> 10 Item10 random_dif -0.0040 0.1710 0.9816       1     -0.3392 0.3313
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1286 0.1591 0.4188  1.0000     -0.1832 0.4405
#> 2   Item2 random_dif -0.1691 0.1549 0.2750  1.0000     -0.4728 0.1345
#> 3   Item3 random_dif -0.1542 0.1661 0.3533  1.0000     -0.4797 0.1714
#> 4   Item4 random_dif -0.2655 0.1561 0.0889  0.8888     -0.5714 0.0404
#> 5   Item5 random_dif  0.0821 0.1635 0.6157  1.0000     -0.2384 0.4025
#> 6   Item6 random_dif -0.1252 0.1637 0.4443  1.0000     -0.4459 0.1956
#> 7   Item7 random_dif  0.2759 0.1616 0.0879  0.8789     -0.0410 0.5927
#> 8   Item8 random_dif -0.1070 0.1606 0.5052  1.0000     -0.4217 0.2077
#> 9   Item9 random_dif  0.1739 0.1603 0.2782  1.0000     -0.1404 0.4881
#> 10 Item10 random_dif  0.1766 0.1593 0.2676  1.0000     -0.1357 0.4889
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2238 0.1622 0.1677  1.0000     -0.5416 0.0941
#> 2   Item2 random_dif  0.1103 0.1698 0.5159  1.0000     -0.2225 0.4431
#> 3   Item3 random_dif  0.0241 0.1662 0.8849  1.0000     -0.3016 0.3498
#> 4   Item4 random_dif  0.0397 0.1687 0.8139  1.0000     -0.2910 0.3704
#> 5   Item5 random_dif -0.2846 0.1705 0.0952  0.9519     -0.6188 0.0497
#> 6   Item6 random_dif -0.1277 0.1713 0.4558  1.0000     -0.4634 0.2080
#> 7   Item7 random_dif  0.2184 0.1559 0.1613  1.0000     -0.0872 0.5239
#> 8   Item8 random_dif  0.0622 0.1742 0.7211  1.0000     -0.2792 0.4035
#> 9   Item9 random_dif  0.1186 0.1690 0.4828  1.0000     -0.2126 0.4497
#> 10 Item10 random_dif  0.0090 0.1803 0.9603  1.0000     -0.3444 0.3624
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0256 0.1672 0.8781       1     -0.3534 0.3021
#> 2   Item2 random_dif -0.0673 0.1837 0.7140       1     -0.4273 0.2927
#> 3   Item3 random_dif  0.0913 0.1887 0.6284       1     -0.2785 0.4612
#> 4   Item4 random_dif  0.1864 0.1667 0.2635       1     -0.1403 0.5132
#> 5   Item5 random_dif -0.0021 0.1795 0.9908       1     -0.3538 0.3497
#> 6   Item6 random_dif -0.0640 0.1722 0.7101       1     -0.4015 0.2735
#> 7   Item7 random_dif  0.2383 0.1749 0.1730       1     -0.1045 0.5812
#> 8   Item8 random_dif -0.1486 0.1623 0.3598       1     -0.4666 0.1694
#> 9   Item9 random_dif -0.1494 0.1712 0.3830       1     -0.4850 0.1862
#> 10 Item10 random_dif -0.0109 0.1784 0.9511       1     -0.3606 0.3387
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0455 0.1712 0.7906  1.0000     -0.2901 0.3810
#> 2   Item2 random_dif  0.0856 0.1700 0.6146  1.0000     -0.2477 0.4189
#> 3   Item3 random_dif -0.0036 0.1715 0.9833  1.0000     -0.3397 0.3326
#> 4   Item4 random_dif -0.0537 0.1679 0.7490  1.0000     -0.3829 0.2754
#> 5   Item5 random_dif -0.2326 0.1626 0.1525  1.0000     -0.5513 0.0860
#> 6   Item6 random_dif  0.2023 0.1755 0.2491  1.0000     -0.1417 0.5462
#> 7   Item7 random_dif  0.0571 0.1756 0.7451  1.0000     -0.2871 0.4012
#> 8   Item8 random_dif -0.2872 0.1584 0.0698  0.6978     -0.5975 0.0232
#> 9   Item9 random_dif  0.1765 0.1701 0.2994  1.0000     -0.1568 0.5098
#> 10 Item10 random_dif  0.0545 0.1694 0.7474  1.0000     -0.2774 0.3865
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0975 0.1638 0.5519  1.0000     -0.4186 0.2236
#> 2   Item2 random_dif -0.0718 0.1698 0.6722  1.0000     -0.4047 0.2610
#> 3   Item3 random_dif  0.0825 0.1664 0.6202  1.0000     -0.2437 0.4086
#> 4   Item4 random_dif  0.0241 0.1716 0.8883  1.0000     -0.3123 0.3605
#> 5   Item5 random_dif  0.2976 0.1509 0.0486  0.4865      0.0018 0.5934
#> 6   Item6 random_dif  0.1577 0.1667 0.3440  1.0000     -0.1690 0.4844
#> 7   Item7 random_dif -0.1875 0.1586 0.2370  1.0000     -0.4983 0.1233
#> 8   Item8 random_dif  0.0204 0.1800 0.9098  1.0000     -0.3325 0.3733
#> 9   Item9 random_dif -0.2040 0.1546 0.1869  1.0000     -0.5071 0.0990
#> 10 Item10 random_dif -0.0232 0.1610 0.8856  1.0000     -0.3388 0.2925
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0937 0.1658 0.5719  1.0000     -0.4186 0.2312
#> 2   Item2 random_dif -0.1961 0.1567 0.2107  1.0000     -0.5033 0.1110
#> 3   Item3 random_dif  0.1092 0.1713 0.5239  1.0000     -0.2266 0.4449
#> 4   Item4 random_dif -0.1565 0.1671 0.3492  1.0000     -0.4840 0.1711
#> 5   Item5 random_dif  0.0799 0.1718 0.6420  1.0000     -0.2568 0.4166
#> 6   Item6 random_dif -0.1357 0.1669 0.4162  1.0000     -0.4627 0.1914
#> 7   Item7 random_dif -0.2588 0.1734 0.1355  1.0000     -0.5986 0.0810
#> 8   Item8 random_dif  0.2777 0.1574 0.0777  0.7773     -0.0308 0.5862
#> 9   Item9 random_dif  0.2133 0.1651 0.1965  1.0000     -0.1104 0.5369
#> 10 Item10 random_dif  0.1636 0.1718 0.3410  1.0000     -0.1731 0.5002
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1926 0.1624 0.2358       1     -0.1258 0.5109
#> 2   Item2 random_dif -0.1898 0.1662 0.2534       1     -0.5155 0.1359
#> 3   Item3 random_dif -0.1189 0.1651 0.4715       1     -0.4425 0.2047
#> 4   Item4 random_dif -0.0674 0.1750 0.7001       1     -0.4104 0.2756
#> 5   Item5 random_dif -0.0989 0.1679 0.5558       1     -0.4281 0.2302
#> 6   Item6 random_dif  0.0652 0.1735 0.7072       1     -0.2749 0.4053
#> 7   Item7 random_dif  0.1456 0.1725 0.3987       1     -0.1926 0.4837
#> 8   Item8 random_dif  0.0341 0.1716 0.8424       1     -0.3022 0.3705
#> 9   Item9 random_dif  0.0512 0.1770 0.7723       1     -0.2958 0.3982
#> 10 Item10 random_dif -0.0125 0.1684 0.9409       1     -0.3425 0.3176
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1780 0.1602 0.2664       1     -0.1359 0.4919
#> 2   Item2 random_dif -0.0332 0.1698 0.8451       1     -0.3659 0.2996
#> 3   Item3 random_dif  0.0137 0.1637 0.9331       1     -0.3071 0.3345
#> 4   Item4 random_dif  0.0776 0.1683 0.6445       1     -0.2522 0.4075
#> 5   Item5 random_dif -0.0562 0.1683 0.7386       1     -0.3861 0.2738
#> 6   Item6 random_dif -0.0536 0.1681 0.7497       1     -0.3831 0.2758
#> 7   Item7 random_dif -0.0836 0.1637 0.6096       1     -0.4044 0.2372
#> 8   Item8 random_dif -0.1056 0.1763 0.5491       1     -0.4512 0.2400
#> 9   Item9 random_dif  0.0872 0.1737 0.6157       1     -0.2532 0.4276
#> 10 Item10 random_dif -0.0454 0.1699 0.7894       1     -0.3784 0.2876
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0333 0.1652 0.8401       1     -0.3571 0.2905
#> 2   Item2 random_dif  0.0270 0.1736 0.8763       1     -0.3133 0.3674
#> 3   Item3 random_dif  0.0889 0.1696 0.6000       1     -0.2434 0.4213
#> 4   Item4 random_dif -0.0084 0.1802 0.9628       1     -0.3617 0.3449
#> 5   Item5 random_dif  0.1626 0.1764 0.3567       1     -0.1832 0.5084
#> 6   Item6 random_dif  0.0325 0.1709 0.8492       1     -0.3025 0.3675
#> 7   Item7 random_dif -0.0440 0.1674 0.7929       1     -0.3721 0.2842
#> 8   Item8 random_dif  0.0307 0.1742 0.8604       1     -0.3108 0.3721
#> 9   Item9 random_dif -0.1061 0.1656 0.5216       1     -0.4308 0.2185
#> 10 Item10 random_dif -0.1389 0.1652 0.4004       1     -0.4626 0.1849
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2740 0.1581 0.0832  0.8317     -0.0360 0.5840
#> 2   Item2 random_dif -0.0959 0.1734 0.5802  1.0000     -0.4357 0.2440
#> 3   Item3 random_dif  0.0621 0.1631 0.7035  1.0000     -0.2576 0.3817
#> 4   Item4 random_dif -0.0581 0.1754 0.7404  1.0000     -0.4018 0.2856
#> 5   Item5 random_dif  0.1892 0.1669 0.2570  1.0000     -0.1379 0.5163
#> 6   Item6 random_dif -0.0111 0.1722 0.9486  1.0000     -0.3487 0.3265
#> 7   Item7 random_dif  0.0221 0.1780 0.9010  1.0000     -0.3267 0.3709
#> 8   Item8 random_dif -0.2332 0.1658 0.1597  1.0000     -0.5581 0.0918
#> 9   Item9 random_dif  0.0956 0.1682 0.5697  1.0000     -0.2340 0.4252
#> 10 Item10 random_dif -0.2667 0.1620 0.0998  0.9979     -0.5842 0.0509
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0783 0.1723 0.6497  1.0000     -0.2594 0.4160
#> 2   Item2 random_dif -0.0368 0.1794 0.8377  1.0000     -0.3884 0.3149
#> 3   Item3 random_dif -0.0669 0.1717 0.6968  1.0000     -0.4034 0.2696
#> 4   Item4 random_dif -0.0714 0.1776 0.6876  1.0000     -0.4195 0.2767
#> 5   Item5 random_dif  0.0303 0.1703 0.8587  1.0000     -0.3034 0.3640
#> 6   Item6 random_dif  0.0000 0.1775 1.0000  1.0000     -0.3478 0.3478
#> 7   Item7 random_dif  0.0730 0.1732 0.6734  1.0000     -0.2665 0.4125
#> 8   Item8 random_dif -0.1178 0.1661 0.4784  1.0000     -0.4433 0.2078
#> 9   Item9 random_dif  0.3243 0.1666 0.0516  0.5159     -0.0022 0.6508
#> 10 Item10 random_dif -0.2598 0.1858 0.1621  1.0000     -0.6239 0.1044
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0051 0.1685 0.9761       1     -0.3354 0.3253
#> 2   Item2 random_dif  0.0461 0.1692 0.7852       1     -0.2855 0.3777
#> 3   Item3 random_dif  0.0033 0.1643 0.9840       1     -0.3188 0.3254
#> 4   Item4 random_dif  0.0870 0.1634 0.5942       1     -0.2331 0.4072
#> 5   Item5 random_dif -0.1208 0.1612 0.4539       1     -0.4368 0.1953
#> 6   Item6 random_dif -0.1582 0.1632 0.3322       1     -0.4781 0.1616
#> 7   Item7 random_dif  0.1627 0.1639 0.3209       1     -0.1585 0.4840
#> 8   Item8 random_dif  0.0035 0.1731 0.9840       1     -0.3358 0.3428
#> 9   Item9 random_dif -0.1649 0.1684 0.3273       1     -0.4949 0.1651
#> 10 Item10 random_dif  0.1391 0.1631 0.3937       1     -0.1805 0.4586
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.2076 0.1636 0.2045  1.0000     -0.5283  0.1131
#> 2   Item2 random_dif  0.0968 0.1611 0.5481  1.0000     -0.2190  0.4125
#> 3   Item3 random_dif -0.0732 0.1729 0.6721  1.0000     -0.4120  0.2657
#> 4   Item4 random_dif -0.0016 0.1683 0.9925  1.0000     -0.3314  0.3282
#> 5   Item5 random_dif  0.2257 0.1579 0.1528  1.0000     -0.0837  0.5351
#> 6   Item6 random_dif -0.0793 0.1714 0.6437  1.0000     -0.4151  0.2566
#> 7   Item7 random_dif  0.1889 0.1671 0.2583  1.0000     -0.1386  0.5164
#> 8   Item8 random_dif  0.1339 0.1697 0.4301  1.0000     -0.1987  0.4664
#> 9   Item9 random_dif -0.3855 0.1395 0.0057  0.0572   . -0.6589 -0.1121
#> 10 Item10 random_dif  0.1246 0.1677 0.4574  1.0000     -0.2040  0.4532
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0794 0.1659 0.6322       1     -0.4046 0.2458
#> 2   Item2 random_dif  0.1679 0.1557 0.2807       1     -0.1372 0.4730
#> 3   Item3 random_dif  0.0915 0.1656 0.5807       1     -0.2332 0.4161
#> 4   Item4 random_dif  0.1155 0.1686 0.4932       1     -0.2149 0.4459
#> 5   Item5 random_dif  0.0312 0.1686 0.8532       1     -0.2992 0.3616
#> 6   Item6 random_dif  0.1131 0.1622 0.4857       1     -0.2049 0.4311
#> 7   Item7 random_dif -0.1691 0.1717 0.3247       1     -0.5055 0.1674
#> 8   Item8 random_dif  0.0108 0.1683 0.9491       1     -0.3191 0.3406
#> 9   Item9 random_dif -0.0701 0.1671 0.6746       1     -0.3977 0.2574
#> 10 Item10 random_dif -0.2626 0.1604 0.1015       1     -0.5770 0.0517
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1057 0.1739 0.5434   1.000     -0.2352  0.4466
#> 2   Item2 random_dif  0.1032 0.1614 0.5225   1.000     -0.2131  0.4195
#> 3   Item3 random_dif -0.3382 0.1569 0.0311   0.311     -0.6457 -0.0307
#> 4   Item4 random_dif  0.0315 0.1709 0.8535   1.000     -0.3033  0.3664
#> 5   Item5 random_dif -0.0957 0.1650 0.5621   1.000     -0.4191  0.2278
#> 6   Item6 random_dif  0.0672 0.1744 0.7001   1.000     -0.2747  0.4091
#> 7   Item7 random_dif  0.0088 0.1692 0.9584   1.000     -0.3228  0.3404
#> 8   Item8 random_dif -0.0281 0.1825 0.8777   1.000     -0.3857  0.3296
#> 9   Item9 random_dif  0.1090 0.1730 0.5286   1.000     -0.2300  0.4480
#> 10 Item10 random_dif  0.0549 0.1674 0.7427   1.000     -0.2732  0.3831
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2590 0.1543 0.0933  0.9330     -0.5615 0.0435
#> 2   Item2 random_dif  0.0458 0.1801 0.7993  1.0000     -0.3072 0.3988
#> 3   Item3 random_dif  0.0754 0.1653 0.6482  1.0000     -0.2486 0.3995
#> 4   Item4 random_dif  0.2308 0.1589 0.1464  1.0000     -0.0807 0.5422
#> 5   Item5 random_dif -0.0155 0.1711 0.9276  1.0000     -0.3510 0.3199
#> 6   Item6 random_dif -0.1376 0.1664 0.4083  1.0000     -0.4639 0.1886
#> 7   Item7 random_dif -0.1510 0.1686 0.3704  1.0000     -0.4815 0.1794
#> 8   Item8 random_dif  0.3626 0.1489 0.0148  0.1485      0.0709 0.6544
#> 9   Item9 random_dif -0.2441 0.1628 0.1337  1.0000     -0.5632 0.0749
#> 10 Item10 random_dif  0.0764 0.1758 0.6640  1.0000     -0.2681 0.4208
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1381 0.1968 0.4828  1.0000     -0.5237 0.2475
#> 2   Item2 random_dif -0.0307 0.1667 0.8538  1.0000     -0.3574 0.2959
#> 3   Item3 random_dif  0.1150 0.1641 0.4834  1.0000     -0.2066 0.4365
#> 4   Item4 random_dif -0.1645 0.1707 0.3351  1.0000     -0.4991 0.1700
#> 5   Item5 random_dif -0.0201 0.1724 0.9071  1.0000     -0.3579 0.3177
#> 6   Item6 random_dif  0.2954 0.1585 0.0624  0.6241     -0.0153 0.6061
#> 7   Item7 random_dif -0.2138 0.1583 0.1770  1.0000     -0.5241 0.0966
#> 8   Item8 random_dif  0.0415 0.1750 0.8125  1.0000     -0.3014 0.3844
#> 9   Item9 random_dif -0.1123 0.1691 0.5066  1.0000     -0.4437 0.2191
#> 10 Item10 random_dif  0.1986 0.1633 0.2240  1.0000     -0.1215 0.5187
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0511 0.1783 0.7744       1     -0.4005 0.2983
#> 2   Item2 random_dif -0.0204 0.1691 0.9040       1     -0.3519 0.3111
#> 3   Item3 random_dif -0.1611 0.1692 0.3410       1     -0.4926 0.1705
#> 4   Item4 random_dif  0.0870 0.1670 0.6025       1     -0.2403 0.4142
#> 5   Item5 random_dif -0.0293 0.1764 0.8681       1     -0.3750 0.3164
#> 6   Item6 random_dif -0.0761 0.1728 0.6596       1     -0.4147 0.2625
#> 7   Item7 random_dif -0.0849 0.1641 0.6049       1     -0.4065 0.2367
#> 8   Item8 random_dif  0.1978 0.1753 0.2592       1     -0.1458 0.5413
#> 9   Item9 random_dif  0.0049 0.1685 0.9769       1     -0.3254 0.3351
#> 10 Item10 random_dif  0.1477 0.1673 0.3776       1     -0.1803 0.4756
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.2377 0.1604 0.1384  1.0000     -0.5522  0.0767
#> 2   Item2 random_dif  0.0365 0.1642 0.8240  1.0000     -0.2853  0.3583
#> 3   Item3 random_dif  0.2207 0.1573 0.1607  1.0000     -0.0876  0.5290
#> 4   Item4 random_dif -0.0843 0.1742 0.6284  1.0000     -0.4256  0.2571
#> 5   Item5 random_dif  0.3165 0.1535 0.0392  0.3924      0.0156  0.6173
#> 6   Item6 random_dif -0.1815 0.1621 0.2628  1.0000     -0.4991  0.1362
#> 7   Item7 random_dif  0.0812 0.1738 0.6404  1.0000     -0.2594  0.4218
#> 8   Item8 random_dif -0.3626 0.1495 0.0153  0.1528     -0.6557 -0.0696
#> 9   Item9 random_dif  0.1839 0.1619 0.2559  1.0000     -0.1333  0.5012
#> 10 Item10 random_dif -0.0084 0.1729 0.9611  1.0000     -0.3473  0.3305
cutoff_res$item_cutoffs
#>      Item  gamma_low gamma_high
#> 1   Item1 -0.4190476  0.3918869
#> 2   Item2 -0.2969502  0.3341646
#> 3   Item3 -0.4864595  0.4043210
#> 4   Item4 -0.3954457  0.3377265
#> 5   Item5 -0.3787375  0.3684211
#> 6   Item6 -0.3541153  0.3467337
#> 7   Item7 -0.2727273  0.3906899
#> 8   Item8 -0.3248998  0.3626374
#> 9   Item9 -0.4188101  0.3724247
#> 10 Item10 -0.3540313  0.3863216
# }
```
