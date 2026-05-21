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
RMpgDIFcutoff(
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

1.  Resamples person parameters (thetas) with replacement from ML
    estimates.

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

Supports both **dichotomous** data (via
[`eRm::RM()`](https://rdrr.io/pkg/eRm/man/RM.html) and
[`psychotools::rrm()`](https://rdrr.io/pkg/psychotools/man/rrm.html))
and **polytomous** data (via
[`eRm::PCM()`](https://rdrr.io/pkg/eRm/man/PCM.html) and an internal
partial credit score simulator).

Parallel processing is provided by the `mirai` package (optional).
Install it with `install.packages("mirai")` to enable parallelisation.

The `iarm` package must be installed (it is in Suggests, not Imports).

## References

Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
(1998). Differential item functioning in the Danish translation of the
SF-36. *Journal of Clinical Epidemiology*, 51(11), 1189–1202.

Henninger, M., Radek, J., Sengewald, M.-A., & Strobl, C. (2024). Partial
credit trees meet the partial gamma coefficient for quantifying DIF and
DSF in polytomous items. OSF Preprints.
[doi:10.31234/osf.io/47sah](https://doi.org/10.31234/osf.io/47sah)

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
cutoff_res <- RMpgDIFcutoff(sim_data, dif_var = dif_sex,
                                  iterations = 100, parallel = FALSE,
                                  seed = 42)
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0247 0.1819 0.8921  1.0000     -0.3811 0.3318
#> 2   Item2 random_dif -0.1753 0.1628 0.2816  1.0000     -0.4944 0.1438
#> 3   Item3 random_dif  0.1813 0.1680 0.2806  1.0000     -0.1480 0.5105
#> 4   Item4 random_dif -0.1566 0.1810 0.3870  1.0000     -0.5114 0.1982
#> 5   Item5 random_dif -0.0747 0.1724 0.6645  1.0000     -0.4126 0.2631
#> 6   Item6 random_dif  0.3052 0.1605 0.0572  0.5723     -0.0094 0.6197
#> 7   Item7 random_dif  0.1742 0.1672 0.2973  1.0000     -0.1534 0.5019
#> 8   Item8 random_dif -0.1226 0.1733 0.4793  1.0000     -0.4623 0.2171
#> 9   Item9 random_dif -0.0496 0.1747 0.7766  1.0000     -0.3919 0.2928
#> 10 Item10 random_dif -0.0775 0.1717 0.6517  1.0000     -0.4140 0.2590
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1070 0.1758 0.5427  1.0000     -0.4514  0.2375
#> 2   Item2 random_dif -0.0604 0.1760 0.7314  1.0000     -0.4054  0.2845
#> 3   Item3 random_dif -0.4021 0.1522 0.0082  0.0824   . -0.7004 -0.1038
#> 4   Item4 random_dif  0.5154 0.1493 0.0006  0.0056  **  0.2228  0.8080
#> 5   Item5 random_dif -0.0520 0.1697 0.7593  1.0000     -0.3845  0.2806
#> 6   Item6 random_dif -0.2083 0.1772 0.2397  1.0000     -0.5556  0.1390
#> 7   Item7 random_dif  0.4915 0.1313 0.0002  0.0018  **  0.2342  0.7488
#> 8   Item8 random_dif  0.1761 0.1669 0.2913  1.0000     -0.1510  0.5033
#> 9   Item9 random_dif -0.2184 0.1624 0.1785  1.0000     -0.5366  0.0998
#> 10 Item10 random_dif -0.0769 0.1710 0.6529  1.0000     -0.4121  0.2583
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.4376 0.1380 0.0015  0.0153   * -0.7082 -0.1670
#> 2   Item2 random_dif -0.0025 0.1644 0.9879  1.0000     -0.3248  0.3198
#> 3   Item3 random_dif -0.1741 0.1667 0.2963  1.0000     -0.5009  0.1526
#> 4   Item4 random_dif  0.0300 0.1725 0.8621  1.0000     -0.3082  0.3681
#> 5   Item5 random_dif  0.0761 0.1699 0.6542  1.0000     -0.2569  0.4091
#> 6   Item6 random_dif -0.1062 0.1731 0.5395  1.0000     -0.4454  0.2330
#> 7   Item7 random_dif  0.1807 0.1640 0.2706  1.0000     -0.1407  0.5021
#> 8   Item8 random_dif  0.1832 0.1647 0.2659  1.0000     -0.1396  0.5059
#> 9   Item9 random_dif  0.0907 0.1646 0.5817  1.0000     -0.2319  0.4132
#> 10 Item10 random_dif  0.1641 0.1619 0.3107  1.0000     -0.1532  0.4815
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0567 0.1751 0.7460       1     -0.4000 0.2865
#> 2   Item2 random_dif  0.1634 0.1598 0.3064       1     -0.1497 0.4765
#> 3   Item3 random_dif -0.0317 0.1772 0.8578       1     -0.3791 0.3156
#> 4   Item4 random_dif  0.0104 0.1843 0.9549       1     -0.3507 0.3716
#> 5   Item5 random_dif -0.0106 0.1719 0.9510       1     -0.3475 0.3264
#> 6   Item6 random_dif -0.0610 0.1720 0.7227       1     -0.3981 0.2760
#> 7   Item7 random_dif -0.0802 0.1732 0.6431       1     -0.4197 0.2592
#> 8   Item8 random_dif  0.0347 0.1777 0.8452       1     -0.3136 0.3831
#> 9   Item9 random_dif -0.0221 0.1643 0.8930       1     -0.3441 0.2999
#> 10 Item10 random_dif  0.0230 0.1657 0.8896       1     -0.3018 0.3478
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1098 0.1678 0.5127  1.0000     -0.2190 0.4387
#> 2   Item2 random_dif -0.1709 0.1619 0.2912  1.0000     -0.4881 0.1464
#> 3   Item3 random_dif  0.2038 0.1591 0.2004  1.0000     -0.1081 0.5157
#> 4   Item4 random_dif  0.0894 0.1663 0.5907  1.0000     -0.2364 0.4153
#> 5   Item5 random_dif -0.2580 0.1592 0.1050  1.0000     -0.5700 0.0539
#> 6   Item6 random_dif -0.1075 0.1681 0.5224  1.0000     -0.4369 0.2219
#> 7   Item7 random_dif  0.0921 0.1654 0.5775  1.0000     -0.2320 0.4162
#> 8   Item8 random_dif -0.1183 0.1716 0.4905  1.0000     -0.4547 0.2180
#> 9   Item9 random_dif  0.2685 0.1525 0.0782  0.7825     -0.0303 0.5674
#> 10 Item10 random_dif -0.1285 0.1599 0.4214  1.0000     -0.4419 0.1848
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1590 0.1870 0.3952  1.0000     -0.2075 0.5254
#> 2   Item2 random_dif -0.0126 0.1752 0.9426  1.0000     -0.3559 0.3307
#> 3   Item3 random_dif  0.0301 0.1741 0.8626  1.0000     -0.3111 0.3714
#> 4   Item4 random_dif -0.3333 0.1722 0.0529  0.5287     -0.6708 0.0041
#> 5   Item5 random_dif  0.1214 0.1736 0.4845  1.0000     -0.2189 0.4616
#> 6   Item6 random_dif -0.2497 0.1694 0.1405  1.0000     -0.5816 0.0823
#> 7   Item7 random_dif  0.1584 0.1730 0.3599  1.0000     -0.1807 0.4974
#> 8   Item8 random_dif  0.1603 0.1698 0.3451  1.0000     -0.1724 0.4930
#> 9   Item9 random_dif -0.1350 0.1731 0.4355  1.0000     -0.4742 0.2042
#> 10 Item10 random_dif  0.0996 0.1693 0.5563  1.0000     -0.2322 0.4314
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0365 0.1728 0.8330  1.0000     -0.3752 0.3023
#> 2   Item2 random_dif  0.2139 0.1717 0.2127  1.0000     -0.1225 0.5504
#> 3   Item3 random_dif -0.0547 0.1763 0.7563  1.0000     -0.4003 0.2909
#> 4   Item4 random_dif  0.0673 0.1740 0.6989  1.0000     -0.2737 0.4083
#> 5   Item5 random_dif  0.2909 0.1579 0.0655  0.6548     -0.0186 0.6005
#> 6   Item6 random_dif  0.0853 0.1735 0.6231  1.0000     -0.2548 0.4254
#> 7   Item7 random_dif -0.1724 0.1637 0.2922  1.0000     -0.4933 0.1484
#> 8   Item8 random_dif -0.1903 0.1597 0.2332  1.0000     -0.5033 0.1226
#> 9   Item9 random_dif -0.1172 0.1662 0.4806  1.0000     -0.4430 0.2085
#> 10 Item10 random_dif -0.0458 0.1636 0.7795  1.0000     -0.3665 0.2749
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0634 0.1884 0.7364       1     -0.4327 0.3058
#> 2   Item2 random_dif -0.0308 0.1715 0.8573       1     -0.3669 0.3052
#> 3   Item3 random_dif  0.2064 0.1732 0.2334       1     -0.1331 0.5459
#> 4   Item4 random_dif  0.1513 0.1646 0.3578       1     -0.1712 0.4739
#> 5   Item5 random_dif -0.1060 0.1681 0.5282       1     -0.4354 0.2234
#> 6   Item6 random_dif -0.1235 0.1671 0.4597       1     -0.4510 0.2040
#> 7   Item7 random_dif -0.1433 0.1746 0.4119       1     -0.4856 0.1990
#> 8   Item8 random_dif  0.1952 0.1688 0.2475       1     -0.1357 0.5260
#> 9   Item9 random_dif -0.1922 0.1664 0.2481       1     -0.5183 0.1339
#> 10 Item10 random_dif  0.0944 0.1672 0.5725       1     -0.2334 0.4222
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0031 0.1788 0.9862       1     -0.3534 0.3473
#> 2   Item2 random_dif  0.0266 0.1655 0.8722       1     -0.2978 0.3510
#> 3   Item3 random_dif  0.1272 0.1721 0.4600       1     -0.2102 0.4645
#> 4   Item4 random_dif -0.0218 0.1629 0.8934       1     -0.3411 0.2974
#> 5   Item5 random_dif -0.1776 0.1686 0.2923       1     -0.5081 0.1529
#> 6   Item6 random_dif  0.2483 0.1583 0.1167       1     -0.0619 0.5586
#> 7   Item7 random_dif  0.0000 0.1695 1.0000       1     -0.3321 0.3321
#> 8   Item8 random_dif -0.1682 0.1642 0.3057       1     -0.4901 0.1537
#> 9   Item9 random_dif  0.0517 0.1684 0.7586       1     -0.2782 0.3817
#> 10 Item10 random_dif -0.0709 0.1647 0.6670       1     -0.3937 0.2519
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0791 0.1731 0.6476  1.0000     -0.2602 0.4185
#> 2   Item2 random_dif -0.0627 0.1627 0.7002  1.0000     -0.3816 0.2563
#> 3   Item3 random_dif -0.1643 0.1569 0.2949  1.0000     -0.4717 0.1431
#> 4   Item4 random_dif  0.3130 0.1519 0.0393  0.3929      0.0154 0.6106
#> 5   Item5 random_dif -0.0400 0.1651 0.8086  1.0000     -0.3636 0.2836
#> 6   Item6 random_dif -0.0115 0.1702 0.9460  1.0000     -0.3452 0.3221
#> 7   Item7 random_dif -0.1975 0.1698 0.2448  1.0000     -0.5304 0.1353
#> 8   Item8 random_dif -0.0060 0.1741 0.9726  1.0000     -0.3472 0.3353
#> 9   Item9 random_dif  0.0759 0.1576 0.6299  1.0000     -0.2330 0.3849
#> 10 Item10 random_dif -0.0094 0.1674 0.9550  1.0000     -0.3376 0.3187
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0205 0.1816 0.9102  1.0000     -0.3764 0.3354
#> 2   Item2 random_dif  0.3299 0.1509 0.0289  0.2887      0.0340 0.6257
#> 3   Item3 random_dif  0.0519 0.1693 0.7591  1.0000     -0.2798 0.3836
#> 4   Item4 random_dif  0.1699 0.1811 0.3480  1.0000     -0.1850 0.5248
#> 5   Item5 random_dif -0.0056 0.1719 0.9740  1.0000     -0.3425 0.3313
#> 6   Item6 random_dif -0.2142 0.1667 0.1988  1.0000     -0.5409 0.1125
#> 7   Item7 random_dif  0.0208 0.1692 0.9021  1.0000     -0.3107 0.3523
#> 8   Item8 random_dif  0.0934 0.1756 0.5947  1.0000     -0.2507 0.4375
#> 9   Item9 random_dif -0.3107 0.1728 0.0722  0.7224     -0.6494 0.0281
#> 10 Item10 random_dif -0.1421 0.1593 0.3723  1.0000     -0.4544 0.1701
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1672 0.1671 0.3173  1.0000     -0.4947 0.1604
#> 2   Item2 random_dif -0.2917 0.1494 0.0509  0.5086     -0.5845 0.0011
#> 3   Item3 random_dif  0.0493 0.1771 0.7808  1.0000     -0.2978 0.3964
#> 4   Item4 random_dif  0.2619 0.1659 0.1143  1.0000     -0.0632 0.5871
#> 5   Item5 random_dif  0.2822 0.1587 0.0753  0.7534     -0.0288 0.5932
#> 6   Item6 random_dif  0.1796 0.1655 0.2778  1.0000     -0.1447 0.5038
#> 7   Item7 random_dif  0.0720 0.1591 0.6510  1.0000     -0.2399 0.3839
#> 8   Item8 random_dif -0.1948 0.1637 0.2341  1.0000     -0.5157 0.1261
#> 9   Item9 random_dif  0.0938 0.1732 0.5879  1.0000     -0.2456 0.4333
#> 10 Item10 random_dif -0.2231 0.1589 0.1602  1.0000     -0.5345 0.0883
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0628 0.1823 0.7304  1.0000     -0.2945 0.4202
#> 2   Item2 random_dif  0.0125 0.1740 0.9426  1.0000     -0.3286 0.3536
#> 3   Item3 random_dif  0.3577 0.1544 0.0205  0.2052      0.0551 0.6603
#> 4   Item4 random_dif -0.0724 0.1717 0.6732  1.0000     -0.4089 0.2641
#> 5   Item5 random_dif -0.1650 0.1668 0.3227  1.0000     -0.4919 0.1620
#> 6   Item6 random_dif -0.1189 0.1635 0.4671  1.0000     -0.4393 0.2015
#> 7   Item7 random_dif -0.1777 0.1690 0.2929  1.0000     -0.5089 0.1535
#> 8   Item8 random_dif -0.1923 0.1657 0.2460  1.0000     -0.5170 0.1325
#> 9   Item9 random_dif  0.2607 0.1592 0.1015  1.0000     -0.0514 0.5728
#> 10 Item10 random_dif  0.0376 0.1635 0.8183  1.0000     -0.2829 0.3580
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2092 0.1841 0.2560   1.000     -0.5700 0.1517
#> 2   Item2 random_dif -0.1074 0.1682 0.5231   1.000     -0.4370 0.2222
#> 3   Item3 random_dif -0.0638 0.1677 0.7037   1.000     -0.3924 0.2649
#> 4   Item4 random_dif  0.2630 0.1620 0.1044   1.000     -0.0544 0.5804
#> 5   Item5 random_dif -0.1304 0.1676 0.4365   1.000     -0.4590 0.1981
#> 6   Item6 random_dif  0.3994 0.1548 0.0099   0.099   .  0.0959 0.7028
#> 7   Item7 random_dif  0.1756 0.1659 0.2898   1.000     -0.1495 0.5007
#> 8   Item8 random_dif -0.1844 0.1639 0.2606   1.000     -0.5057 0.1369
#> 9   Item9 random_dif -0.0046 0.1736 0.9787   1.000     -0.3449 0.3356
#> 10 Item10 random_dif -0.1429 0.1634 0.3820   1.000     -0.4631 0.1774
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1009 0.1642 0.5390  1.0000     -0.2210 0.4228
#> 2   Item2 random_dif  0.1056 0.1702 0.5350  1.0000     -0.2280 0.4392
#> 3   Item3 random_dif  0.0413 0.1767 0.8151  1.0000     -0.3049 0.3876
#> 4   Item4 random_dif -0.1636 0.1767 0.3543  1.0000     -0.5099 0.1826
#> 5   Item5 random_dif -0.1502 0.1694 0.3754  1.0000     -0.4822 0.1818
#> 6   Item6 random_dif -0.0165 0.1851 0.9292  1.0000     -0.3793 0.3464
#> 7   Item7 random_dif  0.0509 0.1780 0.7749  1.0000     -0.2979 0.3997
#> 8   Item8 random_dif  0.2763 0.1572 0.0788  0.7875     -0.0317 0.5843
#> 9   Item9 random_dif -0.0602 0.1693 0.7220  1.0000     -0.3921 0.2716
#> 10 Item10 random_dif -0.2036 0.1619 0.2086  1.0000     -0.5210 0.1137
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2168 0.1670 0.1942       1     -0.1105 0.5442
#> 2   Item2 random_dif -0.2345 0.1634 0.1512       1     -0.5547 0.0858
#> 3   Item3 random_dif  0.0554 0.1662 0.7390       1     -0.2704 0.3812
#> 4   Item4 random_dif  0.0131 0.1697 0.9385       1     -0.3194 0.3456
#> 5   Item5 random_dif  0.0590 0.1622 0.7160       1     -0.2588 0.3768
#> 6   Item6 random_dif  0.1939 0.1644 0.2381       1     -0.1282 0.5161
#> 7   Item7 random_dif  0.0696 0.1783 0.6963       1     -0.2799 0.4192
#> 8   Item8 random_dif -0.2123 0.1729 0.2194       1     -0.5511 0.1265
#> 9   Item9 random_dif -0.1971 0.1668 0.2376       1     -0.5241 0.1299
#> 10 Item10 random_dif  0.0316 0.1599 0.8432       1     -0.2818 0.3450
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0919 0.1695 0.5879  1.0000     -0.4241 0.2404
#> 2   Item2 random_dif  0.1242 0.1673 0.4577  1.0000     -0.2036 0.4520
#> 3   Item3 random_dif  0.4552 0.1426 0.0014  0.0141   *  0.1757 0.7348
#> 4   Item4 random_dif -0.0465 0.1627 0.7750  1.0000     -0.3653 0.2723
#> 5   Item5 random_dif -0.0165 0.1615 0.9186  1.0000     -0.3330 0.3000
#> 6   Item6 random_dif -0.0409 0.1760 0.8164  1.0000     -0.3859 0.3042
#> 7   Item7 random_dif -0.1041 0.1665 0.5316  1.0000     -0.4305 0.2222
#> 8   Item8 random_dif -0.1076 0.1651 0.5147  1.0000     -0.4312 0.2161
#> 9   Item9 random_dif -0.0046 0.1617 0.9772  1.0000     -0.3216 0.3124
#> 10 Item10 random_dif -0.1494 0.1564 0.3394  1.0000     -0.4560 0.1571
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0935 0.1674 0.5766  1.0000     -0.4216 0.2347
#> 2   Item2 random_dif -0.0630 0.1709 0.7126  1.0000     -0.3980 0.2721
#> 3   Item3 random_dif  0.0078 0.1626 0.9620  1.0000     -0.3110 0.3265
#> 4   Item4 random_dif -0.2368 0.1611 0.1416  1.0000     -0.5524 0.0789
#> 5   Item5 random_dif  0.0187 0.1695 0.9124  1.0000     -0.3135 0.3508
#> 6   Item6 random_dif  0.0349 0.1784 0.8448  1.0000     -0.3147 0.3846
#> 7   Item7 random_dif -0.0260 0.1736 0.8809  1.0000     -0.3664 0.3143
#> 8   Item8 random_dif  0.3450 0.1436 0.0162  0.1624      0.0637 0.6264
#> 9   Item9 random_dif -0.0165 0.1694 0.9225  1.0000     -0.3484 0.3155
#> 10 Item10 random_dif -0.0093 0.1630 0.9544  1.0000     -0.3287 0.3101
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0265 0.1781 0.8815       1     -0.3756 0.3225
#> 2   Item2 random_dif  0.0190 0.1608 0.9059       1     -0.2961 0.3341
#> 3   Item3 random_dif -0.0645 0.1667 0.6990       1     -0.3913 0.2623
#> 4   Item4 random_dif  0.0628 0.1761 0.7212       1     -0.2823 0.4080
#> 5   Item5 random_dif  0.0483 0.1659 0.7709       1     -0.2769 0.3736
#> 6   Item6 random_dif -0.2052 0.1663 0.2173       1     -0.5312 0.1208
#> 7   Item7 random_dif  0.0486 0.1717 0.7770       1     -0.2879 0.3851
#> 8   Item8 random_dif  0.1052 0.1627 0.5179       1     -0.2137 0.4241
#> 9   Item9 random_dif -0.1400 0.1651 0.3963       1     -0.4636 0.1835
#> 10 Item10 random_dif  0.1373 0.1617 0.3958       1     -0.1796 0.4542
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1248 0.1799 0.4877  1.0000     -0.2277  0.4773
#> 2   Item2 random_dif  0.0970 0.1660 0.5591  1.0000     -0.2283  0.4222
#> 3   Item3 random_dif -0.3028 0.1505 0.0442  0.4421     -0.5978 -0.0078
#> 4   Item4 random_dif -0.0687 0.1763 0.6968  1.0000     -0.4144  0.2769
#> 5   Item5 random_dif -0.1178 0.1661 0.4782  1.0000     -0.4434  0.2078
#> 6   Item6 random_dif -0.2612 0.1617 0.1063  1.0000     -0.5781  0.0557
#> 7   Item7 random_dif -0.2430 0.1665 0.1445  1.0000     -0.5694  0.0834
#> 8   Item8 random_dif  0.2778 0.1596 0.0818  0.8182     -0.0351  0.5906
#> 9   Item9 random_dif  0.2071 0.1562 0.1849  1.0000     -0.0991  0.5133
#> 10 Item10 random_dif  0.2691 0.1529 0.0785  0.7848     -0.0306  0.5688
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1449 0.1810 0.4233       1     -0.2098 0.4996
#> 2   Item2 random_dif  0.0283 0.1655 0.8644       1     -0.2962 0.3527
#> 3   Item3 random_dif  0.0460 0.1625 0.7769       1     -0.2724 0.3645
#> 4   Item4 random_dif -0.1098 0.1671 0.5113       1     -0.4374 0.2178
#> 5   Item5 random_dif -0.0244 0.1628 0.8809       1     -0.3434 0.2946
#> 6   Item6 random_dif -0.0645 0.1608 0.6882       1     -0.3797 0.2506
#> 7   Item7 random_dif  0.1258 0.1721 0.4649       1     -0.2116 0.4632
#> 8   Item8 random_dif -0.1182 0.1628 0.4679       1     -0.4372 0.2009
#> 9   Item9 random_dif  0.0096 0.1835 0.9582       1     -0.3500 0.3692
#> 10 Item10 random_dif  0.0136 0.1637 0.9337       1     -0.3072 0.3344
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1801 0.1683 0.2844  1.0000     -0.1496  0.5099
#> 2   Item2 random_dif -0.0638 0.1752 0.7159  1.0000     -0.4072  0.2796
#> 3   Item3 random_dif -0.3828 0.1608 0.0173  0.1728     -0.6979 -0.0677
#> 4   Item4 random_dif  0.0291 0.1802 0.8717  1.0000     -0.3241  0.3823
#> 5   Item5 random_dif -0.3149 0.1554 0.0427  0.4271     -0.6195 -0.0104
#> 6   Item6 random_dif  0.2524 0.1613 0.1178  1.0000     -0.0638  0.5686
#> 7   Item7 random_dif  0.3137 0.1498 0.0362  0.3621      0.0202  0.6072
#> 8   Item8 random_dif  0.0100 0.1869 0.9573  1.0000     -0.3563  0.3763
#> 9   Item9 random_dif -0.0028 0.1721 0.9871  1.0000     -0.3402  0.3346
#> 10 Item10 random_dif -0.0341 0.1650 0.8361  1.0000     -0.3575  0.2892
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2764 0.1597 0.0835  0.8348     -0.0366 0.5894
#> 2   Item2 random_dif -0.0221 0.1663 0.8944  1.0000     -0.3480 0.3039
#> 3   Item3 random_dif -0.0182 0.1687 0.9142  1.0000     -0.3487 0.3124
#> 4   Item4 random_dif  0.1581 0.1696 0.3512  1.0000     -0.1743 0.4906
#> 5   Item5 random_dif -0.1106 0.1636 0.4990  1.0000     -0.4313 0.2100
#> 6   Item6 random_dif -0.2847 0.1584 0.0723  0.7228     -0.5952 0.0258
#> 7   Item7 random_dif  0.0515 0.1687 0.7602  1.0000     -0.2791 0.3821
#> 8   Item8 random_dif  0.0336 0.1631 0.8369  1.0000     -0.2862 0.3534
#> 9   Item9 random_dif -0.2436 0.1564 0.1193  1.0000     -0.5501 0.0629
#> 10 Item10 random_dif  0.1652 0.1576 0.2944  1.0000     -0.1436 0.4740
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1543 0.1638 0.3464       1     -0.4754 0.1668
#> 2   Item2 random_dif  0.2215 0.1571 0.1584       1     -0.0863 0.5293
#> 3   Item3 random_dif  0.0763 0.1625 0.6385       1     -0.2421 0.3947
#> 4   Item4 random_dif -0.0196 0.1654 0.9059       1     -0.3436 0.3045
#> 5   Item5 random_dif -0.0212 0.1603 0.8950       1     -0.3354 0.2931
#> 6   Item6 random_dif  0.2187 0.1681 0.1932       1     -0.1107 0.5481
#> 7   Item7 random_dif -0.2400 0.1653 0.1466       1     -0.5640 0.0840
#> 8   Item8 random_dif  0.0662 0.1644 0.6874       1     -0.2561 0.3885
#> 9   Item9 random_dif  0.0102 0.1626 0.9502       1     -0.3085 0.3288
#> 10 Item10 random_dif -0.1736 0.1571 0.2693       1     -0.4815 0.1344
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3345 0.1629 0.0401  0.4005     -0.6538 -0.0152
#> 2   Item2 random_dif -0.0201 0.1652 0.9032  1.0000     -0.3439  0.3037
#> 3   Item3 random_dif  0.0665 0.1663 0.6895  1.0000     -0.2596  0.3925
#> 4   Item4 random_dif -0.0276 0.1734 0.8738  1.0000     -0.3675  0.3123
#> 5   Item5 random_dif  0.2622 0.1528 0.0862  0.8623     -0.0373  0.5618
#> 6   Item6 random_dif -0.2658 0.1602 0.0970  0.9703     -0.5797  0.0481
#> 7   Item7 random_dif  0.0363 0.1642 0.8253  1.0000     -0.2856  0.3582
#> 8   Item8 random_dif  0.1682 0.1798 0.3495  1.0000     -0.1842  0.5205
#> 9   Item9 random_dif  0.2685 0.1531 0.0795  0.7949     -0.0316  0.5686
#> 10 Item10 random_dif -0.2214 0.1606 0.1682  1.0000     -0.5362  0.0935
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0461 0.1731 0.7899  1.0000     -0.3853 0.2931
#> 2   Item2 random_dif  0.0026 0.1678 0.9875  1.0000     -0.3263 0.3316
#> 3   Item3 random_dif -0.0828 0.1748 0.6356  1.0000     -0.4255 0.2598
#> 4   Item4 random_dif  0.0565 0.1732 0.7444  1.0000     -0.2831 0.3960
#> 5   Item5 random_dif -0.2677 0.1550 0.0842  0.8422     -0.5716 0.0362
#> 6   Item6 random_dif  0.0122 0.1923 0.9495  1.0000     -0.3648 0.3891
#> 7   Item7 random_dif  0.0516 0.1801 0.7743  1.0000     -0.3014 0.4047
#> 8   Item8 random_dif  0.0753 0.1723 0.6620  1.0000     -0.2625 0.4131
#> 9   Item9 random_dif  0.1681 0.1797 0.3496  1.0000     -0.1841 0.5202
#> 10 Item10 random_dif  0.0819 0.1663 0.6224  1.0000     -0.2440 0.4077
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3496 0.1729 0.0432  0.4318     -0.6885 -0.0107
#> 2   Item2 random_dif -0.3015 0.1545 0.0510  0.5105     -0.6044  0.0014
#> 3   Item3 random_dif  0.0504 0.1799 0.7793  1.0000     -0.3022  0.4030
#> 4   Item4 random_dif -0.2904 0.1745 0.0961  0.9613     -0.6325  0.0517
#> 5   Item5 random_dif -0.0333 0.1739 0.8480  1.0000     -0.3742  0.3075
#> 6   Item6 random_dif  0.2268 0.1612 0.1593  1.0000     -0.0891  0.5427
#> 7   Item7 random_dif  0.3037 0.1737 0.0803  0.8034     -0.0367  0.6442
#> 8   Item8 random_dif  0.3134 0.1547 0.0428  0.4279      0.0102  0.6165
#> 9   Item9 random_dif  0.1514 0.1638 0.3555  1.0000     -0.1697  0.4724
#> 10 Item10 random_dif -0.1634 0.1708 0.3388  1.0000     -0.4981  0.1714
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1111 0.1666 0.5048   1.000     -0.4376 0.2154
#> 2   Item2 random_dif  0.1725 0.1601 0.2813   1.000     -0.1413 0.4864
#> 3   Item3 random_dif  0.3249 0.1652 0.0492   0.492      0.0011 0.6487
#> 4   Item4 random_dif  0.1345 0.1740 0.4395   1.000     -0.2065 0.4756
#> 5   Item5 random_dif -0.0185 0.1729 0.9150   1.000     -0.3573 0.3204
#> 6   Item6 random_dif -0.0684 0.1676 0.6831   1.000     -0.3969 0.2600
#> 7   Item7 random_dif  0.0700 0.1771 0.6924   1.000     -0.2770 0.4171
#> 8   Item8 random_dif -0.0641 0.1748 0.7137   1.000     -0.4068 0.2785
#> 9   Item9 random_dif -0.2468 0.1675 0.1406   1.000     -0.5751 0.0815
#> 10 Item10 random_dif -0.1804 0.1617 0.2644   1.000     -0.4973 0.1364
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1869 0.1621 0.2490  1.0000     -0.5047 0.1309
#> 2   Item2 random_dif  0.0827 0.1563 0.5966  1.0000     -0.2236 0.3891
#> 3   Item3 random_dif -0.1130 0.1679 0.5009  1.0000     -0.4420 0.2160
#> 4   Item4 random_dif -0.2261 0.1519 0.1364  1.0000     -0.5238 0.0715
#> 5   Item5 random_dif -0.0197 0.1706 0.9080  1.0000     -0.3540 0.3146
#> 6   Item6 random_dif -0.0349 0.1654 0.8327  1.0000     -0.3591 0.2892
#> 7   Item7 random_dif  0.4409 0.1372 0.0013  0.0131   *  0.1721 0.7098
#> 8   Item8 random_dif -0.0055 0.1656 0.9735  1.0000     -0.3301 0.3191
#> 9   Item9 random_dif  0.0350 0.1628 0.8298  1.0000     -0.2841 0.3541
#> 10 Item10 random_dif -0.0025 0.1600 0.9873  1.0000     -0.3162 0.3111
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1950 0.1763 0.2686       1     -0.5405 0.1505
#> 2   Item2 random_dif  0.0942 0.1604 0.5569       1     -0.2202 0.4087
#> 3   Item3 random_dif  0.1660 0.1644 0.3128       1     -0.1563 0.4883
#> 4   Item4 random_dif  0.0223 0.1638 0.8918       1     -0.2987 0.3433
#> 5   Item5 random_dif  0.1710 0.1605 0.2868       1     -0.1436 0.4856
#> 6   Item6 random_dif -0.2277 0.1670 0.1728       1     -0.5550 0.0996
#> 7   Item7 random_dif  0.0819 0.1590 0.6064       1     -0.2297 0.3936
#> 8   Item8 random_dif -0.0525 0.1721 0.7604       1     -0.3897 0.2848
#> 9   Item9 random_dif -0.0081 0.1657 0.9611       1     -0.3328 0.3166
#> 10 Item10 random_dif -0.1155 0.1713 0.4999       1     -0.4513 0.2202
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0167 0.1859 0.9286       1     -0.3810 0.3477
#> 2   Item2 random_dif -0.2637 0.1625 0.1045       1     -0.5822 0.0547
#> 3   Item3 random_dif  0.1022 0.1809 0.5719       1     -0.2522 0.4567
#> 4   Item4 random_dif  0.1145 0.1662 0.4909       1     -0.2113 0.4403
#> 5   Item5 random_dif -0.0560 0.1634 0.7319       1     -0.3762 0.2642
#> 6   Item6 random_dif -0.0302 0.1746 0.8626       1     -0.3724 0.3119
#> 7   Item7 random_dif  0.1364 0.1664 0.4124       1     -0.1897 0.4624
#> 8   Item8 random_dif -0.1550 0.1766 0.3801       1     -0.5012 0.1911
#> 9   Item9 random_dif  0.0316 0.1693 0.8520       1     -0.3002 0.3634
#> 10 Item10 random_dif  0.1088 0.1676 0.5162       1     -0.2197 0.4373
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0565 0.1784 0.7515  1.0000     -0.4062 0.2932
#> 2   Item2 random_dif  0.2318 0.1577 0.1417  1.0000     -0.0774 0.5410
#> 3   Item3 random_dif -0.2073 0.1755 0.2376  1.0000     -0.5513 0.1367
#> 4   Item4 random_dif -0.0704 0.1631 0.6659  1.0000     -0.3901 0.2492
#> 5   Item5 random_dif -0.0646 0.1662 0.6977  1.0000     -0.3904 0.2612
#> 6   Item6 random_dif -0.1451 0.1670 0.3850  1.0000     -0.4725 0.1823
#> 7   Item7 random_dif  0.1783 0.1632 0.2747  1.0000     -0.1416 0.4982
#> 8   Item8 random_dif -0.1129 0.1739 0.5161  1.0000     -0.4537 0.2279
#> 9   Item9 random_dif -0.0696 0.1642 0.6717  1.0000     -0.3913 0.2522
#> 10 Item10 random_dif  0.2679 0.1555 0.0848  0.8484     -0.0368 0.5726
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0943 0.1764 0.5929       1     -0.2514 0.4399
#> 2   Item2 random_dif -0.1765 0.1573 0.2618       1     -0.4847 0.1318
#> 3   Item3 random_dif -0.0874 0.1766 0.6207       1     -0.4336 0.2587
#> 4   Item4 random_dif  0.0373 0.1785 0.8347       1     -0.3127 0.3872
#> 5   Item5 random_dif  0.1741 0.1669 0.2966       1     -0.1529 0.5012
#> 6   Item6 random_dif  0.0115 0.1710 0.9466       1     -0.3238 0.3467
#> 7   Item7 random_dif -0.0617 0.1782 0.7294       1     -0.4110 0.2877
#> 8   Item8 random_dif -0.0593 0.1689 0.7254       1     -0.3903 0.2717
#> 9   Item9 random_dif  0.1690 0.1662 0.3091       1     -0.1567 0.4947
#> 10 Item10 random_dif -0.0874 0.1653 0.5970       1     -0.4115 0.2366
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1782 0.1694 0.2930       1     -0.5102 0.1539
#> 2   Item2 random_dif -0.0679 0.1708 0.6908       1     -0.4027 0.2668
#> 3   Item3 random_dif -0.0074 0.1780 0.9667       1     -0.3562 0.3414
#> 4   Item4 random_dif -0.1391 0.1780 0.4346       1     -0.4878 0.2097
#> 5   Item5 random_dif  0.0761 0.1688 0.6522       1     -0.2548 0.4070
#> 6   Item6 random_dif  0.0252 0.1761 0.8861       1     -0.3199 0.3703
#> 7   Item7 random_dif  0.0931 0.1818 0.6085       1     -0.2632 0.4494
#> 8   Item8 random_dif  0.2504 0.1764 0.1558       1     -0.0953 0.5961
#> 9   Item9 random_dif  0.0090 0.1719 0.9582       1     -0.3278 0.3459
#> 10 Item10 random_dif -0.0350 0.1671 0.8340       1     -0.3626 0.2926
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0867 0.1664 0.6025  1.0000     -0.2395 0.4129
#> 2   Item2 random_dif  0.4002 0.1403 0.0043  0.0432   *  0.1253 0.6751
#> 3   Item3 random_dif  0.1125 0.1648 0.4946  1.0000     -0.2104 0.4355
#> 4   Item4 random_dif -0.1796 0.1673 0.2831  1.0000     -0.5076 0.1484
#> 5   Item5 random_dif -0.3043 0.1597 0.0566  0.5665     -0.6173 0.0086
#> 6   Item6 random_dif  0.2563 0.1566 0.1017  1.0000     -0.0506 0.5633
#> 7   Item7 random_dif -0.0343 0.1692 0.8393  1.0000     -0.3659 0.2973
#> 8   Item8 random_dif -0.0867 0.1605 0.5888  1.0000     -0.4013 0.2278
#> 9   Item9 random_dif -0.0923 0.1651 0.5760  1.0000     -0.4158 0.2312
#> 10 Item10 random_dif -0.2010 0.1585 0.2047  1.0000     -0.5116 0.1096
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0395 0.1784 0.8247       1     -0.3892 0.3102
#> 2   Item2 random_dif  0.1157 0.1714 0.4997       1     -0.2203 0.4517
#> 3   Item3 random_dif  0.1187 0.1674 0.4782       1     -0.2093 0.4467
#> 4   Item4 random_dif  0.0570 0.1755 0.7455       1     -0.2870 0.4010
#> 5   Item5 random_dif  0.0522 0.1703 0.7593       1     -0.2816 0.3860
#> 6   Item6 random_dif -0.0490 0.1870 0.7934       1     -0.4154 0.3175
#> 7   Item7 random_dif -0.0767 0.1730 0.6575       1     -0.4158 0.2624
#> 8   Item8 random_dif  0.0220 0.1765 0.9007       1     -0.3239 0.3680
#> 9   Item9 random_dif  0.0111 0.1741 0.9490       1     -0.3301 0.3524
#> 10 Item10 random_dif -0.2118 0.1624 0.1922       1     -0.5301 0.1065
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1528 0.1675 0.3616  1.0000     -0.1755 0.4811
#> 2   Item2 random_dif -0.1848 0.1620 0.2542  1.0000     -0.5023 0.1328
#> 3   Item3 random_dif -0.1878 0.1641 0.2524  1.0000     -0.5095 0.1338
#> 4   Item4 random_dif  0.2673 0.1716 0.1194  1.0000     -0.0691 0.6037
#> 5   Item5 random_dif -0.0128 0.1678 0.9392  1.0000     -0.3418 0.3162
#> 6   Item6 random_dif  0.0793 0.1720 0.6448  1.0000     -0.2578 0.4163
#> 7   Item7 random_dif -0.0015 0.1763 0.9934  1.0000     -0.3469 0.3440
#> 8   Item8 random_dif  0.2278 0.1774 0.1990  1.0000     -0.1198 0.5755
#> 9   Item9 random_dif -0.2747 0.1507 0.0684  0.6837     -0.5702 0.0207
#> 10 Item10 random_dif  0.0261 0.1671 0.8759  1.0000     -0.3014 0.3536
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1910 0.1814 0.2923  1.0000     -0.5465 0.1645
#> 2   Item2 random_dif -0.0493 0.1631 0.7625  1.0000     -0.3690 0.2704
#> 3   Item3 random_dif -0.2070 0.1633 0.2049  1.0000     -0.5270 0.1130
#> 4   Item4 random_dif -0.1634 0.1675 0.3293  1.0000     -0.4918 0.1649
#> 5   Item5 random_dif  0.0888 0.1681 0.5974  1.0000     -0.2407 0.4182
#> 6   Item6 random_dif -0.1834 0.1740 0.2918  1.0000     -0.5245 0.1576
#> 7   Item7 random_dif  0.1202 0.1716 0.4835  1.0000     -0.2161 0.4565
#> 8   Item8 random_dif  0.3296 0.1592 0.0384  0.3837      0.0177 0.6416
#> 9   Item9 random_dif  0.0000 0.1619 1.0000  1.0000     -0.3172 0.3172
#> 10 Item10 random_dif  0.2217 0.1588 0.1627  1.0000     -0.0895 0.5328
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0208 0.1752 0.9056  1.0000     -0.3641 0.3225
#> 2   Item2 random_dif -0.1210 0.1656 0.4650  1.0000     -0.4456 0.2036
#> 3   Item3 random_dif  0.2938 0.1585 0.0637  0.6374     -0.0168 0.6044
#> 4   Item4 random_dif -0.1782 0.1737 0.3051  1.0000     -0.5187 0.1623
#> 5   Item5 random_dif -0.1734 0.1675 0.3006  1.0000     -0.5018 0.1549
#> 6   Item6 random_dif  0.1959 0.1713 0.2526  1.0000     -0.1397 0.5316
#> 7   Item7 random_dif  0.0032 0.1769 0.9857  1.0000     -0.3436 0.3499
#> 8   Item8 random_dif -0.1653 0.1663 0.3204  1.0000     -0.4912 0.1607
#> 9   Item9 random_dif  0.0851 0.1665 0.6092  1.0000     -0.2412 0.4114
#> 10 Item10 random_dif  0.0568 0.1599 0.7225  1.0000     -0.2567 0.3703
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2432 0.1638 0.1378  1.0000     -0.5643 0.0780
#> 2   Item2 random_dif  0.1377 0.1616 0.3942  1.0000     -0.1790 0.4545
#> 3   Item3 random_dif -0.1194 0.1584 0.4510  1.0000     -0.4299 0.1911
#> 4   Item4 random_dif -0.0185 0.1781 0.9172  1.0000     -0.3676 0.3306
#> 5   Item5 random_dif -0.1466 0.1606 0.3614  1.0000     -0.4614 0.1682
#> 6   Item6 random_dif  0.1013 0.1640 0.5369  1.0000     -0.2202 0.4228
#> 7   Item7 random_dif -0.0904 0.1717 0.5984  1.0000     -0.4269 0.2461
#> 8   Item8 random_dif -0.0116 0.1636 0.9435  1.0000     -0.3321 0.3090
#> 9   Item9 random_dif  0.2739 0.1522 0.0719  0.7195     -0.0244 0.5722
#> 10 Item10 random_dif  0.0818 0.1573 0.6031  1.0000     -0.2266 0.3902
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2353 0.1592 0.1393       1     -0.0767 0.5473
#> 2   Item2 random_dif -0.1971 0.1622 0.2241       1     -0.5150 0.1207
#> 3   Item3 random_dif -0.0076 0.1701 0.9645       1     -0.3409 0.3258
#> 4   Item4 random_dif  0.1523 0.1674 0.3628       1     -0.1758 0.4805
#> 5   Item5 random_dif  0.1616 0.1587 0.3086       1     -0.1494 0.4726
#> 6   Item6 random_dif -0.1932 0.1660 0.2445       1     -0.5185 0.1321
#> 7   Item7 random_dif -0.0292 0.1675 0.8614       1     -0.3575 0.2991
#> 8   Item8 random_dif -0.1272 0.1621 0.4325       1     -0.4449 0.1905
#> 9   Item9 random_dif  0.1105 0.1613 0.4934       1     -0.2057 0.4267
#> 10 Item10 random_dif -0.0972 0.1650 0.5559       1     -0.4206 0.2263
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3303 0.1642 0.0443  0.4429     -0.6522 -0.0084
#> 2   Item2 random_dif  0.0075 0.1673 0.9642  1.0000     -0.3204  0.3354
#> 3   Item3 random_dif -0.0569 0.1690 0.7363  1.0000     -0.3882  0.2744
#> 4   Item4 random_dif  0.1238 0.1731 0.4746  1.0000     -0.2155  0.4630
#> 5   Item5 random_dif -0.1638 0.1560 0.2938  1.0000     -0.4695  0.1419
#> 6   Item6 random_dif  0.2670 0.1483 0.0719  0.7186     -0.0237  0.5577
#> 7   Item7 random_dif  0.0744 0.1617 0.6453  1.0000     -0.2425  0.3914
#> 8   Item8 random_dif -0.0426 0.1675 0.7994  1.0000     -0.3708  0.2857
#> 9   Item9 random_dif  0.0041 0.1711 0.9808  1.0000     -0.3313  0.3395
#> 10 Item10 random_dif  0.0619 0.1615 0.7015  1.0000     -0.2546  0.3785
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2040 0.1687 0.2267       1     -0.1267 0.5347
#> 2   Item2 random_dif -0.0243 0.1682 0.8854       1     -0.3539 0.3054
#> 3   Item3 random_dif  0.0361 0.1720 0.8336       1     -0.3010 0.3733
#> 4   Item4 random_dif -0.2066 0.1636 0.2065       1     -0.5272 0.1140
#> 5   Item5 random_dif -0.0550 0.1698 0.7458       1     -0.3879 0.2778
#> 6   Item6 random_dif  0.1685 0.1856 0.3642       1     -0.1954 0.5323
#> 7   Item7 random_dif -0.0171 0.1700 0.9200       1     -0.3503 0.3161
#> 8   Item8 random_dif -0.2718 0.1678 0.1052       1     -0.6007 0.0570
#> 9   Item9 random_dif  0.1731 0.1594 0.2777       1     -0.1394 0.4856
#> 10 Item10 random_dif -0.0286 0.1668 0.8640       1     -0.3555 0.2984
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1610 0.1675 0.3364       1     -0.1673 0.4893
#> 2   Item2 random_dif -0.0397 0.1627 0.8073       1     -0.3586 0.2792
#> 3   Item3 random_dif  0.1012 0.1616 0.5310       1     -0.2154 0.4179
#> 4   Item4 random_dif  0.0878 0.1651 0.5950       1     -0.2358 0.4114
#> 5   Item5 random_dif -0.0993 0.1600 0.5351       1     -0.4129 0.2144
#> 6   Item6 random_dif -0.0994 0.1645 0.5458       1     -0.4217 0.2230
#> 7   Item7 random_dif -0.0279 0.1806 0.8772       1     -0.3819 0.3261
#> 8   Item8 random_dif -0.1517 0.1551 0.3282       1     -0.4556 0.1523
#> 9   Item9 random_dif  0.2358 0.1574 0.1342       1     -0.0727 0.5443
#> 10 Item10 random_dif -0.1788 0.1647 0.2777       1     -0.5016 0.1440
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1837 0.1754 0.2950       1     -0.1601 0.5276
#> 2   Item2 random_dif  0.1169 0.1626 0.4724       1     -0.2019 0.4357
#> 3   Item3 random_dif -0.0282 0.1700 0.8684       1     -0.3614 0.3050
#> 4   Item4 random_dif -0.0970 0.1690 0.5659       1     -0.4283 0.2342
#> 5   Item5 random_dif -0.1683 0.1634 0.3029       1     -0.4885 0.1519
#> 6   Item6 random_dif  0.0573 0.1679 0.7327       1     -0.2717 0.3863
#> 7   Item7 random_dif -0.1102 0.1747 0.5281       1     -0.4526 0.2322
#> 8   Item8 random_dif  0.2211 0.1610 0.1696       1     -0.0944 0.5367
#> 9   Item9 random_dif  0.0335 0.1678 0.8419       1     -0.2954 0.3623
#> 10 Item10 random_dif -0.1920 0.1607 0.2319       1     -0.5069 0.1228
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0563 0.1688 0.7387       1     -0.3870 0.2744
#> 2   Item2 random_dif  0.2030 0.1585 0.2003       1     -0.1076 0.5136
#> 3   Item3 random_dif -0.0681 0.1696 0.6882       1     -0.4005 0.2644
#> 4   Item4 random_dif  0.1517 0.1741 0.3835       1     -0.1895 0.4929
#> 5   Item5 random_dif  0.1550 0.1559 0.3201       1     -0.1506 0.4606
#> 6   Item6 random_dif -0.0373 0.1661 0.8222       1     -0.3629 0.2883
#> 7   Item7 random_dif  0.0109 0.1702 0.9488       1     -0.3227 0.3446
#> 8   Item8 random_dif -0.1180 0.1750 0.5002       1     -0.4611 0.2251
#> 9   Item9 random_dif -0.1881 0.1611 0.2428       1     -0.5038 0.1275
#> 10 Item10 random_dif -0.0689 0.1633 0.6731       1     -0.3891 0.2512
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1003 0.1821 0.5817  1.0000     -0.2566 0.4572
#> 2   Item2 random_dif -0.1297 0.1555 0.4045  1.0000     -0.4345 0.1752
#> 3   Item3 random_dif  0.0076 0.1662 0.9637  1.0000     -0.3183 0.3334
#> 4   Item4 random_dif -0.2271 0.1617 0.1601  1.0000     -0.5440 0.0897
#> 5   Item5 random_dif -0.1241 0.1629 0.4463  1.0000     -0.4433 0.1952
#> 6   Item6 random_dif  0.4207 0.1439 0.0035  0.0345   *  0.1387 0.7027
#> 7   Item7 random_dif -0.0820 0.1664 0.6220  1.0000     -0.4081 0.2440
#> 8   Item8 random_dif  0.0833 0.1677 0.6193  1.0000     -0.2454 0.4121
#> 9   Item9 random_dif -0.0337 0.1684 0.8412  1.0000     -0.3638 0.2963
#> 10 Item10 random_dif  0.0373 0.1615 0.8174  1.0000     -0.2793 0.3539
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0873 0.1669 0.6010  1.0000     -0.4144 0.2398
#> 2   Item2 random_dif  0.1111 0.1642 0.4986  1.0000     -0.2107 0.4329
#> 3   Item3 random_dif  0.0027 0.1730 0.9876  1.0000     -0.3363 0.3417
#> 4   Item4 random_dif  0.0462 0.1709 0.7869  1.0000     -0.2887 0.3811
#> 5   Item5 random_dif -0.0083 0.1616 0.9591  1.0000     -0.3250 0.3084
#> 6   Item6 random_dif  0.1960 0.1698 0.2484  1.0000     -0.1368 0.5288
#> 7   Item7 random_dif -0.1083 0.1687 0.5210  1.0000     -0.4390 0.2224
#> 8   Item8 random_dif -0.1027 0.1747 0.5567  1.0000     -0.4452 0.2398
#> 9   Item9 random_dif -0.2959 0.1537 0.0542  0.5423     -0.5973 0.0054
#> 10 Item10 random_dif  0.2268 0.1558 0.1456  1.0000     -0.0787 0.5322
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0042 0.1700 0.9802  1.0000     -0.3290 0.3374
#> 2   Item2 random_dif -0.0999 0.1565 0.5233  1.0000     -0.4066 0.2068
#> 3   Item3 random_dif -0.0587 0.1697 0.7297  1.0000     -0.3913 0.2740
#> 4   Item4 random_dif -0.2759 0.1616 0.0878  0.8784     -0.5927 0.0409
#> 5   Item5 random_dif  0.2996 0.1594 0.0602  0.6018     -0.0128 0.6120
#> 6   Item6 random_dif  0.0129 0.1723 0.9401  1.0000     -0.3247 0.3506
#> 7   Item7 random_dif  0.1631 0.1646 0.3216  1.0000     -0.1594 0.4856
#> 8   Item8 random_dif -0.0458 0.1788 0.7978  1.0000     -0.3963 0.3047
#> 9   Item9 random_dif -0.2058 0.1578 0.1923  1.0000     -0.5152 0.1036
#> 10 Item10 random_dif  0.1768 0.1587 0.2655  1.0000     -0.1343 0.4879
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0158 0.1975 0.9361  1.0000     -0.4029  0.3713
#> 2   Item2 random_dif  0.1297 0.1711 0.4486  1.0000     -0.2058  0.4651
#> 3   Item3 random_dif  0.1039 0.1690 0.5386  1.0000     -0.2273  0.4352
#> 4   Item4 random_dif  0.1576 0.1643 0.3373  1.0000     -0.1644  0.4796
#> 5   Item5 random_dif -0.0511 0.1715 0.7659  1.0000     -0.3872  0.2851
#> 6   Item6 random_dif -0.2094 0.1731 0.2262  1.0000     -0.5486  0.1298
#> 7   Item7 random_dif  0.0460 0.1769 0.7947  1.0000     -0.3007  0.3928
#> 8   Item8 random_dif -0.0493 0.1670 0.7676  1.0000     -0.3766  0.2779
#> 9   Item9 random_dif  0.2267 0.1649 0.1693  1.0000     -0.0965  0.5499
#> 10 Item10 random_dif -0.3282 0.1519 0.0307  0.3073     -0.6259 -0.0305
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1830 0.1676 0.2748  1.0000     -0.1454  0.5114
#> 2   Item2 random_dif  0.0741 0.1651 0.6534  1.0000     -0.2494  0.3977
#> 3   Item3 random_dif -0.2093 0.1723 0.2245  1.0000     -0.5471  0.1284
#> 4   Item4 random_dif -0.0861 0.1646 0.6008  1.0000     -0.4086  0.2364
#> 5   Item5 random_dif  0.0853 0.1705 0.6170  1.0000     -0.2489  0.4194
#> 6   Item6 random_dif  0.1101 0.1746 0.5283  1.0000     -0.2322  0.4524
#> 7   Item7 random_dif  0.0495 0.1650 0.7643  1.0000     -0.2739  0.3729
#> 8   Item8 random_dif  0.1596 0.1616 0.3233  1.0000     -0.1572  0.4764
#> 9   Item9 random_dif -0.3358 0.1467 0.0221  0.2206     -0.6234 -0.0483
#> 10 Item10 random_dif -0.0113 0.1652 0.9456  1.0000     -0.3351  0.3125
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1747 0.1668 0.2948  1.0000     -0.5015  0.1521
#> 2   Item2 random_dif -0.0364 0.1662 0.8267  1.0000     -0.3622  0.2894
#> 3   Item3 random_dif -0.1520 0.1618 0.3477  1.0000     -0.4692  0.1652
#> 4   Item4 random_dif -0.2232 0.1600 0.1630  1.0000     -0.5367  0.0904
#> 5   Item5 random_dif  0.0474 0.1679 0.7778  1.0000     -0.2817  0.3764
#> 6   Item6 random_dif -0.3409 0.1537 0.0265  0.2655     -0.6421 -0.0397
#> 7   Item7 random_dif  0.1298 0.1759 0.4606  1.0000     -0.2150  0.4746
#> 8   Item8 random_dif  0.4524 0.1330 0.0007  0.0067  **  0.1918  0.7131
#> 9   Item9 random_dif -0.0143 0.1655 0.9313  1.0000     -0.3387  0.3102
#> 10 Item10 random_dif  0.2766 0.1498 0.0647  0.6472     -0.0169  0.5702
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1872 0.1630 0.2508       1     -0.1323 0.5066
#> 2   Item2 random_dif -0.0223 0.1654 0.8929       1     -0.3465 0.3019
#> 3   Item3 random_dif  0.1866 0.1591 0.2410       1     -0.1253 0.4984
#> 4   Item4 random_dif  0.0054 0.1671 0.9744       1     -0.3222 0.3330
#> 5   Item5 random_dif  0.2290 0.1575 0.1458       1     -0.0796 0.5376
#> 6   Item6 random_dif -0.2633 0.1701 0.1215       1     -0.5966 0.0700
#> 7   Item7 random_dif -0.0118 0.1660 0.9432       1     -0.3371 0.3134
#> 8   Item8 random_dif -0.1497 0.1611 0.3526       1     -0.4654 0.1659
#> 9   Item9 random_dif  0.0222 0.1738 0.8986       1     -0.3185 0.3628
#> 10 Item10 random_dif -0.2284 0.1686 0.1755       1     -0.5588 0.1020
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1406 0.1745 0.4201  1.0000     -0.2013  0.4826
#> 2   Item2 random_dif -0.2224 0.1611 0.1674  1.0000     -0.5380  0.0933
#> 3   Item3 random_dif  0.1169 0.1755 0.5053  1.0000     -0.2270  0.4609
#> 4   Item4 random_dif  0.3758 0.1586 0.0178  0.1781      0.0650  0.6867
#> 5   Item5 random_dif  0.0857 0.1663 0.6062  1.0000     -0.2401  0.4116
#> 6   Item6 random_dif -0.3645 0.1551 0.0188  0.1877     -0.6684 -0.0605
#> 7   Item7 random_dif  0.0634 0.1670 0.7044  1.0000     -0.2640  0.3907
#> 8   Item8 random_dif  0.0142 0.1767 0.9357  1.0000     -0.3320  0.3605
#> 9   Item9 random_dif -0.1324 0.1730 0.4439  1.0000     -0.4714  0.2066
#> 10 Item10 random_dif -0.0444 0.1617 0.7834  1.0000     -0.3613  0.2725
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.3175 0.1648 0.0540  0.5401     -0.6405 0.0055
#> 2   Item2 random_dif  0.1058 0.1659 0.5238  1.0000     -0.2194 0.4310
#> 3   Item3 random_dif  0.1313 0.1742 0.4510  1.0000     -0.2101 0.4727
#> 4   Item4 random_dif  0.0123 0.1743 0.9435  1.0000     -0.3293 0.3540
#> 5   Item5 random_dif  0.0600 0.1660 0.7178  1.0000     -0.2653 0.3853
#> 6   Item6 random_dif -0.0889 0.1680 0.5965  1.0000     -0.4181 0.2403
#> 7   Item7 random_dif -0.2020 0.1752 0.2488  1.0000     -0.5454 0.1413
#> 8   Item8 random_dif  0.2271 0.1567 0.1472  1.0000     -0.0800 0.5342
#> 9   Item9 random_dif  0.0691 0.1614 0.6685  1.0000     -0.2472 0.3855
#> 10 Item10 random_dif -0.0556 0.1648 0.7356  1.0000     -0.3785 0.2673
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0310 0.1772 0.8612       1     -0.3163 0.3782
#> 2   Item2 random_dif  0.0688 0.1734 0.6916       1     -0.2711 0.4087
#> 3   Item3 random_dif -0.1371 0.1697 0.4192       1     -0.4696 0.1955
#> 4   Item4 random_dif -0.0559 0.1807 0.7571       1     -0.4100 0.2982
#> 5   Item5 random_dif  0.1146 0.1634 0.4830       1     -0.2056 0.4348
#> 6   Item6 random_dif -0.1524 0.1662 0.3593       1     -0.4782 0.1734
#> 7   Item7 random_dif -0.0514 0.1704 0.7630       1     -0.3853 0.2826
#> 8   Item8 random_dif  0.0604 0.1831 0.7414       1     -0.2984 0.4193
#> 9   Item9 random_dif  0.1272 0.1600 0.4268       1     -0.1865 0.4408
#> 10 Item10 random_dif -0.0158 0.1796 0.9298       1     -0.3679 0.3362
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.3812 0.1545 0.0136  0.1364      0.0783 0.6840
#> 2   Item2 random_dif -0.2151 0.1558 0.1673  1.0000     -0.5205 0.0902
#> 3   Item3 random_dif -0.0043 0.1784 0.9808  1.0000     -0.3540 0.3454
#> 4   Item4 random_dif  0.0746 0.1661 0.6532  1.0000     -0.2509 0.4001
#> 5   Item5 random_dif -0.1864 0.1747 0.2860  1.0000     -0.5288 0.1560
#> 6   Item6 random_dif -0.0172 0.1714 0.9200  1.0000     -0.3532 0.3187
#> 7   Item7 random_dif  0.2127 0.1647 0.1966  1.0000     -0.1101 0.5354
#> 8   Item8 random_dif -0.0190 0.1708 0.9113  1.0000     -0.3537 0.3157
#> 9   Item9 random_dif  0.0697 0.1670 0.6764  1.0000     -0.2577 0.3971
#> 10 Item10 random_dif -0.2641 0.1569 0.0925  0.9248     -0.5717 0.0436
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0598 0.2065 0.7720  1.0000     -0.4646 0.3449
#> 2   Item2 random_dif -0.1278 0.1660 0.4413  1.0000     -0.4532 0.1975
#> 3   Item3 random_dif -0.2299 0.1668 0.1681  1.0000     -0.5568 0.0970
#> 4   Item4 random_dif -0.0505 0.1708 0.7674  1.0000     -0.3852 0.2842
#> 5   Item5 random_dif  0.1327 0.1656 0.4230  1.0000     -0.1918 0.4572
#> 6   Item6 random_dif  0.0016 0.1791 0.9929  1.0000     -0.3495 0.3527
#> 7   Item7 random_dif  0.0456 0.1724 0.7914  1.0000     -0.2923 0.3835
#> 8   Item8 random_dif -0.2096 0.1694 0.2160  1.0000     -0.5416 0.1224
#> 9   Item9 random_dif  0.1657 0.1700 0.3299  1.0000     -0.1676 0.4989
#> 10 Item10 random_dif  0.3195 0.1595 0.0451  0.4513      0.0069 0.6321
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0630 0.1761 0.7207  1.0000     -0.2822  0.4081
#> 2   Item2 random_dif -0.0814 0.1671 0.6259  1.0000     -0.4089  0.2460
#> 3   Item3 random_dif  0.1022 0.1674 0.5418  1.0000     -0.2260  0.4303
#> 4   Item4 random_dif  0.1612 0.1676 0.3361  1.0000     -0.1673  0.4897
#> 5   Item5 random_dif  0.1063 0.1677 0.5261  1.0000     -0.2224  0.4350
#> 6   Item6 random_dif  0.0547 0.1671 0.7433  1.0000     -0.2728  0.3823
#> 7   Item7 random_dif -0.3638 0.1664 0.0288  0.2881     -0.6899 -0.0376
#> 8   Item8 random_dif -0.0634 0.1718 0.7122  1.0000     -0.4000  0.2733
#> 9   Item9 random_dif -0.1481 0.1714 0.3877  1.0000     -0.4841  0.1879
#> 10 Item10 random_dif  0.1205 0.1604 0.4526  1.0000     -0.1939  0.4350
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0518 0.1724 0.7637       1     -0.3897 0.2860
#> 2   Item2 random_dif -0.0785 0.1616 0.6272       1     -0.3951 0.2382
#> 3   Item3 random_dif  0.2409 0.1671 0.1496       1     -0.0867 0.5685
#> 4   Item4 random_dif  0.2319 0.1660 0.1624       1     -0.0935 0.5573
#> 5   Item5 random_dif -0.1285 0.1712 0.4529       1     -0.4641 0.2071
#> 6   Item6 random_dif -0.0725 0.1726 0.6745       1     -0.4108 0.2658
#> 7   Item7 random_dif -0.0611 0.1725 0.7230       1     -0.3991 0.2769
#> 8   Item8 random_dif -0.0966 0.1757 0.5825       1     -0.4410 0.2478
#> 9   Item9 random_dif  0.1680 0.1736 0.3332       1     -0.1722 0.5082
#> 10 Item10 random_dif -0.1197 0.1646 0.4670       1     -0.4423 0.2029
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0098 0.1672 0.9533  1.0000     -0.3180 0.3376
#> 2   Item2 random_dif -0.0384 0.1664 0.8177  1.0000     -0.3645 0.2878
#> 3   Item3 random_dif  0.1455 0.1618 0.3685  1.0000     -0.1716 0.4626
#> 4   Item4 random_dif  0.1429 0.1784 0.4232  1.0000     -0.2068 0.4925
#> 5   Item5 random_dif  0.0783 0.1725 0.6499  1.0000     -0.2598 0.4163
#> 6   Item6 random_dif -0.2899 0.1614 0.0725  0.7251     -0.6064 0.0265
#> 7   Item7 random_dif  0.0287 0.1696 0.8655  1.0000     -0.3037 0.3612
#> 8   Item8 random_dif -0.0014 0.1636 0.9932  1.0000     -0.3220 0.3192
#> 9   Item9 random_dif  0.0298 0.1659 0.8574  1.0000     -0.2953 0.3549
#> 10 Item10 random_dif -0.1111 0.1635 0.4967  1.0000     -0.4315 0.2093
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2727 0.1748 0.1187       1     -0.6153 0.0698
#> 2   Item2 random_dif  0.0557 0.1806 0.7577       1     -0.2983 0.4098
#> 3   Item3 random_dif -0.1408 0.1701 0.4078       1     -0.4742 0.1926
#> 4   Item4 random_dif -0.2149 0.1654 0.1939       1     -0.5391 0.1093
#> 5   Item5 random_dif  0.0293 0.1755 0.8675       1     -0.3146 0.3732
#> 6   Item6 random_dif  0.1649 0.1747 0.3451       1     -0.1775 0.5073
#> 7   Item7 random_dif  0.1416 0.1724 0.4113       1     -0.1962 0.4794
#> 8   Item8 random_dif  0.0529 0.1707 0.7566       1     -0.2816 0.3874
#> 9   Item9 random_dif -0.1117 0.1643 0.4967       1     -0.4338 0.2104
#> 10 Item10 random_dif  0.2585 0.1610 0.1084       1     -0.0571 0.5742
#>      Item        Var   gamma     se pvalue padj.BH  sig   lower   upper
#> 1   Item1 random_dif  0.1208 0.1855 0.5149  1.0000      -0.2428  0.4845
#> 2   Item2 random_dif  0.1317 0.1703 0.4395  1.0000      -0.2021  0.4655
#> 3   Item3 random_dif -0.0759 0.1723 0.6597  1.0000      -0.4136  0.2619
#> 4   Item4 random_dif  0.3556 0.1614 0.0276  0.2762       0.0392  0.6719
#> 5   Item5 random_dif  0.0917 0.1675 0.5844  1.0000      -0.2367  0.4200
#> 6   Item6 random_dif -0.0478 0.1784 0.7889  1.0000      -0.3975  0.3020
#> 7   Item7 random_dif -0.0236 0.1754 0.8930  1.0000      -0.3675  0.3203
#> 8   Item8 random_dif -0.0260 0.1770 0.8831  1.0000      -0.3730  0.3209
#> 9   Item9 random_dif -0.5414 0.1223 0.0000  0.0001  *** -0.7811 -0.3017
#> 10 Item10 random_dif  0.1219 0.1732 0.4816  1.0000      -0.2176  0.4613
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1383 0.1737 0.4260  1.0000     -0.4787  0.2022
#> 2   Item2 random_dif  0.0652 0.1687 0.6990  1.0000     -0.2654  0.3958
#> 3   Item3 random_dif -0.1589 0.1803 0.3782  1.0000     -0.5122  0.1945
#> 4   Item4 random_dif  0.2318 0.1672 0.1656  1.0000     -0.0959  0.5594
#> 5   Item5 random_dif -0.1137 0.1637 0.4874  1.0000     -0.4345  0.2071
#> 6   Item6 random_dif  0.3438 0.1614 0.0332  0.3316      0.0275  0.6600
#> 7   Item7 random_dif  0.0317 0.1694 0.8514  1.0000     -0.3003  0.3638
#> 8   Item8 random_dif -0.3184 0.1615 0.0486  0.4860     -0.6349 -0.0020
#> 9   Item9 random_dif -0.0086 0.1736 0.9606  1.0000     -0.3488  0.3316
#> 10 Item10 random_dif  0.0429 0.1665 0.7966  1.0000     -0.2835  0.3693
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.5133 0.1442 0.0004  0.0037  ** -0.7959 -0.2307
#> 2   Item2 random_dif -0.0506 0.1638 0.7572  1.0000     -0.3716  0.2703
#> 3   Item3 random_dif -0.0486 0.1770 0.7837  1.0000     -0.3956  0.2984
#> 4   Item4 random_dif -0.1995 0.1554 0.1991  1.0000     -0.5041  0.1050
#> 5   Item5 random_dif  0.1240 0.1709 0.4680  1.0000     -0.2109  0.4589
#> 6   Item6 random_dif  0.2722 0.1824 0.1356  1.0000     -0.0853  0.6296
#> 7   Item7 random_dif  0.0740 0.1669 0.6575  1.0000     -0.2531  0.4010
#> 8   Item8 random_dif  0.0804 0.1608 0.6172  1.0000     -0.2348  0.3956
#> 9   Item9 random_dif  0.1079 0.1666 0.5173  1.0000     -0.2187  0.4345
#> 10 Item10 random_dif  0.1377 0.1570 0.3803  1.0000     -0.1700  0.4454
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2366 0.1672 0.1572  1.0000     -0.5643 0.0912
#> 2   Item2 random_dif  0.3201 0.1569 0.0414  0.4137      0.0125 0.6277
#> 3   Item3 random_dif  0.0546 0.1673 0.7440  1.0000     -0.2733 0.3826
#> 4   Item4 random_dif  0.0679 0.1704 0.6901  1.0000     -0.2660 0.4019
#> 5   Item5 random_dif -0.2831 0.1715 0.0988  0.9879     -0.6192 0.0530
#> 6   Item6 random_dif  0.0875 0.1657 0.5977  1.0000     -0.2374 0.4123
#> 7   Item7 random_dif -0.0165 0.1684 0.9219  1.0000     -0.3467 0.3136
#> 8   Item8 random_dif -0.0759 0.1698 0.6550  1.0000     -0.4086 0.2569
#> 9   Item9 random_dif  0.0264 0.1710 0.8774  1.0000     -0.3088 0.3616
#> 10 Item10 random_dif  0.0050 0.1690 0.9763  1.0000     -0.3263 0.3363
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1230 0.1727 0.4763  1.0000     -0.2155 0.4614
#> 2   Item2 random_dif -0.0475 0.1697 0.7796  1.0000     -0.3802 0.2852
#> 3   Item3 random_dif -0.2373 0.1702 0.1631  1.0000     -0.5708 0.0962
#> 4   Item4 random_dif  0.1526 0.1721 0.3752  1.0000     -0.1847 0.4899
#> 5   Item5 random_dif -0.2341 0.1695 0.1673  1.0000     -0.5663 0.0981
#> 6   Item6 random_dif -0.2978 0.1680 0.0763  0.7626     -0.6271 0.0314
#> 7   Item7 random_dif  0.4264 0.1434 0.0029  0.0293   *  0.1454 0.7074
#> 8   Item8 random_dif  0.1963 0.1684 0.2439  1.0000     -0.1338 0.5264
#> 9   Item9 random_dif -0.0454 0.1671 0.7859  1.0000     -0.3729 0.2821
#> 10 Item10 random_dif -0.0695 0.1650 0.6736  1.0000     -0.3929 0.2539
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.2118 0.1667 0.2040   1.000     -0.1150  0.5385
#> 2   Item2 random_dif  0.2490 0.1542 0.1065   1.000     -0.0533  0.5513
#> 3   Item3 random_dif  0.0712 0.1711 0.6772   1.000     -0.2641  0.4065
#> 4   Item4 random_dif -0.1108 0.1664 0.5055   1.000     -0.4369  0.2153
#> 5   Item5 random_dif  0.0698 0.1700 0.6814   1.000     -0.2633  0.4029
#> 6   Item6 random_dif  0.0715 0.1682 0.6706   1.000     -0.2581  0.4011
#> 7   Item7 random_dif  0.0925 0.1696 0.5853   1.000     -0.2399  0.4249
#> 8   Item8 random_dif -0.1335 0.1723 0.4383   1.000     -0.4712  0.2041
#> 9   Item9 random_dif -0.1392 0.1683 0.4081   1.000     -0.4691  0.1906
#> 10 Item10 random_dif -0.3596 0.1490 0.0158   0.158     -0.6516 -0.0676
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1570 0.1770 0.3753  1.0000     -0.5039 0.1900
#> 2   Item2 random_dif  0.3760 0.1503 0.0124  0.1235      0.0814 0.6706
#> 3   Item3 random_dif -0.0854 0.1700 0.6153  1.0000     -0.4185 0.2477
#> 4   Item4 random_dif  0.1120 0.1764 0.5253  1.0000     -0.2336 0.4577
#> 5   Item5 random_dif  0.3029 0.1608 0.0596  0.5964     -0.0123 0.6181
#> 6   Item6 random_dif -0.0668 0.1683 0.6913  1.0000     -0.3967 0.2631
#> 7   Item7 random_dif -0.2282 0.1724 0.1856  1.0000     -0.5662 0.1097
#> 8   Item8 random_dif -0.0533 0.1831 0.7711  1.0000     -0.4121 0.3056
#> 9   Item9 random_dif -0.1501 0.1692 0.3750  1.0000     -0.4818 0.1816
#> 10 Item10 random_dif -0.0872 0.1680 0.6039  1.0000     -0.4165 0.2421
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1373 0.1744 0.4311  1.0000     -0.2045 0.4792
#> 2   Item2 random_dif  0.0896 0.1734 0.6054  1.0000     -0.2503 0.4295
#> 3   Item3 random_dif -0.1243 0.1683 0.4601  1.0000     -0.4543 0.2056
#> 4   Item4 random_dif  0.1546 0.1601 0.3343  1.0000     -0.1592 0.4685
#> 5   Item5 random_dif -0.0230 0.1611 0.8863  1.0000     -0.3387 0.2927
#> 6   Item6 random_dif  0.0120 0.1639 0.9416  1.0000     -0.3093 0.3333
#> 7   Item7 random_dif -0.1075 0.1709 0.5295  1.0000     -0.4425 0.2275
#> 8   Item8 random_dif  0.0015 0.1763 0.9933  1.0000     -0.3440 0.3470
#> 9   Item9 random_dif -0.2634 0.1556 0.0905  0.9048     -0.5684 0.0416
#> 10 Item10 random_dif  0.1363 0.1588 0.3906  1.0000     -0.1749 0.4475
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.5057 0.1341 0.0002  0.0016  ** -0.7686 -0.2429
#> 2   Item2 random_dif  0.0267 0.1634 0.8701  1.0000     -0.2935  0.3469
#> 3   Item3 random_dif -0.0149 0.1706 0.9305  1.0000     -0.3493  0.3195
#> 4   Item4 random_dif -0.1337 0.1681 0.4262  1.0000     -0.4632  0.1957
#> 5   Item5 random_dif  0.2284 0.1567 0.1451  1.0000     -0.0788  0.5355
#> 6   Item6 random_dif  0.0299 0.1772 0.8659  1.0000     -0.3175  0.3773
#> 7   Item7 random_dif  0.1833 0.1634 0.2621  1.0000     -0.1370  0.5036
#> 8   Item8 random_dif  0.1718 0.1692 0.3100  1.0000     -0.1599  0.5034
#> 9   Item9 random_dif  0.0247 0.1645 0.8808  1.0000     -0.2978  0.3472
#> 10 Item10 random_dif -0.0066 0.1590 0.9668  1.0000     -0.3183  0.3050
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1880 0.1739 0.2796       1     -0.5288 0.1528
#> 2   Item2 random_dif -0.0433 0.1709 0.7998       1     -0.3783 0.2917
#> 3   Item3 random_dif  0.2423 0.1617 0.1339       1     -0.0745 0.5592
#> 4   Item4 random_dif -0.1652 0.1846 0.3709       1     -0.5269 0.1966
#> 5   Item5 random_dif  0.0213 0.1647 0.8969       1     -0.3014 0.3441
#> 6   Item6 random_dif -0.0969 0.1636 0.5537       1     -0.4175 0.2237
#> 7   Item7 random_dif -0.0491 0.1713 0.7745       1     -0.3848 0.2867
#> 8   Item8 random_dif  0.0292 0.1697 0.8635       1     -0.3035 0.3619
#> 9   Item9 random_dif  0.1587 0.1699 0.3503       1     -0.1743 0.4917
#> 10 Item10 random_dif  0.0392 0.1673 0.8148       1     -0.2886 0.3670
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3509 0.1741 0.0439  0.4387     -0.6921 -0.0096
#> 2   Item2 random_dif  0.1145 0.1656 0.4892  1.0000     -0.2100  0.4391
#> 3   Item3 random_dif -0.1908 0.1729 0.2698  1.0000     -0.5296  0.1481
#> 4   Item4 random_dif  0.0867 0.1631 0.5950  1.0000     -0.2329  0.4062
#> 5   Item5 random_dif  0.0680 0.1626 0.6756  1.0000     -0.2506  0.3866
#> 6   Item6 random_dif -0.1679 0.1732 0.3324  1.0000     -0.5074  0.1716
#> 7   Item7 random_dif  0.1203 0.1670 0.4711  1.0000     -0.2069  0.4476
#> 8   Item8 random_dif  0.1085 0.1692 0.5212  1.0000     -0.2231  0.4402
#> 9   Item9 random_dif -0.1149 0.1655 0.4873  1.0000     -0.4392  0.2094
#> 10 Item10 random_dif  0.2000 0.1572 0.2034  1.0000     -0.1082  0.5082
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0269 0.1697 0.8739  1.0000     -0.3596  0.3057
#> 2   Item2 random_dif  0.1658 0.1595 0.2987  1.0000     -0.1469  0.4784
#> 3   Item3 random_dif -0.0859 0.1715 0.6165  1.0000     -0.4220  0.2502
#> 4   Item4 random_dif -0.3274 0.1562 0.0360  0.3602     -0.6335 -0.0214
#> 5   Item5 random_dif -0.1801 0.1565 0.2498  1.0000     -0.4868  0.1266
#> 6   Item6 random_dif -0.0730 0.1781 0.6819  1.0000     -0.4221  0.2760
#> 7   Item7 random_dif -0.1234 0.1635 0.4506  1.0000     -0.4439  0.1971
#> 8   Item8 random_dif -0.0041 0.1662 0.9801  1.0000     -0.3299  0.3216
#> 9   Item9 random_dif  0.1333 0.1665 0.4232  1.0000     -0.1930  0.4596
#> 10 Item10 random_dif  0.4553 0.1324 0.0006  0.0059  **  0.1958  0.7148
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1638 0.1896 0.3878  1.0000     -0.2079 0.5354
#> 2   Item2 random_dif -0.0496 0.1675 0.7672  1.0000     -0.3778 0.2787
#> 3   Item3 random_dif  0.0452 0.1667 0.7862  1.0000     -0.2814 0.3719
#> 4   Item4 random_dif -0.2463 0.1529 0.1071  1.0000     -0.5459 0.0533
#> 5   Item5 random_dif -0.3006 0.1557 0.0535  0.5351     -0.6058 0.0045
#> 6   Item6 random_dif  0.0792 0.1756 0.6518  1.0000     -0.2649 0.4233
#> 7   Item7 random_dif  0.0194 0.1781 0.9134  1.0000     -0.3296 0.3684
#> 8   Item8 random_dif  0.2104 0.1575 0.1817  1.0000     -0.0984 0.5192
#> 9   Item9 random_dif  0.0540 0.1734 0.7556  1.0000     -0.2859 0.3939
#> 10 Item10 random_dif  0.1042 0.1636 0.5244  1.0000     -0.2165 0.4249
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1052 0.1716 0.5399       1     -0.2312 0.4415
#> 2   Item2 random_dif -0.0419 0.1612 0.7949       1     -0.3578 0.2740
#> 3   Item3 random_dif -0.0718 0.1739 0.6796       1     -0.4127 0.2691
#> 4   Item4 random_dif  0.1572 0.1769 0.3742       1     -0.1895 0.5039
#> 5   Item5 random_dif  0.0463 0.1592 0.7713       1     -0.2657 0.3583
#> 6   Item6 random_dif -0.1323 0.1812 0.4652       1     -0.4874 0.2228
#> 7   Item7 random_dif -0.1653 0.1636 0.3123       1     -0.4859 0.1553
#> 8   Item8 random_dif -0.0652 0.1773 0.7131       1     -0.4128 0.2824
#> 9   Item9 random_dif  0.2060 0.1639 0.2087       1     -0.1152 0.5273
#> 10 Item10 random_dif -0.0326 0.1643 0.8427       1     -0.3547 0.2895
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0621 0.1749 0.7226       1     -0.2807 0.4048
#> 2   Item2 random_dif  0.2079 0.1631 0.2026       1     -0.1119 0.5276
#> 3   Item3 random_dif -0.1218 0.1656 0.4623       1     -0.4464 0.2029
#> 4   Item4 random_dif  0.2358 0.1552 0.1287       1     -0.0684 0.5401
#> 5   Item5 random_dif -0.0087 0.1661 0.9581       1     -0.3342 0.3168
#> 6   Item6 random_dif  0.0199 0.1737 0.9089       1     -0.3205 0.3602
#> 7   Item7 random_dif -0.1375 0.1718 0.4234       1     -0.4743 0.1992
#> 8   Item8 random_dif -0.1613 0.1670 0.3341       1     -0.4886 0.1660
#> 9   Item9 random_dif -0.1838 0.1588 0.2472       1     -0.4950 0.1275
#> 10 Item10 random_dif  0.0810 0.1641 0.6215       1     -0.2406 0.4026
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1545 0.1678 0.3574  1.0000     -0.4835 0.1745
#> 2   Item2 random_dif -0.0179 0.1833 0.9222  1.0000     -0.3772 0.3414
#> 3   Item3 random_dif  0.1162 0.1685 0.4906  1.0000     -0.2141 0.4464
#> 4   Item4 random_dif  0.3289 0.1569 0.0361  0.3608      0.0213 0.6364
#> 5   Item5 random_dif  0.0986 0.1679 0.5571  1.0000     -0.2305 0.4277
#> 6   Item6 random_dif  0.0441 0.1738 0.7997  1.0000     -0.2965 0.3847
#> 7   Item7 random_dif -0.0316 0.1665 0.8494  1.0000     -0.3579 0.2946
#> 8   Item8 random_dif -0.2030 0.1591 0.2020  1.0000     -0.5149 0.1089
#> 9   Item9 random_dif -0.0821 0.1693 0.6276  1.0000     -0.4140 0.2498
#> 10 Item10 random_dif -0.0907 0.1707 0.5953  1.0000     -0.4252 0.2438
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0523 0.1794 0.7706  1.0000     -0.2993  0.4040
#> 2   Item2 random_dif -0.2301 0.1720 0.1810  1.0000     -0.5673  0.1070
#> 3   Item3 random_dif  0.0977 0.1606 0.5427  1.0000     -0.2170  0.4124
#> 4   Item4 random_dif  0.2891 0.1587 0.0685  0.6847     -0.0219  0.6000
#> 5   Item5 random_dif -0.0180 0.1700 0.9158  1.0000     -0.3512  0.3152
#> 6   Item6 random_dif -0.0654 0.1692 0.6992  1.0000     -0.3970  0.2663
#> 7   Item7 random_dif -0.4129 0.1402 0.0032  0.0322   * -0.6876 -0.1382
#> 8   Item8 random_dif  0.0623 0.1637 0.7037  1.0000     -0.2586  0.3831
#> 9   Item9 random_dif  0.0494 0.1616 0.7599  1.0000     -0.2673  0.3661
#> 10 Item10 random_dif  0.1598 0.1557 0.3048  1.0000     -0.1454  0.4650
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0538 0.1760 0.7597       1     -0.2912 0.3988
#> 2   Item2 random_dif -0.1193 0.1695 0.4817       1     -0.4515 0.2130
#> 3   Item3 random_dif -0.0678 0.1629 0.6773       1     -0.3871 0.2515
#> 4   Item4 random_dif  0.0396 0.1760 0.8220       1     -0.3054 0.3846
#> 5   Item5 random_dif  0.2373 0.1688 0.1597       1     -0.0935 0.5680
#> 6   Item6 random_dif -0.0513 0.1730 0.7669       1     -0.3904 0.2878
#> 7   Item7 random_dif -0.0900 0.1694 0.5951       1     -0.4222 0.2421
#> 8   Item8 random_dif  0.0316 0.1765 0.8578       1     -0.3142 0.3775
#> 9   Item9 random_dif  0.0687 0.1627 0.6729       1     -0.2502 0.3876
#> 10 Item10 random_dif -0.0758 0.1640 0.6441       1     -0.3972 0.2457
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1032 0.1794 0.5651       1     -0.4549 0.2485
#> 2   Item2 random_dif -0.0416 0.1732 0.8101       1     -0.3810 0.2978
#> 3   Item3 random_dif -0.0776 0.1725 0.6529       1     -0.4157 0.2606
#> 4   Item4 random_dif  0.0347 0.1749 0.8428       1     -0.3082 0.3776
#> 5   Item5 random_dif  0.2169 0.1559 0.1641       1     -0.0886 0.5225
#> 6   Item6 random_dif  0.0890 0.1787 0.6187       1     -0.2614 0.4393
#> 7   Item7 random_dif -0.2069 0.1593 0.1938       1     -0.5191 0.1052
#> 8   Item8 random_dif -0.0525 0.1716 0.7598       1     -0.3889 0.2839
#> 9   Item9 random_dif -0.0025 0.1611 0.9875       1     -0.3183 0.3133
#> 10 Item10 random_dif  0.1106 0.1571 0.4815       1     -0.1973 0.4185
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0710 0.1694 0.6750  1.0000     -0.2610  0.4031
#> 2   Item2 random_dif -0.0431 0.1625 0.7911  1.0000     -0.3616  0.2755
#> 3   Item3 random_dif -0.1008 0.1719 0.5575  1.0000     -0.4377  0.2361
#> 4   Item4 random_dif -0.3797 0.1525 0.0128  0.1275     -0.6785 -0.0809
#> 5   Item5 random_dif  0.0511 0.1629 0.7539  1.0000     -0.2682  0.3703
#> 6   Item6 random_dif -0.0952 0.1699 0.5755  1.0000     -0.4283  0.2379
#> 7   Item7 random_dif  0.1995 0.1614 0.2166  1.0000     -0.1169  0.5159
#> 8   Item8 random_dif -0.1215 0.1637 0.4580  1.0000     -0.4424  0.1994
#> 9   Item9 random_dif  0.0966 0.1631 0.5535  1.0000     -0.2230  0.4162
#> 10 Item10 random_dif  0.2995 0.1567 0.0559  0.5595     -0.0076  0.6066
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1545 0.1697 0.3624   1.000     -0.4871  0.1780
#> 2   Item2 random_dif  0.1377 0.1562 0.3781   1.000     -0.1684  0.4438
#> 3   Item3 random_dif  0.1500 0.1642 0.3610   1.000     -0.1718  0.4718
#> 4   Item4 random_dif -0.1254 0.1683 0.4564   1.000     -0.4553  0.2045
#> 5   Item5 random_dif -0.3492 0.1575 0.0266   0.266     -0.6579 -0.0405
#> 6   Item6 random_dif -0.0207 0.1805 0.9088   1.000     -0.3745  0.3331
#> 7   Item7 random_dif  0.2251 0.1571 0.1518   1.000     -0.0827  0.5329
#> 8   Item8 random_dif -0.0251 0.1679 0.8810   1.000     -0.3543  0.3040
#> 9   Item9 random_dif -0.0182 0.1640 0.9117   1.000     -0.3397  0.3033
#> 10 Item10 random_dif  0.1005 0.1615 0.5337   1.000     -0.2160  0.4171
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1406 0.1744 0.4200  1.0000     -0.4824 0.2012
#> 2   Item2 random_dif  0.0242 0.1662 0.8845  1.0000     -0.3016 0.3499
#> 3   Item3 random_dif  0.1927 0.1680 0.2515  1.0000     -0.1366 0.5219
#> 4   Item4 random_dif  0.1908 0.1708 0.2638  1.0000     -0.1439 0.5256
#> 5   Item5 random_dif  0.2122 0.1655 0.1998  1.0000     -0.1122 0.5367
#> 6   Item6 random_dif -0.1960 0.1622 0.2270  1.0000     -0.5140 0.1220
#> 7   Item7 random_dif  0.1734 0.1798 0.3348  1.0000     -0.1789 0.5257
#> 8   Item8 random_dif -0.0163 0.1759 0.9260  1.0000     -0.3610 0.3284
#> 9   Item9 random_dif -0.3086 0.1724 0.0735  0.7347     -0.6466 0.0293
#> 10 Item10 random_dif -0.1410 0.1671 0.3988  1.0000     -0.4684 0.1865
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2433 0.1650 0.1403  1.0000     -0.0801 0.5667
#> 2   Item2 random_dif  0.1167 0.1644 0.4777  1.0000     -0.2055 0.4390
#> 3   Item3 random_dif  0.1450 0.1688 0.3903  1.0000     -0.1858 0.4757
#> 4   Item4 random_dif  0.0812 0.1681 0.6288  1.0000     -0.2482 0.4107
#> 5   Item5 random_dif -0.2579 0.1580 0.1025  1.0000     -0.5675 0.0517
#> 6   Item6 random_dif  0.1036 0.1664 0.5338  1.0000     -0.2226 0.4297
#> 7   Item7 random_dif  0.1330 0.1849 0.4721  1.0000     -0.2295 0.4955
#> 8   Item8 random_dif -0.3098 0.1656 0.0614  0.6137     -0.6344 0.0148
#> 9   Item9 random_dif -0.1474 0.1670 0.3775  1.0000     -0.4746 0.1799
#> 10 Item10 random_dif -0.0637 0.1653 0.6999  1.0000     -0.3877 0.2603
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1758 0.1814 0.3325  1.0000     -0.1797 0.5312
#> 2   Item2 random_dif  0.0437 0.1719 0.7992  1.0000     -0.2932 0.3806
#> 3   Item3 random_dif  0.0437 0.1756 0.8033  1.0000     -0.3004 0.3878
#> 4   Item4 random_dif -0.1622 0.1747 0.3531  1.0000     -0.5045 0.1801
#> 5   Item5 random_dif  0.3353 0.1636 0.0405  0.4045      0.0146 0.6560
#> 6   Item6 random_dif  0.1098 0.1759 0.5325  1.0000     -0.2349 0.4545
#> 7   Item7 random_dif -0.1217 0.1663 0.4644  1.0000     -0.4477 0.2043
#> 8   Item8 random_dif -0.2338 0.1619 0.1487  1.0000     -0.5512 0.0835
#> 9   Item9 random_dif -0.0209 0.1617 0.8970  1.0000     -0.3379 0.2960
#> 10 Item10 random_dif -0.1318 0.1632 0.4193  1.0000     -0.4517 0.1881
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0326 0.1812 0.8574  1.0000     -0.3225 0.3876
#> 2   Item2 random_dif -0.2876 0.1558 0.0649  0.6493     -0.5929 0.0178
#> 3   Item3 random_dif  0.0757 0.1703 0.6568  1.0000     -0.2582 0.4095
#> 4   Item4 random_dif  0.1486 0.1651 0.3678  1.0000     -0.1748 0.4721
#> 5   Item5 random_dif -0.0887 0.1670 0.5953  1.0000     -0.4160 0.2386
#> 6   Item6 random_dif -0.1756 0.1693 0.2996  1.0000     -0.5074 0.1562
#> 7   Item7 random_dif  0.1159 0.1761 0.5104  1.0000     -0.2292 0.4610
#> 8   Item8 random_dif  0.2565 0.1612 0.1115  1.0000     -0.0594 0.5724
#> 9   Item9 random_dif  0.0174 0.1698 0.9186  1.0000     -0.3155 0.3502
#> 10 Item10 random_dif -0.0775 0.1642 0.6368  1.0000     -0.3994 0.2443
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1718 0.1740 0.3234       1     -0.1692 0.5128
#> 2   Item2 random_dif -0.2091 0.1576 0.1846       1     -0.5180 0.0998
#> 3   Item3 random_dif -0.0447 0.1679 0.7899       1     -0.3738 0.2844
#> 4   Item4 random_dif -0.0806 0.1689 0.6331       1     -0.4117 0.2504
#> 5   Item5 random_dif  0.0673 0.1662 0.6856       1     -0.2585 0.3931
#> 6   Item6 random_dif  0.0338 0.1709 0.8430       1     -0.3011 0.3688
#> 7   Item7 random_dif  0.1711 0.1728 0.3222       1     -0.1676 0.5098
#> 8   Item8 random_dif  0.0221 0.1741 0.8991       1     -0.3192 0.3633
#> 9   Item9 random_dif -0.1852 0.1603 0.2481       1     -0.4994 0.1291
#> 10 Item10 random_dif  0.0985 0.1620 0.5432       1     -0.2191 0.4161
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1150 0.1627 0.4794       1     -0.2038 0.4338
#> 2   Item2 random_dif  0.1303 0.1655 0.4310       1     -0.1940 0.4547
#> 3   Item3 random_dif  0.0937 0.1650 0.5703       1     -0.2298 0.4171
#> 4   Item4 random_dif  0.0954 0.1771 0.5899       1     -0.2516 0.4425
#> 5   Item5 random_dif -0.0026 0.1690 0.9879       1     -0.3338 0.3286
#> 6   Item6 random_dif -0.1632 0.1660 0.3255       1     -0.4886 0.1622
#> 7   Item7 random_dif -0.1111 0.1718 0.5178       1     -0.4478 0.2256
#> 8   Item8 random_dif -0.0013 0.1713 0.9937       1     -0.3370 0.3343
#> 9   Item9 random_dif -0.0902 0.1673 0.5899       1     -0.4182 0.2378
#> 10 Item10 random_dif -0.0606 0.1595 0.7040       1     -0.3732 0.2520
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2124 0.1629 0.1922       1     -0.1068 0.5317
#> 2   Item2 random_dif -0.1432 0.1692 0.3973       1     -0.4748 0.1884
#> 3   Item3 random_dif -0.0040 0.1761 0.9820       1     -0.3492 0.3412
#> 4   Item4 random_dif -0.1008 0.1673 0.5469       1     -0.4286 0.2271
#> 5   Item5 random_dif  0.2132 0.1716 0.2141       1     -0.1231 0.5494
#> 6   Item6 random_dif  0.0197 0.1741 0.9098       1     -0.3215 0.3609
#> 7   Item7 random_dif -0.1373 0.1657 0.4075       1     -0.4620 0.1875
#> 8   Item8 random_dif -0.1642 0.1777 0.3554       1     -0.5125 0.1840
#> 9   Item9 random_dif  0.1556 0.1797 0.3865       1     -0.1966 0.5079
#> 10 Item10 random_dif -0.0215 0.1669 0.8975       1     -0.3486 0.3056
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2133 0.1702 0.2101       1     -0.1203 0.5468
#> 2   Item2 random_dif  0.0280 0.1671 0.8668       1     -0.2996 0.3556
#> 3   Item3 random_dif -0.1242 0.1683 0.4606       1     -0.4540 0.2057
#> 4   Item4 random_dif -0.1086 0.1712 0.5260       1     -0.4441 0.2270
#> 5   Item5 random_dif  0.0646 0.1704 0.7048       1     -0.2695 0.3986
#> 6   Item6 random_dif  0.0633 0.1742 0.7166       1     -0.2782 0.4047
#> 7   Item7 random_dif -0.0145 0.1637 0.9295       1     -0.3354 0.3064
#> 8   Item8 random_dif -0.0839 0.1658 0.6130       1     -0.4089 0.2411
#> 9   Item9 random_dif -0.0498 0.1621 0.7587       1     -0.3676 0.2680
#> 10 Item10 random_dif  0.0281 0.1614 0.8620       1     -0.2883 0.3444
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.3222 0.1556 0.0384  0.3842      0.0172 0.6271
#> 2   Item2 random_dif  0.0756 0.1702 0.6569  1.0000     -0.2580 0.4093
#> 3   Item3 random_dif -0.0387 0.1695 0.8192  1.0000     -0.3709 0.2935
#> 4   Item4 random_dif -0.0215 0.1704 0.8995  1.0000     -0.3555 0.3125
#> 5   Item5 random_dif  0.0639 0.1777 0.7189  1.0000     -0.2843 0.4122
#> 6   Item6 random_dif -0.2372 0.1672 0.1560  1.0000     -0.5650 0.0905
#> 7   Item7 random_dif -0.1862 0.1616 0.2490  1.0000     -0.5029 0.1304
#> 8   Item8 random_dif -0.0936 0.1700 0.5820  1.0000     -0.4268 0.2396
#> 9   Item9 random_dif  0.2977 0.1534 0.0523  0.5234     -0.0030 0.5984
#> 10 Item10 random_dif -0.2039 0.1546 0.1871  1.0000     -0.5069 0.0991
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0760 0.1746 0.6632  1.0000     -0.2661 0.4182
#> 2   Item2 random_dif  0.0256 0.1703 0.8806  1.0000     -0.3083 0.3594
#> 3   Item3 random_dif -0.0574 0.1766 0.7450  1.0000     -0.4035 0.2886
#> 4   Item4 random_dif -0.0110 0.1717 0.9487  1.0000     -0.3475 0.3254
#> 5   Item5 random_dif -0.2677 0.1606 0.0956  0.9564     -0.5825 0.0472
#> 6   Item6 random_dif -0.0349 0.1688 0.8365  1.0000     -0.3658 0.2961
#> 7   Item7 random_dif  0.3060 0.1555 0.0491  0.4913      0.0012 0.6108
#> 8   Item8 random_dif  0.0699 0.1772 0.6932  1.0000     -0.2773 0.4171
#> 9   Item9 random_dif -0.2267 0.1637 0.1662  1.0000     -0.5477 0.0942
#> 10 Item10 random_dif  0.1406 0.1715 0.4124  1.0000     -0.1956 0.4767
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0849 0.1677 0.6127  1.0000     -0.4135 0.2437
#> 2   Item2 random_dif -0.1040 0.1655 0.5298  1.0000     -0.4284 0.2204
#> 3   Item3 random_dif  0.0767 0.1562 0.6231  1.0000     -0.2293 0.3828
#> 4   Item4 random_dif  0.0756 0.1606 0.6380  1.0000     -0.2392 0.3903
#> 5   Item5 random_dif  0.2369 0.1568 0.1308  1.0000     -0.0704 0.5443
#> 6   Item6 random_dif -0.2765 0.1636 0.0909  0.9093     -0.5971 0.0441
#> 7   Item7 random_dif  0.1616 0.1667 0.3326  1.0000     -0.1652 0.4884
#> 8   Item8 random_dif  0.1642 0.1683 0.3294  1.0000     -0.1657 0.4941
#> 9   Item9 random_dif -0.0744 0.1606 0.6432  1.0000     -0.3892 0.2404
#> 10 Item10 random_dif -0.1804 0.1593 0.2575  1.0000     -0.4927 0.1319
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0829 0.1766 0.6387  1.0000     -0.2633 0.4292
#> 2   Item2 random_dif  0.2729 0.1463 0.0621  0.6215     -0.0139 0.5597
#> 3   Item3 random_dif  0.0760 0.1623 0.6398  1.0000     -0.2422 0.3941
#> 4   Item4 random_dif  0.1301 0.1653 0.4312  1.0000     -0.1938 0.4540
#> 5   Item5 random_dif  0.0084 0.1674 0.9601  1.0000     -0.3197 0.3365
#> 6   Item6 random_dif -0.1123 0.1672 0.5016  1.0000     -0.4400 0.2153
#> 7   Item7 random_dif -0.1883 0.1587 0.2356  1.0000     -0.4994 0.1229
#> 8   Item8 random_dif  0.2069 0.1537 0.1783  1.0000     -0.0943 0.5081
#> 9   Item9 random_dif -0.2328 0.1520 0.1257  1.0000     -0.5308 0.0652
#> 10 Item10 random_dif -0.2548 0.1515 0.0926  0.9257     -0.5516 0.0421
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1878 0.1612 0.2442  1.0000     -0.1282  0.5038
#> 2   Item2 random_dif -0.0646 0.1638 0.6930  1.0000     -0.3857  0.2564
#> 3   Item3 random_dif -0.4078 0.1531 0.0078  0.0775   . -0.7079 -0.1076
#> 4   Item4 random_dif  0.0755 0.1709 0.6585  1.0000     -0.2594  0.4105
#> 5   Item5 random_dif -0.0959 0.1666 0.5647  1.0000     -0.4224  0.2305
#> 6   Item6 random_dif  0.2084 0.1667 0.2113  1.0000     -0.1184  0.5351
#> 7   Item7 random_dif  0.0498 0.1818 0.7840  1.0000     -0.3064  0.4061
#> 8   Item8 random_dif  0.0294 0.1780 0.8688  1.0000     -0.3195  0.3784
#> 9   Item9 random_dif -0.0322 0.1705 0.8504  1.0000     -0.3664  0.3021
#> 10 Item10 random_dif  0.0217 0.1637 0.8946  1.0000     -0.2992  0.3425
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3645 0.1462 0.0126  0.1262     -0.6510 -0.0781
#> 2   Item2 random_dif -0.0852 0.1711 0.6184  1.0000     -0.4205  0.2501
#> 3   Item3 random_dif -0.0899 0.1689 0.5943  1.0000     -0.4209  0.2410
#> 4   Item4 random_dif  0.3350 0.1494 0.0250  0.2499      0.0421  0.6278
#> 5   Item5 random_dif -0.0027 0.1750 0.9877  1.0000     -0.3457  0.3402
#> 6   Item6 random_dif -0.1239 0.1791 0.4892  1.0000     -0.4750  0.2272
#> 7   Item7 random_dif  0.0000 0.1793 1.0000  1.0000     -0.3515  0.3515
#> 8   Item8 random_dif  0.3615 0.1520 0.0174  0.1741      0.0635  0.6595
#> 9   Item9 random_dif -0.1826 0.1689 0.2797  1.0000     -0.5136  0.1484
#> 10 Item10 random_dif  0.1179 0.1635 0.4708  1.0000     -0.2026  0.4385
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2737 0.1842 0.1372       1     -0.6347 0.0872
#> 2   Item2 random_dif  0.1590 0.1680 0.3441       1     -0.1704 0.4883
#> 3   Item3 random_dif  0.1596 0.1683 0.3429       1     -0.1702 0.4895
#> 4   Item4 random_dif -0.0613 0.1805 0.7343       1     -0.4151 0.2926
#> 5   Item5 random_dif  0.0999 0.1737 0.5653       1     -0.2405 0.4402
#> 6   Item6 random_dif  0.2301 0.1700 0.1759       1     -0.1031 0.5633
#> 7   Item7 random_dif -0.1221 0.1727 0.4796       1     -0.4606 0.2164
#> 8   Item8 random_dif -0.0912 0.1800 0.6123       1     -0.4439 0.2615
#> 9   Item9 random_dif -0.2167 0.1651 0.1892       1     -0.5402 0.1068
#> 10 Item10 random_dif  0.0598 0.1706 0.7258       1     -0.2745 0.3941
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0043 0.1737 0.9803  1.0000     -0.3447  0.3361
#> 2   Item2 random_dif -0.0421 0.1621 0.7949  1.0000     -0.3599  0.2756
#> 3   Item3 random_dif -0.3310 0.1800 0.0659  0.6594     -0.6839  0.0218
#> 4   Item4 random_dif  0.0912 0.1786 0.6097  1.0000     -0.2589  0.4413
#> 5   Item5 random_dif -0.3437 0.1629 0.0349  0.3489     -0.6629 -0.0244
#> 6   Item6 random_dif -0.0454 0.1747 0.7950  1.0000     -0.3878  0.2970
#> 7   Item7 random_dif  0.0026 0.1701 0.9876  1.0000     -0.3308  0.3361
#> 8   Item8 random_dif  0.1977 0.1674 0.2376  1.0000     -0.1304  0.5259
#> 9   Item9 random_dif  0.0743 0.1622 0.6468  1.0000     -0.2435  0.3921
#> 10 Item10 random_dif  0.2535 0.1489 0.0887  0.8874     -0.0384  0.5454
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0346 0.1777 0.8456  1.0000     -0.3137 0.3829
#> 2   Item2 random_dif  0.1127 0.1641 0.4925  1.0000     -0.2090 0.4343
#> 3   Item3 random_dif  0.2051 0.1676 0.2211  1.0000     -0.1234 0.5336
#> 4   Item4 random_dif -0.0912 0.1702 0.5922  1.0000     -0.4248 0.2425
#> 5   Item5 random_dif  0.3625 0.1504 0.0160  0.1596      0.0677 0.6574
#> 6   Item6 random_dif -0.2389 0.1731 0.1675  1.0000     -0.5780 0.1003
#> 7   Item7 random_dif -0.0597 0.1709 0.7267  1.0000     -0.3947 0.2752
#> 8   Item8 random_dif -0.0649 0.1857 0.7265  1.0000     -0.4288 0.2990
#> 9   Item9 random_dif -0.2032 0.1678 0.2259  1.0000     -0.5320 0.1256
#> 10 Item10 random_dif -0.1138 0.1708 0.5053  1.0000     -0.4486 0.2210
cutoff_res$item_cutoffs
#>      Item  gamma_low gamma_high
#> 1   Item1 -0.5132743  0.3516545
#> 2   Item2 -0.3015007  0.3881181
#> 3   Item3 -0.4077519  0.4064600
#> 4   Item4 -0.3796940  0.4456095
#> 5   Item5 -0.3492063  0.3489032
#> 6   Item6 -0.3526854  0.4207077
#> 7   Item7 -0.4129353  0.4662202
#> 8   Item8 -0.3184358  0.4069738
#> 9   Item9 -0.4386159  0.2977099
#> 10 Item10 -0.3595801  0.3873950
# }
```
