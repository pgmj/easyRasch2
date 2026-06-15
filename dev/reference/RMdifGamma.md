# Partial Gamma DIF Analysis

Computes partial gamma coefficients for Differential Item Functioning
(DIF) using
[`iarm::partgam_DIF()`](https://rdrr.io/pkg/iarm/man/partgam_DIF.html).
Each item is tested for association with a single categorical DIF
variable, controlling for the total score.

## Usage

``` r
RMdifGamma(data, dif_var, cutoff = NULL, output = "kable")
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
  high", and "Flagged" are included.

- If `output = "dataframe"`: a data.frame with columns `Item`, `gamma`,
  `se`, `lower`, `upper`, `padj_bh`, `Significance`. When `cutoff` is
  provided, columns `gamma_low`, `gamma_high`, and `flagged` are also
  included.

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

## References

Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
(1998). Differential item functioning in the Danish translation of the
SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189–1202.
[doi:10.1016/S0895-4356(98)00111-5](https://doi.org/10.1016/S0895-4356%2898%2900111-5)

## See also

[`RMdifGammaCutoff`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md)

## Examples

``` r
# \donttest{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)
dif_group <- factor(sample(c("A", "B"), 200, replace = TRUE))

# Default kable output
RMdifGamma(sim_data, dif_group)
#> 
#> 
#> Table: Partial gamma DIF analysis (n = 200 complete cases). Positive gamma indicates higher scores in higher DIF group levels.
#> 
#> |Item   | Partial gamma|    SE| Lower CI| Upper CI| Adj. p-value (BH)|p-value sign. |
#> |:------|-------------:|-----:|--------:|--------:|-----------------:|:-------------|
#> |Item1  |         0.029| 0.161|   -0.287|    0.345|                 1|              |
#> |Item2  |        -0.168| 0.156|   -0.474|    0.138|                 1|              |
#> |Item3  |         0.111| 0.168|   -0.219|    0.441|                 1|              |
#> |Item4  |        -0.030| 0.155|   -0.333|    0.273|                 1|              |
#> |Item5  |         0.076| 0.159|   -0.235|    0.387|                 1|              |
#> |Item6  |         0.042| 0.156|   -0.264|    0.349|                 1|              |
#> |Item7  |        -0.146| 0.160|   -0.460|    0.169|                 1|              |
#> |Item8  |         0.019| 0.162|   -0.298|    0.336|                 1|              |
#> |Item9  |        -0.141| 0.161|   -0.457|    0.175|                 1|              |
#> |Item10 |         0.178| 0.156|   -0.127|    0.483|                 1|              |

# Return as data.frame
RMdifGamma(sim_data, dif_group, output = "dataframe")
#>      Item  gamma    se  lower upper padj_bh Significance
#> 1   Item1  0.029 0.161 -0.287 0.345       1             
#> 2   Item2 -0.168 0.156 -0.474 0.138       1             
#> 3   Item3  0.111 0.168 -0.219 0.441       1             
#> 4   Item4 -0.030 0.155 -0.333 0.273       1             
#> 5   Item5  0.076 0.159 -0.235 0.387       1             
#> 6   Item6  0.042 0.156 -0.264 0.349       1             
#> 7   Item7 -0.146 0.160 -0.460 0.169       1             
#> 8   Item8  0.019 0.162 -0.298 0.336       1             
#> 9   Item9 -0.141 0.161 -0.457 0.175       1             
#> 10 Item10  0.178 0.156 -0.127 0.483       1             
# }
# \donttest{
# Simulation-based cutoffs (100 Monte-Carlo iterations)
cutoff_res <- RMdifGammaCutoff(sim_data, dif_var = dif_group,
                                  iterations = 100, parallel = FALSE,
                                  seed = 42)
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0247 0.1819 0.8921  1.0000     -0.3318 0.3811
#> 2   Item2 random_dif  0.1753 0.1628 0.2816  1.0000     -0.1438 0.4944
#> 3   Item3 random_dif -0.1813 0.1680 0.2806  1.0000     -0.5105 0.1480
#> 4   Item4 random_dif  0.1566 0.1810 0.3870  1.0000     -0.1982 0.5114
#> 5   Item5 random_dif  0.0747 0.1724 0.6645  1.0000     -0.2631 0.4126
#> 6   Item6 random_dif -0.3052 0.1605 0.0572  0.5723     -0.6197 0.0094
#> 7   Item7 random_dif -0.1742 0.1672 0.2973  1.0000     -0.5019 0.1534
#> 8   Item8 random_dif  0.1226 0.1733 0.4793  1.0000     -0.2171 0.4623
#> 9   Item9 random_dif  0.0496 0.1747 0.7766  1.0000     -0.2928 0.3919
#> 10 Item10 random_dif  0.0775 0.1717 0.6517  1.0000     -0.2590 0.4140
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1070 0.1758 0.5427  1.0000     -0.2375  0.4514
#> 2   Item2 random_dif  0.0604 0.1760 0.7314  1.0000     -0.2845  0.4054
#> 3   Item3 random_dif  0.4021 0.1522 0.0082  0.0824   .  0.1038  0.7004
#> 4   Item4 random_dif -0.5154 0.1493 0.0006  0.0056  ** -0.8080 -0.2228
#> 5   Item5 random_dif  0.0520 0.1697 0.7593  1.0000     -0.2806  0.3845
#> 6   Item6 random_dif  0.2083 0.1772 0.2397  1.0000     -0.1390  0.5556
#> 7   Item7 random_dif -0.4915 0.1313 0.0002  0.0018  ** -0.7488 -0.2342
#> 8   Item8 random_dif -0.1761 0.1669 0.2913  1.0000     -0.5033  0.1510
#> 9   Item9 random_dif  0.2184 0.1624 0.1785  1.0000     -0.0998  0.5366
#> 10 Item10 random_dif  0.0769 0.1710 0.6529  1.0000     -0.2583  0.4121
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.4376 0.1380 0.0015  0.0153   *  0.1670 0.7082
#> 2   Item2 random_dif  0.0025 0.1644 0.9879  1.0000     -0.3198 0.3248
#> 3   Item3 random_dif  0.1741 0.1667 0.2963  1.0000     -0.1526 0.5009
#> 4   Item4 random_dif -0.0300 0.1725 0.8621  1.0000     -0.3681 0.3082
#> 5   Item5 random_dif -0.0761 0.1699 0.6542  1.0000     -0.4091 0.2569
#> 6   Item6 random_dif  0.1062 0.1731 0.5395  1.0000     -0.2330 0.4454
#> 7   Item7 random_dif -0.1807 0.1640 0.2706  1.0000     -0.5021 0.1407
#> 8   Item8 random_dif -0.1832 0.1647 0.2659  1.0000     -0.5059 0.1396
#> 9   Item9 random_dif -0.0907 0.1646 0.5817  1.0000     -0.4132 0.2319
#> 10 Item10 random_dif -0.1641 0.1619 0.3107  1.0000     -0.4815 0.1532
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0567 0.1751 0.7460       1     -0.2865 0.4000
#> 2   Item2 random_dif -0.1634 0.1598 0.3064       1     -0.4765 0.1497
#> 3   Item3 random_dif  0.0317 0.1772 0.8578       1     -0.3156 0.3791
#> 4   Item4 random_dif -0.0104 0.1843 0.9549       1     -0.3716 0.3507
#> 5   Item5 random_dif  0.0106 0.1719 0.9510       1     -0.3264 0.3475
#> 6   Item6 random_dif  0.0610 0.1720 0.7227       1     -0.2760 0.3981
#> 7   Item7 random_dif  0.0802 0.1732 0.6431       1     -0.2592 0.4197
#> 8   Item8 random_dif -0.0347 0.1777 0.8452       1     -0.3831 0.3136
#> 9   Item9 random_dif  0.0221 0.1643 0.8930       1     -0.2999 0.3441
#> 10 Item10 random_dif -0.0230 0.1657 0.8896       1     -0.3478 0.3018
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1098 0.1678 0.5127  1.0000     -0.4387 0.2190
#> 2   Item2 random_dif  0.1709 0.1619 0.2912  1.0000     -0.1464 0.4881
#> 3   Item3 random_dif -0.2038 0.1591 0.2004  1.0000     -0.5157 0.1081
#> 4   Item4 random_dif -0.0894 0.1663 0.5907  1.0000     -0.4153 0.2364
#> 5   Item5 random_dif  0.2580 0.1592 0.1050  1.0000     -0.0539 0.5700
#> 6   Item6 random_dif  0.1075 0.1681 0.5224  1.0000     -0.2219 0.4369
#> 7   Item7 random_dif -0.0921 0.1654 0.5775  1.0000     -0.4162 0.2320
#> 8   Item8 random_dif  0.1183 0.1716 0.4905  1.0000     -0.2180 0.4547
#> 9   Item9 random_dif -0.2685 0.1525 0.0782  0.7825     -0.5674 0.0303
#> 10 Item10 random_dif  0.1285 0.1599 0.4214  1.0000     -0.1848 0.4419
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1590 0.1870 0.3952  1.0000     -0.5254 0.2075
#> 2   Item2 random_dif  0.0126 0.1752 0.9426  1.0000     -0.3307 0.3559
#> 3   Item3 random_dif -0.0301 0.1741 0.8626  1.0000     -0.3714 0.3111
#> 4   Item4 random_dif  0.3333 0.1722 0.0529  0.5287     -0.0041 0.6708
#> 5   Item5 random_dif -0.1214 0.1736 0.4845  1.0000     -0.4616 0.2189
#> 6   Item6 random_dif  0.2497 0.1694 0.1405  1.0000     -0.0823 0.5816
#> 7   Item7 random_dif -0.1584 0.1730 0.3599  1.0000     -0.4974 0.1807
#> 8   Item8 random_dif -0.1603 0.1698 0.3451  1.0000     -0.4930 0.1724
#> 9   Item9 random_dif  0.1350 0.1731 0.4355  1.0000     -0.2042 0.4742
#> 10 Item10 random_dif -0.0996 0.1693 0.5563  1.0000     -0.4314 0.2322
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0365 0.1728 0.8330  1.0000     -0.3023 0.3752
#> 2   Item2 random_dif -0.2139 0.1717 0.2127  1.0000     -0.5504 0.1225
#> 3   Item3 random_dif  0.0547 0.1763 0.7563  1.0000     -0.2909 0.4003
#> 4   Item4 random_dif -0.0673 0.1740 0.6989  1.0000     -0.4083 0.2737
#> 5   Item5 random_dif -0.2909 0.1579 0.0655  0.6548     -0.6005 0.0186
#> 6   Item6 random_dif -0.0853 0.1735 0.6231  1.0000     -0.4254 0.2548
#> 7   Item7 random_dif  0.1724 0.1637 0.2922  1.0000     -0.1484 0.4933
#> 8   Item8 random_dif  0.1903 0.1597 0.2332  1.0000     -0.1226 0.5033
#> 9   Item9 random_dif  0.1172 0.1662 0.4806  1.0000     -0.2085 0.4430
#> 10 Item10 random_dif  0.0458 0.1636 0.7795  1.0000     -0.2749 0.3665
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0634 0.1884 0.7364       1     -0.3058 0.4327
#> 2   Item2 random_dif  0.0308 0.1715 0.8573       1     -0.3052 0.3669
#> 3   Item3 random_dif -0.2064 0.1732 0.2334       1     -0.5459 0.1331
#> 4   Item4 random_dif -0.1513 0.1646 0.3578       1     -0.4739 0.1712
#> 5   Item5 random_dif  0.1060 0.1681 0.5282       1     -0.2234 0.4354
#> 6   Item6 random_dif  0.1235 0.1671 0.4597       1     -0.2040 0.4510
#> 7   Item7 random_dif  0.1433 0.1746 0.4119       1     -0.1990 0.4856
#> 8   Item8 random_dif -0.1952 0.1688 0.2475       1     -0.5260 0.1357
#> 9   Item9 random_dif  0.1922 0.1664 0.2481       1     -0.1339 0.5183
#> 10 Item10 random_dif -0.0944 0.1672 0.5725       1     -0.4222 0.2334
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0031 0.1788 0.9862       1     -0.3473 0.3534
#> 2   Item2 random_dif -0.0266 0.1655 0.8722       1     -0.3510 0.2978
#> 3   Item3 random_dif -0.1272 0.1721 0.4600       1     -0.4645 0.2102
#> 4   Item4 random_dif  0.0218 0.1629 0.8934       1     -0.2974 0.3411
#> 5   Item5 random_dif  0.1776 0.1686 0.2923       1     -0.1529 0.5081
#> 6   Item6 random_dif -0.2483 0.1583 0.1167       1     -0.5586 0.0619
#> 7   Item7 random_dif  0.0000 0.1695 1.0000       1     -0.3321 0.3321
#> 8   Item8 random_dif  0.1682 0.1642 0.3057       1     -0.1537 0.4901
#> 9   Item9 random_dif -0.0517 0.1684 0.7586       1     -0.3817 0.2782
#> 10 Item10 random_dif  0.0709 0.1647 0.6670       1     -0.2519 0.3937
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0791 0.1731 0.6476  1.0000     -0.4185  0.2602
#> 2   Item2 random_dif  0.0627 0.1627 0.7002  1.0000     -0.2563  0.3816
#> 3   Item3 random_dif  0.1643 0.1569 0.2949  1.0000     -0.1431  0.4717
#> 4   Item4 random_dif -0.3130 0.1519 0.0393  0.3929     -0.6106 -0.0154
#> 5   Item5 random_dif  0.0400 0.1651 0.8086  1.0000     -0.2836  0.3636
#> 6   Item6 random_dif  0.0115 0.1702 0.9460  1.0000     -0.3221  0.3452
#> 7   Item7 random_dif  0.1975 0.1698 0.2448  1.0000     -0.1353  0.5304
#> 8   Item8 random_dif  0.0060 0.1741 0.9726  1.0000     -0.3353  0.3472
#> 9   Item9 random_dif -0.0759 0.1576 0.6299  1.0000     -0.3849  0.2330
#> 10 Item10 random_dif  0.0094 0.1674 0.9550  1.0000     -0.3187  0.3376
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0205 0.1816 0.9102  1.0000     -0.3354  0.3764
#> 2   Item2 random_dif -0.3299 0.1509 0.0289  0.2887     -0.6257 -0.0340
#> 3   Item3 random_dif -0.0519 0.1693 0.7591  1.0000     -0.3836  0.2798
#> 4   Item4 random_dif -0.1699 0.1811 0.3480  1.0000     -0.5248  0.1850
#> 5   Item5 random_dif  0.0056 0.1719 0.9740  1.0000     -0.3313  0.3425
#> 6   Item6 random_dif  0.2142 0.1667 0.1988  1.0000     -0.1125  0.5409
#> 7   Item7 random_dif -0.0208 0.1692 0.9021  1.0000     -0.3523  0.3107
#> 8   Item8 random_dif -0.0934 0.1756 0.5947  1.0000     -0.4375  0.2507
#> 9   Item9 random_dif  0.3107 0.1728 0.0722  0.7224     -0.0281  0.6494
#> 10 Item10 random_dif  0.1421 0.1593 0.3723  1.0000     -0.1701  0.4544
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1672 0.1671 0.3173  1.0000     -0.1604 0.4947
#> 2   Item2 random_dif  0.2917 0.1494 0.0509  0.5086     -0.0011 0.5845
#> 3   Item3 random_dif -0.0493 0.1771 0.7808  1.0000     -0.3964 0.2978
#> 4   Item4 random_dif -0.2619 0.1659 0.1143  1.0000     -0.5871 0.0632
#> 5   Item5 random_dif -0.2822 0.1587 0.0753  0.7534     -0.5932 0.0288
#> 6   Item6 random_dif -0.1796 0.1655 0.2778  1.0000     -0.5038 0.1447
#> 7   Item7 random_dif -0.0720 0.1591 0.6510  1.0000     -0.3839 0.2399
#> 8   Item8 random_dif  0.1948 0.1637 0.2341  1.0000     -0.1261 0.5157
#> 9   Item9 random_dif -0.0938 0.1732 0.5879  1.0000     -0.4333 0.2456
#> 10 Item10 random_dif  0.2231 0.1589 0.1602  1.0000     -0.0883 0.5345
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0628 0.1823 0.7304  1.0000     -0.4202  0.2945
#> 2   Item2 random_dif -0.0125 0.1740 0.9426  1.0000     -0.3536  0.3286
#> 3   Item3 random_dif -0.3577 0.1544 0.0205  0.2052     -0.6603 -0.0551
#> 4   Item4 random_dif  0.0724 0.1717 0.6732  1.0000     -0.2641  0.4089
#> 5   Item5 random_dif  0.1650 0.1668 0.3227  1.0000     -0.1620  0.4919
#> 6   Item6 random_dif  0.1189 0.1635 0.4671  1.0000     -0.2015  0.4393
#> 7   Item7 random_dif  0.1777 0.1690 0.2929  1.0000     -0.1535  0.5089
#> 8   Item8 random_dif  0.1923 0.1657 0.2460  1.0000     -0.1325  0.5170
#> 9   Item9 random_dif -0.2607 0.1592 0.1015  1.0000     -0.5728  0.0514
#> 10 Item10 random_dif -0.0376 0.1635 0.8183  1.0000     -0.3580  0.2829
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.2092 0.1841 0.2560   1.000     -0.1517  0.5700
#> 2   Item2 random_dif  0.1074 0.1682 0.5231   1.000     -0.2222  0.4370
#> 3   Item3 random_dif  0.0638 0.1677 0.7037   1.000     -0.2649  0.3924
#> 4   Item4 random_dif -0.2630 0.1620 0.1044   1.000     -0.5804  0.0544
#> 5   Item5 random_dif  0.1304 0.1676 0.4365   1.000     -0.1981  0.4590
#> 6   Item6 random_dif -0.3994 0.1548 0.0099   0.099   . -0.7028 -0.0959
#> 7   Item7 random_dif -0.1756 0.1659 0.2898   1.000     -0.5007  0.1495
#> 8   Item8 random_dif  0.1844 0.1639 0.2606   1.000     -0.1369  0.5057
#> 9   Item9 random_dif  0.0046 0.1736 0.9787   1.000     -0.3356  0.3449
#> 10 Item10 random_dif  0.1429 0.1634 0.3820   1.000     -0.1774  0.4631
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1009 0.1642 0.5390  1.0000     -0.4228 0.2210
#> 2   Item2 random_dif -0.1056 0.1702 0.5350  1.0000     -0.4392 0.2280
#> 3   Item3 random_dif -0.0413 0.1767 0.8151  1.0000     -0.3876 0.3049
#> 4   Item4 random_dif  0.1636 0.1767 0.3543  1.0000     -0.1826 0.5099
#> 5   Item5 random_dif  0.1502 0.1694 0.3754  1.0000     -0.1818 0.4822
#> 6   Item6 random_dif  0.0165 0.1851 0.9292  1.0000     -0.3464 0.3793
#> 7   Item7 random_dif -0.0509 0.1780 0.7749  1.0000     -0.3997 0.2979
#> 8   Item8 random_dif -0.2763 0.1572 0.0788  0.7875     -0.5843 0.0317
#> 9   Item9 random_dif  0.0602 0.1693 0.7220  1.0000     -0.2716 0.3921
#> 10 Item10 random_dif  0.2036 0.1619 0.2086  1.0000     -0.1137 0.5210
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2168 0.1670 0.1942       1     -0.5442 0.1105
#> 2   Item2 random_dif  0.2345 0.1634 0.1512       1     -0.0858 0.5547
#> 3   Item3 random_dif -0.0554 0.1662 0.7390       1     -0.3812 0.2704
#> 4   Item4 random_dif -0.0131 0.1697 0.9385       1     -0.3456 0.3194
#> 5   Item5 random_dif -0.0590 0.1622 0.7160       1     -0.3768 0.2588
#> 6   Item6 random_dif -0.1939 0.1644 0.2381       1     -0.5161 0.1282
#> 7   Item7 random_dif -0.0696 0.1783 0.6963       1     -0.4192 0.2799
#> 8   Item8 random_dif  0.2123 0.1729 0.2194       1     -0.1265 0.5511
#> 9   Item9 random_dif  0.1971 0.1668 0.2376       1     -0.1299 0.5241
#> 10 Item10 random_dif -0.0316 0.1599 0.8432       1     -0.3450 0.2818
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0919 0.1695 0.5879  1.0000     -0.2404  0.4241
#> 2   Item2 random_dif -0.1242 0.1673 0.4577  1.0000     -0.4520  0.2036
#> 3   Item3 random_dif -0.4552 0.1426 0.0014  0.0141   * -0.7348 -0.1757
#> 4   Item4 random_dif  0.0465 0.1627 0.7750  1.0000     -0.2723  0.3653
#> 5   Item5 random_dif  0.0165 0.1615 0.9186  1.0000     -0.3000  0.3330
#> 6   Item6 random_dif  0.0409 0.1760 0.8164  1.0000     -0.3042  0.3859
#> 7   Item7 random_dif  0.1041 0.1665 0.5316  1.0000     -0.2222  0.4305
#> 8   Item8 random_dif  0.1076 0.1651 0.5147  1.0000     -0.2161  0.4312
#> 9   Item9 random_dif  0.0046 0.1617 0.9772  1.0000     -0.3124  0.3216
#> 10 Item10 random_dif  0.1494 0.1564 0.3394  1.0000     -0.1571  0.4560
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0935 0.1674 0.5766  1.0000     -0.2347  0.4216
#> 2   Item2 random_dif  0.0630 0.1709 0.7126  1.0000     -0.2721  0.3980
#> 3   Item3 random_dif -0.0078 0.1626 0.9620  1.0000     -0.3265  0.3110
#> 4   Item4 random_dif  0.2368 0.1611 0.1416  1.0000     -0.0789  0.5524
#> 5   Item5 random_dif -0.0187 0.1695 0.9124  1.0000     -0.3508  0.3135
#> 6   Item6 random_dif -0.0349 0.1784 0.8448  1.0000     -0.3846  0.3147
#> 7   Item7 random_dif  0.0260 0.1736 0.8809  1.0000     -0.3143  0.3664
#> 8   Item8 random_dif -0.3450 0.1436 0.0162  0.1624     -0.6264 -0.0637
#> 9   Item9 random_dif  0.0165 0.1694 0.9225  1.0000     -0.3155  0.3484
#> 10 Item10 random_dif  0.0093 0.1630 0.9544  1.0000     -0.3101  0.3287
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0265 0.1781 0.8815       1     -0.3225 0.3756
#> 2   Item2 random_dif -0.0190 0.1608 0.9059       1     -0.3341 0.2961
#> 3   Item3 random_dif  0.0645 0.1667 0.6990       1     -0.2623 0.3913
#> 4   Item4 random_dif -0.0628 0.1761 0.7212       1     -0.4080 0.2823
#> 5   Item5 random_dif -0.0483 0.1659 0.7709       1     -0.3736 0.2769
#> 6   Item6 random_dif  0.2052 0.1663 0.2173       1     -0.1208 0.5312
#> 7   Item7 random_dif -0.0486 0.1717 0.7770       1     -0.3851 0.2879
#> 8   Item8 random_dif -0.1052 0.1627 0.5179       1     -0.4241 0.2137
#> 9   Item9 random_dif  0.1400 0.1651 0.3963       1     -0.1835 0.4636
#> 10 Item10 random_dif -0.1373 0.1617 0.3958       1     -0.4542 0.1796
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1248 0.1799 0.4877  1.0000     -0.4773 0.2277
#> 2   Item2 random_dif -0.0970 0.1660 0.5591  1.0000     -0.4222 0.2283
#> 3   Item3 random_dif  0.3028 0.1505 0.0442  0.4421      0.0078 0.5978
#> 4   Item4 random_dif  0.0687 0.1763 0.6968  1.0000     -0.2769 0.4144
#> 5   Item5 random_dif  0.1178 0.1661 0.4782  1.0000     -0.2078 0.4434
#> 6   Item6 random_dif  0.2612 0.1617 0.1063  1.0000     -0.0557 0.5781
#> 7   Item7 random_dif  0.2430 0.1665 0.1445  1.0000     -0.0834 0.5694
#> 8   Item8 random_dif -0.2778 0.1596 0.0818  0.8182     -0.5906 0.0351
#> 9   Item9 random_dif -0.2071 0.1562 0.1849  1.0000     -0.5133 0.0991
#> 10 Item10 random_dif -0.2691 0.1529 0.0785  0.7848     -0.5688 0.0306
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1449 0.1810 0.4233       1     -0.4996 0.2098
#> 2   Item2 random_dif -0.0283 0.1655 0.8644       1     -0.3527 0.2962
#> 3   Item3 random_dif -0.0460 0.1625 0.7769       1     -0.3645 0.2724
#> 4   Item4 random_dif  0.1098 0.1671 0.5113       1     -0.2178 0.4374
#> 5   Item5 random_dif  0.0244 0.1628 0.8809       1     -0.2946 0.3434
#> 6   Item6 random_dif  0.0645 0.1608 0.6882       1     -0.2506 0.3797
#> 7   Item7 random_dif -0.1258 0.1721 0.4649       1     -0.4632 0.2116
#> 8   Item8 random_dif  0.1182 0.1628 0.4679       1     -0.2009 0.4372
#> 9   Item9 random_dif -0.0096 0.1835 0.9582       1     -0.3692 0.3500
#> 10 Item10 random_dif -0.0136 0.1637 0.9337       1     -0.3344 0.3072
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1801 0.1683 0.2844  1.0000     -0.5099  0.1496
#> 2   Item2 random_dif  0.0638 0.1752 0.7159  1.0000     -0.2796  0.4072
#> 3   Item3 random_dif  0.3828 0.1608 0.0173  0.1728      0.0677  0.6979
#> 4   Item4 random_dif -0.0291 0.1802 0.8717  1.0000     -0.3823  0.3241
#> 5   Item5 random_dif  0.3149 0.1554 0.0427  0.4271      0.0104  0.6195
#> 6   Item6 random_dif -0.2524 0.1613 0.1178  1.0000     -0.5686  0.0638
#> 7   Item7 random_dif -0.3137 0.1498 0.0362  0.3621     -0.6072 -0.0202
#> 8   Item8 random_dif -0.0100 0.1869 0.9573  1.0000     -0.3763  0.3563
#> 9   Item9 random_dif  0.0028 0.1721 0.9871  1.0000     -0.3346  0.3402
#> 10 Item10 random_dif  0.0341 0.1650 0.8361  1.0000     -0.2892  0.3575
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2764 0.1597 0.0835  0.8348     -0.5894 0.0366
#> 2   Item2 random_dif  0.0221 0.1663 0.8944  1.0000     -0.3039 0.3480
#> 3   Item3 random_dif  0.0182 0.1687 0.9142  1.0000     -0.3124 0.3487
#> 4   Item4 random_dif -0.1581 0.1696 0.3512  1.0000     -0.4906 0.1743
#> 5   Item5 random_dif  0.1106 0.1636 0.4990  1.0000     -0.2100 0.4313
#> 6   Item6 random_dif  0.2847 0.1584 0.0723  0.7228     -0.0258 0.5952
#> 7   Item7 random_dif -0.0515 0.1687 0.7602  1.0000     -0.3821 0.2791
#> 8   Item8 random_dif -0.0336 0.1631 0.8369  1.0000     -0.3534 0.2862
#> 9   Item9 random_dif  0.2436 0.1564 0.1193  1.0000     -0.0629 0.5501
#> 10 Item10 random_dif -0.1652 0.1576 0.2944  1.0000     -0.4740 0.1436
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1543 0.1638 0.3464       1     -0.1668 0.4754
#> 2   Item2 random_dif -0.2215 0.1571 0.1584       1     -0.5293 0.0863
#> 3   Item3 random_dif -0.0763 0.1625 0.6385       1     -0.3947 0.2421
#> 4   Item4 random_dif  0.0196 0.1654 0.9059       1     -0.3045 0.3436
#> 5   Item5 random_dif  0.0212 0.1603 0.8950       1     -0.2931 0.3354
#> 6   Item6 random_dif -0.2187 0.1681 0.1932       1     -0.5481 0.1107
#> 7   Item7 random_dif  0.2400 0.1653 0.1466       1     -0.0840 0.5640
#> 8   Item8 random_dif -0.0662 0.1644 0.6874       1     -0.3885 0.2561
#> 9   Item9 random_dif -0.0102 0.1626 0.9502       1     -0.3288 0.3085
#> 10 Item10 random_dif  0.1736 0.1571 0.2693       1     -0.1344 0.4815
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.3345 0.1629 0.0401  0.4005      0.0152 0.6538
#> 2   Item2 random_dif  0.0201 0.1652 0.9032  1.0000     -0.3037 0.3439
#> 3   Item3 random_dif -0.0665 0.1663 0.6895  1.0000     -0.3925 0.2596
#> 4   Item4 random_dif  0.0276 0.1734 0.8738  1.0000     -0.3123 0.3675
#> 5   Item5 random_dif -0.2622 0.1528 0.0862  0.8623     -0.5618 0.0373
#> 6   Item6 random_dif  0.2658 0.1602 0.0970  0.9703     -0.0481 0.5797
#> 7   Item7 random_dif -0.0363 0.1642 0.8253  1.0000     -0.3582 0.2856
#> 8   Item8 random_dif -0.1682 0.1798 0.3495  1.0000     -0.5205 0.1842
#> 9   Item9 random_dif -0.2685 0.1531 0.0795  0.7949     -0.5686 0.0316
#> 10 Item10 random_dif  0.2214 0.1606 0.1682  1.0000     -0.0935 0.5362
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0461 0.1731 0.7899  1.0000     -0.2931 0.3853
#> 2   Item2 random_dif -0.0026 0.1678 0.9875  1.0000     -0.3316 0.3263
#> 3   Item3 random_dif  0.0828 0.1748 0.6356  1.0000     -0.2598 0.4255
#> 4   Item4 random_dif -0.0565 0.1732 0.7444  1.0000     -0.3960 0.2831
#> 5   Item5 random_dif  0.2677 0.1550 0.0842  0.8422     -0.0362 0.5716
#> 6   Item6 random_dif -0.0122 0.1923 0.9495  1.0000     -0.3891 0.3648
#> 7   Item7 random_dif -0.0516 0.1801 0.7743  1.0000     -0.4047 0.3014
#> 8   Item8 random_dif -0.0753 0.1723 0.6620  1.0000     -0.4131 0.2625
#> 9   Item9 random_dif -0.1681 0.1797 0.3496  1.0000     -0.5202 0.1841
#> 10 Item10 random_dif -0.0819 0.1663 0.6224  1.0000     -0.4077 0.2440
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.3496 0.1729 0.0432  0.4318      0.0107  0.6885
#> 2   Item2 random_dif  0.3015 0.1545 0.0510  0.5105     -0.0014  0.6044
#> 3   Item3 random_dif -0.0504 0.1799 0.7793  1.0000     -0.4030  0.3022
#> 4   Item4 random_dif  0.2904 0.1745 0.0961  0.9613     -0.0517  0.6325
#> 5   Item5 random_dif  0.0333 0.1739 0.8480  1.0000     -0.3075  0.3742
#> 6   Item6 random_dif -0.2268 0.1612 0.1593  1.0000     -0.5427  0.0891
#> 7   Item7 random_dif -0.3037 0.1737 0.0803  0.8034     -0.6442  0.0367
#> 8   Item8 random_dif -0.3134 0.1547 0.0428  0.4279     -0.6165 -0.0102
#> 9   Item9 random_dif -0.1514 0.1638 0.3555  1.0000     -0.4724  0.1697
#> 10 Item10 random_dif  0.1634 0.1708 0.3388  1.0000     -0.1714  0.4981
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1111 0.1666 0.5048   1.000     -0.2154  0.4376
#> 2   Item2 random_dif -0.1725 0.1601 0.2813   1.000     -0.4864  0.1413
#> 3   Item3 random_dif -0.3249 0.1652 0.0492   0.492     -0.6487 -0.0011
#> 4   Item4 random_dif -0.1345 0.1740 0.4395   1.000     -0.4756  0.2065
#> 5   Item5 random_dif  0.0185 0.1729 0.9150   1.000     -0.3204  0.3573
#> 6   Item6 random_dif  0.0684 0.1676 0.6831   1.000     -0.2600  0.3969
#> 7   Item7 random_dif -0.0700 0.1771 0.6924   1.000     -0.4171  0.2770
#> 8   Item8 random_dif  0.0641 0.1748 0.7137   1.000     -0.2785  0.4068
#> 9   Item9 random_dif  0.2468 0.1675 0.1406   1.000     -0.0815  0.5751
#> 10 Item10 random_dif  0.1804 0.1617 0.2644   1.000     -0.1364  0.4973
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1869 0.1621 0.2490  1.0000     -0.1309  0.5047
#> 2   Item2 random_dif -0.0827 0.1563 0.5966  1.0000     -0.3891  0.2236
#> 3   Item3 random_dif  0.1130 0.1679 0.5009  1.0000     -0.2160  0.4420
#> 4   Item4 random_dif  0.2261 0.1519 0.1364  1.0000     -0.0715  0.5238
#> 5   Item5 random_dif  0.0197 0.1706 0.9080  1.0000     -0.3146  0.3540
#> 6   Item6 random_dif  0.0349 0.1654 0.8327  1.0000     -0.2892  0.3591
#> 7   Item7 random_dif -0.4409 0.1372 0.0013  0.0131   * -0.7098 -0.1721
#> 8   Item8 random_dif  0.0055 0.1656 0.9735  1.0000     -0.3191  0.3301
#> 9   Item9 random_dif -0.0350 0.1628 0.8298  1.0000     -0.3541  0.2841
#> 10 Item10 random_dif  0.0025 0.1600 0.9873  1.0000     -0.3111  0.3162
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1950 0.1763 0.2686       1     -0.1505 0.5405
#> 2   Item2 random_dif -0.0942 0.1604 0.5569       1     -0.4087 0.2202
#> 3   Item3 random_dif -0.1660 0.1644 0.3128       1     -0.4883 0.1563
#> 4   Item4 random_dif -0.0223 0.1638 0.8918       1     -0.3433 0.2987
#> 5   Item5 random_dif -0.1710 0.1605 0.2868       1     -0.4856 0.1436
#> 6   Item6 random_dif  0.2277 0.1670 0.1728       1     -0.0996 0.5550
#> 7   Item7 random_dif -0.0819 0.1590 0.6064       1     -0.3936 0.2297
#> 8   Item8 random_dif  0.0525 0.1721 0.7604       1     -0.2848 0.3897
#> 9   Item9 random_dif  0.0081 0.1657 0.9611       1     -0.3166 0.3328
#> 10 Item10 random_dif  0.1155 0.1713 0.4999       1     -0.2202 0.4513
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0167 0.1859 0.9286       1     -0.3477 0.3810
#> 2   Item2 random_dif  0.2637 0.1625 0.1045       1     -0.0547 0.5822
#> 3   Item3 random_dif -0.1022 0.1809 0.5719       1     -0.4567 0.2522
#> 4   Item4 random_dif -0.1145 0.1662 0.4909       1     -0.4403 0.2113
#> 5   Item5 random_dif  0.0560 0.1634 0.7319       1     -0.2642 0.3762
#> 6   Item6 random_dif  0.0302 0.1746 0.8626       1     -0.3119 0.3724
#> 7   Item7 random_dif -0.1364 0.1664 0.4124       1     -0.4624 0.1897
#> 8   Item8 random_dif  0.1550 0.1766 0.3801       1     -0.1911 0.5012
#> 9   Item9 random_dif -0.0316 0.1693 0.8520       1     -0.3634 0.3002
#> 10 Item10 random_dif -0.1088 0.1676 0.5162       1     -0.4373 0.2197
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0565 0.1784 0.7515  1.0000     -0.2932 0.4062
#> 2   Item2 random_dif -0.2318 0.1577 0.1417  1.0000     -0.5410 0.0774
#> 3   Item3 random_dif  0.2073 0.1755 0.2376  1.0000     -0.1367 0.5513
#> 4   Item4 random_dif  0.0704 0.1631 0.6659  1.0000     -0.2492 0.3901
#> 5   Item5 random_dif  0.0646 0.1662 0.6977  1.0000     -0.2612 0.3904
#> 6   Item6 random_dif  0.1451 0.1670 0.3850  1.0000     -0.1823 0.4725
#> 7   Item7 random_dif -0.1783 0.1632 0.2747  1.0000     -0.4982 0.1416
#> 8   Item8 random_dif  0.1129 0.1739 0.5161  1.0000     -0.2279 0.4537
#> 9   Item9 random_dif  0.0696 0.1642 0.6717  1.0000     -0.2522 0.3913
#> 10 Item10 random_dif -0.2679 0.1555 0.0848  0.8484     -0.5726 0.0368
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0943 0.1764 0.5929       1     -0.4399 0.2514
#> 2   Item2 random_dif  0.1765 0.1573 0.2618       1     -0.1318 0.4847
#> 3   Item3 random_dif  0.0874 0.1766 0.6207       1     -0.2587 0.4336
#> 4   Item4 random_dif -0.0373 0.1785 0.8347       1     -0.3872 0.3127
#> 5   Item5 random_dif -0.1741 0.1669 0.2966       1     -0.5012 0.1529
#> 6   Item6 random_dif -0.0115 0.1710 0.9466       1     -0.3467 0.3238
#> 7   Item7 random_dif  0.0617 0.1782 0.7294       1     -0.2877 0.4110
#> 8   Item8 random_dif  0.0593 0.1689 0.7254       1     -0.2717 0.3903
#> 9   Item9 random_dif -0.1690 0.1662 0.3091       1     -0.4947 0.1567
#> 10 Item10 random_dif  0.0874 0.1653 0.5970       1     -0.2366 0.4115
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1782 0.1694 0.2930       1     -0.1539 0.5102
#> 2   Item2 random_dif  0.0679 0.1708 0.6908       1     -0.2668 0.4027
#> 3   Item3 random_dif  0.0074 0.1780 0.9667       1     -0.3414 0.3562
#> 4   Item4 random_dif  0.1391 0.1780 0.4346       1     -0.2097 0.4878
#> 5   Item5 random_dif -0.0761 0.1688 0.6522       1     -0.4070 0.2548
#> 6   Item6 random_dif -0.0252 0.1761 0.8861       1     -0.3703 0.3199
#> 7   Item7 random_dif -0.0931 0.1818 0.6085       1     -0.4494 0.2632
#> 8   Item8 random_dif -0.2504 0.1764 0.1558       1     -0.5961 0.0953
#> 9   Item9 random_dif -0.0090 0.1719 0.9582       1     -0.3459 0.3278
#> 10 Item10 random_dif  0.0350 0.1671 0.8340       1     -0.2926 0.3626
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0867 0.1664 0.6025  1.0000     -0.4129  0.2395
#> 2   Item2 random_dif -0.4002 0.1403 0.0043  0.0432   * -0.6751 -0.1253
#> 3   Item3 random_dif -0.1125 0.1648 0.4946  1.0000     -0.4355  0.2104
#> 4   Item4 random_dif  0.1796 0.1673 0.2831  1.0000     -0.1484  0.5076
#> 5   Item5 random_dif  0.3043 0.1597 0.0566  0.5665     -0.0086  0.6173
#> 6   Item6 random_dif -0.2563 0.1566 0.1017  1.0000     -0.5633  0.0506
#> 7   Item7 random_dif  0.0343 0.1692 0.8393  1.0000     -0.2973  0.3659
#> 8   Item8 random_dif  0.0867 0.1605 0.5888  1.0000     -0.2278  0.4013
#> 9   Item9 random_dif  0.0923 0.1651 0.5760  1.0000     -0.2312  0.4158
#> 10 Item10 random_dif  0.2010 0.1585 0.2047  1.0000     -0.1096  0.5116
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0395 0.1784 0.8247       1     -0.3102 0.3892
#> 2   Item2 random_dif -0.1157 0.1714 0.4997       1     -0.4517 0.2203
#> 3   Item3 random_dif -0.1187 0.1674 0.4782       1     -0.4467 0.2093
#> 4   Item4 random_dif -0.0570 0.1755 0.7455       1     -0.4010 0.2870
#> 5   Item5 random_dif -0.0522 0.1703 0.7593       1     -0.3860 0.2816
#> 6   Item6 random_dif  0.0490 0.1870 0.7934       1     -0.3175 0.4154
#> 7   Item7 random_dif  0.0767 0.1730 0.6575       1     -0.2624 0.4158
#> 8   Item8 random_dif -0.0220 0.1765 0.9007       1     -0.3680 0.3239
#> 9   Item9 random_dif -0.0111 0.1741 0.9490       1     -0.3524 0.3301
#> 10 Item10 random_dif  0.2118 0.1624 0.1922       1     -0.1065 0.5301
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1528 0.1675 0.3616  1.0000     -0.4811 0.1755
#> 2   Item2 random_dif  0.1848 0.1620 0.2542  1.0000     -0.1328 0.5023
#> 3   Item3 random_dif  0.1878 0.1641 0.2524  1.0000     -0.1338 0.5095
#> 4   Item4 random_dif -0.2673 0.1716 0.1194  1.0000     -0.6037 0.0691
#> 5   Item5 random_dif  0.0128 0.1678 0.9392  1.0000     -0.3162 0.3418
#> 6   Item6 random_dif -0.0793 0.1720 0.6448  1.0000     -0.4163 0.2578
#> 7   Item7 random_dif  0.0015 0.1763 0.9934  1.0000     -0.3440 0.3469
#> 8   Item8 random_dif -0.2278 0.1774 0.1990  1.0000     -0.5755 0.1198
#> 9   Item9 random_dif  0.2747 0.1507 0.0684  0.6837     -0.0207 0.5702
#> 10 Item10 random_dif -0.0261 0.1671 0.8759  1.0000     -0.3536 0.3014
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1910 0.1814 0.2923  1.0000     -0.1645  0.5465
#> 2   Item2 random_dif  0.0493 0.1631 0.7625  1.0000     -0.2704  0.3690
#> 3   Item3 random_dif  0.2070 0.1633 0.2049  1.0000     -0.1130  0.5270
#> 4   Item4 random_dif  0.1634 0.1675 0.3293  1.0000     -0.1649  0.4918
#> 5   Item5 random_dif -0.0888 0.1681 0.5974  1.0000     -0.4182  0.2407
#> 6   Item6 random_dif  0.1834 0.1740 0.2918  1.0000     -0.1576  0.5245
#> 7   Item7 random_dif -0.1202 0.1716 0.4835  1.0000     -0.4565  0.2161
#> 8   Item8 random_dif -0.3296 0.1592 0.0384  0.3837     -0.6416 -0.0177
#> 9   Item9 random_dif  0.0000 0.1619 1.0000  1.0000     -0.3172  0.3172
#> 10 Item10 random_dif -0.2217 0.1588 0.1627  1.0000     -0.5328  0.0895
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0208 0.1752 0.9056  1.0000     -0.3225 0.3641
#> 2   Item2 random_dif  0.1210 0.1656 0.4650  1.0000     -0.2036 0.4456
#> 3   Item3 random_dif -0.2938 0.1585 0.0637  0.6374     -0.6044 0.0168
#> 4   Item4 random_dif  0.1782 0.1737 0.3051  1.0000     -0.1623 0.5187
#> 5   Item5 random_dif  0.1734 0.1675 0.3006  1.0000     -0.1549 0.5018
#> 6   Item6 random_dif -0.1959 0.1713 0.2526  1.0000     -0.5316 0.1397
#> 7   Item7 random_dif -0.0032 0.1769 0.9857  1.0000     -0.3499 0.3436
#> 8   Item8 random_dif  0.1653 0.1663 0.3204  1.0000     -0.1607 0.4912
#> 9   Item9 random_dif -0.0851 0.1665 0.6092  1.0000     -0.4114 0.2412
#> 10 Item10 random_dif -0.0568 0.1599 0.7225  1.0000     -0.3703 0.2567
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2432 0.1638 0.1378  1.0000     -0.0780 0.5643
#> 2   Item2 random_dif -0.1377 0.1616 0.3942  1.0000     -0.4545 0.1790
#> 3   Item3 random_dif  0.1194 0.1584 0.4510  1.0000     -0.1911 0.4299
#> 4   Item4 random_dif  0.0185 0.1781 0.9172  1.0000     -0.3306 0.3676
#> 5   Item5 random_dif  0.1466 0.1606 0.3614  1.0000     -0.1682 0.4614
#> 6   Item6 random_dif -0.1013 0.1640 0.5369  1.0000     -0.4228 0.2202
#> 7   Item7 random_dif  0.0904 0.1717 0.5984  1.0000     -0.2461 0.4269
#> 8   Item8 random_dif  0.0116 0.1636 0.9435  1.0000     -0.3090 0.3321
#> 9   Item9 random_dif -0.2739 0.1522 0.0719  0.7195     -0.5722 0.0244
#> 10 Item10 random_dif -0.0818 0.1573 0.6031  1.0000     -0.3902 0.2266
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2353 0.1592 0.1393       1     -0.5473 0.0767
#> 2   Item2 random_dif  0.1971 0.1622 0.2241       1     -0.1207 0.5150
#> 3   Item3 random_dif  0.0076 0.1701 0.9645       1     -0.3258 0.3409
#> 4   Item4 random_dif -0.1523 0.1674 0.3628       1     -0.4805 0.1758
#> 5   Item5 random_dif -0.1616 0.1587 0.3086       1     -0.4726 0.1494
#> 6   Item6 random_dif  0.1932 0.1660 0.2445       1     -0.1321 0.5185
#> 7   Item7 random_dif  0.0292 0.1675 0.8614       1     -0.2991 0.3575
#> 8   Item8 random_dif  0.1272 0.1621 0.4325       1     -0.1905 0.4449
#> 9   Item9 random_dif -0.1105 0.1613 0.4934       1     -0.4267 0.2057
#> 10 Item10 random_dif  0.0972 0.1650 0.5559       1     -0.2263 0.4206
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.3303 0.1642 0.0443  0.4429      0.0084 0.6522
#> 2   Item2 random_dif -0.0075 0.1673 0.9642  1.0000     -0.3354 0.3204
#> 3   Item3 random_dif  0.0569 0.1690 0.7363  1.0000     -0.2744 0.3882
#> 4   Item4 random_dif -0.1238 0.1731 0.4746  1.0000     -0.4630 0.2155
#> 5   Item5 random_dif  0.1638 0.1560 0.2938  1.0000     -0.1419 0.4695
#> 6   Item6 random_dif -0.2670 0.1483 0.0719  0.7186     -0.5577 0.0237
#> 7   Item7 random_dif -0.0744 0.1617 0.6453  1.0000     -0.3914 0.2425
#> 8   Item8 random_dif  0.0426 0.1675 0.7994  1.0000     -0.2857 0.3708
#> 9   Item9 random_dif -0.0041 0.1711 0.9808  1.0000     -0.3395 0.3313
#> 10 Item10 random_dif -0.0619 0.1615 0.7015  1.0000     -0.3785 0.2546
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2040 0.1687 0.2267       1     -0.5347 0.1267
#> 2   Item2 random_dif  0.0243 0.1682 0.8854       1     -0.3054 0.3539
#> 3   Item3 random_dif -0.0361 0.1720 0.8336       1     -0.3733 0.3010
#> 4   Item4 random_dif  0.2066 0.1636 0.2065       1     -0.1140 0.5272
#> 5   Item5 random_dif  0.0550 0.1698 0.7458       1     -0.2778 0.3879
#> 6   Item6 random_dif -0.1685 0.1856 0.3642       1     -0.5323 0.1954
#> 7   Item7 random_dif  0.0171 0.1700 0.9200       1     -0.3161 0.3503
#> 8   Item8 random_dif  0.2718 0.1678 0.1052       1     -0.0570 0.6007
#> 9   Item9 random_dif -0.1731 0.1594 0.2777       1     -0.4856 0.1394
#> 10 Item10 random_dif  0.0286 0.1668 0.8640       1     -0.2984 0.3555
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1610 0.1675 0.3364       1     -0.4893 0.1673
#> 2   Item2 random_dif  0.0397 0.1627 0.8073       1     -0.2792 0.3586
#> 3   Item3 random_dif -0.1012 0.1616 0.5310       1     -0.4179 0.2154
#> 4   Item4 random_dif -0.0878 0.1651 0.5950       1     -0.4114 0.2358
#> 5   Item5 random_dif  0.0993 0.1600 0.5351       1     -0.2144 0.4129
#> 6   Item6 random_dif  0.0994 0.1645 0.5458       1     -0.2230 0.4217
#> 7   Item7 random_dif  0.0279 0.1806 0.8772       1     -0.3261 0.3819
#> 8   Item8 random_dif  0.1517 0.1551 0.3282       1     -0.1523 0.4556
#> 9   Item9 random_dif -0.2358 0.1574 0.1342       1     -0.5443 0.0727
#> 10 Item10 random_dif  0.1788 0.1647 0.2777       1     -0.1440 0.5016
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1837 0.1754 0.2950       1     -0.5276 0.1601
#> 2   Item2 random_dif -0.1169 0.1626 0.4724       1     -0.4357 0.2019
#> 3   Item3 random_dif  0.0282 0.1700 0.8684       1     -0.3050 0.3614
#> 4   Item4 random_dif  0.0970 0.1690 0.5659       1     -0.2342 0.4283
#> 5   Item5 random_dif  0.1683 0.1634 0.3029       1     -0.1519 0.4885
#> 6   Item6 random_dif -0.0573 0.1679 0.7327       1     -0.3863 0.2717
#> 7   Item7 random_dif  0.1102 0.1747 0.5281       1     -0.2322 0.4526
#> 8   Item8 random_dif -0.2211 0.1610 0.1696       1     -0.5367 0.0944
#> 9   Item9 random_dif -0.0335 0.1678 0.8419       1     -0.3623 0.2954
#> 10 Item10 random_dif  0.1920 0.1607 0.2319       1     -0.1228 0.5069
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0563 0.1688 0.7387       1     -0.2744 0.3870
#> 2   Item2 random_dif -0.2030 0.1585 0.2003       1     -0.5136 0.1076
#> 3   Item3 random_dif  0.0681 0.1696 0.6882       1     -0.2644 0.4005
#> 4   Item4 random_dif -0.1517 0.1741 0.3835       1     -0.4929 0.1895
#> 5   Item5 random_dif -0.1550 0.1559 0.3201       1     -0.4606 0.1506
#> 6   Item6 random_dif  0.0373 0.1661 0.8222       1     -0.2883 0.3629
#> 7   Item7 random_dif -0.0109 0.1702 0.9488       1     -0.3446 0.3227
#> 8   Item8 random_dif  0.1180 0.1750 0.5002       1     -0.2251 0.4611
#> 9   Item9 random_dif  0.1881 0.1611 0.2428       1     -0.1275 0.5038
#> 10 Item10 random_dif  0.0689 0.1633 0.6731       1     -0.2512 0.3891
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1003 0.1821 0.5817  1.0000     -0.4572  0.2566
#> 2   Item2 random_dif  0.1297 0.1555 0.4045  1.0000     -0.1752  0.4345
#> 3   Item3 random_dif -0.0076 0.1662 0.9637  1.0000     -0.3334  0.3183
#> 4   Item4 random_dif  0.2271 0.1617 0.1601  1.0000     -0.0897  0.5440
#> 5   Item5 random_dif  0.1241 0.1629 0.4463  1.0000     -0.1952  0.4433
#> 6   Item6 random_dif -0.4207 0.1439 0.0035  0.0345   * -0.7027 -0.1387
#> 7   Item7 random_dif  0.0820 0.1664 0.6220  1.0000     -0.2440  0.4081
#> 8   Item8 random_dif -0.0833 0.1677 0.6193  1.0000     -0.4121  0.2454
#> 9   Item9 random_dif  0.0337 0.1684 0.8412  1.0000     -0.2963  0.3638
#> 10 Item10 random_dif -0.0373 0.1615 0.8174  1.0000     -0.3539  0.2793
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0873 0.1669 0.6010  1.0000     -0.2398 0.4144
#> 2   Item2 random_dif -0.1111 0.1642 0.4986  1.0000     -0.4329 0.2107
#> 3   Item3 random_dif -0.0027 0.1730 0.9876  1.0000     -0.3417 0.3363
#> 4   Item4 random_dif -0.0462 0.1709 0.7869  1.0000     -0.3811 0.2887
#> 5   Item5 random_dif  0.0083 0.1616 0.9591  1.0000     -0.3084 0.3250
#> 6   Item6 random_dif -0.1960 0.1698 0.2484  1.0000     -0.5288 0.1368
#> 7   Item7 random_dif  0.1083 0.1687 0.5210  1.0000     -0.2224 0.4390
#> 8   Item8 random_dif  0.1027 0.1747 0.5567  1.0000     -0.2398 0.4452
#> 9   Item9 random_dif  0.2959 0.1537 0.0542  0.5423     -0.0054 0.5973
#> 10 Item10 random_dif -0.2268 0.1558 0.1456  1.0000     -0.5322 0.0787
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0042 0.1700 0.9802  1.0000     -0.3374 0.3290
#> 2   Item2 random_dif  0.0999 0.1565 0.5233  1.0000     -0.2068 0.4066
#> 3   Item3 random_dif  0.0587 0.1697 0.7297  1.0000     -0.2740 0.3913
#> 4   Item4 random_dif  0.2759 0.1616 0.0878  0.8784     -0.0409 0.5927
#> 5   Item5 random_dif -0.2996 0.1594 0.0602  0.6018     -0.6120 0.0128
#> 6   Item6 random_dif -0.0129 0.1723 0.9401  1.0000     -0.3506 0.3247
#> 7   Item7 random_dif -0.1631 0.1646 0.3216  1.0000     -0.4856 0.1594
#> 8   Item8 random_dif  0.0458 0.1788 0.7978  1.0000     -0.3047 0.3963
#> 9   Item9 random_dif  0.2058 0.1578 0.1923  1.0000     -0.1036 0.5152
#> 10 Item10 random_dif -0.1768 0.1587 0.2655  1.0000     -0.4879 0.1343
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0158 0.1975 0.9361  1.0000     -0.3713 0.4029
#> 2   Item2 random_dif -0.1297 0.1711 0.4486  1.0000     -0.4651 0.2058
#> 3   Item3 random_dif -0.1039 0.1690 0.5386  1.0000     -0.4352 0.2273
#> 4   Item4 random_dif -0.1576 0.1643 0.3373  1.0000     -0.4796 0.1644
#> 5   Item5 random_dif  0.0511 0.1715 0.7659  1.0000     -0.2851 0.3872
#> 6   Item6 random_dif  0.2094 0.1731 0.2262  1.0000     -0.1298 0.5486
#> 7   Item7 random_dif -0.0460 0.1769 0.7947  1.0000     -0.3928 0.3007
#> 8   Item8 random_dif  0.0493 0.1670 0.7676  1.0000     -0.2779 0.3766
#> 9   Item9 random_dif -0.2267 0.1649 0.1693  1.0000     -0.5499 0.0965
#> 10 Item10 random_dif  0.3282 0.1519 0.0307  0.3073      0.0305 0.6259
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1830 0.1676 0.2748  1.0000     -0.5114 0.1454
#> 2   Item2 random_dif -0.0741 0.1651 0.6534  1.0000     -0.3977 0.2494
#> 3   Item3 random_dif  0.2093 0.1723 0.2245  1.0000     -0.1284 0.5471
#> 4   Item4 random_dif  0.0861 0.1646 0.6008  1.0000     -0.2364 0.4086
#> 5   Item5 random_dif -0.0853 0.1705 0.6170  1.0000     -0.4194 0.2489
#> 6   Item6 random_dif -0.1101 0.1746 0.5283  1.0000     -0.4524 0.2322
#> 7   Item7 random_dif -0.0495 0.1650 0.7643  1.0000     -0.3729 0.2739
#> 8   Item8 random_dif -0.1596 0.1616 0.3233  1.0000     -0.4764 0.1572
#> 9   Item9 random_dif  0.3358 0.1467 0.0221  0.2206      0.0483 0.6234
#> 10 Item10 random_dif  0.0113 0.1652 0.9456  1.0000     -0.3125 0.3351
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1747 0.1668 0.2948  1.0000     -0.1521  0.5015
#> 2   Item2 random_dif  0.0364 0.1662 0.8267  1.0000     -0.2894  0.3622
#> 3   Item3 random_dif  0.1520 0.1618 0.3477  1.0000     -0.1652  0.4692
#> 4   Item4 random_dif  0.2232 0.1600 0.1630  1.0000     -0.0904  0.5367
#> 5   Item5 random_dif -0.0474 0.1679 0.7778  1.0000     -0.3764  0.2817
#> 6   Item6 random_dif  0.3409 0.1537 0.0265  0.2655      0.0397  0.6421
#> 7   Item7 random_dif -0.1298 0.1759 0.4606  1.0000     -0.4746  0.2150
#> 8   Item8 random_dif -0.4524 0.1330 0.0007  0.0067  ** -0.7131 -0.1918
#> 9   Item9 random_dif  0.0143 0.1655 0.9313  1.0000     -0.3102  0.3387
#> 10 Item10 random_dif -0.2766 0.1498 0.0647  0.6472     -0.5702  0.0169
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1872 0.1630 0.2508       1     -0.5066 0.1323
#> 2   Item2 random_dif  0.0223 0.1654 0.8929       1     -0.3019 0.3465
#> 3   Item3 random_dif -0.1866 0.1591 0.2410       1     -0.4984 0.1253
#> 4   Item4 random_dif -0.0054 0.1671 0.9744       1     -0.3330 0.3222
#> 5   Item5 random_dif -0.2290 0.1575 0.1458       1     -0.5376 0.0796
#> 6   Item6 random_dif  0.2633 0.1701 0.1215       1     -0.0700 0.5966
#> 7   Item7 random_dif  0.0118 0.1660 0.9432       1     -0.3134 0.3371
#> 8   Item8 random_dif  0.1497 0.1611 0.3526       1     -0.1659 0.4654
#> 9   Item9 random_dif -0.0222 0.1738 0.8986       1     -0.3628 0.3185
#> 10 Item10 random_dif  0.2284 0.1686 0.1755       1     -0.1020 0.5588
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1406 0.1745 0.4201  1.0000     -0.4826  0.2013
#> 2   Item2 random_dif  0.2224 0.1611 0.1674  1.0000     -0.0933  0.5380
#> 3   Item3 random_dif -0.1169 0.1755 0.5053  1.0000     -0.4609  0.2270
#> 4   Item4 random_dif -0.3758 0.1586 0.0178  0.1781     -0.6867 -0.0650
#> 5   Item5 random_dif -0.0857 0.1663 0.6062  1.0000     -0.4116  0.2401
#> 6   Item6 random_dif  0.3645 0.1551 0.0188  0.1877      0.0605  0.6684
#> 7   Item7 random_dif -0.0634 0.1670 0.7044  1.0000     -0.3907  0.2640
#> 8   Item8 random_dif -0.0142 0.1767 0.9357  1.0000     -0.3605  0.3320
#> 9   Item9 random_dif  0.1324 0.1730 0.4439  1.0000     -0.2066  0.4714
#> 10 Item10 random_dif  0.0444 0.1617 0.7834  1.0000     -0.2725  0.3613
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.3175 0.1648 0.0540  0.5401     -0.0055 0.6405
#> 2   Item2 random_dif -0.1058 0.1659 0.5238  1.0000     -0.4310 0.2194
#> 3   Item3 random_dif -0.1313 0.1742 0.4510  1.0000     -0.4727 0.2101
#> 4   Item4 random_dif -0.0123 0.1743 0.9435  1.0000     -0.3540 0.3293
#> 5   Item5 random_dif -0.0600 0.1660 0.7178  1.0000     -0.3853 0.2653
#> 6   Item6 random_dif  0.0889 0.1680 0.5965  1.0000     -0.2403 0.4181
#> 7   Item7 random_dif  0.2020 0.1752 0.2488  1.0000     -0.1413 0.5454
#> 8   Item8 random_dif -0.2271 0.1567 0.1472  1.0000     -0.5342 0.0800
#> 9   Item9 random_dif -0.0691 0.1614 0.6685  1.0000     -0.3855 0.2472
#> 10 Item10 random_dif  0.0556 0.1648 0.7356  1.0000     -0.2673 0.3785
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0310 0.1772 0.8612       1     -0.3782 0.3163
#> 2   Item2 random_dif -0.0688 0.1734 0.6916       1     -0.4087 0.2711
#> 3   Item3 random_dif  0.1371 0.1697 0.4192       1     -0.1955 0.4696
#> 4   Item4 random_dif  0.0559 0.1807 0.7571       1     -0.2982 0.4100
#> 5   Item5 random_dif -0.1146 0.1634 0.4830       1     -0.4348 0.2056
#> 6   Item6 random_dif  0.1524 0.1662 0.3593       1     -0.1734 0.4782
#> 7   Item7 random_dif  0.0514 0.1704 0.7630       1     -0.2826 0.3853
#> 8   Item8 random_dif -0.0604 0.1831 0.7414       1     -0.4193 0.2984
#> 9   Item9 random_dif -0.1272 0.1600 0.4268       1     -0.4408 0.1865
#> 10 Item10 random_dif  0.0158 0.1796 0.9298       1     -0.3362 0.3679
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3812 0.1545 0.0136  0.1364     -0.6840 -0.0783
#> 2   Item2 random_dif  0.2151 0.1558 0.1673  1.0000     -0.0902  0.5205
#> 3   Item3 random_dif  0.0043 0.1784 0.9808  1.0000     -0.3454  0.3540
#> 4   Item4 random_dif -0.0746 0.1661 0.6532  1.0000     -0.4001  0.2509
#> 5   Item5 random_dif  0.1864 0.1747 0.2860  1.0000     -0.1560  0.5288
#> 6   Item6 random_dif  0.0172 0.1714 0.9200  1.0000     -0.3187  0.3532
#> 7   Item7 random_dif -0.2127 0.1647 0.1966  1.0000     -0.5354  0.1101
#> 8   Item8 random_dif  0.0190 0.1708 0.9113  1.0000     -0.3157  0.3537
#> 9   Item9 random_dif -0.0697 0.1670 0.6764  1.0000     -0.3971  0.2577
#> 10 Item10 random_dif  0.2641 0.1569 0.0925  0.9248     -0.0436  0.5717
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0598 0.2065 0.7720  1.0000     -0.3449  0.4646
#> 2   Item2 random_dif  0.1278 0.1660 0.4413  1.0000     -0.1975  0.4532
#> 3   Item3 random_dif  0.2299 0.1668 0.1681  1.0000     -0.0970  0.5568
#> 4   Item4 random_dif  0.0505 0.1708 0.7674  1.0000     -0.2842  0.3852
#> 5   Item5 random_dif -0.1327 0.1656 0.4230  1.0000     -0.4572  0.1918
#> 6   Item6 random_dif -0.0016 0.1791 0.9929  1.0000     -0.3527  0.3495
#> 7   Item7 random_dif -0.0456 0.1724 0.7914  1.0000     -0.3835  0.2923
#> 8   Item8 random_dif  0.2096 0.1694 0.2160  1.0000     -0.1224  0.5416
#> 9   Item9 random_dif -0.1657 0.1700 0.3299  1.0000     -0.4989  0.1676
#> 10 Item10 random_dif -0.3195 0.1595 0.0451  0.4513     -0.6321 -0.0069
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0630 0.1761 0.7207  1.0000     -0.4081 0.2822
#> 2   Item2 random_dif  0.0814 0.1671 0.6259  1.0000     -0.2460 0.4089
#> 3   Item3 random_dif -0.1022 0.1674 0.5418  1.0000     -0.4303 0.2260
#> 4   Item4 random_dif -0.1612 0.1676 0.3361  1.0000     -0.4897 0.1673
#> 5   Item5 random_dif -0.1063 0.1677 0.5261  1.0000     -0.4350 0.2224
#> 6   Item6 random_dif -0.0547 0.1671 0.7433  1.0000     -0.3823 0.2728
#> 7   Item7 random_dif  0.3638 0.1664 0.0288  0.2881      0.0376 0.6899
#> 8   Item8 random_dif  0.0634 0.1718 0.7122  1.0000     -0.2733 0.4000
#> 9   Item9 random_dif  0.1481 0.1714 0.3877  1.0000     -0.1879 0.4841
#> 10 Item10 random_dif -0.1205 0.1604 0.4526  1.0000     -0.4350 0.1939
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0518 0.1724 0.7637       1     -0.2860 0.3897
#> 2   Item2 random_dif  0.0785 0.1616 0.6272       1     -0.2382 0.3951
#> 3   Item3 random_dif -0.2409 0.1671 0.1496       1     -0.5685 0.0867
#> 4   Item4 random_dif -0.2319 0.1660 0.1624       1     -0.5573 0.0935
#> 5   Item5 random_dif  0.1285 0.1712 0.4529       1     -0.2071 0.4641
#> 6   Item6 random_dif  0.0725 0.1726 0.6745       1     -0.2658 0.4108
#> 7   Item7 random_dif  0.0611 0.1725 0.7230       1     -0.2769 0.3991
#> 8   Item8 random_dif  0.0966 0.1757 0.5825       1     -0.2478 0.4410
#> 9   Item9 random_dif -0.1680 0.1736 0.3332       1     -0.5082 0.1722
#> 10 Item10 random_dif  0.1197 0.1646 0.4670       1     -0.2029 0.4423
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0098 0.1672 0.9533  1.0000     -0.3376 0.3180
#> 2   Item2 random_dif  0.0384 0.1664 0.8177  1.0000     -0.2878 0.3645
#> 3   Item3 random_dif -0.1455 0.1618 0.3685  1.0000     -0.4626 0.1716
#> 4   Item4 random_dif -0.1429 0.1784 0.4232  1.0000     -0.4925 0.2068
#> 5   Item5 random_dif -0.0783 0.1725 0.6499  1.0000     -0.4163 0.2598
#> 6   Item6 random_dif  0.2899 0.1614 0.0725  0.7251     -0.0265 0.6064
#> 7   Item7 random_dif -0.0287 0.1696 0.8655  1.0000     -0.3612 0.3037
#> 8   Item8 random_dif  0.0014 0.1636 0.9932  1.0000     -0.3192 0.3220
#> 9   Item9 random_dif -0.0298 0.1659 0.8574  1.0000     -0.3549 0.2953
#> 10 Item10 random_dif  0.1111 0.1635 0.4967  1.0000     -0.2093 0.4315
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2727 0.1748 0.1187       1     -0.0698 0.6153
#> 2   Item2 random_dif -0.0557 0.1806 0.7577       1     -0.4098 0.2983
#> 3   Item3 random_dif  0.1408 0.1701 0.4078       1     -0.1926 0.4742
#> 4   Item4 random_dif  0.2149 0.1654 0.1939       1     -0.1093 0.5391
#> 5   Item5 random_dif -0.0293 0.1755 0.8675       1     -0.3732 0.3146
#> 6   Item6 random_dif -0.1649 0.1747 0.3451       1     -0.5073 0.1775
#> 7   Item7 random_dif -0.1416 0.1724 0.4113       1     -0.4794 0.1962
#> 8   Item8 random_dif -0.0529 0.1707 0.7566       1     -0.3874 0.2816
#> 9   Item9 random_dif  0.1117 0.1643 0.4967       1     -0.2104 0.4338
#> 10 Item10 random_dif -0.2585 0.1610 0.1084       1     -0.5742 0.0571
#>      Item        Var   gamma     se pvalue padj.BH  sig   lower   upper
#> 1   Item1 random_dif -0.1208 0.1855 0.5149  1.0000      -0.4845  0.2428
#> 2   Item2 random_dif -0.1317 0.1703 0.4395  1.0000      -0.4655  0.2021
#> 3   Item3 random_dif  0.0759 0.1723 0.6597  1.0000      -0.2619  0.4136
#> 4   Item4 random_dif -0.3556 0.1614 0.0276  0.2762      -0.6719 -0.0392
#> 5   Item5 random_dif -0.0917 0.1675 0.5844  1.0000      -0.4200  0.2367
#> 6   Item6 random_dif  0.0478 0.1784 0.7889  1.0000      -0.3020  0.3975
#> 7   Item7 random_dif  0.0236 0.1754 0.8930  1.0000      -0.3203  0.3675
#> 8   Item8 random_dif  0.0260 0.1770 0.8831  1.0000      -0.3209  0.3730
#> 9   Item9 random_dif  0.5414 0.1223 0.0000  0.0001  ***  0.3017  0.7811
#> 10 Item10 random_dif -0.1219 0.1732 0.4816  1.0000      -0.4613  0.2176
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1383 0.1737 0.4260  1.0000     -0.2022  0.4787
#> 2   Item2 random_dif -0.0652 0.1687 0.6990  1.0000     -0.3958  0.2654
#> 3   Item3 random_dif  0.1589 0.1803 0.3782  1.0000     -0.1945  0.5122
#> 4   Item4 random_dif -0.2318 0.1672 0.1656  1.0000     -0.5594  0.0959
#> 5   Item5 random_dif  0.1137 0.1637 0.4874  1.0000     -0.2071  0.4345
#> 6   Item6 random_dif -0.3438 0.1614 0.0332  0.3316     -0.6600 -0.0275
#> 7   Item7 random_dif -0.0317 0.1694 0.8514  1.0000     -0.3638  0.3003
#> 8   Item8 random_dif  0.3184 0.1615 0.0486  0.4860      0.0020  0.6349
#> 9   Item9 random_dif  0.0086 0.1736 0.9606  1.0000     -0.3316  0.3488
#> 10 Item10 random_dif -0.0429 0.1665 0.7966  1.0000     -0.3693  0.2835
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.5133 0.1442 0.0004  0.0037  **  0.2307 0.7959
#> 2   Item2 random_dif  0.0506 0.1638 0.7572  1.0000     -0.2703 0.3716
#> 3   Item3 random_dif  0.0486 0.1770 0.7837  1.0000     -0.2984 0.3956
#> 4   Item4 random_dif  0.1995 0.1554 0.1991  1.0000     -0.1050 0.5041
#> 5   Item5 random_dif -0.1240 0.1709 0.4680  1.0000     -0.4589 0.2109
#> 6   Item6 random_dif -0.2722 0.1824 0.1356  1.0000     -0.6296 0.0853
#> 7   Item7 random_dif -0.0740 0.1669 0.6575  1.0000     -0.4010 0.2531
#> 8   Item8 random_dif -0.0804 0.1608 0.6172  1.0000     -0.3956 0.2348
#> 9   Item9 random_dif -0.1079 0.1666 0.5173  1.0000     -0.4345 0.2187
#> 10 Item10 random_dif -0.1377 0.1570 0.3803  1.0000     -0.4454 0.1700
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.2366 0.1672 0.1572  1.0000     -0.0912  0.5643
#> 2   Item2 random_dif -0.3201 0.1569 0.0414  0.4137     -0.6277 -0.0125
#> 3   Item3 random_dif -0.0546 0.1673 0.7440  1.0000     -0.3826  0.2733
#> 4   Item4 random_dif -0.0679 0.1704 0.6901  1.0000     -0.4019  0.2660
#> 5   Item5 random_dif  0.2831 0.1715 0.0988  0.9879     -0.0530  0.6192
#> 6   Item6 random_dif -0.0875 0.1657 0.5977  1.0000     -0.4123  0.2374
#> 7   Item7 random_dif  0.0165 0.1684 0.9219  1.0000     -0.3136  0.3467
#> 8   Item8 random_dif  0.0759 0.1698 0.6550  1.0000     -0.2569  0.4086
#> 9   Item9 random_dif -0.0264 0.1710 0.8774  1.0000     -0.3616  0.3088
#> 10 Item10 random_dif -0.0050 0.1690 0.9763  1.0000     -0.3363  0.3263
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1230 0.1727 0.4763  1.0000     -0.4614  0.2155
#> 2   Item2 random_dif  0.0475 0.1697 0.7796  1.0000     -0.2852  0.3802
#> 3   Item3 random_dif  0.2373 0.1702 0.1631  1.0000     -0.0962  0.5708
#> 4   Item4 random_dif -0.1526 0.1721 0.3752  1.0000     -0.4899  0.1847
#> 5   Item5 random_dif  0.2341 0.1695 0.1673  1.0000     -0.0981  0.5663
#> 6   Item6 random_dif  0.2978 0.1680 0.0763  0.7626     -0.0314  0.6271
#> 7   Item7 random_dif -0.4264 0.1434 0.0029  0.0293   * -0.7074 -0.1454
#> 8   Item8 random_dif -0.1963 0.1684 0.2439  1.0000     -0.5264  0.1338
#> 9   Item9 random_dif  0.0454 0.1671 0.7859  1.0000     -0.2821  0.3729
#> 10 Item10 random_dif  0.0695 0.1650 0.6736  1.0000     -0.2539  0.3929
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2118 0.1667 0.2040   1.000     -0.5385 0.1150
#> 2   Item2 random_dif -0.2490 0.1542 0.1065   1.000     -0.5513 0.0533
#> 3   Item3 random_dif -0.0712 0.1711 0.6772   1.000     -0.4065 0.2641
#> 4   Item4 random_dif  0.1108 0.1664 0.5055   1.000     -0.2153 0.4369
#> 5   Item5 random_dif -0.0698 0.1700 0.6814   1.000     -0.4029 0.2633
#> 6   Item6 random_dif -0.0715 0.1682 0.6706   1.000     -0.4011 0.2581
#> 7   Item7 random_dif -0.0925 0.1696 0.5853   1.000     -0.4249 0.2399
#> 8   Item8 random_dif  0.1335 0.1723 0.4383   1.000     -0.2041 0.4712
#> 9   Item9 random_dif  0.1392 0.1683 0.4081   1.000     -0.1906 0.4691
#> 10 Item10 random_dif  0.3596 0.1490 0.0158   0.158      0.0676 0.6516
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1570 0.1770 0.3753  1.0000     -0.1900  0.5039
#> 2   Item2 random_dif -0.3760 0.1503 0.0124  0.1235     -0.6706 -0.0814
#> 3   Item3 random_dif  0.0854 0.1700 0.6153  1.0000     -0.2477  0.4185
#> 4   Item4 random_dif -0.1120 0.1764 0.5253  1.0000     -0.4577  0.2336
#> 5   Item5 random_dif -0.3029 0.1608 0.0596  0.5964     -0.6181  0.0123
#> 6   Item6 random_dif  0.0668 0.1683 0.6913  1.0000     -0.2631  0.3967
#> 7   Item7 random_dif  0.2282 0.1724 0.1856  1.0000     -0.1097  0.5662
#> 8   Item8 random_dif  0.0533 0.1831 0.7711  1.0000     -0.3056  0.4121
#> 9   Item9 random_dif  0.1501 0.1692 0.3750  1.0000     -0.1816  0.4818
#> 10 Item10 random_dif  0.0872 0.1680 0.6039  1.0000     -0.2421  0.4165
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1373 0.1744 0.4311  1.0000     -0.4792 0.2045
#> 2   Item2 random_dif -0.0896 0.1734 0.6054  1.0000     -0.4295 0.2503
#> 3   Item3 random_dif  0.1243 0.1683 0.4601  1.0000     -0.2056 0.4543
#> 4   Item4 random_dif -0.1546 0.1601 0.3343  1.0000     -0.4685 0.1592
#> 5   Item5 random_dif  0.0230 0.1611 0.8863  1.0000     -0.2927 0.3387
#> 6   Item6 random_dif -0.0120 0.1639 0.9416  1.0000     -0.3333 0.3093
#> 7   Item7 random_dif  0.1075 0.1709 0.5295  1.0000     -0.2275 0.4425
#> 8   Item8 random_dif -0.0015 0.1763 0.9933  1.0000     -0.3470 0.3440
#> 9   Item9 random_dif  0.2634 0.1556 0.0905  0.9048     -0.0416 0.5684
#> 10 Item10 random_dif -0.1363 0.1588 0.3906  1.0000     -0.4475 0.1749
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.5057 0.1341 0.0002  0.0016  **  0.2429 0.7686
#> 2   Item2 random_dif -0.0267 0.1634 0.8701  1.0000     -0.3469 0.2935
#> 3   Item3 random_dif  0.0149 0.1706 0.9305  1.0000     -0.3195 0.3493
#> 4   Item4 random_dif  0.1337 0.1681 0.4262  1.0000     -0.1957 0.4632
#> 5   Item5 random_dif -0.2284 0.1567 0.1451  1.0000     -0.5355 0.0788
#> 6   Item6 random_dif -0.0299 0.1772 0.8659  1.0000     -0.3773 0.3175
#> 7   Item7 random_dif -0.1833 0.1634 0.2621  1.0000     -0.5036 0.1370
#> 8   Item8 random_dif -0.1718 0.1692 0.3100  1.0000     -0.5034 0.1599
#> 9   Item9 random_dif -0.0247 0.1645 0.8808  1.0000     -0.3472 0.2978
#> 10 Item10 random_dif  0.0066 0.1590 0.9668  1.0000     -0.3050 0.3183
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1880 0.1739 0.2796       1     -0.1528 0.5288
#> 2   Item2 random_dif  0.0433 0.1709 0.7998       1     -0.2917 0.3783
#> 3   Item3 random_dif -0.2423 0.1617 0.1339       1     -0.5592 0.0745
#> 4   Item4 random_dif  0.1652 0.1846 0.3709       1     -0.1966 0.5269
#> 5   Item5 random_dif -0.0213 0.1647 0.8969       1     -0.3441 0.3014
#> 6   Item6 random_dif  0.0969 0.1636 0.5537       1     -0.2237 0.4175
#> 7   Item7 random_dif  0.0491 0.1713 0.7745       1     -0.2867 0.3848
#> 8   Item8 random_dif -0.0292 0.1697 0.8635       1     -0.3619 0.3035
#> 9   Item9 random_dif -0.1587 0.1699 0.3503       1     -0.4917 0.1743
#> 10 Item10 random_dif -0.0392 0.1673 0.8148       1     -0.3670 0.2886
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.3509 0.1741 0.0439  0.4387      0.0096 0.6921
#> 2   Item2 random_dif -0.1145 0.1656 0.4892  1.0000     -0.4391 0.2100
#> 3   Item3 random_dif  0.1908 0.1729 0.2698  1.0000     -0.1481 0.5296
#> 4   Item4 random_dif -0.0867 0.1631 0.5950  1.0000     -0.4062 0.2329
#> 5   Item5 random_dif -0.0680 0.1626 0.6756  1.0000     -0.3866 0.2506
#> 6   Item6 random_dif  0.1679 0.1732 0.3324  1.0000     -0.1716 0.5074
#> 7   Item7 random_dif -0.1203 0.1670 0.4711  1.0000     -0.4476 0.2069
#> 8   Item8 random_dif -0.1085 0.1692 0.5212  1.0000     -0.4402 0.2231
#> 9   Item9 random_dif  0.1149 0.1655 0.4873  1.0000     -0.2094 0.4392
#> 10 Item10 random_dif -0.2000 0.1572 0.2034  1.0000     -0.5082 0.1082
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.0269 0.1697 0.8739  1.0000     -0.3057  0.3596
#> 2   Item2 random_dif -0.1658 0.1595 0.2987  1.0000     -0.4784  0.1469
#> 3   Item3 random_dif  0.0859 0.1715 0.6165  1.0000     -0.2502  0.4220
#> 4   Item4 random_dif  0.3274 0.1562 0.0360  0.3602      0.0214  0.6335
#> 5   Item5 random_dif  0.1801 0.1565 0.2498  1.0000     -0.1266  0.4868
#> 6   Item6 random_dif  0.0730 0.1781 0.6819  1.0000     -0.2760  0.4221
#> 7   Item7 random_dif  0.1234 0.1635 0.4506  1.0000     -0.1971  0.4439
#> 8   Item8 random_dif  0.0041 0.1662 0.9801  1.0000     -0.3216  0.3299
#> 9   Item9 random_dif -0.1333 0.1665 0.4232  1.0000     -0.4596  0.1930
#> 10 Item10 random_dif -0.4553 0.1324 0.0006  0.0059  ** -0.7148 -0.1958
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1638 0.1896 0.3878  1.0000     -0.5354 0.2079
#> 2   Item2 random_dif  0.0496 0.1675 0.7672  1.0000     -0.2787 0.3778
#> 3   Item3 random_dif -0.0452 0.1667 0.7862  1.0000     -0.3719 0.2814
#> 4   Item4 random_dif  0.2463 0.1529 0.1071  1.0000     -0.0533 0.5459
#> 5   Item5 random_dif  0.3006 0.1557 0.0535  0.5351     -0.0045 0.6058
#> 6   Item6 random_dif -0.0792 0.1756 0.6518  1.0000     -0.4233 0.2649
#> 7   Item7 random_dif -0.0194 0.1781 0.9134  1.0000     -0.3684 0.3296
#> 8   Item8 random_dif -0.2104 0.1575 0.1817  1.0000     -0.5192 0.0984
#> 9   Item9 random_dif -0.0540 0.1734 0.7556  1.0000     -0.3939 0.2859
#> 10 Item10 random_dif -0.1042 0.1636 0.5244  1.0000     -0.4249 0.2165
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1052 0.1716 0.5399       1     -0.4415 0.2312
#> 2   Item2 random_dif  0.0419 0.1612 0.7949       1     -0.2740 0.3578
#> 3   Item3 random_dif  0.0718 0.1739 0.6796       1     -0.2691 0.4127
#> 4   Item4 random_dif -0.1572 0.1769 0.3742       1     -0.5039 0.1895
#> 5   Item5 random_dif -0.0463 0.1592 0.7713       1     -0.3583 0.2657
#> 6   Item6 random_dif  0.1323 0.1812 0.4652       1     -0.2228 0.4874
#> 7   Item7 random_dif  0.1653 0.1636 0.3123       1     -0.1553 0.4859
#> 8   Item8 random_dif  0.0652 0.1773 0.7131       1     -0.2824 0.4128
#> 9   Item9 random_dif -0.2060 0.1639 0.2087       1     -0.5273 0.1152
#> 10 Item10 random_dif  0.0326 0.1643 0.8427       1     -0.2895 0.3547
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0621 0.1749 0.7226       1     -0.4048 0.2807
#> 2   Item2 random_dif -0.2079 0.1631 0.2026       1     -0.5276 0.1119
#> 3   Item3 random_dif  0.1218 0.1656 0.4623       1     -0.2029 0.4464
#> 4   Item4 random_dif -0.2358 0.1552 0.1287       1     -0.5401 0.0684
#> 5   Item5 random_dif  0.0087 0.1661 0.9581       1     -0.3168 0.3342
#> 6   Item6 random_dif -0.0199 0.1737 0.9089       1     -0.3602 0.3205
#> 7   Item7 random_dif  0.1375 0.1718 0.4234       1     -0.1992 0.4743
#> 8   Item8 random_dif  0.1613 0.1670 0.3341       1     -0.1660 0.4886
#> 9   Item9 random_dif  0.1838 0.1588 0.2472       1     -0.1275 0.4950
#> 10 Item10 random_dif -0.0810 0.1641 0.6215       1     -0.4026 0.2406
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.1545 0.1678 0.3574  1.0000     -0.1745  0.4835
#> 2   Item2 random_dif  0.0179 0.1833 0.9222  1.0000     -0.3414  0.3772
#> 3   Item3 random_dif -0.1162 0.1685 0.4906  1.0000     -0.4464  0.2141
#> 4   Item4 random_dif -0.3289 0.1569 0.0361  0.3608     -0.6364 -0.0213
#> 5   Item5 random_dif -0.0986 0.1679 0.5571  1.0000     -0.4277  0.2305
#> 6   Item6 random_dif -0.0441 0.1738 0.7997  1.0000     -0.3847  0.2965
#> 7   Item7 random_dif  0.0316 0.1665 0.8494  1.0000     -0.2946  0.3579
#> 8   Item8 random_dif  0.2030 0.1591 0.2020  1.0000     -0.1089  0.5149
#> 9   Item9 random_dif  0.0821 0.1693 0.6276  1.0000     -0.2498  0.4140
#> 10 Item10 random_dif  0.0907 0.1707 0.5953  1.0000     -0.2438  0.4252
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0523 0.1794 0.7706  1.0000     -0.4040 0.2993
#> 2   Item2 random_dif  0.2301 0.1720 0.1810  1.0000     -0.1070 0.5673
#> 3   Item3 random_dif -0.0977 0.1606 0.5427  1.0000     -0.4124 0.2170
#> 4   Item4 random_dif -0.2891 0.1587 0.0685  0.6847     -0.6000 0.0219
#> 5   Item5 random_dif  0.0180 0.1700 0.9158  1.0000     -0.3152 0.3512
#> 6   Item6 random_dif  0.0654 0.1692 0.6992  1.0000     -0.2663 0.3970
#> 7   Item7 random_dif  0.4129 0.1402 0.0032  0.0322   *  0.1382 0.6876
#> 8   Item8 random_dif -0.0623 0.1637 0.7037  1.0000     -0.3831 0.2586
#> 9   Item9 random_dif -0.0494 0.1616 0.7599  1.0000     -0.3661 0.2673
#> 10 Item10 random_dif -0.1598 0.1557 0.3048  1.0000     -0.4650 0.1454
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0538 0.1760 0.7597       1     -0.3988 0.2912
#> 2   Item2 random_dif  0.1193 0.1695 0.4817       1     -0.2130 0.4515
#> 3   Item3 random_dif  0.0678 0.1629 0.6773       1     -0.2515 0.3871
#> 4   Item4 random_dif -0.0396 0.1760 0.8220       1     -0.3846 0.3054
#> 5   Item5 random_dif -0.2373 0.1688 0.1597       1     -0.5680 0.0935
#> 6   Item6 random_dif  0.0513 0.1730 0.7669       1     -0.2878 0.3904
#> 7   Item7 random_dif  0.0900 0.1694 0.5951       1     -0.2421 0.4222
#> 8   Item8 random_dif -0.0316 0.1765 0.8578       1     -0.3775 0.3142
#> 9   Item9 random_dif -0.0687 0.1627 0.6729       1     -0.3876 0.2502
#> 10 Item10 random_dif  0.0758 0.1640 0.6441       1     -0.2457 0.3972
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1032 0.1794 0.5651       1     -0.2485 0.4549
#> 2   Item2 random_dif  0.0416 0.1732 0.8101       1     -0.2978 0.3810
#> 3   Item3 random_dif  0.0776 0.1725 0.6529       1     -0.2606 0.4157
#> 4   Item4 random_dif -0.0347 0.1749 0.8428       1     -0.3776 0.3082
#> 5   Item5 random_dif -0.2169 0.1559 0.1641       1     -0.5225 0.0886
#> 6   Item6 random_dif -0.0890 0.1787 0.6187       1     -0.4393 0.2614
#> 7   Item7 random_dif  0.2069 0.1593 0.1938       1     -0.1052 0.5191
#> 8   Item8 random_dif  0.0525 0.1716 0.7598       1     -0.2839 0.3889
#> 9   Item9 random_dif  0.0025 0.1611 0.9875       1     -0.3133 0.3183
#> 10 Item10 random_dif -0.1106 0.1571 0.4815       1     -0.4185 0.1973
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0710 0.1694 0.6750  1.0000     -0.4031 0.2610
#> 2   Item2 random_dif  0.0431 0.1625 0.7911  1.0000     -0.2755 0.3616
#> 3   Item3 random_dif  0.1008 0.1719 0.5575  1.0000     -0.2361 0.4377
#> 4   Item4 random_dif  0.3797 0.1525 0.0128  0.1275      0.0809 0.6785
#> 5   Item5 random_dif -0.0511 0.1629 0.7539  1.0000     -0.3703 0.2682
#> 6   Item6 random_dif  0.0952 0.1699 0.5755  1.0000     -0.2379 0.4283
#> 7   Item7 random_dif -0.1995 0.1614 0.2166  1.0000     -0.5159 0.1169
#> 8   Item8 random_dif  0.1215 0.1637 0.4580  1.0000     -0.1994 0.4424
#> 9   Item9 random_dif -0.0966 0.1631 0.5535  1.0000     -0.4162 0.2230
#> 10 Item10 random_dif -0.2995 0.1567 0.0559  0.5595     -0.6066 0.0076
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1545 0.1697 0.3624   1.000     -0.1780 0.4871
#> 2   Item2 random_dif -0.1377 0.1562 0.3781   1.000     -0.4438 0.1684
#> 3   Item3 random_dif -0.1500 0.1642 0.3610   1.000     -0.4718 0.1718
#> 4   Item4 random_dif  0.1254 0.1683 0.4564   1.000     -0.2045 0.4553
#> 5   Item5 random_dif  0.3492 0.1575 0.0266   0.266      0.0405 0.6579
#> 6   Item6 random_dif  0.0207 0.1805 0.9088   1.000     -0.3331 0.3745
#> 7   Item7 random_dif -0.2251 0.1571 0.1518   1.000     -0.5329 0.0827
#> 8   Item8 random_dif  0.0251 0.1679 0.8810   1.000     -0.3040 0.3543
#> 9   Item9 random_dif  0.0182 0.1640 0.9117   1.000     -0.3033 0.3397
#> 10 Item10 random_dif -0.1005 0.1615 0.5337   1.000     -0.4171 0.2160
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.1406 0.1744 0.4200  1.0000     -0.2012 0.4824
#> 2   Item2 random_dif -0.0242 0.1662 0.8845  1.0000     -0.3499 0.3016
#> 3   Item3 random_dif -0.1927 0.1680 0.2515  1.0000     -0.5219 0.1366
#> 4   Item4 random_dif -0.1908 0.1708 0.2638  1.0000     -0.5256 0.1439
#> 5   Item5 random_dif -0.2122 0.1655 0.1998  1.0000     -0.5367 0.1122
#> 6   Item6 random_dif  0.1960 0.1622 0.2270  1.0000     -0.1220 0.5140
#> 7   Item7 random_dif -0.1734 0.1798 0.3348  1.0000     -0.5257 0.1789
#> 8   Item8 random_dif  0.0163 0.1759 0.9260  1.0000     -0.3284 0.3610
#> 9   Item9 random_dif  0.3086 0.1724 0.0735  0.7347     -0.0293 0.6466
#> 10 Item10 random_dif  0.1410 0.1671 0.3988  1.0000     -0.1865 0.4684
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2433 0.1650 0.1403  1.0000     -0.5667 0.0801
#> 2   Item2 random_dif -0.1167 0.1644 0.4777  1.0000     -0.4390 0.2055
#> 3   Item3 random_dif -0.1450 0.1688 0.3903  1.0000     -0.4757 0.1858
#> 4   Item4 random_dif -0.0812 0.1681 0.6288  1.0000     -0.4107 0.2482
#> 5   Item5 random_dif  0.2579 0.1580 0.1025  1.0000     -0.0517 0.5675
#> 6   Item6 random_dif -0.1036 0.1664 0.5338  1.0000     -0.4297 0.2226
#> 7   Item7 random_dif -0.1330 0.1849 0.4721  1.0000     -0.4955 0.2295
#> 8   Item8 random_dif  0.3098 0.1656 0.0614  0.6137     -0.0148 0.6344
#> 9   Item9 random_dif  0.1474 0.1670 0.3775  1.0000     -0.1799 0.4746
#> 10 Item10 random_dif  0.0637 0.1653 0.6999  1.0000     -0.2603 0.3877
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.1758 0.1814 0.3325  1.0000     -0.5312  0.1797
#> 2   Item2 random_dif -0.0437 0.1719 0.7992  1.0000     -0.3806  0.2932
#> 3   Item3 random_dif -0.0437 0.1756 0.8033  1.0000     -0.3878  0.3004
#> 4   Item4 random_dif  0.1622 0.1747 0.3531  1.0000     -0.1801  0.5045
#> 5   Item5 random_dif -0.3353 0.1636 0.0405  0.4045     -0.6560 -0.0146
#> 6   Item6 random_dif -0.1098 0.1759 0.5325  1.0000     -0.4545  0.2349
#> 7   Item7 random_dif  0.1217 0.1663 0.4644  1.0000     -0.2043  0.4477
#> 8   Item8 random_dif  0.2338 0.1619 0.1487  1.0000     -0.0835  0.5512
#> 9   Item9 random_dif  0.0209 0.1617 0.8970  1.0000     -0.2960  0.3379
#> 10 Item10 random_dif  0.1318 0.1632 0.4193  1.0000     -0.1881  0.4517
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0326 0.1812 0.8574  1.0000     -0.3876 0.3225
#> 2   Item2 random_dif  0.2876 0.1558 0.0649  0.6493     -0.0178 0.5929
#> 3   Item3 random_dif -0.0757 0.1703 0.6568  1.0000     -0.4095 0.2582
#> 4   Item4 random_dif -0.1486 0.1651 0.3678  1.0000     -0.4721 0.1748
#> 5   Item5 random_dif  0.0887 0.1670 0.5953  1.0000     -0.2386 0.4160
#> 6   Item6 random_dif  0.1756 0.1693 0.2996  1.0000     -0.1562 0.5074
#> 7   Item7 random_dif -0.1159 0.1761 0.5104  1.0000     -0.4610 0.2292
#> 8   Item8 random_dif -0.2565 0.1612 0.1115  1.0000     -0.5724 0.0594
#> 9   Item9 random_dif -0.0174 0.1698 0.9186  1.0000     -0.3502 0.3155
#> 10 Item10 random_dif  0.0775 0.1642 0.6368  1.0000     -0.2443 0.3994
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1718 0.1740 0.3234       1     -0.5128 0.1692
#> 2   Item2 random_dif  0.2091 0.1576 0.1846       1     -0.0998 0.5180
#> 3   Item3 random_dif  0.0447 0.1679 0.7899       1     -0.2844 0.3738
#> 4   Item4 random_dif  0.0806 0.1689 0.6331       1     -0.2504 0.4117
#> 5   Item5 random_dif -0.0673 0.1662 0.6856       1     -0.3931 0.2585
#> 6   Item6 random_dif -0.0338 0.1709 0.8430       1     -0.3688 0.3011
#> 7   Item7 random_dif -0.1711 0.1728 0.3222       1     -0.5098 0.1676
#> 8   Item8 random_dif -0.0221 0.1741 0.8991       1     -0.3633 0.3192
#> 9   Item9 random_dif  0.1852 0.1603 0.2481       1     -0.1291 0.4994
#> 10 Item10 random_dif -0.0985 0.1620 0.5432       1     -0.4161 0.2191
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1150 0.1627 0.4794       1     -0.4338 0.2038
#> 2   Item2 random_dif -0.1303 0.1655 0.4310       1     -0.4547 0.1940
#> 3   Item3 random_dif -0.0937 0.1650 0.5703       1     -0.4171 0.2298
#> 4   Item4 random_dif -0.0954 0.1771 0.5899       1     -0.4425 0.2516
#> 5   Item5 random_dif  0.0026 0.1690 0.9879       1     -0.3286 0.3338
#> 6   Item6 random_dif  0.1632 0.1660 0.3255       1     -0.1622 0.4886
#> 7   Item7 random_dif  0.1111 0.1718 0.5178       1     -0.2256 0.4478
#> 8   Item8 random_dif  0.0013 0.1713 0.9937       1     -0.3343 0.3370
#> 9   Item9 random_dif  0.0902 0.1673 0.5899       1     -0.2378 0.4182
#> 10 Item10 random_dif  0.0606 0.1595 0.7040       1     -0.2520 0.3732
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2124 0.1629 0.1922       1     -0.5317 0.1068
#> 2   Item2 random_dif  0.1432 0.1692 0.3973       1     -0.1884 0.4748
#> 3   Item3 random_dif  0.0040 0.1761 0.9820       1     -0.3412 0.3492
#> 4   Item4 random_dif  0.1008 0.1673 0.5469       1     -0.2271 0.4286
#> 5   Item5 random_dif -0.2132 0.1716 0.2141       1     -0.5494 0.1231
#> 6   Item6 random_dif -0.0197 0.1741 0.9098       1     -0.3609 0.3215
#> 7   Item7 random_dif  0.1373 0.1657 0.4075       1     -0.1875 0.4620
#> 8   Item8 random_dif  0.1642 0.1777 0.3554       1     -0.1840 0.5125
#> 9   Item9 random_dif -0.1556 0.1797 0.3865       1     -0.5079 0.1966
#> 10 Item10 random_dif  0.0215 0.1669 0.8975       1     -0.3056 0.3486
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.2133 0.1702 0.2101       1     -0.5468 0.1203
#> 2   Item2 random_dif -0.0280 0.1671 0.8668       1     -0.3556 0.2996
#> 3   Item3 random_dif  0.1242 0.1683 0.4606       1     -0.2057 0.4540
#> 4   Item4 random_dif  0.1086 0.1712 0.5260       1     -0.2270 0.4441
#> 5   Item5 random_dif -0.0646 0.1704 0.7048       1     -0.3986 0.2695
#> 6   Item6 random_dif -0.0633 0.1742 0.7166       1     -0.4047 0.2782
#> 7   Item7 random_dif  0.0145 0.1637 0.9295       1     -0.3064 0.3354
#> 8   Item8 random_dif  0.0839 0.1658 0.6130       1     -0.2411 0.4089
#> 9   Item9 random_dif  0.0498 0.1621 0.7587       1     -0.2680 0.3676
#> 10 Item10 random_dif -0.0281 0.1614 0.8620       1     -0.3444 0.2883
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.3222 0.1556 0.0384  0.3842     -0.6271 -0.0172
#> 2   Item2 random_dif -0.0756 0.1702 0.6569  1.0000     -0.4093  0.2580
#> 3   Item3 random_dif  0.0387 0.1695 0.8192  1.0000     -0.2935  0.3709
#> 4   Item4 random_dif  0.0215 0.1704 0.8995  1.0000     -0.3125  0.3555
#> 5   Item5 random_dif -0.0639 0.1777 0.7189  1.0000     -0.4122  0.2843
#> 6   Item6 random_dif  0.2372 0.1672 0.1560  1.0000     -0.0905  0.5650
#> 7   Item7 random_dif  0.1862 0.1616 0.2490  1.0000     -0.1304  0.5029
#> 8   Item8 random_dif  0.0936 0.1700 0.5820  1.0000     -0.2396  0.4268
#> 9   Item9 random_dif -0.2977 0.1534 0.0523  0.5234     -0.5984  0.0030
#> 10 Item10 random_dif  0.2039 0.1546 0.1871  1.0000     -0.0991  0.5069
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0760 0.1746 0.6632  1.0000     -0.4182  0.2661
#> 2   Item2 random_dif -0.0256 0.1703 0.8806  1.0000     -0.3594  0.3083
#> 3   Item3 random_dif  0.0574 0.1766 0.7450  1.0000     -0.2886  0.4035
#> 4   Item4 random_dif  0.0110 0.1717 0.9487  1.0000     -0.3254  0.3475
#> 5   Item5 random_dif  0.2677 0.1606 0.0956  0.9564     -0.0472  0.5825
#> 6   Item6 random_dif  0.0349 0.1688 0.8365  1.0000     -0.2961  0.3658
#> 7   Item7 random_dif -0.3060 0.1555 0.0491  0.4913     -0.6108 -0.0012
#> 8   Item8 random_dif -0.0699 0.1772 0.6932  1.0000     -0.4171  0.2773
#> 9   Item9 random_dif  0.2267 0.1637 0.1662  1.0000     -0.0942  0.5477
#> 10 Item10 random_dif -0.1406 0.1715 0.4124  1.0000     -0.4767  0.1956
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0849 0.1677 0.6127  1.0000     -0.2437 0.4135
#> 2   Item2 random_dif  0.1040 0.1655 0.5298  1.0000     -0.2204 0.4284
#> 3   Item3 random_dif -0.0767 0.1562 0.6231  1.0000     -0.3828 0.2293
#> 4   Item4 random_dif -0.0756 0.1606 0.6380  1.0000     -0.3903 0.2392
#> 5   Item5 random_dif -0.2369 0.1568 0.1308  1.0000     -0.5443 0.0704
#> 6   Item6 random_dif  0.2765 0.1636 0.0909  0.9093     -0.0441 0.5971
#> 7   Item7 random_dif -0.1616 0.1667 0.3326  1.0000     -0.4884 0.1652
#> 8   Item8 random_dif -0.1642 0.1683 0.3294  1.0000     -0.4941 0.1657
#> 9   Item9 random_dif  0.0744 0.1606 0.6432  1.0000     -0.2404 0.3892
#> 10 Item10 random_dif  0.1804 0.1593 0.2575  1.0000     -0.1319 0.4927
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.0829 0.1766 0.6387  1.0000     -0.4292 0.2633
#> 2   Item2 random_dif -0.2729 0.1463 0.0621  0.6215     -0.5597 0.0139
#> 3   Item3 random_dif -0.0760 0.1623 0.6398  1.0000     -0.3941 0.2422
#> 4   Item4 random_dif -0.1301 0.1653 0.4312  1.0000     -0.4540 0.1938
#> 5   Item5 random_dif -0.0084 0.1674 0.9601  1.0000     -0.3365 0.3197
#> 6   Item6 random_dif  0.1123 0.1672 0.5016  1.0000     -0.2153 0.4400
#> 7   Item7 random_dif  0.1883 0.1587 0.2356  1.0000     -0.1229 0.4994
#> 8   Item8 random_dif -0.2069 0.1537 0.1783  1.0000     -0.5081 0.0943
#> 9   Item9 random_dif  0.2328 0.1520 0.1257  1.0000     -0.0652 0.5308
#> 10 Item10 random_dif  0.2548 0.1515 0.0926  0.9257     -0.0421 0.5516
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif -0.1878 0.1612 0.2442  1.0000     -0.5038 0.1282
#> 2   Item2 random_dif  0.0646 0.1638 0.6930  1.0000     -0.2564 0.3857
#> 3   Item3 random_dif  0.4078 0.1531 0.0078  0.0775   .  0.1076 0.7079
#> 4   Item4 random_dif -0.0755 0.1709 0.6585  1.0000     -0.4105 0.2594
#> 5   Item5 random_dif  0.0959 0.1666 0.5647  1.0000     -0.2305 0.4224
#> 6   Item6 random_dif -0.2084 0.1667 0.2113  1.0000     -0.5351 0.1184
#> 7   Item7 random_dif -0.0498 0.1818 0.7840  1.0000     -0.4061 0.3064
#> 8   Item8 random_dif -0.0294 0.1780 0.8688  1.0000     -0.3784 0.3195
#> 9   Item9 random_dif  0.0322 0.1705 0.8504  1.0000     -0.3021 0.3664
#> 10 Item10 random_dif -0.0217 0.1637 0.8946  1.0000     -0.3425 0.2992
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif  0.3645 0.1462 0.0126  0.1262      0.0781  0.6510
#> 2   Item2 random_dif  0.0852 0.1711 0.6184  1.0000     -0.2501  0.4205
#> 3   Item3 random_dif  0.0899 0.1689 0.5943  1.0000     -0.2410  0.4209
#> 4   Item4 random_dif -0.3350 0.1494 0.0250  0.2499     -0.6278 -0.0421
#> 5   Item5 random_dif  0.0027 0.1750 0.9877  1.0000     -0.3402  0.3457
#> 6   Item6 random_dif  0.1239 0.1791 0.4892  1.0000     -0.2272  0.4750
#> 7   Item7 random_dif  0.0000 0.1793 1.0000  1.0000     -0.3515  0.3515
#> 8   Item8 random_dif -0.3615 0.1520 0.0174  0.1741     -0.6595 -0.0635
#> 9   Item9 random_dif  0.1826 0.1689 0.2797  1.0000     -0.1484  0.5136
#> 10 Item10 random_dif -0.1179 0.1635 0.4708  1.0000     -0.4385  0.2026
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.2737 0.1842 0.1372       1     -0.0872 0.6347
#> 2   Item2 random_dif -0.1590 0.1680 0.3441       1     -0.4883 0.1704
#> 3   Item3 random_dif -0.1596 0.1683 0.3429       1     -0.4895 0.1702
#> 4   Item4 random_dif  0.0613 0.1805 0.7343       1     -0.2926 0.4151
#> 5   Item5 random_dif -0.0999 0.1737 0.5653       1     -0.4402 0.2405
#> 6   Item6 random_dif -0.2301 0.1700 0.1759       1     -0.5633 0.1031
#> 7   Item7 random_dif  0.1221 0.1727 0.4796       1     -0.2164 0.4606
#> 8   Item8 random_dif  0.0912 0.1800 0.6123       1     -0.2615 0.4439
#> 9   Item9 random_dif  0.2167 0.1651 0.1892       1     -0.1068 0.5402
#> 10 Item10 random_dif -0.0598 0.1706 0.7258       1     -0.3941 0.2745
#>      Item        Var   gamma     se pvalue padj.BH sig   lower  upper
#> 1   Item1 random_dif  0.0043 0.1737 0.9803  1.0000     -0.3361 0.3447
#> 2   Item2 random_dif  0.0421 0.1621 0.7949  1.0000     -0.2756 0.3599
#> 3   Item3 random_dif  0.3310 0.1800 0.0659  0.6594     -0.0218 0.6839
#> 4   Item4 random_dif -0.0912 0.1786 0.6097  1.0000     -0.4413 0.2589
#> 5   Item5 random_dif  0.3437 0.1629 0.0349  0.3489      0.0244 0.6629
#> 6   Item6 random_dif  0.0454 0.1747 0.7950  1.0000     -0.2970 0.3878
#> 7   Item7 random_dif -0.0026 0.1701 0.9876  1.0000     -0.3361 0.3308
#> 8   Item8 random_dif -0.1977 0.1674 0.2376  1.0000     -0.5259 0.1304
#> 9   Item9 random_dif -0.0743 0.1622 0.6468  1.0000     -0.3921 0.2435
#> 10 Item10 random_dif -0.2535 0.1489 0.0887  0.8874     -0.5454 0.0384
#>      Item        Var   gamma     se pvalue padj.BH sig   lower   upper
#> 1   Item1 random_dif -0.0346 0.1777 0.8456  1.0000     -0.3829  0.3137
#> 2   Item2 random_dif -0.1127 0.1641 0.4925  1.0000     -0.4343  0.2090
#> 3   Item3 random_dif -0.2051 0.1676 0.2211  1.0000     -0.5336  0.1234
#> 4   Item4 random_dif  0.0912 0.1702 0.5922  1.0000     -0.2425  0.4248
#> 5   Item5 random_dif -0.3625 0.1504 0.0160  0.1596     -0.6574 -0.0677
#> 6   Item6 random_dif  0.2389 0.1731 0.1675  1.0000     -0.1003  0.5780
#> 7   Item7 random_dif  0.0597 0.1709 0.7267  1.0000     -0.2752  0.3947
#> 8   Item8 random_dif  0.0649 0.1857 0.7265  1.0000     -0.2990  0.4288
#> 9   Item9 random_dif  0.2032 0.1678 0.2259  1.0000     -0.1256  0.5320
#> 10 Item10 random_dif  0.1138 0.1708 0.5053  1.0000     -0.2210  0.4486
RMdifGamma(sim_data, dif_group, cutoff = cutoff_res)
#> 
#> 
#> Table: Partial gamma DIF analysis (n = 200 complete cases). Cutoff values based on 100 simulation iterations (99% HDCI).
#> 
#> |Item   | Partial gamma|    SE| Lower CI| Upper CI| Adj. p-value (BH)|p-value sign. | Gamma low| Gamma high|Flagged |
#> |:------|-------------:|-----:|--------:|--------:|-----------------:|:-------------|---------:|----------:|:-------|
#> |Item1  |         0.029| 0.161|   -0.287|    0.345|                 1|              |    -0.352|      0.513|FALSE   |
#> |Item2  |        -0.168| 0.156|   -0.474|    0.138|                 1|              |    -0.388|      0.302|FALSE   |
#> |Item3  |         0.111| 0.168|   -0.219|    0.441|                 1|              |    -0.406|      0.408|FALSE   |
#> |Item4  |        -0.030| 0.155|   -0.333|    0.273|                 1|              |    -0.446|      0.380|FALSE   |
#> |Item5  |         0.076| 0.159|   -0.235|    0.387|                 1|              |    -0.349|      0.349|FALSE   |
#> |Item6  |         0.042| 0.156|   -0.264|    0.349|                 1|              |    -0.421|      0.353|FALSE   |
#> |Item7  |        -0.146| 0.160|   -0.460|    0.169|                 1|              |    -0.466|      0.413|FALSE   |
#> |Item8  |         0.019| 0.162|   -0.298|    0.336|                 1|              |    -0.407|      0.318|FALSE   |
#> |Item9  |        -0.141| 0.161|   -0.457|    0.175|                 1|              |    -0.298|      0.439|FALSE   |
#> |Item10 |         0.178| 0.156|   -0.127|    0.483|                 1|              |    -0.387|      0.360|FALSE   |
# }
```
