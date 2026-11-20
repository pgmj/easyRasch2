# Item and threshold hierarchy plot

Plots the item and threshold locations of a partial credit model on the
logit scale. Items are sorted by location (the mean of their
thresholds). Each item's location is a black diamond and its thresholds
are coloured dots with confidence intervals.

## Usage

``` r
RMitemHierarchy(
  data,
  show_numbers = TRUE,
  sem_multiplier = 1.405,
  item_labels = NULL,
  output = c("ggplot", "dataframe")
)
```

## Arguments

- data:

  A data.frame or matrix of polytomous item responses, items scored from
  0, with at least one item above 1. One column per item, one row per
  person. Dichotomous data give an error.

- show_numbers:

  Logical. When `TRUE` (default), prints the numerical item location and
  threshold values next to the points on the plot. When `FALSE`, only
  the threshold labels (`T1`, `T2`, ...) are shown.

- sem_multiplier:

  Numeric multiplier for the threshold SE used to draw the error bars.
  Default `1.405` (84% CI). Common alternatives: `1.96` (95% CI),
  `2.576` (99% CI).

- item_labels:

  Optional character vector of length `ncol(data)` providing descriptive
  labels for items. Default `NULL` uses the column names of `data`.
  Labels are appended to the column names on the y-axis (`name - label`)
  and wrapped at 36 characters via base-R
  [`strwrap()`](https://rdrr.io/r/base/strwrap.html).

- output:

  `"ggplot"` (default) for the figure, or `"dataframe"` for its data,
  one row per item and threshold.

## Value

Either a `ggplot` (default) or a data.frame with columns `Item`,
`ItemLabel`, `Threshold`, `ThresholdLocation`, `ThresholdSE`, and
`ItemLocation` (the per-item mean of the centred thresholds).

## Details

The thresholds are centred on the mean of all thresholds, so the dashed
line at 0 marks the average threshold of the scale.

The intervals are 84% by default (`sem_multiplier = 1.405`), following
Payton, Greenstone and Schenker (2003). Two 84% intervals that do not
overlap correspond roughly to a two-sample test at \\\alpha = 0.05\\.
Use `sem_multiplier = 1.96` for 95% intervals.

**Polytomous only.** A dichotomous item has a single threshold, its
location, so the plot would show nothing more. For dichotomous data use
[`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md)
or
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md).

**Centring.** The CML thresholds
([`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html))
are shifted so that their mean is zero, and each item's location is the
mean of its centred thresholds. The dashed line marks this zero.

## References

Payton, M. E., Greenstone, M. H., & Schenker, N. (2003). Overlapping
confidence intervals or standard error intervals: What do they mean in
terms of statistical significance? *Journal of Insect Science, 3*(34),
1-6. [doi:10.1093/jis/3.1.34](https://doi.org/10.1093/jis/3.1.34)

## See also

[`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md)
for the thresholds together with the person locations,
[`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md)
for the estimates as a table,
[`RMitemCatProb()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatProb.md),
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md).

## Examples

``` r
# \donttest{
if (requireNamespace("eRm", quietly = TRUE)) {
  data("pcmdat2", package = "eRm")
  RMitemHierarchy(pcmdat2)

  # 95% CI instead of 84%
  RMitemHierarchy(pcmdat2, sem_multiplier = 1.96)

  # Underlying data.frame
  RMitemHierarchy(pcmdat2, output = "dataframe")
}
#>   Item ItemLabel Threshold ThresholdLocation ThresholdSE ItemLocation
#> 1   I1        I1        T1        -0.4581328   0.1545659    0.6840543
#> 2   I1        I1        T2         1.8262414   0.1815586    0.6840543
#> 3   I2        I2        T1         0.2175018   0.1496752    0.9295949
#> 4   I2        I2        T2         1.6416881   0.1907856    0.9295949
#> 5   I3        I3        T1        -2.6684648   0.3556982   -1.2127668
#> 6   I3        I3        T2         0.2429311   0.1419154   -1.2127668
#> 7   I4        I4        T1        -1.3370901   0.2126350   -0.4008824
#> 8   I4        I4        T2         0.5353253   0.1476954   -0.4008824
# }
```
