# Item-Threshold Hierarchy Plot

Visualises item and threshold locations on the logit scale for a Partial
Credit Model. Items are sorted by their (mean-threshold) location; each
item shows its location as a black diamond and its individual thresholds
as coloured dots with confidence-interval error bars.

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

  A data.frame or matrix of polytomous item responses (non-negative
  integers, 0-based, max value \> 1). One column per item, one row per
  person.

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

  One of `"ggplot"` (default) – a faceted hierarchy figure – or
  `"dataframe"` – the long-format underlying data with one row per (item
  × threshold).

## Value

Either a `ggplot` (default) or a data.frame with columns `Item`,
`ItemLabel`, `Threshold`, `ThresholdLocation`, `ThresholdSE`, and
`ItemLocation` (the per-item mean of the centred thresholds).

## Details

Threshold locations are centred at the grand mean of all thresholds
across items, so the dashed reference line at 0 represents the mean
threshold location across the scale.

Confidence intervals around thresholds are 84% by default
(`sem_multiplier = 1.405`), following Payton, Greenstone, & Schenker
(2003) – non-overlap of 84% intervals approximately corresponds to a
two-sample significance test at \\\alpha = 0.05\\. Use
`sem_multiplier = 1.96` for 95% intervals.

**Polytomous only.** Dichotomous items have a single threshold that
coincides with the item location; the hierarchy plot is visually
degenerate in that case. For dichotomous data use
[`RMtargeting()`](https://pgmj.github.io/easyRasch2/reference/RMtargeting.md)
or
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/reference/RMscoreSE.md)
instead.

**Centring convention.** The PCM thresholds returned by
[`eRm::thresholds()`](https://rdrr.io/pkg/eRm/man/thresholds.html) are
shifted so that their grand mean is zero; each item's location is then
the mean of its centred thresholds. The dashed horizontal reference line
on the plot marks this zero – i.e., the average threshold across all
items.

## References

Payton, M. E., Greenstone, M. H., & Schenker, N. (2003). Overlapping
confidence intervals or standard error intervals: What do they mean in
terms of statistical significance? *Journal of Insect Science, 3*(34),
1-6. [doi:10.1093/jis/3.1.34](https://doi.org/10.1093/jis/3.1.34)

## See also

[`RMtargeting()`](https://pgmj.github.io/easyRasch2/reference/RMtargeting.md),
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/reference/RMscoreSE.md),
[`RMciccPlot()`](https://pgmj.github.io/easyRasch2/reference/RMciccPlot.md)

## Examples

``` r
# \donttest{
data("pcmdat2", package = "eRm")
RMitemHierarchy(pcmdat2)


# 95% CI instead of 84%
RMitemHierarchy(pcmdat2, sem_multiplier = 1.96)


# Underlying data.frame
RMitemHierarchy(pcmdat2, output = "dataframe")
#>   Item ItemLabel Threshold ThresholdLocation ThresholdSE ItemLocation
#> 1   I1        I1        T1        -0.4581338   0.1519290    0.6840541
#> 2   I1        I1        T2         1.8262421   0.2013173    0.6840541
#> 3   I2        I2        T1         0.2174992   0.1484292    0.9295923
#> 4   I2        I2        T2         1.6416854   0.2110598    0.9295923
#> 5   I3        I3        T1        -2.6684535   0.3364599   -1.2127610
#> 6   I3        I3        T2         0.2429314   0.1582221   -1.2127610
#> 7   I4        I4        T1        -1.3370918   0.2026607   -0.4008854
#> 8   I4        I4        T2         0.5353210   0.1652776   -0.4008854
# }
```
