# Residuals of the joint subscore table for the Martin-Löf test

Standardised residuals for each cell of the joint table of subscores,
compared with its distribution under unidimensionality (Christensen,
Bjorner, Kreiner, & Petersen, 2002, eq. 13).
[`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md)
tests whether a partition departs from unidimensionality, and this table
shows where.

## Usage

``` r
RMdimMartinLofResiduals(
  data,
  partition,
  output = c("kable", "dataframe", "ggplot"),
  flag_threshold = NULL,
  color_by = c("residual", "n"),
  color_limits = NULL,
  min_expected = NULL
)
```

## Arguments

- data:

  A data.frame or matrix of item responses, items scored from 0. Only
  complete cases on the items in `partition` are used. Missing values on
  other items are ignored.

- partition:

  As in
  [`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md):
  a list of item-name or index vectors, or a vector of length
  `ncol(data)` with each item's subscale. Each subscale must have at
  least 2 items.

- output:

  Character. `"kable"` (default) for a two-way residual table when D = 2
  (a long table when D \> 2), `"dataframe"` for the long data.frame, or
  `"ggplot"` for a heatmap (D = 2, or D = 3 with a panel per value of
  `t3`, an error for D \> 3).

- flag_threshold:

  Numeric or `NULL`. If set, cells with `|residual| > flag_threshold`
  are marked, in bold in the kable and in a `Flagged` column of the
  data.frame. Default `NULL`, no marks: a fixed value is not a
  calibrated test, and with many cells some large residuals are expected
  under unidimensionality.
  [`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md)
  is the test.

- color_by:

  Character. What the tile colour shows with `output = "ggplot"`.
  `"residual"` (default), a diverging red-white-blue scale centred at 0.
  `"n"`, a sequential blue scale of the observed cell count, which shows
  whether large residuals come from sparse cells. In both, the residual
  is printed in each cell.

- color_limits:

  Numeric vector of length 2, or `NULL`. Limits of the colour scale with
  `output = "ggplot"`. Default `c(-5, 5)` with `color_by = "residual"`,
  since residuals in cells with tiny expected counts can be very large
  and would compress the rest of the scale. The cell labels show the
  values without the limits. `NULL` with `color_by = "n"`, which uses
  the range of the data.

- min_expected:

  Numeric or `NULL`. If set, cells with an expected count below
  `min_expected` get the residual `NA` (grey in the heatmap, never
  flagged). Default `NULL`. A value of 1 or 5 is the analogue of the
  Cochran rule for sparse cells in chi-square tests, and removes
  residuals whose normal approximation is unreliable.

## Value

- `output = "kable"`: a `knitr_kable` object. For D = 2, a wide table
  with `t1` in rows, `t2` in columns and the standardised residual in
  the cells (bold if beyond `flag_threshold`, a dash if `NA`). For D \>
  2, one row per cell.

- `output = "dataframe"`: a data.frame with one row per cell and the
  columns `t1`, ..., `tD`, `total`, `observed`, `expected` and
  `residual`, plus `Flagged` when `flag_threshold` is set.

- `output = "ggplot"`: a `geom_tile()` heatmap (D = 2 or 3 only).

## Details

For each cell of the joint subscore table (indexed by \\(t_1, \ldots,
t_D)\\), the conditional probability under H0 given the total score \\t
= \sum_d t_d\\ is \$\$p(t_1, \ldots, t_D \mid t) = \prod_d
\gamma^{(d)}\_{t_d} / \gamma_t,\$\$ the expected count is \\e = n_t
\cdot p\\, and the residual is \\(o - e) / \sqrt{n_t \cdot p \cdot (1 -
p)}\\. CML estimates from the unidimensional model are used for the
\\\gamma\\-functions.

Reading the table (D = 2): under unidimensionality, the residuals should
show no pattern and be roughly N(0, 1). Multidimensionality with
positively correlated dimensions typically shows up as **positive**
residuals at the corners of each antidiagonal (high `t1` + low `t2`, low
`t1` + high `t2`) and **negative** residuals near the antidiagonal
centre (matched subscores). Negatively correlated dimensions show
positive residuals at the table corners (high/low and low/high) and
negative residuals at high/high and low/low. See Christensen et al.
(2002, section 7) for a worked example.

Cells where the total score has no observed cases (`n_t = 0`) are
uninformative and are dropped from the output.

Complete cases only, as in
[`RMdimMartinLof`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md),
and rows with `NA` are dropped only after items outside `partition` have
been removed. The reported `n = X of Y respondents` counts respondents
complete on the partitioned items against the raw input rows.

Expected counts and residuals agree to numerical precision with the
`pml` SAS macro (Christensen, 2004), kindly shared by Karl Bang
Christensen. The one deliberate difference is at the extreme totals
`t = 0` and `t = max`, where the conditional probability is 1 and the
residual is undefined: the macro prints 0, this function returns `NA`.

## References

Christensen, K. B. (2004). *pml: A SAS macro for testing
unidimensionality in polytomous Rasch models* (Technical note). National
Institute of Occupational Health, Denmark, and Department of
Biostatistics, University of Copenhagen.

Christensen, K. B., Bjorner, J. B., Kreiner, S., & Petersen, J. H.
(2002). Testing unidimensionality in polytomous Rasch models.
*Psychometrika, 67*(4), 563-574.
[doi:10.1007/BF02295131](https://doi.org/10.1007/BF02295131)

## See also

[`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md)

## Examples

``` r
# \donttest{
set.seed(1)
dat <- as.data.frame(matrix(sample(0:1, 400 * 8, replace = TRUE),
                            nrow = 400, ncol = 8))
colnames(dat) <- paste0("I", 1:8)

# Wide kable table for D = 2
RMdimMartinLofResiduals(dat,
                     partition = list(c("I1","I2","I3","I4"),
                                      c("I5","I6","I7","I8")))
#> Warning: NaNs produced
#> 
#> 
#> Table: Standardised residuals (Christensen et al. 2002, eq. 13). Rows = subscale 1 score, columns = subscale 2 score. -- = uncomputable. n = 400 respondents.
#> 
#> |t1\t2 |0     |1     |2     |3     |4     |
#> |:-----|:-----|:-----|:-----|:-----|:-----|
#> |0     |      |0.95  |2.38  |-0.32 |-1.31 |
#> |1     |-0.95 |-2.17 |-0.12 |0.32  |0.22  |
#> |2     |0.22  |0.11  |-0.41 |0.10  |0.84  |
#> |3     |0.35  |0.22  |-0.46 |-0.66 |-1.67 |
#> |4     |1.22  |0.47  |-0.06 |1.67  |--    |

# Heatmap
if (requireNamespace("ggplot2", quietly = TRUE)) {
  RMdimMartinLofResiduals(dat,
                       partition = c(1,1,1,1,2,2,2,2),
                       output = "ggplot")
}
#> Warning: NaNs produced


# Underlying data.frame for custom analysis
df <- RMdimMartinLofResiduals(dat,
                           partition = c(1,1,1,1,2,2,2,2),
                           flag_threshold = 2,
                           output = "dataframe")
#> Warning: NaNs produced
df[df$Flagged, ]
#>   t1 t2 total observed expected  residual Flagged
#> 3  0  2     2       18 11.02774  2.378139    TRUE
#> 4  1  1     2       21 28.59561 -2.170933    TRUE
# }
```
