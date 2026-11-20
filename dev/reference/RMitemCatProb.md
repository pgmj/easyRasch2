# Item category probability curves

Plots the model's probability of each response category against the
latent variable \\\theta\\, one panel per item. The item parameters are
estimated by conditional maximum likelihood with
[`psychotools::pcmodel()`](https://rdrr.io/pkg/psychotools/man/pcmodel.html)
(a dichotomous item is a 2-category partial credit model). The curves
are coloured from low to high categories with a viridis palette.
Comparable to
[`eRm::plotICC()`](https://rdrr.io/pkg/eRm/man/plotICC.html) and the
trace plots of mirt, with optional labels for items and categories.

## Usage

``` r
RMitemCatProb(
  data,
  item_labels = NULL,
  category_labels = NULL,
  theta_range = NULL,
  n_points = 200L,
  viridis_option = "D",
  viridis_end = 0.95,
  facet_ncol = NULL,
  label_wrap = 25L,
  line_width = 0.9,
  font = "sans",
  output = c("ggplot", "dataframe"),
  label_curves = c("legend", "path"),
  item = NULL,
  text_size = 4
)
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0. Missing values (`NA`) are allowed, since CML uses the
  items each respondent answered.

- item_labels:

  Optional character vector of descriptive item labels (facet strip
  titles). Must be the same length as `ncol(data)`. If `NULL` (the
  default), column names are used.

- category_labels:

  Optional character vector of labels for the response categories
  (legend). Must be the same length as the number of categories spanning
  from 0 to the maximum observed value. If `NULL` (the default), numeric
  category values are used.

- theta_range:

  Numeric length 2. Range of the latent trait \\\theta\\ (logits)
  plotted along the x-axis. When `NULL` (default) the range is
  `c(-4, 4)` for `label_curves = "legend"` and `c(-5, 5)` for
  `label_curves = "path"`, which gives geomtextpath more room to place
  the labels at the peaks of the curves. Pass `c(lo, hi)` to override.

- n_points:

  Integer. Number of evenly-spaced \\\theta\\ values at which the curves
  are evaluated. Default `200L`.

- viridis_option:

  Character. Viridis palette identifier. One of `"A"` through `"H"`.
  Default `"D"` (viridis green).

- viridis_end:

  Numeric in (0, 1\]. Upper end of the viridis palette range. Lower
  values keep the palette in its mid-tones and avoid the very bright
  yellow at `1.0`. Default `0.95`.

- facet_ncol:

  Optional integer. Number of columns in the facet layout. Default
  `NULL`
  ([`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  auto-layout).

- label_wrap:

  Integer. Characters per line for facet-strip label wrapping. Default
  `25L`.

- line_width:

  Numeric. Line width for the probability curves. Default `0.9`.

- font:

  Character. Font family for all text. Default `"sans"`.

- output:

  Character. Either `"ggplot"` (default) for a `ggplot` object, or
  `"dataframe"` for the underlying long-format probability table.

- label_curves:

  Character. How response categories are identified. `"legend"`
  (default) draws lines and a colour legend with one entry per category.
  `"path"` uses
  [`geomtextpath::geom_textpath()`](https://allancameron.github.io/geomtextpath/reference/geom_textpath.html)
  to write each category's label along its curve, as in classic trace
  plots, without a legend. It works for one item at a time, since narrow
  panels leave too little room for the labels. The model is still fitted
  on all items in `data`, as CML needs, and `item` selects the item to
  plot. Needs the optional geomtextpath package.

- item:

  Character or integer. Used only when `label_curves = "path"`. Either
  an item name (matching a column in `data`) or a column index,
  identifying which item's category curves to plot. Required for path
  mode.

- text_size:

  Numeric. Used only when `label_curves = "path"`. Size of the path
  labels in mm. Default `4`.

## Value

- If `output = "ggplot"`: a
  [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object, one facet per item.

- If `output = "dataframe"`: a long-format data.frame with columns
  `Item` (factor in column order), `Category` (integer), `Theta` and
  `Probability`, one row per item, category and grid value of theta.

## Details

For each polytomous item *i* with response categories \\0, 1, \ldots,
K_i\\ and threshold parameters \\\delta\_{i,1}, \ldots,
\delta\_{i,K_i}\\ (CML, `psychotools`), the PCM category probability is

\$\$P(X_i = k \mid \theta) = \frac{\exp(\sum\_{j=1}^{k} (\theta -
\delta\_{i,j}))}{\sum\_{k'=0}^{K_i} \exp(\sum\_{j=1}^{k'} (\theta -
\delta\_{i,j}))},\$\$

with the empty sum (when \\k = 0\\) taken as zero. A dichotomous item
has a single threshold, its difficulty, which gives the usual logistic
item characteristic curve.

The colours follow the category number, so low categories are at one end
of the palette and high categories at the other. With `category_labels`,
the legend shows those labels (for example "Never", "Sometimes",
"Often").

An item with fewer categories than the others (for example one
three-category item in a four-category scale) shows only its own
categories. The y-axis still spans \\\[0, 1\]\\.

## See also

[`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md)
for the categories' ranges together with the person locations,
[`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemHierarchy.md)
for the threshold locations,
[`RMitemCatFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatFit.md)
for a test of whether the categories work as scored,
[`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md)
for conditional ICCs.

## Examples

``` r
# \donttest{
if (requireNamespace("eRm", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  data(pcmdat2, package = "eRm")

  # Default plot
  RMitemCatProb(pcmdat2)

  # Custom item and category labels
  RMitemCatProb(
    pcmdat2,
    item_labels     = c("Mood", "Sleep", "Appetite", "Energy"),
    category_labels = c("Never", "Sometimes", "Often")
  )

  # Underlying probability data
  df <- RMitemCatProb(pcmdat2, output = "dataframe")
  head(df)

  # One item with the labels written along the curves. The model is
  # still fitted on all four items, and `item` picks the one to plot.
  if (requireNamespace("geomtextpath", quietly = TRUE)) {
    RMitemCatProb(
      pcmdat2,
      category_labels = c("Never", "Sometimes", "Often"),
      label_curves    = "path",
      item            = "I1"
    )
  }
}

# }
```
