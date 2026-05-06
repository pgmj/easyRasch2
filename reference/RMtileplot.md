# Tile Plot of Item Response Distributions

Creates a tile (heat map) plot showing the distribution of responses
across all items and response categories. Each cell displays the count
(or percentage) of responses, with optional conditional highlighting for
cells with low counts. Optional faceting by a grouping variable is
provided for inspecting subgroup response distributions before DIF
analyses — particularly useful for spotting empty categories or
under-represented subgroups before fitting Rasch models per group.

## Usage

``` r
RMtileplot(
  data,
  group = NULL,
  cutoff = 10,
  highlight = TRUE,
  percent = FALSE,
  text_color = "orange",
  item_labels = NULL,
  category_labels = NULL,
  group_labels = NULL,
  facet_ncol = NULL,
  output = c("ggplot", "dataframe")
)
```

## Arguments

- data:

  A data.frame in wide format containing only the item response columns.
  Each column is one item, each row is one person. All columns must be
  numeric (integer-valued). Response categories may be coded starting
  from 0 or 1. Do not include person IDs, grouping variables, or other
  non-item columns — supply the grouping variable separately via
  `group`.

- group:

  Optional vector of length `nrow(data)` (factor, character, or numeric)
  defining a grouping variable. When provided, the plot is faceted by
  group and counts / percentages are computed within each group. Default
  `NULL` (no faceting). Persons with `NA` group are excluded.

- cutoff:

  Integer. Cells with counts below this value are highlighted (when
  `highlight = TRUE`). Default `10`.

- highlight:

  Logical. If `TRUE` (default), cell labels with counts below `cutoff`
  are displayed in red. This includes empty cells (`n = 0`), useful for
  identifying gaps in the response distribution.

- percent:

  Logical. If `TRUE`, cell labels show percentages instead of raw
  counts. Percentages are computed within item (and within group, when
  `group` is supplied). Default `FALSE`.

- text_color:

  Character. Colour for non-highlighted cell labels. Default `"orange"`.

- item_labels:

  Optional character vector of descriptive labels for the items
  (y-axis), same length as `ncol(data)`. Default `NULL` uses the column
  names.

- category_labels:

  Optional character vector of labels for the response categories
  (x-axis), same length as the number of categories spanning `min` to
  `max` observed value. Default `NULL` uses the numeric category values.

- group_labels:

  Optional character vector of length `nlevels(as.factor(group))` to
  override the displayed facet labels. Order corresponds to
  `levels(as.factor(group))`.

- facet_ncol:

  Integer or `NULL`. Number of columns in the facet grid when `group` is
  supplied. Default `NULL` (ggplot2 chooses).

- output:

  Character. `"ggplot"` (default) returns the plot; `"dataframe"`
  returns the underlying per-cell counts (one row per item × category,
  plus group when supplied).

## Value

Either a `ggplot` object or a data.frame, depending on `output`.

## Details

Adapted from `easyRaschBayes::plot_tile()` and extended with the `group`
parameter for faceted display.

Items are placed on the y-axis (in the same order as the columns of
`data`, top to bottom) and response categories on the x-axis. Cell
shading represents the count of responses (darker = more responses).
Categories with zero responses are explicitly shown (`n = 0`), which
helps identify gaps in the response distribution — one of the primary
purposes of the plot, especially before DIF analyses where
under-represented categories within a subgroup can break model fitting
on that subgroup.

When `group` is supplied, percentages and the highlight cutoff are
applied within each group, so a cell labelled "5" in the group-A facet
contains the count for group A only.

## Examples

``` r
if (FALSE) { # \dontrun{
library(eRm)

# Basic tile plot
RMtileplot(pcmdat2)

# With percentages
RMtileplot(pcmdat2, percent = TRUE)

# Faceted by an external grouping variable
set.seed(1)
grp <- sample(c("A", "B"), nrow(pcmdat2), replace = TRUE)
RMtileplot(pcmdat2, group = grp)

# With custom labels and tighter cutoff
RMtileplot(pcmdat2,
           group = grp,
           group_labels = c("Female", "Male"),
           cutoff = 5,
           facet_ncol = 2)

# Underlying counts as a data.frame
RMtileplot(pcmdat2, group = grp, output = "dataframe")
} # }
```
