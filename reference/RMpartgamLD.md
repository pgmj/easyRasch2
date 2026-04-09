# Partial Gamma Local Dependence Analysis

Computes partial gamma coefficients for Local Dependence (LD) assessment
using
[`iarm::partgam_LD()`](https://rdrr.io/pkg/iarm/man/partgam_LD.html).
Each pair of items is tested for residual association, controlling for
the rest score (total score minus one of the items in the pair).

## Usage

``` r
RMpartgamLD(data, cutoff = NULL, output = "kable")
```

## Arguments

- data:

  A data.frame or matrix of item responses. Items must be scored
  starting at 0 (non-negative integers). Missing values (`NA`) are
  allowed, but at least one complete case must exist.

- cutoff:

  Optional. Default `NULL` (no cutoff applied). Can be:

  - The return value of
    [`RMpgLDcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgLDcutoff.md)
    (a list with `$pair_cutoffs`): the data.frame is extracted
    automatically and simulation metadata is included in the kable
    caption.

  - The `$pair_cutoffs` data.frame from
    [`RMpgLDcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgLDcutoff.md)
    directly: must have columns `Item1`, `Item2`, `gamma_low`,
    `gamma_high`. When provided, adds columns `Gamma_low`, `Gamma_high`,
    and `Flagged` (logical; `TRUE` when the observed partial gamma falls
    outside the credible range) to the result.

- output:

  Character string controlling the return value. Either `"kable"`
  (default) for a formatted
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) table, or
  `"dataframe"` for the underlying data.frame.

## Value

- If `output = "kable"`: a `knitr_kable` object with columns "Item 1",
  "Item 2", "Partial gamma", and "Adjusted p-value (BH)". When `cutoff`
  is provided, additional columns "Gamma low", "Gamma high", and
  "Flagged" are included. Two tables are returned (one per rest-score
  direction), combined into a single output.

- If `output = "dataframe"`: a list of two data.frames (one per
  rest-score direction) with columns `Item1`, `Item2`, `gamma`,
  `padj_bh`. When `cutoff` is provided, columns `gamma_low`,
  `gamma_high`, and `flagged` are also included.

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

## References

Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) *Rasch Models in
Health*. Iste and Wiley (2013), pp. 133–135.

## See also

[`RMpgLDcutoff`](https://pgmj.github.io/easyRasch2/reference/RMpgLDcutoff.md),
[`RMpgLDplot`](https://pgmj.github.io/easyRasch2/reference/RMpgLDplot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Default kable output
RMpartgamLD(sim_data)

# Return as data.frame list
RMpartgamLD(sim_data, output = "dataframe")

# With simulation-based cutoffs
cutoff_res <- RMpgLDcutoff(sim_data, iterations = 100, parallel = FALSE,
                           seed = 42)
RMpartgamLD(sim_data, cutoff = cutoff_res)
} # }
```
