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
# \donttest{
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Default kable output
RMpartgamLD(sim_data)
#>  [1] "Table: Partial gamma LD analysis (n = 200 complete cases). Positive gamma indicates positive local dependence between items.\n\nDirection 1: rest score = total - Item2 \n\n Table: Partial gamma LD analysis (n = 200 complete cases). Positive gamma indicates positive local dependence between items.\n\nDirection 2: rest score = total - Item1"
#>  [2] " \n\n "                                                                                                                                                                                                                                                                                                                                              
#>  [3] "|Item 1 |Item 2 | Partial gamma| Adjusted p-value (BH)| \n\n |Item 1 |Item 2 | Partial gamma| Adjusted p-value (BH)|"                                                                                                                                                                                                                                
#>  [4] "|:------|:------|-------------:|---------------------:| \n\n |:------|:------|-------------:|---------------------:|"                                                                                                                                                                                                                                
#>  [5] "|Item1  |Item2  |         0.320|                 1.000| \n\n |Item2  |Item1  |         0.300|                 1.000|"                                                                                                                                                                                                                                
#>  [6] "|Item1  |Item3  |        -0.143|                 1.000| \n\n |Item3  |Item1  |        -0.158|                 1.000|"                                                                                                                                                                                                                                
#>  [7] "|Item1  |Item4  |        -0.069|                 1.000| \n\n |Item3  |Item2  |        -0.230|                 1.000|"                                                                                                                                                                                                                                
#>  [8] "|Item1  |Item5  |        -0.153|                 1.000| \n\n |Item4  |Item1  |        -0.141|                 1.000|"                                                                                                                                                                                                                                
#>  [9] "|Item1  |Item6  |        -0.031|                 1.000| \n\n |Item4  |Item2  |        -0.263|                 1.000|"                                                                                                                                                                                                                                
#> [10] "|Item1  |Item7  |        -0.163|                 1.000| \n\n |Item4  |Item3  |        -0.042|                 1.000|"                                                                                                                                                                                                                                
#> [11] "|Item1  |Item8  |         0.191|                 1.000| \n\n |Item5  |Item1  |        -0.152|                 1.000|"                                                                                                                                                                                                                                
#> [12] "|Item1  |Item9  |         0.124|                 1.000| \n\n |Item5  |Item2  |         0.039|                 1.000|"                                                                                                                                                                                                                                
#> [13] "|Item1  |Item10 |        -0.039|                 1.000| \n\n |Item5  |Item3  |         0.201|                 1.000|"                                                                                                                                                                                                                                
#> [14] "|Item2  |Item3  |        -0.206|                 1.000| \n\n |Item5  |Item4  |         0.021|                 1.000|"                                                                                                                                                                                                                                
#> [15] "|Item2  |Item4  |        -0.202|                 1.000| \n\n |Item6  |Item1  |        -0.102|                 1.000|"                                                                                                                                                                                                                                
#> [16] "|Item2  |Item5  |         0.058|                 1.000| \n\n |Item6  |Item2  |        -0.063|                 1.000|"                                                                                                                                                                                                                                
#> [17] "|Item2  |Item6  |         0.025|                 1.000| \n\n |Item6  |Item3  |        -0.208|                 1.000|"                                                                                                                                                                                                                                
#> [18] "|Item2  |Item7  |        -0.230|                 1.000| \n\n |Item6  |Item4  |         0.206|                 1.000|"                                                                                                                                                                                                                                
#> [19] "|Item2  |Item8  |         0.150|                 1.000| \n\n |Item6  |Item5  |        -0.201|                 1.000|"                                                                                                                                                                                                                                
#> [20] "|Item2  |Item9  |         0.154|                 1.000| \n\n |Item7  |Item1  |        -0.116|                 1.000|"                                                                                                                                                                                                                                
#> [21] "|Item2  |Item10 |         0.024|                 1.000| \n\n |Item7  |Item2  |        -0.195|                 1.000|"                                                                                                                                                                                                                                
#> [22] "|Item3  |Item4  |         0.021|                 1.000| \n\n |Item7  |Item3  |         0.211|                 1.000|"                                                                                                                                                                                                                                
#> [23] "|Item3  |Item5  |         0.190|                 1.000| \n\n |Item7  |Item4  |         0.428|                 0.126|"                                                                                                                                                                                                                                
#> [24] "|Item3  |Item6  |        -0.143|                 1.000| \n\n |Item7  |Item5  |         0.218|                 1.000|"                                                                                                                                                                                                                                
#> [25] "|Item3  |Item7  |         0.086|                 1.000| \n\n |Item7  |Item6  |         0.093|                 1.000|"                                                                                                                                                                                                                                
#> [26] "|Item3  |Item8  |        -0.042|                 1.000| \n\n |Item8  |Item1  |         0.237|                 1.000|"                                                                                                                                                                                                                                
#> [27] "|Item3  |Item9  |        -0.074|                 1.000| \n\n |Item8  |Item2  |         0.220|                 1.000|"                                                                                                                                                                                                                                
#> [28] "|Item3  |Item10 |         0.221|                 1.000| \n\n |Item8  |Item3  |         0.024|                 1.000|"                                                                                                                                                                                                                                
#> [29] "|Item4  |Item5  |        -0.024|                 1.000| \n\n |Item8  |Item4  |        -0.319|                 1.000|"                                                                                                                                                                                                                                
#> [30] "|Item4  |Item6  |         0.224|                 1.000| \n\n |Item8  |Item5  |         0.152|                 1.000|"                                                                                                                                                                                                                                
#> [31] "|Item4  |Item7  |         0.333|                 1.000| \n\n |Item8  |Item6  |        -0.022|                 1.000|"                                                                                                                                                                                                                                
#> [32] "|Item4  |Item8  |        -0.406|                 0.305| \n\n |Item8  |Item7  |        -0.222|                 1.000|"                                                                                                                                                                                                                                
#> [33] "|Item4  |Item9  |        -0.153|                 1.000| \n\n |Item9  |Item1  |         0.193|                 1.000|"                                                                                                                                                                                                                                
#> [34] "|Item4  |Item10 |        -0.084|                 1.000| \n\n |Item9  |Item2  |         0.199|                 1.000|"                                                                                                                                                                                                                                
#> [35] "|Item5  |Item6  |        -0.136|                 1.000| \n\n |Item9  |Item3  |        -0.008|                 1.000|"                                                                                                                                                                                                                                
#> [36] "|Item5  |Item7  |         0.157|                 1.000| \n\n |Item9  |Item4  |        -0.071|                 1.000|"                                                                                                                                                                                                                                
#> [37] "|Item5  |Item8  |         0.130|                 1.000| \n\n |Item9  |Item5  |        -0.256|                 1.000|"                                                                                                                                                                                                                                
#> [38] "|Item5  |Item9  |        -0.348|                 1.000| \n\n |Item9  |Item6  |         0.218|                 1.000|"                                                                                                                                                                                                                                
#> [39] "|Item5  |Item10 |         0.014|                 1.000| \n\n |Item9  |Item7  |        -0.001|                 1.000|"                                                                                                                                                                                                                                
#> [40] "|Item6  |Item7  |         0.004|                 1.000| \n\n |Item9  |Item8  |         0.357|                 1.000|"                                                                                                                                                                                                                                
#> [41] "|Item6  |Item8  |        -0.133|                 1.000| \n\n |Item10 |Item1  |        -0.071|                 1.000|"                                                                                                                                                                                                                                
#> [42] "|Item6  |Item9  |         0.062|                 1.000| \n\n |Item10 |Item2  |        -0.023|                 1.000|"                                                                                                                                                                                                                                
#> [43] "|Item6  |Item10 |        -0.243|                 1.000| \n\n |Item10 |Item3  |         0.176|                 1.000|"                                                                                                                                                                                                                                
#> [44] "|Item7  |Item8  |        -0.198|                 1.000| \n\n |Item10 |Item4  |        -0.088|                 1.000|"                                                                                                                                                                                                                                
#> [45] "|Item7  |Item9  |        -0.005|                 1.000| \n\n |Item10 |Item5  |        -0.029|                 1.000|"                                                                                                                                                                                                                                
#> [46] "|Item7  |Item10 |         0.021|                 1.000| \n\n |Item10 |Item6  |        -0.235|                 1.000|"                                                                                                                                                                                                                                
#> [47] "|Item8  |Item9  |         0.301|                 1.000| \n\n |Item10 |Item7  |        -0.040|                 1.000|"                                                                                                                                                                                                                                
#> [48] "|Item8  |Item10 |         0.045|                 1.000| \n\n |Item10 |Item8  |        -0.047|                 1.000|"                                                                                                                                                                                                                                
#> [49] "|Item9  |Item10 |         0.036|                 1.000| \n\n |Item10 |Item9  |        -0.069|                 1.000|"                                                                                                                                                                                                                                
#> attr(,"class")
#> [1] "knit_asis"
#> attr(,"knit_cacheable")
#> [1] NA

# Return as data.frame list
RMpartgamLD(sim_data, output = "dataframe")
#> [[1]]
#>    Item1  Item2  gamma padj_bh
#> 1  Item1  Item2  0.320   1.000
#> 2  Item1  Item3 -0.143   1.000
#> 3  Item1  Item4 -0.069   1.000
#> 4  Item1  Item5 -0.153   1.000
#> 5  Item1  Item6 -0.031   1.000
#> 6  Item1  Item7 -0.163   1.000
#> 7  Item1  Item8  0.191   1.000
#> 8  Item1  Item9  0.124   1.000
#> 9  Item1 Item10 -0.039   1.000
#> 10 Item2  Item3 -0.206   1.000
#> 11 Item2  Item4 -0.202   1.000
#> 12 Item2  Item5  0.058   1.000
#> 13 Item2  Item6  0.025   1.000
#> 14 Item2  Item7 -0.230   1.000
#> 15 Item2  Item8  0.150   1.000
#> 16 Item2  Item9  0.154   1.000
#> 17 Item2 Item10  0.024   1.000
#> 18 Item3  Item4  0.021   1.000
#> 19 Item3  Item5  0.190   1.000
#> 20 Item3  Item6 -0.143   1.000
#> 21 Item3  Item7  0.086   1.000
#> 22 Item3  Item8 -0.042   1.000
#> 23 Item3  Item9 -0.074   1.000
#> 24 Item3 Item10  0.221   1.000
#> 25 Item4  Item5 -0.024   1.000
#> 26 Item4  Item6  0.224   1.000
#> 27 Item4  Item7  0.333   1.000
#> 28 Item4  Item8 -0.406   0.305
#> 29 Item4  Item9 -0.153   1.000
#> 30 Item4 Item10 -0.084   1.000
#> 31 Item5  Item6 -0.136   1.000
#> 32 Item5  Item7  0.157   1.000
#> 33 Item5  Item8  0.130   1.000
#> 34 Item5  Item9 -0.348   1.000
#> 35 Item5 Item10  0.014   1.000
#> 36 Item6  Item7  0.004   1.000
#> 37 Item6  Item8 -0.133   1.000
#> 38 Item6  Item9  0.062   1.000
#> 39 Item6 Item10 -0.243   1.000
#> 40 Item7  Item8 -0.198   1.000
#> 41 Item7  Item9 -0.005   1.000
#> 42 Item7 Item10  0.021   1.000
#> 43 Item8  Item9  0.301   1.000
#> 44 Item8 Item10  0.045   1.000
#> 45 Item9 Item10  0.036   1.000
#> 
#> [[2]]
#>     Item1 Item2  gamma padj_bh
#> 1   Item2 Item1  0.300   1.000
#> 2   Item3 Item1 -0.158   1.000
#> 3   Item3 Item2 -0.230   1.000
#> 4   Item4 Item1 -0.141   1.000
#> 5   Item4 Item2 -0.263   1.000
#> 6   Item4 Item3 -0.042   1.000
#> 7   Item5 Item1 -0.152   1.000
#> 8   Item5 Item2  0.039   1.000
#> 9   Item5 Item3  0.201   1.000
#> 10  Item5 Item4  0.021   1.000
#> 11  Item6 Item1 -0.102   1.000
#> 12  Item6 Item2 -0.063   1.000
#> 13  Item6 Item3 -0.208   1.000
#> 14  Item6 Item4  0.206   1.000
#> 15  Item6 Item5 -0.201   1.000
#> 16  Item7 Item1 -0.116   1.000
#> 17  Item7 Item2 -0.195   1.000
#> 18  Item7 Item3  0.211   1.000
#> 19  Item7 Item4  0.428   0.126
#> 20  Item7 Item5  0.218   1.000
#> 21  Item7 Item6  0.093   1.000
#> 22  Item8 Item1  0.237   1.000
#> 23  Item8 Item2  0.220   1.000
#> 24  Item8 Item3  0.024   1.000
#> 25  Item8 Item4 -0.319   1.000
#> 26  Item8 Item5  0.152   1.000
#> 27  Item8 Item6 -0.022   1.000
#> 28  Item8 Item7 -0.222   1.000
#> 29  Item9 Item1  0.193   1.000
#> 30  Item9 Item2  0.199   1.000
#> 31  Item9 Item3 -0.008   1.000
#> 32  Item9 Item4 -0.071   1.000
#> 33  Item9 Item5 -0.256   1.000
#> 34  Item9 Item6  0.218   1.000
#> 35  Item9 Item7 -0.001   1.000
#> 36  Item9 Item8  0.357   1.000
#> 37 Item10 Item1 -0.071   1.000
#> 38 Item10 Item2 -0.023   1.000
#> 39 Item10 Item3  0.176   1.000
#> 40 Item10 Item4 -0.088   1.000
#> 41 Item10 Item5 -0.029   1.000
#> 42 Item10 Item6 -0.235   1.000
#> 43 Item10 Item7 -0.040   1.000
#> 44 Item10 Item8 -0.047   1.000
#> 45 Item10 Item9 -0.069   1.000
#> 
# }
# \donttest{
# Simulation-based cutoffs (slow): 100+ Monte-Carlo iterations
cutoff_res <- RMpgLDcutoff(sim_data, iterations = 100, parallel = FALSE,
                           seed = 42)
RMpartgamLD(sim_data, cutoff = cutoff_res)
#>  [1] "Table: Partial gamma LD analysis (n = 200 complete cases). Cutoff values based on 100 simulation iterations (99% HDCI).\n\nDirection 1: rest score = total - Item2 \n\n Table: Partial gamma LD analysis (n = 200 complete cases). Cutoff values based on 100 simulation iterations (99% HDCI).\n\nDirection 2: rest score = total - Item1"
#>  [2] " \n\n "                                                                                                                                                                                                                                                                                                                                    
#>  [3] "|Item 1 |Item 2 | Partial gamma| Adjusted p-value (BH)| Gamma low| Gamma high|Flagged | \n\n |Item 1 |Item 2 | Partial gamma| Adjusted p-value (BH)| Gamma low| Gamma high|Flagged |"                                                                                                                                                      
#>  [4] "|:------|:------|-------------:|---------------------:|---------:|----------:|:-------| \n\n |:------|:------|-------------:|---------------------:|---------:|----------:|:-------|"                                                                                                                                                      
#>  [5] "|Item1  |Item2  |         0.320|                 1.000|    -0.431|      0.431|FALSE   | \n\n |Item2  |Item1  |         0.300|                 1.000|    -0.431|      0.431|FALSE   |"                                                                                                                                                      
#>  [6] "|Item1  |Item3  |        -0.143|                 1.000|    -0.358|      0.325|FALSE   | \n\n |Item3  |Item1  |        -0.158|                 1.000|    -0.358|      0.325|FALSE   |"                                                                                                                                                      
#>  [7] "|Item1  |Item4  |        -0.069|                 1.000|    -0.509|      0.430|FALSE   | \n\n |Item3  |Item2  |        -0.230|                 1.000|    -0.420|      0.399|FALSE   |"                                                                                                                                                      
#>  [8] "|Item1  |Item5  |        -0.153|                 1.000|    -0.395|      0.307|FALSE   | \n\n |Item4  |Item1  |        -0.141|                 1.000|    -0.509|      0.430|FALSE   |"                                                                                                                                                      
#>  [9] "|Item1  |Item6  |        -0.031|                 1.000|    -0.464|      0.524|FALSE   | \n\n |Item4  |Item2  |        -0.263|                 1.000|    -0.364|      0.353|FALSE   |"                                                                                                                                                      
#> [10] "|Item1  |Item7  |        -0.163|                 1.000|    -0.473|      0.415|FALSE   | \n\n |Item4  |Item3  |        -0.042|                 1.000|    -0.393|      0.411|FALSE   |"                                                                                                                                                      
#> [11] "|Item1  |Item8  |         0.191|                 1.000|    -0.502|      0.339|FALSE   | \n\n |Item5  |Item1  |        -0.152|                 1.000|    -0.395|      0.307|FALSE   |"                                                                                                                                                      
#> [12] "|Item1  |Item9  |         0.124|                 1.000|    -0.525|      0.382|FALSE   | \n\n |Item5  |Item2  |         0.039|                 1.000|    -0.447|      0.488|FALSE   |"                                                                                                                                                      
#> [13] "|Item1  |Item10 |        -0.039|                 1.000|    -0.376|      0.370|FALSE   | \n\n |Item5  |Item3  |         0.201|                 1.000|    -0.414|      0.330|FALSE   |"                                                                                                                                                      
#> [14] "|Item2  |Item3  |        -0.206|                 1.000|    -0.420|      0.399|FALSE   | \n\n |Item5  |Item4  |         0.021|                 1.000|    -0.390|      0.369|FALSE   |"                                                                                                                                                      
#> [15] "|Item2  |Item4  |        -0.202|                 1.000|    -0.364|      0.353|FALSE   | \n\n |Item6  |Item1  |        -0.102|                 1.000|    -0.464|      0.524|FALSE   |"                                                                                                                                                      
#> [16] "|Item2  |Item5  |         0.058|                 1.000|    -0.447|      0.488|FALSE   | \n\n |Item6  |Item2  |        -0.063|                 1.000|    -0.373|      0.373|FALSE   |"                                                                                                                                                      
#> [17] "|Item2  |Item6  |         0.025|                 1.000|    -0.373|      0.373|FALSE   | \n\n |Item6  |Item3  |        -0.208|                 1.000|    -0.382|      0.430|FALSE   |"                                                                                                                                                      
#> [18] "|Item2  |Item7  |        -0.230|                 1.000|    -0.438|      0.409|FALSE   | \n\n |Item6  |Item4  |         0.206|                 1.000|    -0.405|      0.410|FALSE   |"                                                                                                                                                      
#> [19] "|Item2  |Item8  |         0.150|                 1.000|    -0.402|      0.459|FALSE   | \n\n |Item6  |Item5  |        -0.201|                 1.000|    -0.355|      0.489|FALSE   |"                                                                                                                                                      
#> [20] "|Item2  |Item9  |         0.154|                 1.000|    -0.339|      0.394|FALSE   | \n\n |Item7  |Item1  |        -0.116|                 1.000|    -0.473|      0.415|FALSE   |"                                                                                                                                                      
#> [21] "|Item2  |Item10 |         0.024|                 1.000|    -0.354|      0.309|FALSE   | \n\n |Item7  |Item2  |        -0.195|                 1.000|    -0.438|      0.409|FALSE   |"                                                                                                                                                      
#> [22] "|Item3  |Item4  |         0.021|                 1.000|    -0.393|      0.411|FALSE   | \n\n |Item7  |Item3  |         0.211|                 1.000|    -0.435|      0.384|FALSE   |"                                                                                                                                                      
#> [23] "|Item3  |Item5  |         0.190|                 1.000|    -0.414|      0.330|FALSE   | \n\n |Item7  |Item4  |         0.428|                 0.126|    -0.466|      0.357|TRUE    |"                                                                                                                                                      
#> [24] "|Item3  |Item6  |        -0.143|                 1.000|    -0.382|      0.430|FALSE   | \n\n |Item7  |Item5  |         0.218|                 1.000|    -0.412|      0.392|FALSE   |"                                                                                                                                                      
#> [25] "|Item3  |Item7  |         0.086|                 1.000|    -0.435|      0.384|FALSE   | \n\n |Item7  |Item6  |         0.093|                 1.000|    -0.398|      0.465|FALSE   |"                                                                                                                                                      
#> [26] "|Item3  |Item8  |        -0.042|                 1.000|    -0.411|      0.377|FALSE   | \n\n |Item8  |Item1  |         0.237|                 1.000|    -0.502|      0.339|FALSE   |"                                                                                                                                                      
#> [27] "|Item3  |Item9  |        -0.074|                 1.000|    -0.385|      0.384|FALSE   | \n\n |Item8  |Item2  |         0.220|                 1.000|    -0.402|      0.459|FALSE   |"                                                                                                                                                      
#> [28] "|Item3  |Item10 |         0.221|                 1.000|    -0.375|      0.425|FALSE   | \n\n |Item8  |Item3  |         0.024|                 1.000|    -0.411|      0.377|FALSE   |"                                                                                                                                                      
#> [29] "|Item4  |Item5  |        -0.024|                 1.000|    -0.390|      0.369|FALSE   | \n\n |Item8  |Item4  |        -0.319|                 1.000|    -0.496|      0.471|FALSE   |"                                                                                                                                                      
#> [30] "|Item4  |Item6  |         0.224|                 1.000|    -0.405|      0.410|FALSE   | \n\n |Item8  |Item5  |         0.152|                 1.000|    -0.436|      0.448|FALSE   |"                                                                                                                                                      
#> [31] "|Item4  |Item7  |         0.333|                 1.000|    -0.466|      0.357|FALSE   | \n\n |Item8  |Item6  |        -0.022|                 1.000|    -0.376|      0.393|FALSE   |"                                                                                                                                                      
#> [32] "|Item4  |Item8  |        -0.406|                 0.305|    -0.496|      0.471|FALSE   | \n\n |Item8  |Item7  |        -0.222|                 1.000|    -0.351|      0.497|FALSE   |"                                                                                                                                                      
#> [33] "|Item4  |Item9  |        -0.153|                 1.000|    -0.441|      0.467|FALSE   | \n\n |Item9  |Item1  |         0.193|                 1.000|    -0.525|      0.382|FALSE   |"                                                                                                                                                      
#> [34] "|Item4  |Item10 |        -0.084|                 1.000|    -0.330|      0.296|FALSE   | \n\n |Item9  |Item2  |         0.199|                 1.000|    -0.339|      0.394|FALSE   |"                                                                                                                                                      
#> [35] "|Item5  |Item6  |        -0.136|                 1.000|    -0.355|      0.489|FALSE   | \n\n |Item9  |Item3  |        -0.008|                 1.000|    -0.385|      0.384|FALSE   |"                                                                                                                                                      
#> [36] "|Item5  |Item7  |         0.157|                 1.000|    -0.412|      0.392|FALSE   | \n\n |Item9  |Item4  |        -0.071|                 1.000|    -0.441|      0.467|FALSE   |"                                                                                                                                                      
#> [37] "|Item5  |Item8  |         0.130|                 1.000|    -0.436|      0.448|FALSE   | \n\n |Item9  |Item5  |        -0.256|                 1.000|    -0.325|      0.342|FALSE   |"                                                                                                                                                      
#> [38] "|Item5  |Item9  |        -0.348|                 1.000|    -0.325|      0.342|TRUE    | \n\n |Item9  |Item6  |         0.218|                 1.000|    -0.367|      0.417|FALSE   |"                                                                                                                                                      
#> [39] "|Item5  |Item10 |         0.014|                 1.000|    -0.407|      0.344|FALSE   | \n\n |Item9  |Item7  |        -0.001|                 1.000|    -0.336|      0.480|FALSE   |"                                                                                                                                                      
#> [40] "|Item6  |Item7  |         0.004|                 1.000|    -0.398|      0.465|FALSE   | \n\n |Item9  |Item8  |         0.357|                 1.000|    -0.445|      0.431|FALSE   |"                                                                                                                                                      
#> [41] "|Item6  |Item8  |        -0.133|                 1.000|    -0.376|      0.393|FALSE   | \n\n |Item10 |Item1  |        -0.071|                 1.000|    -0.376|      0.370|FALSE   |"                                                                                                                                                      
#> [42] "|Item6  |Item9  |         0.062|                 1.000|    -0.367|      0.417|FALSE   | \n\n |Item10 |Item2  |        -0.023|                 1.000|    -0.354|      0.309|FALSE   |"                                                                                                                                                      
#> [43] "|Item6  |Item10 |        -0.243|                 1.000|    -0.410|      0.451|FALSE   | \n\n |Item10 |Item3  |         0.176|                 1.000|    -0.375|      0.425|FALSE   |"                                                                                                                                                      
#> [44] "|Item7  |Item8  |        -0.198|                 1.000|    -0.351|      0.497|FALSE   | \n\n |Item10 |Item4  |        -0.088|                 1.000|    -0.330|      0.296|FALSE   |"                                                                                                                                                      
#> [45] "|Item7  |Item9  |        -0.005|                 1.000|    -0.336|      0.480|FALSE   | \n\n |Item10 |Item5  |        -0.029|                 1.000|    -0.407|      0.344|FALSE   |"                                                                                                                                                      
#> [46] "|Item7  |Item10 |         0.021|                 1.000|    -0.353|      0.384|FALSE   | \n\n |Item10 |Item6  |        -0.235|                 1.000|    -0.410|      0.451|FALSE   |"                                                                                                                                                      
#> [47] "|Item8  |Item9  |         0.301|                 1.000|    -0.445|      0.431|FALSE   | \n\n |Item10 |Item7  |        -0.040|                 1.000|    -0.353|      0.384|FALSE   |"                                                                                                                                                      
#> [48] "|Item8  |Item10 |         0.045|                 1.000|    -0.340|      0.471|FALSE   | \n\n |Item10 |Item8  |        -0.047|                 1.000|    -0.340|      0.471|FALSE   |"                                                                                                                                                      
#> [49] "|Item9  |Item10 |         0.036|                 1.000|    -0.381|      0.411|FALSE   | \n\n |Item10 |Item9  |        -0.069|                 1.000|    -0.381|      0.411|FALSE   |"                                                                                                                                                      
#> attr(,"class")
#> [1] "knit_asis"
#> attr(,"knit_cacheable")
#> [1] NA
# }
```
