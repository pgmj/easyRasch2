# easyRasch2

<!-- badges: start -->
<!-- badges: end -->

`easyRasch2` is a CRAN-targeted R package that provides functions for Rasch
measurement theory analysis workflows. It is the successor to
[easyRasch](https://github.com/pgmj/easyRasch), offering a lightweight,
CRAN-ready structure with proper namespacing and minimal dependencies.

## Key design principles

- **Estimation**: Conditional Maximum Likelihood (CML) via `eRm`,
  `psychotools`, and `iarm` for item parameters; Weighted Likelihood
  Estimation (WLE) for person parameters; `mirt` (MML) only where CML
  alternatives do not exist (e.g., Q3 residual correlations).
- **Output**: `knitr::kable()` for tables (Quarto-friendly) and `ggplot2`
  for figures, with `"dataframe"` output options for further processing.
- **Naming**: Functions use the `RM` prefix (e.g., `RMlocdepQ3()`).

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("pgmj/easyRasch2")
```

## Example

```r
library(easyRasch2)

# Simulate binary item response data
set.seed(42)
sim_data <- as.data.frame(
  matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
)
colnames(sim_data) <- paste0("Item", 1:10)

# Q3 residual correlations (dynamic cut-off = mean Q3 + 0.2)
RMlocdepQ3(sim_data, cutoff = 0.2)

# Return the underlying data.frame instead
q3_df <- RMlocdepQ3(sim_data, cutoff = 0.2, output = "dataframe")
```

## License

GPL (>= 3)
