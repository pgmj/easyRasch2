# Format a human-readable label for the gamma cutoff method

Format a human-readable label for the gamma cutoff method

## Usage

``` r
.format_gamma_cutoff_method_label(cutoff_method, hdci_width)
```

## Arguments

- cutoff_method:

  Character. `"hdci"`, `"quantile"`, or `NULL`.

- hdci_width:

  Numeric or `NULL`. HDCI width (e.g., `0.99`).

## Value

A character label, or `NULL` if the method is unknown/unset.
