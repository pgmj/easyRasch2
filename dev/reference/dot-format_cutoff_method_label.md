# Format a human-readable label for the cutoff method

Format a human-readable label for the cutoff method

## Usage

``` r
.format_cutoff_method_label(cutoff_method, hdci_width)
```

## Arguments

- cutoff_method:

  Character. `"hdci"`, `"quantile"`, or `NULL`.

- hdci_width:

  Numeric or `NULL`. HDCI width (e.g., `0.999`).

## Value

A character label, or `NULL` if the method is unknown/unset.
