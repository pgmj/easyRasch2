# Run a single infit simulation iteration

Run a single infit simulation iteration

## Usage

``` r
run_single_infit_sim(seed, data_list)
```

## Arguments

- seed:

  Integer seed for reproducibility.

- data_list:

  List produced inside
  [`RMinfitcutoff()`](https://pgmj.github.io/easyRasch2/reference/RMinfitcutoff.md).

## Value

A data.frame with columns `Item`, `InfitMSQ`, `OutfitMSQ`, or a
character string on failure.
