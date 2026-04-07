# Run a single partial gamma DIF simulation iteration

Run a single partial gamma DIF simulation iteration

## Usage

``` r
run_single_partgam_sim(seed, data_list)
```

## Arguments

- seed:

  Integer seed for reproducibility.

- data_list:

  List produced inside
  [`RMpgDIFcutoff()`](https://pgmj.github.io/easyRasch2/reference/RMpgDIFcutoff.md).

## Value

A data.frame with columns `Item` and `gamma`, or a character string on
failure.
