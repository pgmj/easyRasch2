# Run a single partial gamma LD simulation iteration

Run a single partial gamma LD simulation iteration

## Usage

``` r
run_single_partgam_LD_sim(seed, data_list)
```

## Arguments

- seed:

  Integer seed for reproducibility.

- data_list:

  List produced inside
  [`RMpgLDcutoff()`](https://pgmj.github.io/easyRasch2/reference/RMpgLDcutoff.md).

## Value

A data.frame with columns `Item1`, `Item2`, and `gamma`, or a character
string on failure.
