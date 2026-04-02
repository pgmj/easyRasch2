# Run a single Q3 simulation iteration

Run a single Q3 simulation iteration

## Usage

``` r
run_single_q3_sim(seed, data_list)
```

## Arguments

- seed:

  Integer seed for reproducibility.

- data_list:

  List produced inside
  [`RMlocdepQ3cutoff()`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3cutoff.md).

## Value

A list with `mean` and `max` Q3, or a character string on failure.
