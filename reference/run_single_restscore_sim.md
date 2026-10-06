# Run a single item-restscore simulation iteration

Run a single item-restscore simulation iteration

## Usage

``` r
run_single_restscore_sim(seed, data_list)
```

## Arguments

- seed:

  Integer seed for reproducibility.

- data_list:

  List produced inside
  [`RMitemRestscoreCutoff()`](https://pgmj.github.io/easyRasch2/reference/RMitemRestscoreCutoff.md).

## Value

A data.frame with columns `Item`, `Observed`, `Expected`, `Difference`,
or a character string on failure.
