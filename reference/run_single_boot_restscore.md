# Run a single item-restscore bootstrap iteration

Run a single item-restscore bootstrap iteration

## Usage

``` r
run_single_boot_restscore(seed, data_list)
```

## Arguments

- seed:

  Integer seed for reproducibility.

- data_list:

  List produced inside
  [`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/reference/RMitemRestscoreBoot.md)
  containing `data`, `samplesize`, `is_polytomous`, `item_names`.

## Value

A data.frame with columns `Item`, `item_restscore`, `diff`, `diff_abs`,
or a character string on failure.
