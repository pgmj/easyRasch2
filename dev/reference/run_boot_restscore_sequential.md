# Run item-restscore bootstrap iterations sequentially

Run item-restscore bootstrap iterations sequentially

## Usage

``` r
run_boot_restscore_sequential(
  iterations,
  boot_seeds,
  boot_data_list,
  verbose = FALSE
)
```

## Arguments

- iterations:

  Number of iterations.

- boot_seeds:

  Integer vector of per-iteration seeds.

- boot_data_list:

  List of data passed to each worker.

- verbose:

  Show progress bar.

## Value

List of raw results (one element per iteration).
