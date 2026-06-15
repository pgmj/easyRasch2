# Run item-restscore bootstrap iterations in parallel using mirai

Run item-restscore bootstrap iterations in parallel using mirai

## Usage

``` r
run_boot_restscore_parallel(
  iterations,
  boot_seeds,
  boot_data_list,
  n_cores,
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

- n_cores:

  Number of mirai daemons.

- verbose:

  Show progress bar.

## Value

List of raw results (one element per iteration).
