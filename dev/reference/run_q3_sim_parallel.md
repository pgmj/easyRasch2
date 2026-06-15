# Run Q3 simulations in parallel using mirai

Run Q3 simulations in parallel using mirai

## Usage

``` r
run_q3_sim_parallel(
  iterations,
  sim_seeds,
  sim_data_list,
  n_cores,
  verbose = FALSE
)
```

## Arguments

- iterations:

  Number of iterations.

- sim_seeds:

  Integer vector of per-iteration seeds.

- sim_data_list:

  List of data passed to each worker.

- n_cores:

  Number of mirai daemons.

- verbose:

  Show progress bar.

## Value

List of raw results (one element per iteration).
