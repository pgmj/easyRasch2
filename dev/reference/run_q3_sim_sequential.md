# Run Q3 simulations sequentially

Run Q3 simulations sequentially

## Usage

``` r
run_q3_sim_sequential(iterations, sim_seeds, sim_data_list, verbose = FALSE)
```

## Arguments

- iterations:

  Number of iterations.

- sim_seeds:

  Integer vector of per-iteration seeds.

- sim_data_list:

  List of data passed to each worker.

- verbose:

  Show progress bar.

## Value

List of raw results (one element per iteration).
