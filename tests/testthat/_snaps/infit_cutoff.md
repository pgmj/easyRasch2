# errors when no complete cases remain after na.omit

    Code
      RMitemInfitCutoff(d, iterations = 3, parallel = FALSE)
    Condition
      Error:
      ! No complete cases in data. All rows contain at least one NA.

# parallel = TRUE without n_cores or mc.cores warns and falls back

    Code
      res <- RMitemInfitCutoff(make_dich(), iterations = 4, parallel = TRUE, seed = 1,
      cutoff_method = "quantile")
    Condition
      Warning:
      For parallel processing, specify n_cores or set options(mc.cores = N).
      (Use `parallel::detectCores()` to see how many cores are available.)
      Falling back to sequential (single core) processing.

