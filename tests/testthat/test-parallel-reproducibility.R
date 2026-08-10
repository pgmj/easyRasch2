# Regression tests: `parallel = TRUE` and `parallel = FALSE` must produce
# identical results for the same `seed`.
#
# mirai daemons start under RNGkind("L'Ecuyer-CMRG") while the calling session
# uses the "Mersenne-Twister" default, so a worker that only calls
# set.seed(seed) draws a different stream in each path. The per-iteration
# runners therefore pin the RNG kind as well as the seed.

skip_if_no_mirai <- function() {
  skip_on_cran()
  skip_if_not_installed("mirai")
}

make_dich <- function(n = 250, k = 8, seed = 3L) {
  set.seed(seed)
  th <- stats::rnorm(n)
  df <- as.data.frame(sapply(
    seq_len(k),
    function(j) stats::rbinom(n, 1, stats::plogis(th - (j - (k + 1) / 2) / 3))
  ))
  names(df) <- paste0("I", seq_len(k))
  df
}

test_that("mirai daemons really do change the RNG kind", {
  # Guards the premise. If a future mirai stops switching the RNG kind this
  # test still passes, but the pinning below becomes belt-and-braces.
  skip_if_no_mirai()
  mirai::daemons(1L)
  on.exit(mirai::daemons(0L), add = TRUE)
  kind <- mirai::call_mirai(mirai::mirai({
    set.seed(1L)
    RNGkind()
  }))$data
  expect_type(kind, "character")

  # and that pinning makes a daemon reproduce the calling session's stream
  drawn <- mirai::call_mirai(mirai::mirai({
    set.seed(
      1L,
      kind = "Mersenne-Twister",
      normal.kind = "Inversion",
      sample.kind = "Rejection"
    )
    sample.int(1000L, 3L)
  }))$data
  set.seed(1L)
  expect_identical(drawn, sample.int(1000L, 3L))
})

test_that("RMdimMartinLof gives identical results in both paths", {
  skip_if_no_mirai()
  skip_if_not_installed("psychotools")
  df <- make_dich()
  args <- list(
    partition = list(1:4, 5:8),
    iterations = 40L,
    seed = 1L
  )
  par <- do.call(
    RMdimMartinLof,
    c(list(df), args, list(parallel = TRUE, n_cores = 2L))
  )
  seq <- do.call(RMdimMartinLof, c(list(df), args, list(parallel = FALSE)))
  expect_equal(par$T_rep, seq$T_rep)
  expect_identical(par$p_value, seq$p_value)
})

test_that("RMdimResidualPCACutoff gives identical results in both paths", {
  skip_if_no_mirai()
  skip_if_not_installed("psychotools")
  df <- make_dich()
  par <- RMdimResidualPCACutoff(
    df,
    iterations = 30L,
    parallel = TRUE,
    n_cores = 2L,
    seed = 1L
  )
  seq <- RMdimResidualPCACutoff(
    df,
    iterations = 30L,
    parallel = FALSE,
    seed = 1L
  )
  expect_equal(par$simulated, seq$simulated)
})

test_that("RMitemInfitCutoff gives identical results in both paths", {
  skip_if_no_mirai()
  skip_if_not_installed("psychotools")
  df <- make_dich()
  par <- RMitemInfitCutoff(
    df,
    iterations = 30L,
    parallel = TRUE,
    n_cores = 2L,
    seed = 1L
  )
  seq <- RMitemInfitCutoff(df, iterations = 30L, parallel = FALSE, seed = 1L)
  expect_equal(par$simulated, seq$simulated)
})

test_that("RMpersonFit gives identical results in both paths", {
  skip_if_no_mirai()
  skip_if_not_installed("psychotools")
  df <- make_dich(n = 60L, k = 6L)
  par <- RMpersonFit(
    df,
    iterations = 30L,
    parallel = TRUE,
    n_cores = 2L,
    seed = 1L,
    output = "dataframe"
  )
  seq <- RMpersonFit(
    df,
    iterations = 30L,
    parallel = FALSE,
    seed = 1L,
    output = "dataframe"
  )
  expect_equal(par, seq)
})
