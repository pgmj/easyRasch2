# Tests for RMdimMartinLof() and RMdimMartinLofResiduals()

make_polytomous <- function(n = 200, k = 8, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_dichotomous <- function(n = 300, k = 8, seed = 2L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMdimMartinLof errors when partition has fewer than 2 subscales", {
  df <- make_polytomous()
  expect_error(
    RMdimMartinLof(df, partition = rep(1L, ncol(df)),
                iterations = 5L, parallel = FALSE, seed = 1L),
    regexp = ">= 2 subscales"
  )
})

test_that("RMdimMartinLof errors when a subscale has fewer than 2 items", {
  df <- make_polytomous()
  partition <- c(1L, rep(2L, ncol(df) - 1L))
  expect_error(
    RMdimMartinLof(df, partition = partition,
                iterations = 5L, parallel = FALSE, seed = 1L),
    regexp = "at least 2 items"
  )
})

# ---------------------------------------------------------------------
# Output structure (polytomous)
# ---------------------------------------------------------------------
test_that("RMdimMartinLof returns expected fields on polytomous data", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_polytomous()
  res <- RMdimMartinLof(
    df, partition = list(c("I1","I2","I3","I4"), c("I5","I6","I7","I8")),
    iterations = 20L, parallel = FALSE, seed = 1L,
    stopping = "sequential", h = 10L
  )
  expect_type(res, "list")
  expect_true(all(c("T_obs", "p_value", "actual_iterations", "rejected",
                    "partition", "n_subscales", "is_polytomous",
                    "sample_n", "n_items", "stopping", "h",
                    "T_rep", "wle_scores", "wle_correlation") %in% names(res)))
  expect_true(res$is_polytomous)
  expect_equal(res$n_subscales, 2L)
  expect_true(is.finite(res$T_obs) && res$T_obs >= 0)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

# ---------------------------------------------------------------------
# Output structure (dichotomous)
# ---------------------------------------------------------------------
test_that("RMdimMartinLof works on dichotomous data (Rasch model path)", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  res <- RMdimMartinLof(
    df, partition = c(1,1,1,1,2,2,2,2),
    iterations = 20L, parallel = FALSE, seed = 1L,
    stopping = "sequential", h = 10L
  )
  expect_type(res, "list")
  expect_false(res$is_polytomous)
  expect_true(is.finite(res$T_obs))
})

test_that("RMdimMartinLof p_value is reproducible with the same seed", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_polytomous(n = 150)
  r1 <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                    iterations = 10L, parallel = FALSE, seed = 42L,
                    stopping = "sequential", h = 5L)
  r2 <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                    iterations = 10L, parallel = FALSE, seed = 42L,
                    stopping = "sequential", h = 5L)
  expect_equal(r1$T_obs,  r2$T_obs)
  expect_equal(r1$T_rep,  r2$T_rep)
})

# ---------------------------------------------------------------------
# RMdimMartinLofResiduals
# ---------------------------------------------------------------------
test_that("RMdimMartinLofResiduals output = 'dataframe' returns expected columns", {
  skip_if_not_installed("psychotools")
  df <- make_polytomous()
  res <- RMdimMartinLofResiduals(
    df, partition = list(c("I1","I2","I3","I4"), c("I5","I6","I7","I8")),
    output = "dataframe"
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("t1", "t2", "total", "observed", "expected",
                    "residual", "flagged") %in% names(res)))
  expect_type(res$flagged, "logical")
})

test_that("RMdimMartinLofResiduals output = 'kable' returns a kable-like object", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("knitr")
  df  <- make_polytomous()
  out <- RMdimMartinLofResiduals(
    df, partition = c(1,1,1,1,2,2,2,2),
    output = "kable"
  )
  expect_true(inherits(out, "knitr_kable") || inherits(out, "knit_asis"))
})

test_that("RMdimMartinLofResiduals output = 'ggplot' returns a ggplot", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p  <- RMdimMartinLofResiduals(
    df, partition = c(1,1,1,1,2,2,2,2),
    output = "ggplot"
  )
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------
# stopping = "none" Monte Carlo runner (the default; existing tests above
# all use stopping = "sequential")
# ---------------------------------------------------------------------
test_that("RMdimMartinLof stopping = 'none' returns a valid p-value", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  res <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                        iterations = 20L, stopping = "none",
                        parallel = FALSE, seed = 1L)
  expect_equal(res$stopping, "none")
  expect_equal(res$actual_iterations, 20L)
  expect_true(res$p_value > 0 && res$p_value <= 1)
})

test_that("RMdimMartinLof verbose = TRUE runs the sequential progress path", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  suppressMessages(invisible(utils::capture.output(
    res <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                          iterations = 10L, stopping = "none",
                          parallel = FALSE, seed = 1L, verbose = TRUE)
  )))
  expect_true(res$p_value <= 1)
})

test_that("RMdimMartinLof parallel path returns a valid result", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  skip_if_not_installed("mirai")
  df <- make_dichotomous()
  res <- suppressMessages(
    RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                   iterations = 8L, stopping = "none",
                   parallel = TRUE, n_cores = 2L, seed = 1L)
  )
  expect_true(res$p_value > 0 && res$p_value <= 1)
})

test_that("RMdimMartinLof parallel without n_cores/mc.cores warns and falls back", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  skip_if_not_installed("mirai")
  withr::local_options(mc.cores = NULL)
  df <- make_dichotomous()
  expect_warning(
    res <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                          iterations = 6L, stopping = "none",
                          parallel = TRUE, seed = 1L),
    regexp = "Falling back to sequential"
  )
  expect_true(res$p_value <= 1)
})

test_that("RMdimMartinLof errors with fewer than 30 complete cases", {
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()[1:20, ]
  expect_error(
    RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                   iterations = 5L, parallel = FALSE),
    regexp = "at least 30 complete cases"
  )
})

# ---------------------------------------------------------------------
# Partition normalisation branches
# ---------------------------------------------------------------------
test_that("RMdimMartinLof warns and drops items not assigned to a subscale", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()  # 8 items
  expect_warning(
    res <- RMdimMartinLof(df, partition = list(c(1, 2, 3), c(4, 5, 6)),
                          iterations = 10L, parallel = FALSE, seed = 1L),
    regexp = "not assigned"
  )
  expect_equal(res$n_items, 6L)
})

test_that("RMdimMartinLof errors on overlapping partition", {
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  expect_error(
    RMdimMartinLof(df, partition = list(c(1, 2, 3), c(3, 4, 5)),
                   iterations = 5L, parallel = FALSE),
    regexp = "overlapping"
  )
})

test_that("RMdimMartinLof errors on an invalid partition format", {
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  expect_error(
    RMdimMartinLof(df, partition = TRUE, iterations = 5L, parallel = FALSE),
    regexp = "Invalid partition format"
  )
})

# ---------------------------------------------------------------------
# RMdimMartinLofResiduals: option branches
# ---------------------------------------------------------------------
test_that("RMdimMartinLofResiduals min_expected blanks sparse cells", {
  skip_if_not_installed("psychotools")
  df  <- make_polytomous()
  res <- RMdimMartinLofResiduals(df, partition = c(1,1,1,1,2,2,2,2),
                                 output = "dataframe", min_expected = 1)
  expect_true(all(is.na(res$residual[res$expected < 1])))
})

test_that("RMdimMartinLofResiduals ggplot color_by = 'n' returns a ggplot", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p  <- RMdimMartinLofResiduals(df, partition = c(1,1,1,1,2,2,2,2),
                                output = "ggplot", color_by = "n",
                                color_limits = c(0, 50))
  expect_s3_class(p, "ggplot")
})

test_that("RMdimMartinLofResiduals supports D = 3 for kable and ggplot", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("knitr")
  df <- make_polytomous(k = 6)
  k  <- RMdimMartinLofResiduals(df, partition = c(1,1,2,2,3,3), output = "kable")
  expect_true(inherits(k, "knitr_kable") || inherits(k, "knit_asis"))
  skip_if_not_installed("ggplot2")
  p <- RMdimMartinLofResiduals(df, partition = c(1,1,2,2,3,3), output = "ggplot")
  expect_s3_class(p, "ggplot")
})

test_that("RMdimMartinLofResiduals rejects bad min_expected / color_limits / D>3 ggplot", {
  skip_if_not_installed("psychotools")
  df <- make_polytomous()
  expect_error(
    RMdimMartinLofResiduals(df, partition = c(1,1,1,1,2,2,2,2),
                            min_expected = -1),
    regexp = "min_expected"
  )
  expect_error(
    RMdimMartinLofResiduals(df, partition = c(1,1,1,1,2,2,2,2),
                            color_limits = c(5, 1)),
    regexp = "color_limits"
  )
  expect_error(
    RMdimMartinLofResiduals(df, partition = c(1,1,2,2,3,3,4,4),
                            output = "ggplot"),
    regexp = "D = 2 or 3"
  )
})

test_that("p_value_floor reports the smallest attainable MC p-value", {
  set.seed(11)
  dat <- as.data.frame(replicate(6, sample(0:2, 120, replace = TRUE)))
  names(dat) <- paste0("q", 1:6)
  res <- RMdimMartinLof(dat, partition = c(1, 1, 1, 2, 2, 2),
                        iterations = 100, parallel = FALSE, seed = 3)
  expect_equal(res$p_value_floor, 1 / (res$actual_iterations + 1))
  expect_gte(res$p_value, res$p_value_floor)
})

test_that("polytomous sampler draws the exact conditional pattern distribution", {
  # Regression for the 1.0.0 sign-convention bug: the numerator weights and
  # the gamma functions must use the same convention (category k of item j
  # carries weight exp(-par_jk)).
  taus <- list(c(-0.8, 0.4), c(0.2, 0.9), c(-0.3, -0.1))
  params <- lapply(taus, cumsum)
  psi <- lapply(params, function(p) exp(-c(0, p)))
  pats <- expand.grid(x1 = 0:2, x2 = 0:2, x3 = 0:2)
  pats$w <- psi[[1]][pats$x1 + 1] * psi[[2]][pats$x2 + 1] *
    psi[[3]][pats$x3 + 1]
  pats$t <- pats$x1 + pats$x2 + pats$x3
  sub <- pats[pats$t == 3L, ]
  exact <- sub$w / sum(sub$w)

  set.seed(1)
  tabs <- build_ml_gamma_tables(params)
  draws <- t(replicate(20000, sample_pattern_at_score(3L, tabs)))
  key_d <- apply(draws, 1, paste, collapse = "")
  key_e <- paste0(sub$x1, sub$x2, sub$x3)
  freq <- as.numeric(table(factor(key_d, levels = key_e)) / nrow(draws))
  expect_lt(max(abs(freq - exact)), 0.015)

  # and invariant (in distribution) to item order
  set.seed(1)
  tabs_p <- build_ml_gamma_tables(params[c(3, 1, 2)])
  draws_p <- t(replicate(20000, sample_pattern_at_score(3L, tabs_p)))
  key_p <- apply(draws_p[, c(2, 3, 1)], 1, paste, collapse = "")
  freq_p <- as.numeric(table(factor(key_p, levels = key_e)) / nrow(draws_p))
  expect_lt(max(abs(freq_p - exact)), 0.015)
})

test_that("dichotomous sampler draws the exact conditional pattern distribution", {
  # Regression for the 1.1.1 dichotomous-sampler bug: the previous path used
  # sample.int(prob = , replace = FALSE), which is successive sampling
  # (Wallenius) and not the Rasch conditional distribution
  # p(x | t) prop. to prod_i phi_i^{x_i}. This is also the point at which
  # the implementation departs from Christensen & Kreiner (2007, p. 23),
  # whose dichotomous shortcut describes exactly that successive scheme.
  b <- c(-1.5, -0.4, 0.3, 1.6)
  phi <- exp(-b)
  pats <- as.matrix(expand.grid(rep(list(0:1), 4)))
  pats <- pats[rowSums(pats) == 2L, , drop = FALSE]
  w <- apply(pats, 1, function(x) prod(phi^x))
  exact <- w / sum(w)
  key_e <- apply(pats, 1, paste, collapse = "")

  set.seed(1)
  tabs <- build_ml_gamma_tables(b)
  draws <- t(replicate(20000, sample_pattern_at_score(2L, tabs)))
  freq <- as.numeric(
    table(factor(apply(draws, 1, paste, collapse = ""), levels = key_e))
  ) / nrow(draws)
  expect_lt(max(abs(freq - exact)), 0.015)

  # the boundary scores are deterministic
  expect_equal(sample_pattern_at_score(0L, tabs), integer(4))
  expect_equal(sample_pattern_at_score(4L, tabs), rep(1L, 4))
})

test_that("build_ml_gamma_tables matches the psychotools gamma functions", {
  # Up to the per-item rescaling, which is constant within a gamma vector.
  params <- lapply(list(c(-0.8, 0.4), c(0.2, 0.9), c(-0.3, -0.1)), cumsum)
  tabs <- build_ml_gamma_tables(params)
  ref <- compute_esf_gamma(params)
  ratio <- tabs$gamma[[3]] / ref
  expect_equal(max(ratio) - min(ratio), 0, tolerance = 1e-12)
  expect_equal(tabs$m_i, c(2L, 2L, 2L))
  expect_equal(tabs$M_total, 6L)
})

test_that("same seed gives the same p-value regardless of column arrangement", {
  set.seed(9)
  n <- 150
  th <- rnorm(n)
  dat <- as.data.frame(sapply(1:6, function(j)
    findInterval(plogis(th + rnorm(n, 0, 0.8)), c(0.3, 0.6, 0.85))))
  names(dat) <- paste0("w", 1:6)
  s1 <- paste0("w", 1:3); s2 <- paste0("w", 4:6)

  r1 <- RMdimMartinLof(dat, partition = list(s1, s2), iterations = 60,
                       parallel = FALSE, seed = 42)
  r2 <- RMdimMartinLof(dat[, c(s2, s1)], partition = list(s2, s1),
                       iterations = 60, parallel = FALSE, seed = 42)
  expect_identical(r1$p_value, r2$p_value)
  expect_equal(r1$T_obs, r2$T_obs)
  expect_equal(r1$T_rep, r2$T_rep, tolerance = 1e-9)
})

test_that("missingness on unpartitioned items does not drop respondents", {
  df <- make_dichotomous(n = 200, k = 8)
  part <- list(paste0("I", 1:3), paste0("I", 4:6))
  # NA only on I7/I8, which the partition never uses
  df$I7[1:60] <- NA
  df$I8[61:100] <- NA
  expect_equal(sum(stats::complete.cases(df)), 100L)

  res <- suppressWarnings(RMdimMartinLof(
    df, partition = part, iterations = 20L, parallel = FALSE, seed = 1L))
  expect_equal(res$sample_n, 200L)
  expect_equal(res$sample_n_total, 200L)
  expect_false(res$sample_has_na)
  expect_equal(nrow(res$wle_scores), 200L)

  rdf <- suppressWarnings(RMdimMartinLofResiduals(
    df, partition = part, output = "dataframe"))
  expect_equal(sum(rdf$observed), 200L)
})

test_that("missingness on partitioned items is reported the house way", {
  df <- make_dichotomous(n = 200, k = 8)
  part <- list(paste0("I", 1:3), paste0("I", 4:6))
  df$I2[1:20] <- NA

  res <- suppressWarnings(RMdimMartinLof(
    df, partition = part, iterations = 20L, parallel = FALSE, seed = 1L))
  expect_equal(res$sample_n, 180L)
  expect_equal(res$sample_n_total, 200L)
  expect_true(res$sample_has_na)

  k <- suppressWarnings(RMdimMartinLofResiduals(df, partition = part))
  expect_match(
    paste(as.character(k), collapse = " "),
    "n = 180 of 200 respondents \\(complete cases\\)"
  )
})

test_that("the 30-case minimum counts partitioned items only", {
  df <- make_dichotomous(n = 60, k = 8)
  part <- list(paste0("I", 1:3), paste0("I", 4:6))
  df$I8[1:40] <- NA  # only 20 rows complete overall, but 60 on I1-I6
  expect_no_error(suppressWarnings(RMdimMartinLof(
    df, partition = part, iterations = 5L, parallel = FALSE, seed = 1L)))

  df2 <- make_dichotomous(n = 60, k = 8)
  df2$I1[1:40] <- NA
  expect_error(
    suppressWarnings(RMdimMartinLof(
      df2, partition = part, iterations = 5L, parallel = FALSE, seed = 1L)),
    regexp = "30 complete cases"
  )
})

test_that("residual expected counts match brute-force enumeration", {
  set.seed(3)
  n <- 200
  th <- rnorm(n)
  dat <- as.data.frame(sapply(1:4, function(j)
    findInterval(plogis(th + rnorm(n, 0, 0.8)), c(0.3, 0.6, 0.85))))
  names(dat) <- paste0("v", 1:4)
  part <- list(c("v1", "v2"), c("v3", "v4"))
  res_df <- suppressWarnings(RMdimMartinLofResiduals(
    dat, partition = part, output = "dataframe"))

  fp <- extract_ml_sampling_params(dat, TRUE)
  g1 <- compute_esf_gamma(fp[1:2])
  g2 <- compute_esf_gamma(fp[3:4])
  gT <- compute_esf_gamma(fp)
  tot <- rowSums(dat)
  row_pick <- res_df[which.max(res_df$observed), ]
  E_bf <- sum(tot == row_pick$t1 + row_pick$t2) *
    g1[row_pick$t1 + 1] * g2[row_pick$t2 + 1] /
    gT[row_pick$t1 + row_pick$t2 + 1]
  expect_equal(row_pick$expected, as.numeric(E_bf), tolerance = 1e-6)
})
