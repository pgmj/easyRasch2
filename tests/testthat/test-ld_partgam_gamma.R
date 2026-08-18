# Tests for .partgam_ld_gamma(), the vectorised partial gamma used for the
# simulated null in RMlocdepGammaCutoff() and for the observed statistic tested
# by RMlocdepGamma().
#
# Exact agreement with iarm::partgam_LD() is load-bearing: the bootstrap null
# and the observed statistic both come from this path, and the `gamma` column
# users see still comes from iarm, so any data shape where the two diverge
# would silently bias the p-values and contradict the displayed coefficient.

make_poly <- function(n = 300, k = 6, cats = 4L, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(
    matrix(sample(0:(cats - 1L), n * k, replace = TRUE), n, k)
  )
  colnames(df) <- paste0("Item", seq_len(k))
  df
}

# Gammas from iarm for one rest-score direction, keyed by unordered pair.
iarm_gamma <- function(df, direction) {
  out <- NULL
  sink(nullfile())
  on.exit(sink(), add = TRUE)
  raw <- iarm::partgam_LD(as.data.frame(df))
  out <- raw[[direction]]
  stats::setNames(
    as.numeric(out$gamma),
    paste(
      pmin(as.character(out$Item1), as.character(out$Item2)),
      pmax(as.character(out$Item1), as.character(out$Item2)),
      sep = "___"
    )
  )
}

fast_gamma <- function(df, direction) {
  out <- .partgam_ld_gamma(df, direction = direction)
  stats::setNames(
    out$gamma,
    paste(pmin(out$Item1, out$Item2), pmax(out$Item1, out$Item2), sep = "___")
  )
}

expect_matches_iarm <- function(df, label) {
  for (d in 1:2) {
    theirs <- iarm_gamma(df, d)
    ours <- fast_gamma(df, d)
    expect_setequal(names(ours), names(theirs))
    expect_equal(
      ours[names(theirs)],
      theirs,
      tolerance = 0,
      info = paste(label, "direction", d)
    )
  }
}

# ---------------------------------------------------------------------
# Agreement with iarm across data shapes
# ---------------------------------------------------------------------
test_that(".partgam_ld_gamma matches iarm on polytomous data", {
  skip_if_not_installed("iarm")
  expect_matches_iarm(make_poly(n = 300), "polytomous 4 categories")
})

test_that(".partgam_ld_gamma matches iarm on dichotomous data", {
  skip_if_not_installed("iarm")
  expect_matches_iarm(make_poly(n = 300, cats = 2L), "dichotomous")
})

test_that(".partgam_ld_gamma matches iarm when items differ in categories", {
  skip_if_not_installed("iarm")
  df <- make_poly(n = 300)
  df[[1]] <- pmin(df[[1]], 1L)
  df[[2]] <- pmin(df[[2]], 2L)
  expect_matches_iarm(df, "mixed category counts")
})

test_that(".partgam_ld_gamma matches iarm when a middle category is unused", {
  skip_if_not_installed("iarm")
  df <- make_poly(n = 300)
  df[[1]][df[[1]] == 2L] <- 3L
  expect_matches_iarm(df, "unobserved middle category")
})

test_that(".partgam_ld_gamma matches iarm when strata are sparse", {
  skip_if_not_installed("iarm")
  expect_matches_iarm(make_poly(n = 40, seed = 7L), "sparse strata")
})

test_that(".partgam_ld_gamma matches iarm with more categories", {
  skip_if_not_installed("iarm")
  expect_matches_iarm(make_poly(n = 400, cats = 6L), "six categories")
})

test_that(".partgam_ld_gamma matches iarm on the minimum item count", {
  skip_if_not_installed("iarm")
  expect_matches_iarm(make_poly(n = 300, k = 3L), "three items")
})

# ---------------------------------------------------------------------
# Row order and shape
# ---------------------------------------------------------------------
test_that(".partgam_ld_gamma returns iarm's pair order", {
  skip_if_not_installed("iarm")
  df <- make_poly(n = 300, k = 5L)
  sink(nullfile())
  raw <- iarm::partgam_LD(as.data.frame(df))
  sink()
  for (d in 1:2) {
    ours <- .partgam_ld_gamma(df, direction = d)
    expect_identical(ours$Item1, as.character(raw[[d]]$Item1))
    expect_identical(ours$Item2, as.character(raw[[d]]$Item2))
  }
})

test_that(".partgam_ld_gamma returns one row per unordered pair", {
  df <- make_poly(n = 200, k = 7L)
  expect_identical(nrow(.partgam_ld_gamma(df)), 21L)
  expect_identical(nrow(.partgam_ld_gamma(df, direction = 2L)), 21L)
})

# ---------------------------------------------------------------------
# Degenerate data, where iarm cannot provide a reference
# ---------------------------------------------------------------------
# iarm::partgam_LD() errors for the whole data set when an item is constant,
# losing every pair rather than the affected ones, so agreement is checked
# against a brute-force count taken straight from the definition instead.
brute_gamma <- function(x, y, z) {
  conc <- 0
  disc <- 0
  for (s in unique(z)) {
    xi <- x[z == s]
    yi <- y[z == s]
    if (length(xi) < 2L) next
    sgn <- sign(outer(xi, xi, "-")) * sign(outer(yi, yi, "-"))
    conc <- conc + sum(sgn > 0) / 2
    disc <- disc + sum(sgn < 0) / 2
  }
  if (conc + disc == 0) NA_real_ else (conc - disc) / (conc + disc)
}

brute_all <- function(df, direction = 1L) {
  X <- as.matrix(df)
  k <- ncol(X)
  score <- rowSums(X)
  pairs <- do.call(
    rbind,
    lapply(seq_len(k), function(i) {
      js <- if (direction == 1L) {
        seq_len(k)[seq_len(k) > i]
      } else {
        seq_len(k)[seq_len(k) < i]
      }
      if (length(js) == 0L) NULL else cbind(i, js)
    })
  )
  vapply(
    seq_len(nrow(pairs)),
    function(p) {
      i <- pairs[p, 1L]
      j <- pairs[p, 2L]
      brute_gamma(X[, i], X[, j], score - X[, j])
    },
    numeric(1L)
  )
}

test_that("the brute-force reference reproduces .partgam_ld_gamma", {
  # Validates the reference itself on data where iarm also agrees, so the
  # degenerate case below rests on a checked comparator.
  df <- make_poly(n = 120, k = 5L)
  for (d in 1:2) {
    expect_equal(.partgam_ld_gamma(df, d)$gamma, brute_all(df, d), tolerance = 0)
  }
})

test_that(".partgam_ld_gamma handles a constant item that breaks iarm", {
  df <- make_poly(n = 120, k = 5L)
  df[[3]] <- 0L

  skip_if_not_installed("iarm")
  sink(nullfile())
  broke <- inherits(try(iarm::partgam_LD(as.data.frame(df)), silent = TRUE), "try-error")
  sink()
  expect_true(broke)

  for (d in 1:2) {
    out <- .partgam_ld_gamma(df, direction = d)
    expect_equal(out$gamma, brute_all(df, d), tolerance = 0)
    # Only the pairs involving the constant item are undefined.
    involved <- out$Item1 == "Item3" | out$Item2 == "Item3"
    expect_true(all(is.na(out$gamma[involved])))
    expect_true(all(is.finite(out$gamma[!involved])))
  }
})

# ---------------------------------------------------------------------
# Stratum sign homogeneity
# ---------------------------------------------------------------------
# Kreiner's condition: pooling C and D over strata is interpretable as a
# partial correlation only when the stratum associations point the same way.

# Two strata built by hand, so the expected answer is known exactly.
two_strata <- function(second) {
  x <- c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L)
  y <- c(0L, 0L, 1L, 1L, second)
  z <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  G <- outer(1:2, 1:2, function(a, b) as.numeric(b > a))
  .partgam_one(x, y, z, 2L, G, t(G), strata = TRUE)
}

test_that("stratum gammas and weights reconstruct the pooled coefficient", {
  df <- make_poly(n = 300, k = 5L)
  out <- .partgam_ld_gamma(df, strata = TRUE)
  for (p in seq_len(nrow(out))) {
    st <- attr(out, "strata")[[p]]
    ok <- !is.na(st$gamma_k)
    # gamma_k * weight_k is C_k - D_k, so the weighted mean is the pooled value.
    expect_equal(
      sum(st$gamma_k[ok] * st$weight_k[ok]) / sum(st$weight_k[ok]),
      out$gamma[p],
      tolerance = 1e-12
    )
  }
})

test_that("opposite-signed strata are reported as heterogeneous", {
  st <- two_strata(c(1L, 1L, 0L, 0L)) # second stratum perfectly discordant
  s <- .partgam_strata_summary(st)
  expect_identical(s$n_pos, 1L)
  expect_identical(s$n_neg, 1L)
  expect_false(s$homogeneous)
  expect_equal(st$gamma_k, c(1, -1))
  # The two strata carry equal weight and cancel exactly.
  expect_equal(st$gamma, 0)
  expect_true(is.na(s$w_opposing)) # undefined when the pooled value is zero
})

test_that("same-signed strata are reported as homogeneous", {
  st <- two_strata(c(0L, 0L, 1L, 1L)) # second stratum perfectly concordant
  s <- .partgam_strata_summary(st)
  expect_identical(s$n_pos, 2L)
  expect_identical(s$n_neg, 0L)
  expect_true(s$homogeneous)
  expect_equal(st$gamma, 1)
})

test_that("w_opposing is the weight share of the sign-opposing strata", {
  # Stratum 1 concordant with 4 comparable pairs, stratum 2 discordant with 1.
  x <- c(0L, 0L, 1L, 1L, 0L, 1L)
  y <- c(0L, 0L, 1L, 1L, 1L, 0L)
  z <- c(1L, 1L, 1L, 1L, 2L, 2L)
  G <- outer(1:2, 1:2, function(a, b) as.numeric(b > a))
  st <- .partgam_one(x, y, z, 2L, G, t(G), strata = TRUE)
  s <- .partgam_strata_summary(st)
  expect_equal(st$weight_k, c(4, 1))
  expect_false(s$homogeneous)
  expect_equal(s$w_opposing, 1 / 5)
})

test_that("strata with no comparable pair are excluded, not counted as zero", {
  # The second stratum holds observations but one item is constant in it.
  st <- two_strata(c(1L, 1L, 1L, 1L))
  s <- .partgam_strata_summary(st)
  expect_true(is.na(st$gamma_k[2]))
  expect_identical(st$n_k, c(4L, 4L))
  expect_identical(s$n_strata, 1L)
  expect_true(s$homogeneous)
})

test_that("strata = TRUE leaves the coefficient unchanged", {
  df <- make_poly(n = 250, k = 6L)
  for (d in 1:2) {
    expect_equal(
      .partgam_ld_gamma(df, direction = d, strata = TRUE)$gamma,
      .partgam_ld_gamma(df, direction = d)$gamma,
      tolerance = 0
    )
  }
})

test_that("strata = FALSE returns no homogeneity columns", {
  df <- make_poly(n = 200, k = 4L)
  expect_identical(names(.partgam_ld_gamma(df)), c("Item1", "Item2", "gamma"))
  expect_null(attr(.partgam_ld_gamma(df), "strata"))
})

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that(".partgam_ld_gamma rejects missing values", {
  df <- make_poly(n = 100, k = 4L)
  df[1, 1] <- NA
  expect_error(.partgam_ld_gamma(df), regexp = "missing values")
})

test_that(".partgam_ld_gamma rejects a single item", {
  df <- make_poly(n = 100, k = 4L)[, 1, drop = FALSE]
  expect_error(.partgam_ld_gamma(df), regexp = "at least two items")
})

test_that(".partgam_ld_gamma agrees with iarm after complete-case filtering", {
  skip_if_not_installed("iarm")
  df <- make_poly(n = 300, k = 5L)
  df[sample(nrow(df), 30), 2] <- NA
  # iarm filters internally, so the fast path must be given the same subset.
  complete <- df[stats::complete.cases(df), , drop = FALSE]
  for (d in 1:2) {
    expect_equal(
      fast_gamma(complete, d)[names(iarm_gamma(df, d))],
      iarm_gamma(df, d),
      tolerance = 0
    )
  }
})
