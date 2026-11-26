# Tests for RMitemCatFitPlot() with `dif_var` and its internal helpers
# (R/item_cat_fit_dif.R). The validation study is in dev/catfit_dif/.

# PCM data with two (or more) groups and no DIF. The groups differ in mean
# ability, which must not matter.
make_dif_data <- function(n = 300, k = 4, G = 2, seed = 1L) {
  set.seed(seed)
  thr <- lapply(seq_len(k), function(i) c(-1, 0, 1) + (i - (k + 1) / 2) * 0.5)
  grp <- rep(LETTERS[seq_len(G)], length.out = n)
  theta <- stats::rnorm(n, ifelse(grp == "A", -0.4, 0.4))
  X <- as.data.frame(sim_partial_score(thr, theta))
  colnames(X) <- paste0("I", seq_len(k))
  list(X = X, grp = grp)
}

dif_parts <- function(d) {
  prep <- easyRasch2:::.catfit_prepare(d$X)
  f <- droplevels(as.factor(d$grp)[prep$rows])
  list(prep = prep, pr = easyRasch2:::.catdif_prepare(prep, as.integer(f)))
}

test_that(".catfit_prepare() rows map back to the input rows", {
  d <- make_dif_data()
  X <- d$X
  X[c(3, 10), 1] <- NA
  prep <- easyRasch2:::.catfit_prepare(X)
  expect_false(any(prep$rows %in% c(3, 10)))
  expect_equal(unname(as.matrix(X[prep$rows, ])), unname(prep$X),
               ignore_attr = TRUE)
  expect_equal(unname(rowSums(as.matrix(X[prep$rows, ]))), prep$r)
})

test_that(".pav_w() on per-score aggregates equals .pav() on the data", {
  d <- make_dif_data()
  p <- dif_parts(d)
  pr <- p$pr
  cal <- easyRasch2:::.catdif_calibration(pr)
  fast <- cal$fun(pr$grp)
  ref <- unlist(lapply(split(seq_len(nrow(cal$lay)), cal$lay$curve),
                       function(i) {
    li <- cal$lay[i, ]
    rows <- pr$grp == li$group[1L]
    e <- li$event[1L]
    r <- easyRasch2:::.pav(pr$Pev[pr$s_idx[rows], e], pr$Y[rows, e])
    r$cep[match(li$predicted, r$p)]
  }))
  expect_equal(fast, unname(ref), tolerance = 1e-12)
})

test_that("permutations keep each group's count at every score", {
  p <- dif_parts(make_dif_data())
  pr <- p$pr
  S <- length(pr$sc)
  n0 <- easyRasch2:::.catdif_counts(pr$Y, pr$s_idx, pr$grp, pr$G, S)$n
  set.seed(3)
  for (b in 1:20) {
    gp <- easyRasch2:::.catdif_permute(pr$grp, pr$s_idx)
    expect_identical(
      easyRasch2:::.catdif_counts(pr$Y, pr$s_idx, gp, pr$G, S)$n, n0)
  }
})

test_that("a label constant within scores gives bands of zero width", {
  p <- dif_parts(make_dif_data())
  pr <- p$pr
  pr$grp <- as.integer(pr$r %% 2L) + 1L
  set.seed(1)
  for (lay in list(easyRasch2:::.catdif_calibration(pr),
                   easyRasch2:::.catdif_adjacent(pr))) {
    b <- easyRasch2:::.catdif_bands(pr, lay, B = 49L, conf_level = 0.95)
    expect_equal(b$lower, b$value, tolerance = 1e-10)
    expect_equal(b$upper, b$value, tolerance = 1e-10)
    expect_false(any(b$outside))
  }
})

test_that("group values match the non-DIF helpers", {
  d <- make_dif_data()
  p <- dif_parts(d)
  pr <- p$pr
  P <- p$prep$P
  # Adjacent shares per group equal .catplot_adjacent() on the group's rows.
  adj <- easyRasch2:::.catdif_adjacent(pr)
  obs <- adj$fun(pr$grp)
  for (g in 1:2) {
    rows <- pr$grp == g
    ref <- easyRasch2:::.catplot_adjacent(pr$X[rows, , drop = FALSE],
                                          pr$r[rows], P, 0.05)
    ref$pair_j <- as.integer(sub(" vs.*", "", ref$pair))
    mine <- adj$lay[adj$lay$group == g, ]
    mine$obs <- obs[adj$lay$group == g]
    mm <- merge(ref, mine, by.x = c("item", "pair_j", "score"),
                by.y = c("item", "pair_j", "score"))
    expect_equal(nrow(mm), nrow(ref))
    expect_equal(mm$observed, mm$obs)
  }
  # One-group residuals equal .catplot_residuals().
  pr1 <- pr
  pr1$grp <- rep(1L, length(pr$grp))
  pr1$G <- 1L
  bin <- easyRasch2:::.cicc_class_bins(pr1$r, "quantile", 4L, NULL)
  sg <- easyRasch2:::.score_group_labels(pr1$r, bin)
  res <- easyRasch2:::.catdif_residuals(pr1, as.integer(sg))
  ref <- easyRasch2:::.catplot_residuals(pr1$X, pr1$r, P, sg)
  ref <- ref[!is.na(ref$z), ]
  mine <- res$lay
  mine$z <- res$fun(pr1$grp)
  mine$score_group <- levels(sg)[mine$q]
  mm <- merge(ref, mine, by = c("item", "score_group", "category"))
  expect_equal(nrow(mm), nrow(ref))
  expect_equal(mm$z.x, mm$z.y)
})

test_that("RMitemCatFitPlot() with dif_var returns the documented columns", {
  skip_if_not_installed("ggplot2")
  d <- make_dif_data()
  common <- c("item", "group", "perm_mean", "lower", "upper", "outside",
              "iterations_used")
  cols <- list(
    calibration = c("event", "predicted", "cep", "n"),
    adjacent = c("pair", "score", "n_pair", "observed", "expected"),
    residuals = c("score_group", "category", "observed", "expected", "z")
  )
  # Legend n: respondents used (complete, non-extreme total score).
  n_used <- table(d$grp[easyRasch2:::.catfit_prepare(d$X)$rows])
  for (ty in names(cols)) {
    tab <- RMitemCatFitPlot(d$X, type = ty, dif_var = d$grp,
                            iterations = 50, seed = 1, output = "dataframe")
    expect_true(all(c(common, cols[[ty]]) %in% names(tab)), info = ty)
    expect_equal(levels(tab$group),
                 paste0(c("A", "B"), " (n = ", as.vector(n_used), ")"),
                 info = ty)
    expect_true(all(tab$lower <= tab$upper + 1e-12), info = ty)
    expect_true(all(tab$iterations_used == 50L), info = ty)
    p <- RMitemCatFitPlot(d$X, type = ty, dif_var = d$grp, iterations = 50,
                          seed = 1)
    expect_s3_class(p, "ggplot")
  }
})

test_that("dif_var output is reproducible with a seed and defaults to 1000", {
  d <- make_dif_data(n = 200)
  a <- RMitemCatFitPlot(d$X, type = "calibration", dif_var = d$grp, seed = 4,
                        output = "dataframe")
  b <- RMitemCatFitPlot(d$X, type = "calibration", dif_var = d$grp, seed = 4,
                        output = "dataframe")
  expect_identical(a, b)
  expect_true(all(a$iterations_used == 1000L))
})

test_that("dif_var checks and messages", {
  d <- make_dif_data(n = 200)
  expect_error(RMitemCatFitPlot(d$X, type = "restscore", dif_var = d$grp),
               "not available with `dif_var`")
  fit_like <- structure(list(), class = "RMitemCatFit")
  expect_error(RMitemCatFitPlot(d$X, fit = fit_like, dif_var = d$grp),
               "cannot be combined")
  expect_error(RMitemCatFitPlot(d$X, dif_var = rep("A", 200)),
               "at least 2 groups")
  expect_error(RMitemCatFitPlot(d$X, dif_var = rep(1:5, 40)),
               "at most 4")
  expect_error(RMitemCatFitPlot(d$X, dif_var = d$grp[-1]), "same length")
  g <- d$grp
  g[1:5] <- NA
  expect_message(
    tab <- RMitemCatFitPlot(d$X, type = "adjacent", dif_var = g,
                            iterations = 20, seed = 1, output = "dataframe"),
    "5 row\\(s\\) with NA in `dif_var` dropped"
  )
  n_used <- sum(as.integer(sub(".*n = ([0-9]+)\\)", "\\1",
                               levels(tab$group))))
  expect_lte(n_used, 195L)
})

test_that("three groups work and the non-DIF output is unchanged", {
  skip_if_not_installed("ggplot2")
  d <- make_dif_data(n = 300, G = 3)
  tab <- RMitemCatFitPlot(d$X, type = "residuals", dif_var = d$grp,
                          iterations = 30, seed = 1, output = "dataframe")
  expect_equal(nlevels(tab$group), 3L)
  plain <- RMitemCatFitPlot(d$X, type = "adjacent", output = "dataframe")
  expect_named(plain, c("item", "pair", "score", "n_pair", "observed",
                        "expected", "lower", "upper"))
})
