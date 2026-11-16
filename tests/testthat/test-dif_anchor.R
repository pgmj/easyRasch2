# Tests for RMdifAnchor()
#
# Covers:
#   - Input validation (groups, reference, n_seeds, seed, alpha, conf)
#   - Agreement with psychotools::anchortest() for the same anchor
#   - Spacing tests do not depend on the anchor
#   - Sign convention and `reference`
#   - Reproducibility, and the caller's RNG left unchanged
#   - Gini anchor: one untested item
#   - Detection of a planted DIF item
#   - Output types and plots

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
sim_pcm_anchor <- function(theta, tau) {
  sapply(tau, function(t) {
    m <- length(t)
    eta <- outer(theta, 0:m) -
      matrix(c(0, cumsum(t)), length(theta), m + 1L, byrow = TRUE)
    eta <- eta - apply(eta, 1L, max)
    p <- exp(eta)
    p <- p / rowSums(p)
    u <- stats::runif(length(theta))
    rowSums(u > t(apply(p, 1L, cumsum))[, -(m + 1L), drop = FALSE])
  })
}

# 8 items with 3 categories, n per group, DIF of `size` logits in item 3.
make_dif_data <- function(n = 400, size = 0.8, seed = 2L) {
  set.seed(seed)
  tau <- lapply(seq(-1.2, 1.2, length.out = 8), function(l) l + c(-0.7, 0.7))
  tau2 <- tau
  tau2[[3]] <- tau2[[3]] + size
  x <- rbind(sim_pcm_anchor(stats::rnorm(n), tau),
             sim_pcm_anchor(stats::rnorm(n, 0.3), tau2))
  x <- as.data.frame(x)
  names(x) <- paste0("I", seq_len(ncol(x)))
  list(data = x, group = rep(c("A", "B"), each = n))
}

dd <- make_dif_data()

# The forward anchor may message about items near .5 inclusion; quiet here.
rda <- function(...) suppressMessages(RMdifAnchor(...))

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMdifAnchor requires exactly two groups", {
  g3 <- rep(c("A", "B", "C"), length.out = nrow(dd$data))
  expect_error(rda(dd$data, g3), "exactly 2")
  expect_error(rda(dd$data, rep("A", nrow(dd$data))), "exactly 2")
})

test_that("RMdifAnchor validates its arguments", {
  expect_error(rda(dd$data, dd$group[-1]), "same length")
  expect_error(rda(dd$data, dd$group, reference = "Z"), "reference")
  expect_error(rda(dd$data, dd$group, n_seeds = 0), "n_seeds")
  expect_error(rda(dd$data, dd$group, n_seeds = 2.5), "n_seeds")
  expect_error(rda(dd$data, dd$group, seed = "a"), "seed")
  expect_error(rda(dd$data, dd$group, seed = 3e9), "between")
  expect_error(rda(dd$data, dd$group, alpha = 1), "alpha")
  expect_error(rda(dd$data, dd$group, conf = 0), "conf")
  expect_error(rda(dd$data, dd$group, sort = NA), "sort")
})

test_that("RMdifAnchor drops rows with NA in dif_var with a message", {
  g <- dd$group
  g[1:5] <- NA
  expect_message(
    out <- RMdifAnchor(dd$data, g, output = "list", n_seeds = 10),
    "5 row"
  )
  expect_equal(out$settings$N, nrow(dd$data) - 5L)
  expect_equal(out$settings$N_excluded, 5L)
})

test_that("spacing needs items with two or more thresholds", {
  x <- dd$data
  x[] <- lapply(x, function(v) as.integer(v > 0))
  expect_error(
    rda(x, dd$group, level = "spacing", n_seeds = 10),
    "two or more thresholds"
  )
})

# ---------------------------------------------------------------------
# Method
# ---------------------------------------------------------------------
test_that("anchored location tests match psychotools::anchortest()", {
  skip_if_not_installed("multcomp") # anchortest() returns a multcomp object
  sub <- split(dd$data, dd$group)
  f1 <- psychotools::pcmodel(sub$A)
  f2 <- psychotools::pcmodel(sub$B)
  set.seed(1)
  at <- psychotools::anchortest(f1, f2, class = "forward", select = "MTT")
  t1 <- psychotools::threshpar(f1, vcov = TRUE)
  t2 <- psychotools::threshpar(f2, vcov = TRUE)
  m <- lengths(t1)
  res <- easyRasch2:::.dif_anchor_tests(
    unlist(t2) - unlist(t1),
    as.matrix(attr(t1, "vcov")) + as.matrix(attr(t2, "vcov")),
    rep(seq_along(m), m), unlist(lapply(m, seq_len)),
    at$anchor_items, length(m), "item", "none"
  )
  tested <- setdiff(seq_along(m), at$anchor_items[1L])
  # anchortest reports group 1 minus group 2, RMdifAnchor the reverse.
  expect_equal(res$z[tested], -as.numeric(at$final_tests$test$tstat),
               tolerance = 1e-8)
})

test_that("spacing tests do not depend on the anchor", {
  a <- rda(dd$data, dd$group, level = "spacing",
                   output = "dataframe", n_seeds = 10)
  b <- rda(dd$data, dd$group, level = "spacing", anchor = "gini",
                   output = "dataframe")
  expect_equal(a$Chisq, b$Chisq)
  expect_equal(a$p, b$p)
})

test_that("reference flips the sign and leaves p-values unchanged", {
  a <- rda(dd$data, dd$group, output = "dataframe", n_seeds = 20)
  b <- rda(dd$data, dd$group, output = "dataframe", n_seeds = 20,
                   reference = "B")
  expect_equal(a$DIF, -b$DIF)
  expect_equal(a$padj, b$padj)
})

test_that("results are reproducible and leave the caller's RNG unchanged", {
  set.seed(42)
  before <- stats::runif(1)
  set.seed(42)
  a <- rda(dd$data, dd$group, output = "dataframe", n_seeds = 20)
  after <- stats::runif(1)
  expect_identical(before, after)
  b <- rda(dd$data, dd$group, output = "dataframe", n_seeds = 20)
  expect_identical(a, b)
  expect_no_error(
    rda(dd$data, dd$group, output = "dataframe", n_seeds = 5,
                seed = NULL)
  )
})

test_that("the seed gives the same result under any RNGkind, which is restored", {
  a <- rda(dd$data, dd$group, output = "dataframe", n_seeds = 20)
  old <- RNGkind()
  on.exit(RNGkind(old[1L], old[2L], old[3L]), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(7)
  state <- .Random.seed
  b <- rda(dd$data, dd$group, output = "dataframe", n_seeds = 20)
  expect_identical(RNGkind()[1L], "L'Ecuyer-CMRG")
  expect_identical(.Random.seed, state)
  expect_identical(a, b)
})

test_that("the forward anchor reports inclusion frequencies", {
  out <- rda(dd$data, dd$group, output = "dataframe", n_seeds = 20)
  expect_true(all(out$Inclusion >= 0 & out$Inclusion <= 1))
  expect_true(all(out$Anchor == (out$Inclusion >= 0.5)))
  # With two or more anchor items, every item is tested.
  expect_false(anyNA(out$padj))
})

test_that("the Gini anchor is one untested item", {
  out <- rda(dd$data, dd$group, anchor = "gini",
                     output = "dataframe")
  expect_equal(sum(out$Anchor), 1L)
  expect_true(is.na(out$padj[out$Anchor]))
  expect_true(is.na(out$DIF[out$Anchor]))
  expect_false("Inclusion" %in% names(out))
  st <- rda(dd$data, dd$group, anchor = "gini", output = "list")$settings
  expect_equal(st$Untested, out$Item[out$Anchor])
})

test_that("a planted DIF item is flagged and kept out of the anchor", {
  out <- rda(dd$data, dd$group, output = "dataframe", n_seeds = 20)
  expect_true(out$Flagged[out$Item == "I3"])
  expect_false(out$Anchor[out$Item == "I3"])
  expect_gt(out$DIF[out$Item == "I3"], 0)
})

test_that("threshold level has one row per threshold", {
  out <- rda(dd$data, dd$group, level = "threshold",
                     output = "dataframe", n_seeds = 10)
  expect_equal(nrow(out), 16L)
  expect_setequal(unique(out$Threshold), c("c1", "c2"))
})

test_that("threshold rows plot in numeric order past c9 (DA2)", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  tau <- c(list(seq(-1.5, 1.5, length.out = 10)),
           lapply(seq(-1, 1, length.out = 5), function(l) l + c(-0.5, 0.5)))
  x <- as.data.frame(rbind(sim_pcm_anchor(stats::rnorm(1500), tau),
                           sim_pcm_anchor(stats::rnorm(1500), tau)))
  names(x) <- paste0("I", seq_len(ncol(x)))
  g <- rep(c("A", "B"), each = 1500)
  p <- rda(x, g, level = "threshold", output = "ggplot", n_seeds = 5)
  expect_identical(rev(levels(p$data$Row))[c(1, 2, 10)], c("c1", "c2", "c10"))
})

test_that("corrections order as expected", {
  none <- rda(dd$data, dd$group, correction = "none",
                      output = "dataframe", n_seeds = 10)
  holm <- rda(dd$data, dd$group, correction = "holm",
                      output = "dataframe", n_seeds = 10)
  fwer <- rda(dd$data, dd$group, correction = "fwer",
                      output = "dataframe", n_seeds = 10)
  expect_equal(none$padj, none$p)
  expect_true(all(holm$padj >= holm$p - 1e-12))
  expect_true(all(fwer$padj >= fwer$p - 1e-3))
})

# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------
test_that("list output holds the table, settings, distances and criterion", {
  out <- rda(dd$data, dd$group, output = "list", n_seeds = 10)
  expect_named(out, c("table", "settings", "distances", "criterion"))
  expect_identical(out$table,
                   rda(dd$data, dd$group, output = "dataframe", n_seeds = 10))
  st <- out$settings
  expect_equal(nrow(st), 1L)
  expect_equal(c(st$n_A, st$n_B), c(400L, 400L))
  expect_identical(st$Comparison, "B minus A")
  expect_identical(st$Anchor_items,
                   paste(out$table$Item[out$table$Anchor], collapse = ", "))
  expect_true(is.na(st$Untested))
  expect_equal(st$N_seeds, 10L)
  expect_null(out$distances)
  expect_equal(nrow(out$criterion), ncol(dd$data))
  expect_equal(out$criterion$Location_diff, out$table$DIF, tolerance = 1e-10)
  sp <- rda(dd$data, dd$group, level = "spacing", output = "list")
  expect_equal(nrow(sp$distances), 8L)
  expect_named(sp$distances, c("Item", "Distance", "Difference", "SE",
                               "lower", "upper"))
  # nothing in attributes
  expect_setequal(names(attributes(out$table)), c("names", "row.names", "class"))
})
test_that("spacing runs no anchor search and reports none (DA11)", {
  sp <- rda(dd$data, dd$group, level = "spacing", output = "list")
  expect_identical(sp$settings$Anchor, "none (spacing needs no anchor)")
  expect_true(is.na(sp$settings$Anchor_items))
  expect_true(is.na(sp$settings$N_seeds))
  expect_null(sp$criterion)
  # no tie-break message, since no anchor is searched for
  expect_no_message(RMdifAnchor(dd$data, dd$group, level = "spacing",
                                output = "dataframe"))
})

test_that("single-step padj is never below the raw p (DA12)", {
  # Two items with the Gini anchor leave one tested unit.
  out <- rda(dd$data[, c("I2", "I3")], dd$group, anchor = "gini",
             output = "dataframe")
  expect_true(all(out$padj >= out$p, na.rm = TRUE))
  sp <- rda(dd$data, dd$group, level = "spacing", output = "dataframe")
  expect_true(all(sp$padj >= sp$p))
})

test_that("missing cells are blank whether or not a row is flagged (DA14)", {
  skip_if_not_installed("knitr")
  x <- make_dif_data(size = 0)
  k <- rda(x$data, x$group, anchor = "gini", correction = "fwer")
  expect_false(any(grepl("\\bNA\\b", k)))
  k2 <- rda(dd$data, dd$group, anchor = "gini")
  expect_false(any(grepl("\\bNA\\b", k2)))
})

test_that("kable output has a caption naming the anchor", {
  skip_if_not_installed("knitr")
  k <- rda(dd$data, dd$group, n_seeds = 10)
  expect_s3_class(k, "knitr_kable")
  expect_true(any(grepl("iterative forward", k)))
  k2 <- rda(dd$data, dd$group, level = "spacing", sort = TRUE)
  expect_true(any(grepl("No anchor is needed", k2)))
})

test_that("sort orders by the size of the difference", {
  out <- rda(dd$data, dd$group, output = "dataframe", n_seeds = 10,
                     sort = TRUE)
  expect_equal(order(-abs(out$DIF)), seq_len(nrow(out)))
  expect_setequal(names(attributes(out)), c("names", "row.names", "class"))
})

test_that("ggplot outputs build for every level and the criterion plot", {
  skip_if_not_installed("ggplot2")
  for (lv in c("item", "threshold", "spacing")) {
    p <- rda(dd$data, dd$group, level = lv, output = "ggplot",
                     n_seeds = 10)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  p <- rda(dd$data, dd$group, output = "ggplot",
                   plot_type = "criterion", anchor = "gini")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})
