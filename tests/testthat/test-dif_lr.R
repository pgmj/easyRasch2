# Tests for RMdifLR() (Andersen's LR test only, since the redesign of
# 2026-11; per-item and per-threshold DIF moved to RMdifAnchor()).
#
# Covers:
#   - Input validation (data, dif_var, model) and removed arguments
#   - Output: one-row data.frame and kable, complete in its columns
#   - The LR test equals eRm::LRtest()
#   - NA handling, category checks, failed fits
#   - The once-per-session pointer to RMdifAnchor()

# RMdifLR() informs once per session; quiet in tests that do not check it.
rlr <- function(...) suppressMessages(RMdifLR(...))

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
make_dichotomous <- function(n = 80, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(
    sample(0:1, n * k, replace = TRUE),
    nrow = n,
    ncol = k
  ))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_polytomous <- function(n = 80, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(
    sample(0:2, n * k, replace = TRUE),
    nrow = n,
    ncol = k
  ))
  colnames(df) <- paste0("I", seq_len(k))
  df
}


# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMdifLR errors when data has non-zero minimum", {
  df <- make_dichotomous() + 1L # min = 1
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  expect_error(rlr(df, dif_var = grp), regexp = "scored starting at 0")
})

test_that("RMdifLR errors when dif_var length does not match nrow(data)", {
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df) - 1L))
  expect_error(rlr(df, dif_var = grp), regexp = "same length as nrow")
})

test_that("RMdifLR errors when dif_var has fewer than 2 distinct levels", {
  df <- make_dichotomous()
  grp <- factor(rep("A", nrow(df)))
  expect_error(rlr(df, dif_var = grp), regexp = "at least 2 distinct")
})

test_that("RMdifLR errors when model = 'RM' on polytomous data", {
  df <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  expect_error(
    rlr(df, dif_var = grp, model = "RM"),
    regexp = "dichotomous"
  )
})



test_that("removed arguments and plots give an error pointing to RMdifAnchor()", {
  df <- make_polytomous()
  grp <- rep(c("A", "B"), length.out = nrow(df))
  expect_error(rlr(df, grp, level = "threshold"), "RMdifAnchor\\(\\)")
  expect_error(rlr(df, grp, p_value = TRUE, alpha = 0.01), "`p_value`, `alpha`")
  expect_error(rlr(df, grp, output = "ggplot"), "no longer draws plots")
  expect_error(rlr(df, grp, foo = 1), "Unused argument")
  # positional calls written for 1.4.0, where `level` came fourth (DA13)
  expect_error(rlr(df, grp, "auto", "threshold"), "RMdifAnchor\\(\\)")
  expect_error(rlr(df, grp, "PCM", "item", "dataframe"), "removed `level`")
})

test_that("RMdifLR points to RMdifAnchor() once per session", {
  df <- make_polytomous()
  grp <- rep(c("A", "B"), length.out = nrow(df))
  rlang::local_options(rlib_message_verbosity = "default")
  rlang::reset_message_verbosity("easyRasch2_difLR_anchor")
  expect_message(RMdifLR(df, grp, output = "dataframe"), "RMdifAnchor")
  expect_no_message(RMdifLR(df, grp, output = "dataframe"))
})

# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------
test_that("dataframe output is one row with the LR test", {
  df <- make_dichotomous()
  grp <- rep(c("A", "B"), length.out = nrow(df))
  out <- rlr(df, grp, output = "dataframe")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_named(out, c("Chisq", "df", "p", "Model", "Items", "Groups", "N",
                      "N_excluded", "N_incomplete", "n_A", "n_B"))
  expect_identical(out$Model, "RM")
  expect_equal(c(out$n_A, out$n_B), c(40L, 40L))
  expect_equal(out$N_excluded, 0L)
  expect_equal(out$N_incomplete, 0L)
  # nothing is attached as attributes
  expect_setequal(names(attributes(out)), c("names", "row.names", "class"))
  expect_identical(rlr(make_polytomous(), grp, output = "dataframe")$Model, "PCM")
})

test_that("kable is the default and names the groups and test", {
  skip_if_not_installed("knitr")
  df <- make_polytomous()
  grp <- rep(c("A", "B"), length.out = nrow(df))
  k <- rlr(df, grp)
  expect_s3_class(k, "knitr_kable")
  txt <- paste(k, collapse = "\n")
  expect_match(txt, "Andersen's likelihood-ratio test")
  expect_match(txt, "A: 40, B: 40")
  expect_no_match(txt, "e-0")
  # a tiny p-value stays plain
  big <- make_polytomous(n = 400)
  g2 <- rep(c("A", "B"), length.out = nrow(big))
  set.seed(7)
  up <- g2 == "B" & stats::runif(nrow(big)) < 0.9
  big$I1[up] <- pmin(big$I1[up] + 1L, 2L)
  k2 <- paste(rlr(big, g2), collapse = "\n")
  expect_match(k2, "\\|(< 0\\.0001|0\\.[0-9]{4}) *\\|")
  expect_no_match(k2, "e-0")
})

test_that("table p-values have 4 decimals and a lower bound (DA15)", {
  expect_identical(.format_p_table(c(0.04321, 0.00004, 0, 1, NA)),
                   c("0.0432", "< 0.0001", "< 0.0001", "1.0000", ""))
})


# NA handling in dif_var
# ---------------------------------------------------------------------
test_that("RMdifLR drops NA rows in dif_var with a message", {
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  grp[c(1L, 5L, 9L)] <- NA

  expect_message(
    out <- RMdifLR(df, dif_var = grp, output = "dataframe"),
    regexp = "3 row\\(s\\) with NA"
  )
  expect_equal(out$N, nrow(df) - 3L)
  expect_equal(out$N_excluded, 3L)
})

test_that(".format_p follows APA-style thresholds", {
  expect_equal(.format_p(0.0000497), "p < 0.001")
  expect_equal(.format_p(0.001), "p = 0.001")
  expect_equal(.format_p(0.004975), "p = 0.005")
  expect_equal(.format_p(0.71638), "p = 0.716")
  expect_equal(.format_p(NA_real_), "p = NA")
})

test_that("RMdifLR drops all-NA respondents instead of erroring", {
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:2, 80 * 5, replace = TRUE), nrow = 80))
  colnames(df) <- paste0("i", 1:5)
  df[3, ] <- NA
  grp <- factor(rep(c("a", "b"), 40))

  expect_message(k <- RMdifLR(df, dif_var = grp), "no responses dropped")
  expect_s3_class(k, "knitr_kable")
})

# ---------------------------------------------------------------------
# psychotools engine (replaced eRm::LRtest in 1.4.0.9002)
# ---------------------------------------------------------------------
test_that("the LR test equals eRm::LRtest(), with and without missing data", {
  skip_if_not_installed("eRm")
  e <- new.env()
  utils::data(pcmdat2, raschdat1, package = "eRm", envir = e)
  set.seed(4)
  pcm <- as.data.frame(e$pcmdat2)
  pcm_na <- pcm
  pcm_na[matrix(stats::runif(length(as.matrix(pcm))) < 0.05, nrow(pcm))] <- NA
  ras <- as.data.frame(e$raschdat1)[, 1:12]
  cases <- list(
    list(pcm, sample(c("a", "b"), nrow(pcm), TRUE), "PCM"),
    # I3 left out: only 7 of its 14 zeros come from respondents without a
    # zero score, so a 3-way split often leaves a group with none there,
    # which now stops (eRm returned I3.c1 = -13.5, SE 276 in one such split)
    list(pcm[, -3], sample(c("x", "y", "z"), nrow(pcm), TRUE), "PCM"),
    list(pcm_na, sample(c("a", "b"), nrow(pcm), TRUE), "PCM"),
    list(ras, sample(1:2, nrow(ras), TRUE), "RM")
  )
  for (cs in cases) {
    lr <- rlr(cs[[1]], cs[[2]], output = "dataframe")
    fit <- if (cs[[3]] == "PCM") eRm::PCM(cs[[1]]) else eRm::RM(cs[[1]])
    ref <- eRm::LRtest(fit, splitcr = factor(cs[[2]]))
    expect_equal(lr$Chisq, unname(ref$LR), tolerance = 1e-4)
    expect_identical(as.integer(lr$df), as.integer(ref$df))
  }
})

test_that("a category unused in one group stops, naming item, category and group", {
  df <- make_polytomous(n = 120)
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  # top category of I2 unused in group B (allowed in the full data, not
  # across groups), and the lowest category of I4 unused in group A
  df$I2[grp == "B" & df$I2 == 2L] <- 1L
  df$I4[grp == "A" & df$I4 == 0L] <- 1L
  err <- tryCatch(rlr(df, grp, output = "dataframe"), error = function(e) e)
  expect_s3_class(err, "easyRasch2_category_error")
  expect_match(conditionMessage(err), "I2 (category 2 unused in group 'B')", fixed = TRUE)
  expect_match(conditionMessage(err), "I4 (category 0 unused in group 'A')", fixed = TRUE)
  expect_identical(err$items, c("I2", "I4"))
  expect_match(conditionMessage(err),
               "Item parameters cannot be compared across DIF groups", fixed = TRUE)
  # an item with no responses in a group
  df2 <- make_polytomous(n = 120)
  df2$I3[grp == "A"] <- NA
  expect_error(rlr(df2, grp, output = "dataframe"), "I3 (no responses in group 'A')", fixed = TRUE)
})

test_that("a category used in a group only by extreme scorers stops (DL1)", {
  df <- make_dichotomous(n = 120, k = 5)
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  # I1 = 1 in group B only for two perfect scorers
  df$I1[grp == "B"] <- 0L
  b <- which(grp == "B")[1:2]
  df[b, ] <- 1L
  err <- tryCatch(rlr(df, grp, output = "dataframe"), error = function(e) e)
  expect_s3_class(err, "easyRasch2_category_error")
  expect_match(conditionMessage(err),
    "I1 (category 1 used only by respondents with a zero or perfect score in group 'B')",
    fixed = TRUE)
})

test_that("the caption says incomplete responses only when some remain (DL2)", {
  df <- make_polytomous(n = 200, k = 6)
  grp <- rep(c("A", "B"), length.out = nrow(df))
  grp[1:3] <- NA
  k <- suppressMessages(RMdifLR(df, grp))
  expect_false(any(grepl("incomplete responses retained", k)))
  df2 <- df
  df2[5, ] <- NA
  k2 <- suppressMessages(RMdifLR(df2, rep(c("A", "B"), length.out = nrow(df2))))
  expect_false(any(grepl("incomplete responses retained", k2)))
  df3 <- df
  df3[5, 2] <- NA
  k3 <- suppressMessages(RMdifLR(df3, rep(c("A", "B"), length.out = nrow(df3))))
  expect_true(any(grepl("incomplete responses retained", k3)))
})

test_that("edge cases: 2 dichotomous items, too few items, duplicated names", {
  df <- make_dichotomous(n = 200, k = 2)
  grp <- rep(c("A", "B"), length.out = nrow(df))
  # raschmodel() fails in its Hessian with two items; .cml_fit() uses pcmodel()
  out <- rlr(df, grp, output = "dataframe")
  expect_true(is.finite(out$Chisq))
  expect_error(rlr(df[, 1, drop = FALSE], grp), "at least 2 items")
  dup <- cbind(df, df)
  expect_error(rlr(dup, grp), "must be unique. Duplicated: I1, I2")
})


test_that("second review: full-data extreme-only, group maxima (DL11, DL13)", {
  # DL11: category 1 of I1 used only by perfect scorers in the full data
  d <- make_dichotomous(n = 100, k = 6)
  d$I1 <- 0L
  d[1:3, ] <- 1L
  err <- tryCatch(rlr(d, rep(1:2, 50), output = "dataframe"), error = function(e) e)
  expect_s3_class(err, "easyRasch2_category_error")
  expect_match(conditionMessage(err), "zero or perfect total score: I1 (category 1)", fixed = TRUE)
  expect_no_match(conditionMessage(err), "in group")

  # DL13: in group B, I1's top category is unused. Judged against B's own
  # maxima, (1, 2, 2) would look perfect and make I2's category 2 look
  # extreme-only; against the full-data maxima only I1 is named.
  set.seed(9)
  a <- as.data.frame(matrix(sample(0:2, 300, TRUE), ncol = 3))
  b <- data.frame(
    V1 = c(rep(0:1, 30), rep(1L, 3)),
    V2 = c(rep(0:1, each = 30), rep(2L, 3)),
    V3 = c(sample(0:2, 60, TRUE), rep(2L, 3))
  )
  d13 <- rbind(a, b)
  names(d13) <- paste0("I", 1:3)
  g <- c(rep("A", nrow(a)), rep("B", nrow(b)))
  expect_gt(nrow(.extreme_only_categories(d13[g == "B", ])), 0L)
  expect_equal(nrow(.extreme_only_categories(d13[g == "B", ], item_max = c(2, 2, 2))), 0L)
  err <- tryCatch(rlr(d13, g, output = "dataframe"), error = function(e) e)
  expect_identical(err$items, "I1")
})

test_that(".cml_fit(auto) keeps raschmodel for 2 items without a Hessian (DL12)", {
  df <- make_dichotomous(n = 100, k = 2)
  expect_s3_class(.cml_fit(df, model = "auto", hessian = FALSE), "raschmodel")
  expect_s3_class(.cml_fit(df, model = "auto", hessian = TRUE), "pcmodel")
  # with a Hessian, "auto" is always pcmodel() (DL17)
  df10 <- make_dichotomous(n = 100, k = 10)
  expect_s3_class(.cml_fit(df10, model = "auto", hessian = TRUE), "pcmodel")
  expect_s3_class(.cml_fit(df10, model = "auto", hessian = FALSE), "raschmodel")
})

test_that("dichotomous data with respondents answering 1-2 items fit (DL17)", {
  set.seed(1)
  d <- as.data.frame(matrix(stats::rbinom(3000, 1, 0.5), 300))
  names(d) <- paste0("I", 1:10)
  d[1:5, 3:10] <- NA
  out <- rlr(d, rep(1:2, 150), output = "dataframe")
  expect_true(is.finite(out$Chisq))
})


dif_fixture <- function() {
  d <- make_polytomous(n = 300, k = 5, seed = 3L)
  list(d = d, g = rep(c("a", "b"), length.out = nrow(d)))
}

test_that("a failed CML fit names where it failed and points to RMplotTile (DL4)", {
  fx <- dif_fixture()
  local_mocked_bindings(.cml_fit = function(...) stop("boom"))
  expect_error(
    rlr(fx$d, fx$g, output = "dataframe"),
    "CML estimation failed in the full sample (boom). Small or sparse groups",
    fixed = TRUE
  )
  # a failure in a group fit names the group
  n <- 0L
  local_mocked_bindings(.cml_fit = function(d, ...) {
    n <<- n + 1L
    if (n == 2L) stop("group boom")
    psychotools::pcmodel(as.matrix(d), hessian = TRUE)
  })
  expect_error(rlr(fx$d, fx$g, output = "dataframe"),
               "CML estimation failed in group 'a' (group boom)", fixed = TRUE)
})

