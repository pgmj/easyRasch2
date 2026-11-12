# Tests for RMdifLR()
#
# Covers:
#   - Input validation (data, dif_var, cutoff, conf, model)
#   - Auto-detect of model (PCM vs RM) by data shape
#   - Output structures: "dataframe", "kable", "ggplot"
#   - Flagged column / MaxDiff
#   - LR test attribute
#   - NA handling in dif_var
#   - sort = TRUE for kable

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
  expect_error(RMdifLR(df, dif_var = grp), regexp = "scored starting at 0")
})

test_that("RMdifLR errors when dif_var length does not match nrow(data)", {
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df) - 1L))
  expect_error(RMdifLR(df, dif_var = grp), regexp = "same length as nrow")
})

test_that("RMdifLR errors when dif_var has fewer than 2 distinct levels", {
  df <- make_dichotomous()
  grp <- factor(rep("A", nrow(df)))
  expect_error(RMdifLR(df, dif_var = grp), regexp = "at least 2 distinct")
})

test_that("RMdifLR errors when model = 'RM' on polytomous data", {
  df <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  expect_error(
    RMdifLR(df, dif_var = grp, model = "RM"),
    regexp = "dichotomous"
  )
})

test_that("RMdifLR errors on invalid cutoff or conf", {
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  expect_error(RMdifLR(df, dif_var = grp, cutoff = -1), regexp = "non-negative")
  expect_error(
    RMdifLR(df, dif_var = grp, cutoff = c(0.1, 0.2)),
    regexp = "non-negative"
  )
  expect_error(RMdifLR(df, dif_var = grp, conf = 0), regexp = "in \\(0, 1\\)")
  expect_error(RMdifLR(df, dif_var = grp, conf = 1), regexp = "in \\(0, 1\\)")
})

# ---------------------------------------------------------------------
# Output structure (dataframe)
# ---------------------------------------------------------------------
test_that("RMdifLR output = 'dataframe' on dichotomous data picks RM and returns one row per item", {
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))

  out <- RMdifLR(df, dif_var = grp, output = "dataframe")

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), ncol(df)) # one row per item
  expect_true(all(
    c("Item", "A", "B", "All", "MaxDiff", "SE_A", "SE_B", "SE_All") %in%
      names(out)
  ))
  expect_equal(as.character(out$Item), names(df))

  lr <- attr(out, "lr_test")
  expect_type(lr, "list")
  expect_true(all(c("LR", "df", "p_value", "model") %in% names(lr)))
  expect_equal(lr$model, "RM") # auto-picked
})

test_that("RMdifLR level = 'threshold' on polytomous data returns one row per item × threshold", {
  df <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))

  out <- RMdifLR(df, dif_var = grp, level = "threshold", output = "dataframe")

  expect_s3_class(out, "data.frame")
  expect_true(all(
    c("Item", "Threshold", "A", "B", "All", "MaxDiff") %in% names(out)
  ))
  # 0/1/2 data → 2 thresholds per item
  expect_equal(nrow(out), ncol(df) * 2L)
  expect_equal(attr(out, "lr_test")$model, "PCM") # auto-picked
})

# ---------------------------------------------------------------------
# Output structure (kable / ggplot)
# ---------------------------------------------------------------------
test_that("RMdifLR output = 'kable' returns a knitr_kable", {
  skip_if_not_installed("knitr")
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- RMdifLR(df, dif_var = grp, output = "kable")
  expect_s3_class(out, "knitr_kable")
})

test_that("RMdifLR defaults to kable (consistent with other kable/ggplot functions)", {
  skip_if_not_installed("knitr")
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  expect_s3_class(RMdifLR(df, dif_var = grp), "knitr_kable")
})

test_that("RMdifLR output = 'ggplot' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- RMdifLR(df, dif_var = grp, output = "ggplot")
  expect_s3_class(out, "ggplot")
})

# ---------------------------------------------------------------------
# Flagged column / MaxDiff
# ---------------------------------------------------------------------
test_that("Flagged column is logical and absent when cutoff = NULL", {
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))

  with_cutoff <- RMdifLR(df, dif_var = grp, output = "dataframe", cutoff = 0.5)
  expect_true("Flagged" %in% names(with_cutoff))
  expect_type(with_cutoff$Flagged, "logical")
  expect_true(all(with_cutoff$MaxDiff[with_cutoff$Flagged] > 0.5, na.rm = TRUE))

  no_cutoff <- RMdifLR(df, dif_var = grp, output = "dataframe", cutoff = NULL)
  expect_false("Flagged" %in% names(no_cutoff))
})

# ---------------------------------------------------------------------
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
  expect_equal(attr(out, "lr_test")$n_persons, nrow(df) - 3L)
})

# ---------------------------------------------------------------------
# Sorting (kable)
# ---------------------------------------------------------------------
test_that("sort = TRUE in dataframe output sorts by MaxDiff descending", {
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- RMdifLR(df, dif_var = grp, output = "dataframe", sort = TRUE)
  expect_true(all(diff(out$MaxDiff) <= 0))
})

# ---------------------------------------------------------------------
# p-value formatting in captions
# ---------------------------------------------------------------------
test_that("captions format p as plain numbers, never scientific notation", {
  skip_if_not_installed("ggplot2")
  df <- make_dichotomous()
  grp <- rep(c("A", "B"), length.out = nrow(df))

  p <- RMdifLR(df, dif_var = grp, output = "ggplot")
  expect_match(p$labels$caption, "p (=|<) 0?\\.\\d{3}")
  expect_no_match(p$labels$caption, "e-0")

  kbl <- RMdifLR(df, dif_var = grp)
  expect_match(paste(kbl, collapse = "\n"), "p (=|<) 0?\\.\\d{3}")
  expect_no_match(paste(kbl, collapse = "\n"), "e-0")
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
    lr <- attr(RMdifLR(cs[[1]], cs[[2]], output = "dataframe"), "lr_test")
    fit <- if (cs[[3]] == "PCM") eRm::PCM(cs[[1]]) else eRm::RM(cs[[1]])
    ref <- eRm::LRtest(fit, splitcr = factor(cs[[2]]))
    expect_equal(lr$LR, unname(ref$LR), tolerance = 1e-4)
    expect_identical(as.integer(lr$df), as.integer(ref$df))
  }
})

test_that("dichotomous locations match eRm; polytomous match up to a group shift", {
  skip_if_not_installed("eRm")
  e <- new.env()
  utils::data(pcmdat2, raschdat1, package = "eRm", envir = e)
  set.seed(5)
  ras <- as.data.frame(e$raschdat1)[, 1:10]
  g <- factor(sample(c("a", "b"), nrow(ras), TRUE))
  out <- RMdifLR(ras, g, output = "dataframe")
  ref <- eRm::LRtest(eRm::RM(ras), splitcr = g)
  expect_equal(out$a, unname(-ref$fitobj[[1]]$betapar), tolerance = 1e-3)
  expect_equal(out$SE_a, unname(ref$fitobj[[1]]$se.beta), tolerance = 1e-3)

  pcm <- as.data.frame(e$pcmdat2)
  g <- factor(sample(c("a", "b"), nrow(pcm), TRUE))
  thr <- RMdifLR(pcm, g, level = "threshold", output = "dataframe")
  ref <- eRm::LRtest(eRm::PCM(pcm), splitcr = g)
  er <- unname(eRm::thresholds(ref$fitobj[[1]])$threshpar)
  # psychotools makes thresholds sum to zero within each group, eRm makes
  # the item-category parameters do so: a constant shift per group.
  d <- er - thr$a
  expect_lt(max(d) - min(d), 1e-3)
})

test_that("item locations and SEs equal psychotools::itempar() for equal category counts", {
  # With equal counts, sum-zero thresholds and sum-zero item means are the
  # same identification, so locations and SEs must agree exactly.
  df <- make_polytomous(n = 200)
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- RMdifLR(df, grp, output = "dataframe")
  ip <- psychotools::itempar(psychotools::pcmodel(df), vcov = TRUE)
  expect_equal(out$All, unname(as.numeric(ip)), tolerance = 1e-8)
  expect_equal(out$SE_All, unname(sqrt(diag(stats::vcov(ip)))), tolerance = 1e-8)
})

test_that("a category unused in one group stops, naming item, category and group", {
  df <- make_polytomous(n = 120)
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  # top category of I2 unused in group B (allowed in the full data, not
  # across groups), and the lowest category of I4 unused in group A
  df$I2[grp == "B" & df$I2 == 2L] <- 1L
  df$I4[grp == "A" & df$I4 == 0L] <- 1L
  err <- tryCatch(RMdifLR(df, grp, output = "dataframe"), error = function(e) e)
  expect_s3_class(err, "easyRasch2_category_error")
  expect_match(conditionMessage(err), "I2 (category 2 unused in group 'B')", fixed = TRUE)
  expect_match(conditionMessage(err), "I4 (category 0 unused in group 'A')", fixed = TRUE)
  expect_identical(err$items, c("I2", "I4"))
  # an item with no responses in a group
  df2 <- make_polytomous(n = 120)
  df2$I3[grp == "A"] <- NA
  expect_error(RMdifLR(df2, grp, output = "dataframe"), "I3 (no responses in group 'A')", fixed = TRUE)
})

test_that("with mixed category counts, item locations equal itempar() up to a constant", {
  # itempar() makes the item means sum to zero, RMdifLR() the thresholds:
  # with different counts the two differ by a constant (SEs then depend on
  # the identification, so only locations are compared).
  set.seed(7)
  df <- data.frame(
    I1 = sample(0:1, 300, TRUE), I2 = sample(0:1, 300, TRUE),
    I3 = sample(0:3, 300, TRUE), I4 = sample(0:3, 300, TRUE),
    I5 = sample(0:2, 300, TRUE)
  )
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- RMdifLR(df, grp, output = "dataframe")
  ip <- as.numeric(psychotools::itempar(psychotools::pcmodel(df)))
  d <- out$All - ip
  expect_lt(max(d) - min(d), 1e-8)
})

test_that("a category used in a group only by extreme scorers stops (DL1)", {
  df <- make_dichotomous(n = 120, k = 5)
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  # I1 = 1 in group B only for two perfect scorers
  df$I1[grp == "B"] <- 0L
  b <- which(grp == "B")[1:2]
  df[b, ] <- 1L
  err <- tryCatch(RMdifLR(df, grp, output = "dataframe"), error = function(e) e)
  expect_s3_class(err, "easyRasch2_category_error")
  expect_match(conditionMessage(err),
    "I1 (category 1 used only by respondents with a zero or perfect score in group 'B')",
    fixed = TRUE)
})

test_that("thresholds sort numerically, c10 after c9 (DL3)", {
  set.seed(8)
  df <- data.frame(I1 = rep(0:10, length.out = 440), I2 = rep(0:10, length.out = 440),
                   I3 = sample(0:2, 440, TRUE))
  df$I1 <- sample(df$I1); df$I2 <- sample(df$I2)
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- suppressWarnings(RMdifLR(df, grp, level = "threshold", output = "dataframe"))
  expect_identical(out$Threshold[out$Item == "I1"], paste0("c", 1:10))
  skip_if_not_installed("ggplot2")
  p <- suppressWarnings(RMdifLR(df, grp, level = "threshold", output = "ggplot"))
  expect_identical(levels(p$data$Threshold)[1:10], paste0("c", 1:10))
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

test_that("edge cases: 2 dichotomous items, too few items, names, reserved levels", {
  df <- make_dichotomous(n = 200, k = 2)
  grp <- rep(c("A", "B"), length.out = nrow(df))
  # raschmodel() fails in its Hessian with two items; .cml_fit() uses pcmodel()
  out <- RMdifLR(df, grp, output = "dataframe")
  expect_equal(nrow(out), 2L)
  expect_true(all(is.finite(out$SE_All)))
  expect_error(RMdifLR(df[, 1, drop = FALSE], grp), "at least 2 items")
  dup <- cbind(df, df)
  expect_error(RMdifLR(dup, grp), "must be unique. Duplicated: I1, I2")
  df5 <- make_dichotomous(n = 200)
  # a level named like an output column would overwrite it silently
  expect_error(RMdifLR(df5, rep(c("All", "x"), 100)), "'All' would clash")
  expect_error(RMdifLR(df5, rep(c("Item", "SE_x"), 100)), "'Item', 'SE_x' would clash")
})

test_that("second review: empty level, full-data extreme-only, group maxima (DL10, DL11, DL13)", {
  df5 <- make_dichotomous(n = 200)
  expect_error(RMdifLR(df5, rep(c("", "x"), 100)), "'' would clash")

  # DL11: category 1 of I1 used only by perfect scorers in the full data
  d <- make_dichotomous(n = 100, k = 6)
  d$I1 <- 0L
  d[1:3, ] <- 1L
  err <- tryCatch(RMdifLR(d, rep(1:2, 50), output = "dataframe"), error = function(e) e)
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
  err <- tryCatch(RMdifLR(d13, g, output = "dataframe"), error = function(e) e)
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
  out <- RMdifLR(d, rep(1:2, 150), output = "dataframe")
  expect_equal(nrow(out), 10L)
  expect_true(all(is.finite(out$SE_All)))
})

# ---------------------------------------------------------------------
# Remaining test gaps from the dif_lr review (2026-11-11)
# ---------------------------------------------------------------------
dif_fixture <- function() {
  d <- make_polytomous(n = 300, k = 5, seed = 3L)
  list(d = d, g = rep(c("a", "b"), length.out = nrow(d)))
}

test_that("the kable bolds flagged rows only, and says so in the caption", {
  skip_if_not_installed("knitr")
  fx <- dif_fixture()
  df <- RMdifLR(fx$d, fx$g, output = "dataframe")
  # a cutoff between the MaxDiff values gives flagged and unflagged rows
  cut <- stats::median(df$MaxDiff)
  flagged <- as.character(df$Item[df$MaxDiff > cut])
  kept <- as.character(df$Item[df$MaxDiff <= cut])
  expect_true(length(flagged) > 0L && length(kept) > 0L)
  k <- RMdifLR(fx$d, fx$g, cutoff = cut)
  rows <- as.character(k)
  row_of <- function(it) rows[grepl(paste0("^\\|", it, "\\s*\\|"), rows)]
  for (it in flagged) {
    r <- row_of(it)
    expect_length(r, 1L)
    expect_match(r, "**yes**", fixed = TRUE)
    # every data cell is bold, the item name is not
    expect_false(grepl("**I", r, fixed = TRUE))
    cells <- trimws(strsplit(r, "|", fixed = TRUE)[[1L]])
    cells <- cells[nzchar(cells)][-1L]
    expect_true(all(startsWith(cells, "**") & endsWith(cells, "**")))
  }
  for (it in kept) {
    r <- row_of(it)
    expect_false(grepl("**", r, fixed = TRUE))
    expect_match(r, "|no ", fixed = TRUE)
  }
  cap <- function(x) grep("^Table:", as.character(x), value = TRUE)
  expect_match(cap(k), paste0("**Bold** = MaxDiff > ", cut, " logits."),
               fixed = TRUE)
  # no cutoff: no Flagged column, no bold, no bold note
  k0 <- RMdifLR(fx$d, fx$g, cutoff = NULL)
  expect_false(any(grepl("**", as.character(k0), fixed = TRUE)))
  expect_length(cap(k0), 1L)
  expect_false(any(grepl("Flagged", as.character(k0))))
  expect_no_match(cap(k0), "Bold")
  # sort = TRUE orders the kable rows by MaxDiff
  ks <- as.character(RMdifLR(fx$d, fx$g, cutoff = NULL, sort = TRUE))
  body <- ks[grepl("^\\|I[0-9]", ks)]
  expect_identical(sub("^\\|(I[0-9]+).*", "\\1", body),
                   as.character(df$Item[order(-df$MaxDiff)]))
})

test_that("the item ggplot draws group locations with item-SE bars of the conf width", {
  skip_if_not_installed("ggplot2")
  fx <- dif_fixture()
  df <- RMdifLR(fx$d, fx$g, output = "dataframe")
  p <- RMdifLR(fx$d, fx$g, output = "ggplot")
  pd <- p$data
  expect_setequal(as.character(unique(pd$DIFgroup)), c("a", "b"))
  for (g in c("a", "b")) {
    sub_g <- pd[pd$DIFgroup == g, ]
    sub_g <- sub_g[match(as.character(df$Item), as.character(sub_g$Item)), ]
    expect_equal(sub_g$Location, df[[g]], tolerance = 1e-10)
    expect_equal(sub_g$SE, df[[paste0("SE_", g)]], tolerance = 1e-10)
  }
  half <- function(pl) {
    eb <- ggplot2::layer_data(pl, 3L) # geom_errorbar
    sort((eb$ymax - eb$ymin) / 2)
  }
  z95 <- stats::qnorm(0.975)
  expect_equal(half(p), sort(pd$SE * z95), tolerance = 1e-10)
  p80 <- RMdifLR(fx$d, fx$g, output = "ggplot", conf = 0.80)
  expect_equal(half(p80) / half(p), rep(stats::qnorm(0.90) / z95, nrow(pd)),
               tolerance = 1e-10)
  expect_match(p80$labels$caption, "80% CI", fixed = TRUE)
  # the whole-sample diamonds sit at the All locations
  dm <- ggplot2::layer_data(p, 4L)
  expect_equal(sort(dm$y), sort(df$All), tolerance = 1e-10)
})

test_that("a failed CML fit names where it failed and points to RMplotTile (DL4)", {
  fx <- dif_fixture()
  local_mocked_bindings(.cml_fit = function(...) stop("boom"))
  expect_error(
    RMdifLR(fx$d, fx$g, output = "dataframe"),
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
  expect_error(RMdifLR(fx$d, fx$g, output = "dataframe"),
               "CML estimation failed in group 'a' (group boom)", fixed = TRUE)
})
