# Direct tests for the shared helpers in R/utils-validation.R.
# Code review 2026-11-11, dev/code-review/utils-validation.md (UV1-UV6).

test_that("validate_response_data accepts valid data and reports each problem (UV1-UV3)", {
  ok <- data.frame(a = c(0, 1, 2, NA), b = c(1, 0, 1, 1))
  expect_invisible(validate_response_data(ok))
  expect_true(validate_response_data(as.matrix(ok)))
  # a matrix without column names is fine
  expect_true(validate_response_data(unname(as.matrix(ok))))
  # logical columns are read as 0/1
  expect_true(validate_response_data(data.frame(a = c(TRUE, FALSE), b = c(0, 1))))

  expect_error(validate_response_data(list(a = 1)), "data.frame or matrix")
  # UV1
  expect_error(validate_response_data(ok[, 1, drop = FALSE]), "at least 2 items")
  dup <- ok
  names(dup) <- c("x", "x")
  expect_error(validate_response_data(dup), "must be unique. Duplicated: x.", fixed = TRUE)
  # UV2: named columns, with the factor hint only for factors
  chr <- data.frame(a = c("0", "1"), b = c(1, 0), stringsAsFactors = FALSE)
  e <- tryCatch(validate_response_data(chr), error = conditionMessage)
  expect_match(e, "Non-numeric column(s): a.", fixed = TRUE)
  expect_no_match(e, "Convert factors")
  fct <- data.frame(a = factor(c("0", "1")), b = c(1, 0))
  expect_error(validate_response_data(fct), "Convert factors")
  expect_error(validate_response_data(matrix(c("0", "1", "1", "0"), 2)), "Non-numeric column\\(s\\): 1, 2")
  # UV3
  expect_error(validate_response_data(data.frame(a = c(0, Inf), b = c(1, 0))), "infinite values")
  # the existing checks still fire
  expect_error(validate_response_data(data.frame(a = c(0, -1), b = c(1, 0))), "negative")
  expect_error(validate_response_data(data.frame(a = c(0, 0.5), b = c(1, 0))), "non-integer")
  expect_error(validate_response_data(data.frame(a = c(1, 2), b = c(1, 2))), "minimum value")
  expect_error(validate_response_data(data.frame(a = c(NA, NA), b = c(NA, NA))), "no non-missing")
})

test_that(".drop_empty_respondents keeps order and row names, matrix and data.frame (UV6)", {
  d <- data.frame(a = c(0, NA, 1, NA), b = c(1, NA, NA, NA), row.names = c("p1", "p2", "p3", "p4"))
  expect_message(out <- .drop_empty_respondents(d), "2 respondent\\(s\\) with no responses dropped")
  expect_identical(rownames(out), c("p1", "p3"))
  expect_identical(out$a, c(0, 1))
  m <- as.matrix(d)
  expect_message(outm <- .drop_empty_respondents(m), "2 respondent")
  expect_true(is.matrix(outm))
  expect_identical(rownames(outm), c("p1", "p3"))
  # nothing to drop: silent and unchanged
  expect_silent(same <- .drop_empty_respondents(d[c(1, 3), ]))
  expect_identical(same, d[c(1, 3), ])
})

test_that(".n_caption covers all four shapes and the noun (UV6)", {
  expect_identical(.n_caption(100, 100), "n = 100 respondents")
  expect_identical(.n_caption(90, 100), "n = 90 of 100 respondents")
  expect_identical(.n_caption(100, 100, "incomplete responses retained"),
                   "n = 100 respondents (incomplete responses retained)")
  expect_identical(.n_caption(90, 100, c("complete cases", "extreme scores excluded")),
                   "n = 90 of 100 respondents (complete cases; extreme scores excluded)")
  expect_identical(.n_caption(50, 60, noun = "respondents assessed"),
                   "n = 50 of 60 respondents assessed")
})

test_that(".mi_completed_to_numeric converts code levels and stops on labels (UV4, UV6)", {
  d <- data.frame(
    a = factor(c("0", "2", "1"), levels = c("0", "1", "2")),
    b = c(1, 0, 1)
  )
  out <- .mi_completed_to_numeric(d)
  # the level values, not the 1-based integer codes
  expect_identical(out$a, c(0, 2, 1))
  expect_identical(out$b, c(1, 0, 1))
  # ordered factors (mice polr) too
  d$a <- factor(c("0", "2", "1"), ordered = TRUE)
  expect_identical(.mi_completed_to_numeric(d)$a, c(0, 2, 1))
  # numeric-only data unchanged
  expect_identical(.mi_completed_to_numeric(data.frame(x = 1:2)), data.frame(x = 1:2))
  lab <- data.frame(a = factor(c("never", "often")), b = c(1, 0),
                    c = factor(c("0", "1")))
  expect_error(.mi_completed_to_numeric(lab), "not response codes: a.", fixed = TRUE)
})

test_that(".validate_filename and .write_output_csv (UV6)", {
  for (bad in list(NULL, 1, c("a.csv", "b.csv"), NA_character_, "")) {
    expect_error(.validate_filename(bad), "single non-empty string")
  }
  expect_invisible(.validate_filename("x.csv"))
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- data.frame(Item = c("I1", "I2"), loc = c(-0.5, 0.5), row.names = c("r1", "r2"))
  expect_message(res <- .write_output_csv(df, f), "Wrote 2 row\\(s\\)")
  expect_identical(res, df)
  back <- utils::read.csv(f)
  expect_identical(names(back), c("Item", "loc"))
  expect_equal(back$loc, df$loc)
  expect_message(.write_output_csv(df, f, row.names = TRUE), "Wrote")
  expect_identical(utils::read.csv(f, row.names = 1)$Item, c("I1", "I2"))
  expect_error(.write_output_csv(df, ""), "single non-empty string")
})

test_that(".all_sims_failed_message names extreme-only categories, else the detail (UV6)", {
  generic <- data.frame(a = c(0, 1, 1, 0), b = c(1, 0, 1, 0))
  expect_identical(.all_sims_failed_message(generic),
                   "All simulation iterations failed. Check your data.")
  expect_identical(.all_sims_failed_message(generic, detail = "Try more."),
                   "All simulation iterations failed. Try more.")
  # category 2 of a used only by a perfect scorer
  ext <- data.frame(a = c(0, 1, 1, 2), b = c(1, 0, 1, 1))
  msg <- .all_sims_failed_message(ext)
  expect_match(msg, "^All simulation iterations failed. Response categories used only")
  expect_match(msg, "a (category 2)", fixed = TRUE)
})

# Review of the UV fixes (FR1-FR6): the item count with fixed item_params, the
# order of checks in RMtargeting, and one up-front check in the MI functions.

make_fixed <- function(n = 60, k = 4, seed = 1L) {
  set.seed(seed)
  thr <- lapply(seq(-1, 1, length.out = k), function(b) b + c(-0.6, 0.6))
  names(thr) <- paste0("I", seq_len(k))
  d1 <- as.data.frame(sim_partial_score(thr, stats::rnorm(n)))
  d2 <- as.data.frame(sim_partial_score(thr, stats::rnorm(n)))
  colnames(d1) <- colnames(d2) <- names(thr)
  list(t1 = d1, t2 = d2, thr = thr)
}

test_that("scoring with item_params needs 2+ items, before and after alignment (FR1, FR2)", {
  d <- make_fixed()
  one <- d$thr[1]
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  # one data column: the validator stops
  expect_error(q(RMpersonParameters(d$t1[, 1, drop = FALSE], item_params = one)),
               "at least 2 items (columns). Found 1.", fixed = TRUE)
  # four data columns, one item in item_params: .align_items() stops
  msg <- "`item_params` must cover at least 2 items in `data`. Found 1."
  expect_error(q(RMpersonParameters(d$t1, item_params = one)), msg, fixed = TRUE)
  expect_error(q(RMreliabilityCurve(d$t1, item_params = one)), msg, fixed = TRUE)
  expect_error(q(RMretestSD(d$t1, d$t2, item_params = one, sim_iter = 20,
                            parallel = FALSE, seed = 1)), msg, fixed = TRUE)
  expect_error(q(RMpersonChange(d$t1, d$t2, item_params = one, parallel = FALSE,
                                seed = 1)), msg, fixed = TRUE)
  # two items still work
  res <- q(RMpersonParameters(d$t1, item_params = d$thr[1:2], output = "dataframe"))
  expect_equal(nrow(res), nrow(d$t1))
})

test_that("RMtargeting gives the validator's message for non-numeric data (FR3)", {
  skip_if_not_installed("ggplot2")
  d <- make_fixed()$t1
  d$I2 <- as.character(d$I2)
  expect_error(RMtargeting(d), "Non-numeric column(s): I2.", fixed = TRUE)
})

test_that("the MI functions check the data once, before the imputation loop (FR4)", {
  skip_if_not_installed("mice")
  d <- make_fixed(n = 80)$t1
  d$I1[1:5] <- NA
  imp <- mice::mice(d, m = 2, printFlag = FALSE, seed = 1)
  imp1 <- imp
  imp1$data <- imp1$data[, 1, drop = FALSE]
  imp1$imp <- imp1$imp[1]
  imp1$where <- imp1$where[, 1, drop = FALSE]
  imp1$nmis <- imp1$nmis[1]
  imp1$method <- imp1$method[1]
  imp1$predictorMatrix <- imp1$predictorMatrix[1, 1, drop = FALSE]
  imp1$visitSequence <- imp1$visitSequence[1]
  for (f in list(
    function() RMitemInfitMI(imp1, output = "dataframe"),
    function() RMitemInfitCutoffMI(imp1, iterations = 4, parallel = FALSE, seed = 1)
  )) {
    expect_no_warning(expect_error(f(), "at least 2 items (columns). Found 1.", fixed = TRUE))
  }
})

# Second review (SR1, SR3, SR4).
test_that("RMreliabilityCurve items = needs 2+ distinct items (SR1, SR3)", {
  d <- make_fixed()$t1
  expect_error(RMreliabilityCurve(d, items = "I1"),
               "`items` must select at least 2 items. Found 1.", fixed = TRUE)
  expect_error(RMreliabilityCurve(d, items = c(1, 1)), "must not repeat an item")
  expect_error(RMreliabilityCurve(d, items = c("I1", "I2", "I1")), "must not repeat an item")
  expect_error(RMreliabilityCurve(d, items = TRUE), "character vector of column names")
})

test_that("duplicated names in item_params stop (SR4)", {
  d <- make_fixed()
  dup <- list(I1 = d$thr$I1, I1 = d$thr$I1)
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  msg <- "Item names (names of `item_params`) must be unique. Duplicated: I1."
  expect_error(q(RMpersonParameters(d$t1, item_params = dup)), msg, fixed = TRUE)
  expect_error(q(RMreliabilityCurve(d$t1, item_params = dup)), msg, fixed = TRUE)
})

# Third review (TR1-TR6).
test_that("item_params pair by name only, else stop (TR1, TR4)", {
  d <- make_fixed()
  t1 <- d$t1[, 1:3]
  thr <- d$thr[1:3]
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  ref <- q(RMpersonParameters(t1, item_params = thr, output = "dataframe"))
  # by name, in any order
  rev <- q(RMpersonParameters(t1, item_params = thr[3:1], output = "dataframe"))
  expect_equal(rev$location, ref$location)
  # no name matches, equal count: used to pair by position, now stops
  pos <- stats::setNames(thr, c("x", "y", "z"))
  expect_error(q(RMpersonParameters(t1, item_params = pos)),
               "not found in `data`: x, y, z.", fixed = TRUE)
  # partial overlap used to pair by position silently
  part <- stats::setNames(thr[c(2, 1, 3)], c("I2", "I1", "zz"))
  expect_error(q(RMpersonParameters(t1, item_params = part)),
               "not found in `data`: zz.", fixed = TRUE)
  # some names empty
  some <- thr
  names(some)[2:3] <- ""
  expect_error(q(RMpersonParameters(t1, item_params = some)),
               "Unnamed: element(s) 2, 3.", fixed = TRUE)
})

test_that("RMreliabilityCurve items = rejects NA and fractional indices (TR2, TR3)", {
  d <- make_fixed()$t1
  expect_error(RMreliabilityCurve(d, items = c(1, NA)), "must not contain NA")
  expect_error(RMreliabilityCurve(d, items = c("I1", NA)), "must not contain NA")
  expect_error(RMreliabilityCurve(d, items = c(1, 1.5)), "whole column indices")
})

test_that("a factor covariate in imputed data gets a hint (TR6)", {
  d <- data.frame(a = c(0, 1), b = c(1, 0), sex = factor(c("f", "m")))
  expect_error(.mi_completed_to_numeric(d), "not items, such as covariates")
})
