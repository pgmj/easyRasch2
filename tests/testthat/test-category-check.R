# The shared category check (.check_response_categories(), via .cml_fit() or
# directly before mirt/eRm fits) must stop every model-fitting function on an
# item with an empty lowest or interior category, or a single category.
# Code review 2026-11-10, dev/code-review/item_parameters.md (step 2).

mid_empty <- function() {
  skip_if_not_installed("eRm")
  e <- new.env()
  utils::data(pcmdat2, package = "eRm", envir = e)
  x <- as.data.frame(e$pcmdat2)
  x[x[, 1] == 1, 1] <- 2L
  x
}

dich_constant <- function() {
  skip_if_not_installed("eRm")
  e <- new.env()
  utils::data(raschdat1, package = "eRm", envir = e)
  x <- as.data.frame(e$raschdat1)[, 1:10]
  x[, 3] <- 1L
  x
}

msg_mid <- "I1 \\(empty category 1\\)"

test_that("CML item-fit and location functions stop on an empty interior category", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggplot2")
  d <- mid_empty()
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  expect_error(q(RMitemInfit(d, output = "dataframe")), msg_mid)
  expect_error(q(RMitemRestscore(d, output = "dataframe")), msg_mid)
  expect_error(q(RMitemHierarchy(d, output = "dataframe")), msg_mid)
  expect_error(q(RMtargeting(d, output = "list")), msg_mid)
  expect_error(
    q(RMitemRestscoreBoot(d, iterations = 5, samplesize = 150,
                          parallel = FALSE, seed = 1)),
    msg_mid
  )
})

test_that("mirt-based paths stop on an empty interior category", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  d <- mid_empty()
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  expect_error(q(RMlocdepQ3(d, output = "dataframe", estimator = "MML")), msg_mid)
  expect_error(q(RMscoreSE(d, method = "EAP", output = "dataframe")), msg_mid)
  expect_error(q(RMtargeting(d, robust = TRUE, output = "list")), msg_mid)
  expect_error(
    q(RMreliability(d, draws = 20, rmu_iter = 3, parallel = FALSE, seed = 1,
                    output = "dataframe")),
    msg_mid
  )
})

test_that("RMdifLR and RMdimMartinLof check the full data", {
  d <- mid_empty()
  g <- factor(rep(c("a", "b"), length.out = nrow(d)))
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  expect_error(q(RMdifLR(d, g, output = "dataframe")), msg_mid)
  dc <- dich_constant()
  expect_error(
    q(RMdimMartinLof(dc, partition = rep(1:2, 5), iterations = 20L,
                     parallel = FALSE, seed = 1)),
    "I3 \\(only category 1\\)"
  )
})

test_that("RMitemInfitMI stops when every imputation lacks the category", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  d <- mid_empty()
  set.seed(2)
  d[sample(nrow(d), 20), 2] <- NA
  imp <- mice::mice(d, m = 2, printFlag = FALSE, seed = 1)
  expect_error(
    suppressWarnings(RMitemInfitMI(imp, output = "dataframe")),
    msg_mid
  )
})

test_that(".notify_dropped_resamples counts resamples by reason (D4a, E7, E8)", {
  cat_err <- structure(
    "Item thresholds cannot be estimated for: ...",
    items = c("Q1 (rev)", "I3")
  )
  res <- c(
    rep(list(list(ok = 1)), 7),
    list(cat_err, cat_err, "PSI not estimable for this resample")
  )
  expect_message(.notify_dropped_resamples(res), "3 of 10 bootstrap resamples were dropped")
  expect_message(.notify_dropped_resamples(res), "The interval uses the remaining 7")
  # names come from the attribute, so names with " (" survive (E8)
  expect_message(.notify_dropped_resamples(res), "Q1 \\(rev\\) in 2, I3 in 2|I3 in 2, Q1 \\(rev\\) in 2")
  expect_message(.notify_dropped_resamples(res), "PSI not estimable")
  expect_message(.notify_dropped_resamples(res, what = "band"), "The band uses")
  # with fewer than 2 left, no "uses the remaining 1" (E7), and nothing that
  # repeats the caller's "not reported" warning (G9)
  few <- c(list(list(ok = 1)), rep(list(cat_err), 4))
  msg <- conditionMessage(rlang::catch_cnd(.notify_dropped_resamples(few)))
  expect_match(msg, "4 of 5 bootstrap resamples were dropped.", fixed = TRUE)
  expect_no_match(msg, "remaining|Too few")
  expect_match(msg, "leaves out resamples that lost", fixed = TRUE)
  expect_silent(.notify_dropped_resamples(rep(list(list(ok = 1)), 5)))
})

test_that("the category error carries the item names", {
  e <- tryCatch(
    .check_response_categories(cbind(a = c(0, 2, 2, 0), b = c(0, 1, 1, 0), c = c(1, 1, 1, 1))),
    error = function(e) e
  )
  expect_s3_class(e, "easyRasch2_category_error")
  expect_equal(e$items, c("a", "c"))
})

test_that("complete-case functions stop when the complete cases lose a category (E2)", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  e <- new.env()
  utils::data(pcmdat2, package = "eRm", envir = e)
  x <- as.data.frame(e$pcmdat2)
  x[x$I1 == 1, "I2"] <- NA # category 1 of I1 occurs only in incomplete rows
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  expect_error(q(RMitemInfit(x, output = "dataframe")), "complete cases only.*I1 \\(empty category 1\\)")
  expect_error(q(RMitemRestscore(x, output = "dataframe")), "complete cases only.*I1 \\(empty category 1\\)")
  expect_error(
    q(RMitemRestscoreBoot(x, iterations = 5, samplesize = 100, parallel = FALSE, seed = 1)),
    "complete cases only"
  )
})

test_that("restscore bootstrap matches iarm rows by item name (E1)", {
  skip_if_not_installed("iarm")
  set.seed(7)
  d <- as.data.frame(matrix(rbinom(300 * 6, 1, 0.5), 300, 6))
  names(d) <- paste0("V", 1:6)
  d$V3 <- 0L
  r <- suppressWarnings(
    run_single_boot_restscore(1, list(data = d, samplesize = 300, item_names = names(d)))
  )
  expect_s3_class(r, "data.frame")
  expect_true(is.na(r$diff[r$Item == "V3"]))
  # Same resample as the worker draws (its seed and RNG kinds).
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  db <- d[sample.int(300, 300, replace = TRUE), , drop = FALSE]
  ref <- suppressWarnings(as.data.frame(iarm::item_restscore(
    psychotools::pcmodel(db[, -3], hessian = FALSE)
  ))[[1L]])
  expect_equal(
    r$diff[r$Item == "V6"],
    as.numeric(ref["V6", 2L]) - as.numeric(ref["V6", 1L])
  )
})

test_that("an empty category in the full data is not blamed on complete cases (F1)", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  d <- mid_empty()
  d[5, 3] <- NA
  err <- tryCatch(suppressWarnings(RMitemInfit(d, output = "dataframe")),
                  error = function(e) conditionMessage(e))
  expect_match(err, "I1 \\(empty category 1\\)")
  expect_false(grepl("complete cases only", err))
})

test_that("RMitemRestscoreBoot gives an item its own base and no NaN (F2)", {
  skip_if_not_installed("iarm")
  set.seed(11)
  d <- as.data.frame(matrix(rbinom(400 * 6, 1, 0.5), 400, 6))
  names(d) <- paste0("V", 1:6)
  d$V3 <- 0L
  d$V3[1] <- 1L # one endorsement: constant in most resamples of 50
  res <- suppressWarnings(RMitemRestscoreBoot(
    d, iterations = 5, samplesize = 50, parallel = FALSE, seed = 1,
    output = "dataframe"
  ))
  expect_false(any(is.nan(res$percent)))
  v3 <- res[res$Item == "V3", ]
  expect_true(all(v3$n_tested < 5L))
  expect_true(all(res$n_tested[res$Item == "V1"] == 5L))
})

test_that("a category error comes without the extreme-only warning first (G5)", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  d <- as.data.frame(matrix(sample(0:2, 300, TRUE), ncol = 3,
                            dimnames = list(NULL, paste0("I", 1:3))))
  # I1: category 1 empty, category 2 used only by a perfect scorer
  d$I1 <- 0L
  d <- rbind(d, data.frame(I1 = 2L, I2 = 2L, I3 = 2L))
  expect_equal(nrow(.extreme_only_categories(d)), 1L)
  calls <- list(
    function() RMitemParameters(d, output = "dataframe"),
    function() RMlocdepQ3(d),
    function() RMdimResidualPCA(d),
    function() RMitemICCPlot(d)
  )
  for (f in calls) {
    expect_no_warning(expect_error(f(), "I1 \\(empty category 1\\)"))
  }
})

test_that("RMitemInfitMI skips imputations that lack a category (G4)", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  d <- mid_empty()
  set.seed(2)
  miss <- sample(nrow(d), 20)
  d[miss, 1] <- NA
  imp <- mice::mice(d, m = 4, printFlag = FALSE, seed = 1)
  # pmm cannot impute the absent category 1; put it into imputations 1-3 by
  # hand, as polr or polyreg could, and leave imputation 4 without it
  imp$imp[[1]][1:5, 1:3] <- 1L
  ws <- character()
  res <- withCallingHandlers(
    RMitemInfitMI(imp, output = "dataframe"),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(res, "data.frame")
  expect_true(any(grepl("1 of 4 imputed datasets lacked a response category", ws, fixed = TRUE)))
  expect_true(any(grepl("I1 in 1", ws, fixed = TRUE)))
  expect_false(any(grepl("failed and were excluded", ws, fixed = TRUE)))
})

test_that(".mml_fit() and .cml_fit() share the category check (G8)", {
  skip_if_not_installed("mirt")
  d <- mid_empty()
  expect_error(.mml_fit(d), msg_mid)
  expect_error(.cml_fit(d), msg_mid)
  # check = FALSE leaves psychotools' own handling (the RMitemRestscoreBoot
  # resamples and Martin-Löf simulated datasets rely on it)
  expect_s3_class(suppressWarnings(.cml_fit(d, check = FALSE)), "pcmodel")
})

test_that("no model fit bypasses .cml_fit() or .mml_fit() (G8)", {
  # Source files are only present when testing from the source tree.
  r_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "R sources not available")
  files <- setdiff(list.files(r_dir, "\\.R$", full.names = TRUE),
                   file.path(r_dir, "utils-theta.R"))
  hits <- unlist(lapply(files, function(f) {
    code <- readLines(f, warn = FALSE)
    code <- gsub('"[^"]*"', '""', code) # drop strings (a "#" inside one)
    code <- sub("#.*$", "", code) # then comments and roxygen
    pat <- "\\b(pcmodel|raschmodel)\\b|\\bmirt::mirt\\b|\\beRm::(PCM|RM|LRtest)\\b"
    i <- grep(pat, code)
    if (length(i)) paste0(basename(f), ":", i) else character()
  }))
  expect_identical(hits, character())
})
