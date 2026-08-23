# Display rounding is a house rule with two halves: `output = "dataframe"`
# returns full precision, and `output = "kable"` rounds for display, three
# decimals for a statistic and four for a Monte-Carlo p-value.
#
# The half that breaks silently is the kable one. Rounding is applied by
# `.round_display(df, digits)`, which rounds only the columns named in
# `digits` and passes everything else through untouched, so a column added to
# a table without a matching entry is printed at full precision. That is what
# happened to `gamma_pair` when it was added to the partial gamma tables in
# 1.1.1.9003: nothing errored, the column simply arrived with seven
# significant digits beside its three-decimal neighbours.
#
# Scanning the rendered output only catches it when the values happen to be
# long, so this instruments `.round_display()` instead and records every
# numeric column that reaches it without a rule, whatever the values look
# like in the test data.

# Columns that are counts rather than statistics, where rounding would be
# meaningless. Keep this list short and say why each one is here.
allowed_unrounded <- c(
  "n" # RMitemRestscoreBoot: iterations in which the item was flagged
)

small_poly <- function(n = 120L, k = 5L, seed = 4L) {
  set.seed(seed)
  th <- stats::rnorm(n)
  df <- as.data.frame(sapply(
    seq_len(k),
    function(j) stats::rbinom(n, 2, stats::plogis(th - (j - (k + 1) / 2) / 3))
  ))
  names(df) <- paste0("q", seq_len(k))
  df
}

test_that("every numeric column shown in a kable has a rounding rule", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")

  seen <- new.env(parent = emptyenv())
  original <- easyRasch2:::.round_display
  local_mocked_bindings(
    .round_display = function(df, digits) {
      numeric_cols <- names(df)[vapply(df, is.numeric, logical(1L))]
      gap <- setdiff(numeric_cols, names(digits))
      for (nm in gap) assign(nm, TRUE, envir = seen)
      original(df, digits)
    },
    .package = "easyRasch2"
  )

  q <- function(expr) suppressWarnings(suppressMessages(expr))
  df <- small_poly()
  grp <- factor(rep(c("a", "b"), length.out = nrow(df)))

  cut_infit <- q(RMitemInfitCutoff(df, iterations = 20, parallel = FALSE, seed = 1))
  cut_q3 <- q(RMlocdepQ3Cutoff(df, iterations = 20, parallel = FALSE, seed = 1))
  cut_gamma <- q(RMlocdepGammaCutoff(df, iterations = 20, parallel = FALSE, seed = 1))
  cut_pca <- q(RMdimResidualPCACutoff(df, iterations = 20, parallel = FALSE, seed = 1))
  cut_dif <- q(RMdifGammaCutoff(df, dif_var = grp, iterations = 20,
                                parallel = FALSE, seed = 1))

  # Both flagging branches of every function that has them, since the two
  # return different column sets.
  invisible(list(
    q(RMitemInfit(df)),
    q(RMitemInfit(df, cutoff = cut_infit)),
    q(RMitemInfit(df, cutoff = cut_infit, p_value = FALSE)),
    q(RMitemInfit(df, cutoff = cut_infit, statistic = "outfit")),
    q(RMlocdepQ3(df, cutoff = cut_q3)),
    q(RMlocdepQ3(df, cutoff = cut_q3, p_value = FALSE)),
    q(RMlocdepGamma(df)),
    q(RMlocdepGamma(df, cutoff = cut_gamma)),
    q(RMlocdepGamma(df, cutoff = cut_gamma, p_value = FALSE)),
    q(RMdifGamma(df, dif_var = grp)),
    q(RMdifGamma(df, dif_var = grp, cutoff = cut_dif, p_value = TRUE)),
    q(RMdimResidualPCA(df, cutoff = cut_pca, p_value = TRUE)),
    q(RMitemRestscore(df)),
    q(RMitemRestscoreBoot(df, iterations = 10, samplesize = 100,
                          parallel = FALSE, seed = 1)),
    q(RMitemParameters(df)),
    q(RMpersonParameters(df)),
    q(RMscoreSE(df, output = "kable")),
    q(RMdifLR(df, dif_var = grp)),
    q(RMpersonFit(df, iterations = 10, seed = 1))
  ))

  expect_setequal(setdiff(ls(seen), allowed_unrounded), character(0))
})

test_that("the guard notices a column with no rounding rule", {
  # Without this the test above passes just as happily when the instrument
  # is broken as when the package is correct.
  seen <- new.env(parent = emptyenv())
  original <- easyRasch2:::.round_display
  instrumented <- function(df, digits) {
    numeric_cols <- names(df)[vapply(df, is.numeric, logical(1L))]
    for (nm in setdiff(numeric_cols, names(digits))) assign(nm, TRUE, envir = seen)
    original(df, digits)
  }
  tbl <- data.frame(item = "q1", stat = 1 / 3, forgotten = 2 / 3)
  out <- instrumented(tbl, c(stat = 3))
  expect_identical(ls(seen), "forgotten")
  expect_equal(out$stat, 0.333)
  expect_equal(out$forgotten, 2 / 3)
})

test_that("kable rounds for display while dataframe keeps full precision", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  df <- small_poly()
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  cut_q3 <- q(RMlocdepQ3Cutoff(df, iterations = 20, parallel = FALSE, seed = 1))

  raw <- q(RMlocdepQ3(df, cutoff = cut_q3, output = "dataframe")$pairs)
  # full precision: at least one value that is not its own 4-decimal rounding
  expect_true(any(raw$Observed != round(raw$Observed, 4)))

  shown <- paste(as.character(q(RMlocdepQ3(df, cutoff = cut_q3))$pairs),
                 collapse = "\n")
  body <- grep("^[|]", strsplit(shown, "\n")[[1]], value = TRUE)
  expect_length(unlist(regmatches(body, gregexpr("[0-9]+[.][0-9]{5,}", body))), 0L)
})
