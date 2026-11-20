# The simulation plots draw the interval stored in the cutoff object, the
# one the tables use (help-page review, findings group A, 2026-11-20).

outer_matches <- function(p, low, high, layer = 2L) {
  b <- ggplot2::ggplot_build(p)$data[[layer]]
  isTRUE(all.equal(sort(b$x), sort(as.numeric(low)))) &&
    isTRUE(all.equal(sort(b$xend), sort(as.numeric(high))))
}

plot_data <- function() phq9[, grep("^q", names(phq9))][, 1:5]

test_that("item plots draw the table's interval", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  skip_if_not_installed("iarm")
  d <- plot_data()

  ic <- RMitemInfitCutoff(d, iterations = 30, parallel = FALSE, seed = 1)
  p <- suppressMessages(RMitemInfitPlot(ic, data = d))
  expect_true(outer_matches(p, ic$item_cutoffs$infit_low,
                            ic$item_cutoffs$infit_high))
  p <- suppressMessages(RMitemInfitPlot(ic, data = d, statistic = "outfit"))
  expect_true(outer_matches(p, ic$item_cutoffs$outfit_low,
                            ic$item_cutoffs$outfit_high))

  rc <- RMitemRestscoreCutoff(d, iterations = 30, parallel = FALSE, seed = 1)
  p <- suppressMessages(RMitemRestscorePlot(rc, data = d))
  expect_true(outer_matches(p, rc$item_cutoffs$diff_low,
                            rc$item_cutoffs$diff_high, layer = 3L))

  set.seed(2)
  grp <- sample(c("a", "b"), nrow(d), replace = TRUE)
  dc <- RMdifGammaCutoff(d, dif_var = grp, iterations = 30,
                         parallel = FALSE, seed = 1)
  p <- suppressMessages(RMdifGammaPlot(dc, data = d, dif_var = grp))
  expect_true(outer_matches(p, dc$item_cutoffs$gamma_low,
                            dc$item_cutoffs$gamma_high))
})

test_that("pair plots draw the table's interval", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  skip_if_not_installed("iarm")
  d <- plot_data()

  qc <- RMlocdepQ3Cutoff(d, iterations = 30, parallel = FALSE, seed = 1)
  p <- suppressMessages(RMlocdepQ3Plot(qc, data = d))$pairs
  expect_true(outer_matches(p, qc$pair_cutoffs$Q3_low,
                            qc$pair_cutoffs$Q3_high))

  gc <- RMlocdepGammaCutoff(d, iterations = 30, parallel = FALSE, seed = 1)
  p <- suppressMessages(RMlocdepGammaPlot(gc, data = d))
  expect_true(outer_matches(p, gc$pair_cutoffs$gamma_low,
                            gc$pair_cutoffs$gamma_high))
})

test_that("RMdimCFAPlot draws the expected range at `percentile`", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  skip_if_not_installed("lavaan")
  d <- plot_data()
  cc <- RMdimCFACutoff(d, iterations = 20, parallel = FALSE, seed = 1)
  p <- suppressMessages(RMdimCFAPlot(cc, data = d))$loadings
  expect_true(outer_matches(p, cc$loading_cutoffs$low,
                            cc$loading_cutoffs$high))
})

test_that("captions name the interval and the table's p-value flagging", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  d <- plot_data()
  qc <- RMlocdepQ3Cutoff(d, iterations = 30, parallel = FALSE, seed = 1,
                         hdci_width = 0.9)
  cap <- suppressMessages(RMlocdepQ3Plot(qc))$pairs$labels$caption
  cap <- gsub("\\s+", " ", cap)
  expect_match(cap, "90% HDCI", fixed = TRUE)
  expect_match(cap, "RMlocdepQ3() flags by adjusted p-values", fixed = TRUE)
})

test_that("RMitemInfitPlot without data honours an explicit `statistic`", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  skip_if_not_installed("iarm")
  d <- plot_data()
  ic <- RMitemInfitCutoff(d, iterations = 30, parallel = FALSE, seed = 1)
  panels <- function(p) {
    length(unique(ggplot2::ggplot_build(p)$data[[1]]$PANEL))
  }
  expect_equal(panels(RMitemInfitPlot(ic)), 2L)
  expect_equal(panels(RMitemInfitPlot(ic, statistic = "both")), 2L)
  p <- RMitemInfitPlot(ic, statistic = "outfit")
  expect_equal(panels(p), 1L)
  expect_equal(unique(p$data$statistic), "OutfitMSQ")
})
