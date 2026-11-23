# Does the rest-score plot read discrimination? One focus item (item 1) with
# a whole-item slope a (generalised PCM, P(X = j) proportional to
# exp(a * (j * theta - sum_{l <= j} delta_l))) among PCM items (a = 1).
# Hypothesis (2026-11-20 conversation): a < 1 gives rest-score SD above
# expected in each category and compressed means (low categories above,
# high categories below expected), a > 1 the reverse. a = 1 is the null
# reference for the outside-the-range rates.
#
# Run from easyRasch2/:
#   Rscript dev/category_conditional_null/R/restscore_discrimination_check.R [reps] [cores]
# Container: reps = 100, cores = 4.

args <- commandArgs(trailingOnly = TRUE)
REPS  <- if (length(args) >= 1) as.integer(args[1]) else 100L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 4L
K <- 10L; M <- 4L; ITER <- 200L

suppressMessages(devtools::load_all(".", quiet = TRUE))
out_dir   <- "dev/category_conditional_null/results"
cells_dir <- file.path(out_dir, "restscore_discr_cells")
dir.create(cells_dir, showWarnings = FALSE, recursive = TRUE)

# Same thresholds as residual_null_check.R. Item locations over [-1, 1],
# so with k = 10 the focus item (item 1) is at -1. Move it to the centre so
# its categories are all well used.
src <- readLines("dev/category_conditional_null/R/residual_null_check.R")
eval(parse(text = src[grep("^true_thr <- ", src):(grep("^sim_pcm <- ", src) - 1L)]))
thr <- true_thr(K, M)
thr[[1]] <- thr[[1]] - mean(thr[[1]])

sim_gpcm <- function(theta, thr, a) {
  vapply(seq_along(thr), function(i) {
    t <- thr[[i]]
    lp <- a[i] * (outer(theta, 0:length(t)) -
      matrix(c(0, cumsum(t)), length(theta), length(t) + 1L, byrow = TRUE))
    p <- exp(lp - apply(lp, 1, max)); p <- p / rowSums(p)
    rowSums(stats::runif(length(theta)) > t(apply(p, 1, cumsum)))
  }, numeric(length(theta)))
}

cells <- expand.grid(a = c(0.5, 1, 2), n = c(300L, 1000L))
t0 <- Sys.time()
for (ci in seq_len(nrow(cells))) {
  cl <- cells[ci, ]
  f <- file.path(cells_dir, sprintf("a%s_n%d.rds", cl$a, cl$n))
  if (file.exists(f)) next
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20261122 + ci)
  res <- parallel::mclapply(seq_len(REPS), function(rp) {
    a <- c(cl$a, rep(1, K - 1L))
    X <- as.data.frame(sim_gpcm(stats::rnorm(cl$n), thr, a))
    names(X) <- paste0("I", seq_len(K))
    fit <- tryCatch(suppressMessages(suppressWarnings(
      RMitemCatFit(X, iterations = ITER, seed = rp, output = "list"))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    d <- RMitemCatFitPlot(X, type = "restscore", fit = fit,
                          output = "dataframe")
    d <- d[as.integer(d$item) == 1L, ]
    data.frame(rep = rp, statistic = d$statistic, category = d$category,
               observed = d$observed, lower = d$lower, upper = d$upper,
               padj_focus = fit$table$padj[1],
               any_other = any(fit$table$padj[-1] < 0.05))
  }, mc.cores = CORES, mc.set.seed = TRUE)
  bad <- vapply(res, function(x) inherits(x, "try-error"), TRUE)
  if (any(bad)) stop("errors in cell ", ci, ": ", res[bad][[1]])
  d <- do.call(rbind, res)
  d$a <- cl$a; d$n <- cl$n
  attr(d, "n_skipped") <- sum(vapply(res, is.null, TRUE))
  saveRDS(d, f)
  message(sprintf("cell %d/%d a=%s n=%d skipped=%d  %.1f min", ci, nrow(cells),
                  cl$a, cl$n, attr(d, "n_skipped"),
                  as.numeric(Sys.time() - t0, units = "mins")))
}

all <- do.call(rbind, lapply(list.files(cells_dir, full.names = TRUE), readRDS))
by_cat <- do.call(rbind, lapply(
  split(all, list(all$a, all$n, all$statistic, all$category), drop = TRUE),
  function(d) data.frame(
    a = d$a[1], n = d$n[1], statistic = d$statistic[1], category = d$category[1],
    reps = nrow(d), mean_diff = mean(d$observed, na.rm = TRUE),
    above = mean(d$observed > d$upper, na.rm = TRUE),
    below = mean(d$observed < d$lower, na.rm = TRUE))))
by_cell <- do.call(rbind, lapply(split(all, list(all$a, all$n), drop = TRUE),
  function(d) {
    r1 <- d[!duplicated(d$rep), ]
    data.frame(a = d$a[1], n = d$n[1], reps = nrow(r1),
               focus_flagged = mean(r1$padj_focus < 0.05),
               other_flagged = mean(r1$any_other))
  }))
out <- list(settings = list(reps = REPS, k = K, m = M, iterations = ITER,
                            thr_focus = thr[[1]]),
            by_cat = by_cat[order(by_cat$a, by_cat$n, by_cat$statistic,
                                  by_cat$category), ],
            by_cell = by_cell)
saveRDS(out, file.path(out_dir, "restscore_discr_summary.rds"))
saveRDS(out, file.path(out_dir, paste0("restscore_discr_summary_", Sys.Date(), ".rds")))
message("done, ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
