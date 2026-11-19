# Step 3c of RMdifImpact (dev/dif-impact-design.md, "Build notes"): rates at
# a sum-score cut with the package's classification rule.
#
# Step 3b validated the rates of a fixed rule with the dev pipeline, whose
# split rule put a focal person at or above a sum-score cut when the
# equated raw score was at or above the cut. RMdifImpact() uses the equated
# score at or above cut - 0.5 (nearest whole score), and at the true
# parameters that rule differs from keep in 8 of the 14 settings, against 2
# for the dev rule (dev/dif_impact_sum_rule_check.R). This step reruns the
# rates with the package itself: every estimate and SE comes from the
# package's internal .dimp_estimate(), the code RMdifImpact() runs, so the
# run also covers the review fixes of 2026-11-19 (one WLE solver, exact
# equating for the rule).
#
# Per replication: DIF data as in step 3b (sim_twin(), true DIF items),
# .dimp_estimate() with the sum-score cut and the step 3 target population.
# Targets: the rates of the fitted rule (rebuilt from the reported lowest
# raw scores) under the true item parameters, in the sample (true latent
# N(0, 1) and N(mu2, 1)) and in the target population. Delta-method
# intervals only, as in the package.
#
# Seeds: one L'Ecuyer-CMRG stream per replication, derived from
# 20261140 + cell, so the results do not depend on the worker count.
#
# Run from the package root. Smoke test (container, 4 workers):
#   DIF_3C_SMOKE=1 LANG=C.UTF-8 LC_ALL=C.UTF-8 Rscript dev/dif_impact_ext_sim_3c.R
# Full run on the host (about 1 to 1.5 h on 10 workers; restarts skip
# finished cells, one file per cell in ext3c_results_cells/):
#   DIF_3C_WORKERS=10 LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 \
#     nohup Rscript dev/dif_impact_ext_sim_3c.R > dev/dif_impact_ext_sim/run_3c.log 2>&1 &
# Optional: DIF_3C_SCALES="sum,theta" adds the theta cut (about twice the
# time). DIF_3C_REPS sets the replications (default 500).

suppressMessages(devtools::load_all(".", quiet = TRUE))  # package internals
source("dev/dif_impact_ext_sim.R")  # grid_scenarios(), sim_twin(), TARGET, ...

TARGET_PKG <- list(n = TARGET$n, mean = TARGET$mu, sd = TARGET$sd)

# Logical rule (keep / split, ref / focal, over raw scores 0..R) from the
# lowest raw score classified at or above the cut.
rule_from_cuts <- function(cuts, R) {
  r <- 0:R
  out <- list()
  for (st in c("keep", "split")) for (g in c("ref", "focal")) {
    out[[st]][[g]] <- r >= cuts[[paste(st, g, sep = ".")]]
  }
  out
}

truth_3c <- function(cell) {
  sc <- make_scale(cell$type, cell$k)
  si <- scenario_items(sc$loc, cell$scenario, cell$position, cell$cut)
  items <- sort(si$j)
  d <- (cell$d * si$sign)[order(si$j)]
  thr <- .dimp_with_shifts(list(ref = sc$thr, focal = sc$thr), items, d)
  cut_raw <- cut_raw_for(cell)
  pt <- .dimp_rate_parts(thr, items, cut_raw)
  list(thr = thr, items = items, cut_raw = cut_raw, pt = pt,
       mu = c(0, cell$mu2), sd = c(1, 1),
       # the package's sum rule at the true parameters
       rule_cuts_sum = .dimp_rule_cuts(.dimp_rule(pt, NULL, "sum", cut_raw)))
}

rep_3c <- function(cell, tr, scales, rep_id) {
  sc <- make_scale(cell$type, cell$k)
  si <- scenario_items(sc$loc, cell$scenario, cell$position, cell$cut)
  dif <- stats::setNames(cell$d * si$sign, names(sc$thr)[si$j])
  tw <- sim_twin(sc$thr, dif, cell$n_ref, cell$n_focal, mu_focal = cell$mu2)
  X <- tw$X_dif
  storage.mode(X) <- "integer"
  focal <- tw$g == "focal"
  do.call(rbind, lapply(scales, function(scale) {
    cut <- if (scale == "sum") tr$cut_raw else cell$cut
    e <- .dimp_estimate(X, focal, tr$items, cut, scale, TARGET_PKG)
    rule <- rule_from_cuts(e$rates$rule, e$R)
    th_c <- if (scale == "sum") tr$pt$theta_c_sum else cell$cut
    tg <- list(
      sample = .dimp_fixed_rates(tr$pt, rule, tr$mu, tr$sd, th_c),
      projection = .dimp_fixed_rates(tr$pt, rule, TARGET$mu, TARGET$sd, th_c)
    )
    rc <- e$rates$rule
    rows <- lapply(c("sample", "projection"), function(ctx) {
      x <- e$rates[[if (ctx == "sample") "sample" else "target"]]
      data.frame(scale = scale, context = ctx, quantity = names(x$est),
                 estimate = unname(x$est), se = unname(x$se),
                 target = unname(tg[[ctx]][names(x$est)]))
    })
    rows[[3L]] <- data.frame(scale = scale, context = "rule",
                             quantity = names(rc), estimate = unname(rc),
                             se = NA_real_,
                             target = if (scale == "sum") {
                               unname(tr$rule_cuts_sum[names(rc)])
                             } else NA_real_)
    out <- do.call(rbind, rows)
    out$rule_differs_ref <- rc[["keep.ref"]] != rc[["split.ref"]]
    out$rule_differs_focal <- rc[["keep.focal"]] != rc[["split.focal"]]
    out
  }))
}

run_cell_3c <- function(cell, reps, seed, workers, scales) {
  tr <- truth_3c(cell)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  streams <- vector("list", reps)
  s <- .Random.seed
  for (r in seq_len(reps)) {
    s <- parallel::nextRNGStream(s)
    streams[[r]] <- s
  }
  out <- parallel::mclapply(seq_len(reps), function(r) {
    assign(".Random.seed", streams[[r]], envir = globalenv())
    n_warn <- 0L
    tryCatch({
      x <- withCallingHandlers(
        rep_3c(cell, tr, scales, r),
        warning = function(w) {
          # e.g. a convergence warning: counted per replication, run goes on
          n_warn <<- n_warn + 1L
          invokeRestart("muffleWarning")
        }
      )
      cbind(rep = r, x, warnings = n_warn)
    }, error = function(e) list(error = conditionMessage(e), rep = r))
  }, mc.cores = workers, mc.preschedule = FALSE)
  err <- vapply(out, function(x) is.list(x) && !is.data.frame(x), TRUE)
  list(cell = cell, truth = tr[c("cut_raw", "rule_cuts_sum", "items")],
       estimates = do.call(rbind, out[!err]),
       errors = if (any(err)) vapply(out[err], `[[`, "", "error") else character(0))
}

# Coverage and SE ratio per quantity, over all replications and, for
# keep - split, over the replications where the group's rule differs.
summarise_3c <- function(res) {
  e <- res$estimates
  e <- e[e$context != "rule", ]
  one <- function(x, subset) {
    err <- x$estimate - x$target
    data.frame(
      scale = x$scale[1], context = x$context[1], quantity = x$quantity[1],
      subset = subset, reps = nrow(x),
      mean_target = mean(x$target), bias = mean(err), sd_err = stats::sd(err),
      se_ratio = mean(x$se) / stats::sd(err),
      coverage = mean(abs(err) <= 1.96 * x$se + 1e-12),
      share_zero = mean(x$estimate == 0 & x$target == 0),
      stringsAsFactors = FALSE)
  }
  key <- paste(e$scale, e$context, e$quantity, sep = "|")
  all_rows <- lapply(split(e, key), one, subset = "all")
  dz <- e[grepl("^diff\\.", e$quantity), ]
  differs <- ifelse(grepl("^diff\\.ref", dz$quantity), dz$rule_differs_ref,
                    dz$rule_differs_focal)
  dz <- dz[differs, ]
  nz_rows <- if (nrow(dz) > 1L) {
    k2 <- paste(dz$scale, dz$context, dz$quantity, sep = "|")
    lapply(split(dz, k2), function(x) {
      if (nrow(x) < 2L) return(NULL)
      one(x, "rule_differs")
    })
  }
  out <- do.call(rbind, c(all_rows, nz_rows))
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  smoke <- identical(Sys.getenv("DIF_3C_SMOKE"), "1")
  workers <- as.integer(Sys.getenv("DIF_3C_WORKERS", "4"))
  reps <- as.integer(Sys.getenv("DIF_3C_REPS", "500"))
  scales <- strsplit(Sys.getenv("DIF_3C_SCALES", "sum"), ",")[[1]]
  stopifnot(all(scales %in% c("sum", "theta")))
  cells <- grid_scenarios()
  rownames(cells) <- NULL
  ids <- seq_len(nrow(cells))
  if (smoke) {
    # k = 5 one item, k = 10 two items same direction (rule moves), opposite
    ids <- c(which(cells$k == 5 & cells$d == 0.8)[1],
             which(cells$scenario == "two_same" & cells$d == 0.8)[1],
             which(cells$scenario == "two_opposite")[1])
    reps <- as.integer(Sys.getenv("DIF_3C_SMOKE_REPS", "8"))
  }
  stem <- if (smoke) "ext3c_smoke" else "ext3c_results"
  cell_dir <- file.path(ext_dir, paste0(stem, "_cells"))
  dir.create(cell_dir, showWarnings = FALSE)
  message(sprintf("%d cells, %d replications, %d workers, scales: %s",
                  length(ids), reps, workers, paste(scales, collapse = ",")))
  t0 <- Sys.time()
  results <- vector("list", length(ids))
  for (ii in seq_along(ids)) {
    i <- ids[ii]
    f_cell <- file.path(cell_dir, sprintf("cell_%02d.rds", i))
    if (file.exists(f_cell)) {
      results[[ii]] <- readRDS(f_cell)
      message(sprintf("cell %d/%d read from %s", i, nrow(cells), f_cell))
      next
    }
    cell <- as.list(cells[i, ])
    res <- run_cell_3c(cell, reps, seed = 20261140L + i, workers = workers,
                       scales = scales)
    res$summary <- if (!is.null(res$estimates)) summarise_3c(res) else NULL
    saveRDS(res, f_cell)
    results[[ii]] <- res
    message(sprintf("cell %d/%d done (%s), errors: %d", i, nrow(cells),
                    format(Sys.time() - t0, digits = 3), length(res$errors)))
    if (length(res$errors)) print(utils::head(unique(res$errors)))
  }
  out <- list(results = results, cells = cells[ids, ], reps = reps,
              scales = scales, target_population = TARGET, workers = workers,
              run_time = Sys.time() - t0, date = Sys.Date(),
              package_version = as.character(utils::packageVersion("easyRasch2")),
              session = utils::sessionInfo())
  saveRDS(out, file.path(ext_dir, paste0(stem, ".rds")))
  saveRDS(out, file.path(ext_dir, sprintf("%s_%s.rds", stem,
                                          format(Sys.Date(), "%Y-%m-%d"))))
  message("Done in ", format(out$run_time))
}
