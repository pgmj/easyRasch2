#' Martin-Lof Test of Unidimensionality
#'
#' Likelihood-ratio test of unidimensionality against an *a priori* specified
#' multidimensional alternative, generalised to polytomous Rasch / partial
#' credit models (Christensen, Bjorner, Kreiner, & Petersen, 2002). The
#' p-value is obtained by parametric-bootstrap (Monte Carlo) sampling under
#' the unidimensional null, following Christensen & Kreiner (2007), because
#' the asymptotic chi-square approximation is biased toward conservatism for
#' realistic sample sizes -- especially with polytomous items, where the
#' degrees of freedom can be very large.
#'
#' This is **not** a routine screening tool. The test requires an *a priori*
#' partition of items into subscales; using it post-hoc on, e.g., the
#' partition suggested by `RMdimResidualPCA()`'s PC1 sign would inflate the
#' Type-I error rate. Both source papers state this explicitly.
#'
#' @param data A data.frame or matrix of item responses (0-based,
#'   non-negative integers). Rows with any `NA` are dropped.
#' @param partition The hypothesised partition of items into subscales. One of:
#'   * a **list** of column-name or column-index vectors, e.g.
#'     `list(c("I1","I2","I3"), c("I4","I5","I6"))`;
#'   * a **vector** of length `ncol(data)` indicating each item's subscale
#'     (factor, character, or integer), e.g. `c(1,1,1,2,2,2)`.
#'   Each subscale must contain at least two items. Subscales must not
#'   overlap; items not assigned to any subscale are dropped with a warning.
#' @param iterations Integer. Maximum number of Monte Carlo iterations
#'   (default `1000`).
#' @param stopping Character. `"none"` (default) runs all iterations.
#'   `"sequential"` uses Besag & Clifford's (1991) sequential rule: stop as
#'   soon as `h` simulated statistics have exceeded the observed value. The
#'   sequential strategy substantially reduces compute time when H0 holds
#'   but cannot be parallelised.
#' @param h Integer. Sequential-stopping count threshold (default `50`).
#'   Ignored when `stopping = "none"`.
#' @param alpha Numeric in (0, 1). Nominal significance level used only for
#'   the `rejected` flag in the result; default `0.05`.
#' @param parallel Logical. Use parallel processing via `mirai` (default
#'   `TRUE`). Ignored when `stopping = "sequential"`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first; if neither is set, falls back
#'   to sequential with a warning.
#' @param verbose Logical. Show a progress bar (default `FALSE`).
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{`T_obs`}{Observed Martin-Lof likelihood-ratio statistic.}
#'   \item{`p_value`}{Monte Carlo p-value with `(n_exceed + 1) / (n + 1)`
#'     correction.}
#'   \item{`actual_iterations`}{Number of successful MC iterations
#'     completed.}
#'   \item{`rejected`}{Logical: is `p_value < alpha`?}
#'   \item{`partition`}{Normalised partition (list of integer indices).}
#'   \item{`n_subscales`}{Number of subscales.}
#'   \item{`is_polytomous`}{Whether a PCM was fitted.}
#'   \item{`sample_n`}{Number of complete cases analysed.}
#'   \item{`n_items`}{Number of items.}
#'   \item{`stopping`}{The stopping strategy used.}
#'   \item{`h`}{The sequential-stopping count, or `NA` for `stopping = "none"`.}
#'   \item{`T_rep`}{Numeric vector of successful MC test statistics.}
#'   \item{`wle_scores`}{data.frame with one row per person and one column
#'     per subscale (`subscale_1_wle`, ..., `subscale_D_wle`), giving Warm's
#'     Weighted Likelihood Estimate of theta from a CML fit on each subscale
#'     alone. Persons whose subscore equals the minimum or maximum on a
#'     subscale produce non-finite WLEs (`Inf` / `-Inf`) and are excluded
#'     from `wle_correlation` pairwise.}
#'   \item{`wle_correlation`}{data.frame of pairwise Pearson correlations
#'     between subscale WLEs, with columns `subscale_a`, `subscale_b`, `r`,
#'     `ci_lower`, `ci_upper` (95% CI from `stats::cor.test`), `p_value`,
#'     and `n` (number of persons with finite WLEs on both subscales). One
#'     row per pair; for D = 2, a single row. Useful as an effect-size
#'     companion to `p_value` -- a rejected test with `r` near 1 indicates a
#'     small effect; `r` clearly below 1 indicates substantive
#'     multidimensionality.}
#' }
#'
#' @details
#' **Test statistic.** With items partitioned into D subscales, total score
#' \eqn{t} and subscores \eqn{(t_1, \ldots, t_D)} (Christensen et al. 2002,
#' eq. 22):
#' \deqn{T = 2\Bigl[\sum_{t_1, \ldots, t_D}
#'   n_{t_1, \ldots, t_D}\log(n_{t_1, \ldots, t_D}/N)
#'   - \sum_t n_t\log(n_t/N)
#'   - \ell_C(\hat{\epsilon}) + \sum_d \ell_C(\hat{\epsilon}^{(d)})\Bigr]}
#' where \eqn{\ell_C} is the conditional log-likelihood and the
#' \eqn{\hat{\epsilon}^{(d)}} are CML estimates on the d-th subscale alone.
#' CML fits use `psychotools::raschmodel()` (RM) or `psychotools::pcmodel()`
#' (PCM) for speed.
#'
#' **Monte Carlo sampling under H0.** Following Christensen & Kreiner (2007):
#' (a) sample N total scores from the empirical score distribution
#' \eqn{n_t/N}; (b) for each sampled score, sample an item-response vector
#' from the conditional distribution \eqn{p(x \mid t, \hat{\epsilon})} given
#' by eq. 4 of the paper. For dichotomous data the fast algorithm of
#' Christensen & Kreiner (2007, p. 23) is used (sample without replacement
#' weighted by item easinesses). For polytomous data the recursive
#' \eqn{\gamma}-function approach is used, with each item's response sampled
#' conditional on the remaining items' joint score distribution (computed via
#' `psychotools::elementary_symmetric_functions()`).
#'
#' Iterations that fail (e.g., simulated dataset has an empty category for an
#' item) are silently dropped.
#'
#' Item parameters are estimated once on the observed data and held fixed
#' across MC iterations. Christensen & Kreiner (2007) use the extended
#' likelihood function (Tjur, 1982) with the empirical score distribution as
#' a non-parametric estimate of the latent distribution, so no distributional
#' assumption about \eqn{\theta} is needed.
#'
#' @references
#' Christensen, K. B., Bjorner, J. B., Kreiner, S., & Petersen, J. H. (2002).
#' Testing unidimensionality in polytomous Rasch models. *Psychometrika,
#' 67*(4), 563-574. \doi{10.1007/BF02295132}
#'
#' Christensen, K. B., & Kreiner, S. (2007). A Monte Carlo approach to
#' unidimensionality testing in polytomous Rasch models. *Applied
#' Psychological Measurement, 31*(1), 20-30.
#' \doi{10.1177/0146621605286204}
#'
#' Besag, J., & Clifford, P. (1991). Sequential Monte Carlo p-values.
#' *Biometrika, 78*(2), 301-304. \doi{10.1093/biomet/78.2.301}
#'
#' @seealso \code{\link{RMdimResidualPCA}}, \code{\link{RMdimResidualPCACutoff}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' # Build 2-dimensional polytomous data: 4 items per subscale, 5 categories
#' n     <- 400
#' theta1 <- rnorm(n)
#' theta2 <- 0.6 * theta1 + sqrt(1 - 0.6^2) * rnorm(n)
#' make_pcm <- function(theta, n_items, taus) {
#'   sapply(seq_len(n_items), function(j) {
#'     # ... toy simulation here
#'     sample(0:4, n, replace = TRUE)
#'   })
#' }
#' dat <- cbind(make_pcm(theta1, 4, NULL), make_pcm(theta2, 4, NULL))
#' colnames(dat) <- paste0("I", 1:8)
#'
#' RMdimMartinLof(dat,
#'             partition = list(c("I1","I2","I3","I4"),
#'                              c("I5","I6","I7","I8")),
#'             iterations = 200, parallel = FALSE, seed = 1)
#'
#' # Sequential stopping: stop as soon as h = 50 simulated statistics exceed
#' # the observed one (cuts compute time under H0).
#' RMdimMartinLof(dat,
#'             partition = c(1,1,1,1,2,2,2,2),
#'             iterations = 1000, stopping = "sequential", h = 50,
#'             seed = 1)
#' }
RMdimMartinLof <- function(data,
                        partition,
                        iterations = 1000L,
                        stopping   = c("none", "sequential"),
                        h          = 50L,
                        alpha      = 0.05,
                        parallel   = TRUE,
                        n_cores    = NULL,
                        verbose    = FALSE,
                        seed       = NULL) {

  stopping <- match.arg(stopping)

  validate_response_data(data)

  data <- stats::na.omit(as.data.frame(data))
  if (nrow(data) < 30L) {
    stop("Need at least 30 complete cases.", call. = FALSE)
  }

  # Normalise partition
  partition_list <- normalize_ml_partition(partition, data)
  D <- length(partition_list)
  if (D < 2L) {
    stop("Partition must specify >= 2 subscales.", call. = FALSE)
  }
  if (any(lengths(partition_list) < 2L)) {
    stop("Each subscale must contain at least 2 items.", call. = FALSE)
  }

  # Restrict data to items used in the partition (drop unassigned items)
  used_idx <- sort(unique(unlist(partition_list)))
  if (length(used_idx) < ncol(data)) {
    data <- data[, used_idx, drop = FALSE]
    # Re-index partition relative to the restricted column set
    remap <- stats::setNames(seq_along(used_idx), used_idx)
    partition_list <- lapply(partition_list, function(idx) unname(remap[as.character(idx)]))
  }

  data_mat      <- as.matrix(data)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L
  N             <- nrow(data)

  # rgl workaround (in case any helper triggers it)
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # Item parameters for sampling under H0 (fixed at observed-data CML estimates)
  sampling_params <- extract_ml_sampling_params(data, is_polytomous)

  # Observed test statistic (also returns log-likelihood components)
  T_obs <- compute_ml_statistic(data, partition_list, is_polytomous)

  # ---- Per-subscale WLE thetas + pairwise correlations ---------------------
  # Diagnostic to interpret rejection / non-rejection: under H0 the WLE
  # estimates from each subscale should be highly correlated; clearly < 1
  # values point to substantive multidimensionality.
  wle_compute <- compute_subscale_wle_and_correlations(
    data, partition_list, is_polytomous
  )
  wle_scores      <- wle_compute$wle_scores
  wle_correlation <- wle_compute$wle_correlation

  # Empirical total-score distribution for sampling
  total_scores <- rowSums(data)
  score_table  <- table(total_scores)
  score_values <- as.integer(names(score_table))
  score_probs  <- as.numeric(score_table) / N

  # Per-iteration seed list
  if (!is.null(seed)) set.seed(seed)
  sim_seeds <- sample.int(.Machine$integer.max, iterations)

  sim_data_list <- list(
    N               = N,
    is_polytomous   = is_polytomous,
    sampling_params = sampling_params,
    score_values    = score_values,
    score_probs     = score_probs,
    item_names      = colnames(data),
    partition_list  = partition_list
  )

  # ---- Run iterations -------------------------------------------------------
  if (stopping == "none") {
    use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)
    if (parallel && !use_parallel) {
      message("Install 'mirai' for parallel processing: install.packages(\"mirai\")")
      message("Running sequentially...")
    }
    if (use_parallel) {
      if (is.null(n_cores)) n_cores <- getOption("mc.cores")
      if (is.null(n_cores)) {
        warning(
          "For parallel processing, specify n_cores or set options(mc.cores = N).\n",
          "Falling back to sequential.",
          call. = FALSE
        )
        use_parallel <- FALSE
      } else {
        n_cores <- min(n_cores, iterations)
      }
    }

    results_raw <- if (use_parallel) {
      run_ml_sim_parallel(iterations, sim_seeds, sim_data_list, n_cores, verbose)
    } else {
      run_ml_sim_sequential(iterations, sim_seeds, sim_data_list, verbose)
    }

    T_rep <- vapply(results_raw,
                    function(x) if (is.numeric(x) && length(x) == 1L) as.numeric(x) else NA_real_,
                    numeric(1L))
    T_rep <- T_rep[is.finite(T_rep)]
    actual_iterations <- length(T_rep)
    if (actual_iterations < 2L) {
      stop("Fewer than 2 successful Monte Carlo iterations; cannot estimate p-value.",
           call. = FALSE)
    }
    n_exceed <- sum(T_rep >= T_obs)
    p_value  <- (n_exceed + 1) / (actual_iterations + 1)

  } else {
    # Sequential stopping (Besag & Clifford 1991)
    if (verbose) pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
    n_exceed <- 0L
    actual_iterations <- 0L
    T_rep <- numeric(0)
    for (i in seq_len(iterations)) {
      result <- run_single_ml_iteration(sim_seeds[i], sim_data_list)
      if (is.numeric(result) && length(result) == 1L && is.finite(result)) {
        actual_iterations <- actual_iterations + 1L
        T_rep <- c(T_rep, result)
        if (result >= T_obs) {
          n_exceed <- n_exceed + 1L
          if (n_exceed >= h) break
        }
      }
      if (verbose) utils::setTxtProgressBar(pb, i)
    }
    if (verbose) { close(pb); message("") }

    if (actual_iterations < 2L) {
      stop("Fewer than 2 successful Monte Carlo iterations; cannot estimate p-value.",
           call. = FALSE)
    }
    p_value <- (n_exceed + 1) / (actual_iterations + 1)
  }

  list(
    T_obs             = T_obs,
    p_value           = p_value,
    actual_iterations = actual_iterations,
    rejected          = p_value < alpha,
    partition         = partition_list,
    n_subscales       = D,
    is_polytomous     = is_polytomous,
    sample_n          = N,
    n_items           = ncol(data),
    stopping          = stopping,
    h                 = if (stopping == "sequential") h else NA_integer_,
    T_rep             = T_rep,
    wle_scores        = wle_scores,
    wle_correlation   = wle_correlation
  )
}

# ===========================================================================
# Internal helpers
# ===========================================================================

#' Normalise a partition argument to a list of integer column indices
#'
#' @keywords internal
#' @noRd
normalize_ml_partition <- function(partition, data) {
  n_items    <- ncol(data)
  item_names <- colnames(data)

  if (is.list(partition)) {
    partition_list <- lapply(partition, function(p) {
      if (is.numeric(p)) {
        idx <- as.integer(p)
        if (any(idx < 1L) || any(idx > n_items)) {
          stop("Partition indices out of range.", call. = FALSE)
        }
        idx
      } else if (is.character(p)) {
        idx <- match(p, item_names)
        if (any(is.na(idx))) {
          missing <- p[is.na(idx)]
          stop("Items not found in data: ", paste(missing, collapse = ", "),
               call. = FALSE)
        }
        idx
      } else {
        stop("Each partition list element must be numeric or character.",
             call. = FALSE)
      }
    })
  } else if (is.factor(partition) || is.character(partition) || is.numeric(partition)) {
    if (length(partition) != n_items) {
      stop("Partition vector must have length equal to ncol(data) (",
           n_items, "); got ", length(partition), ".", call. = FALSE)
    }
    partition_list <- split(seq_len(n_items), as.character(partition))
    partition_list <- partition_list[lengths(partition_list) > 0L]
    names(partition_list) <- NULL
  } else {
    stop("Invalid partition format. Provide a list of column vectors or a ",
         "vector of length ncol(data).", call. = FALSE)
  }

  used_items <- unlist(partition_list)
  if (length(used_items) != length(unique(used_items))) {
    stop("Partition has overlapping items.", call. = FALSE)
  }
  if (length(used_items) < n_items) {
    not_used <- setdiff(seq_len(n_items), used_items)
    warning("Items not assigned to any subscale (dropped): ",
            paste(item_names[not_used], collapse = ", "), call. = FALSE)
  }

  partition_list
}

#' Extract item parameters for conditional sampling under H0
#'
#' For dichotomous data, returns a numeric vector of \eqn{\epsilon_{i,1} =
#' -\beta_i} (item easinesses). For polytomous data, returns a list with one
#' element per item: a numeric vector \eqn{(\epsilon_{i,1}, \ldots,
#' \epsilon_{i,m_i})} with \eqn{\epsilon_{i,x} = -\sum_{k=1}^x \tau_{i,k}}
#' (cumulative negated Andrich thresholds), suitable for direct use with
#' `psychotools::elementary_symmetric_functions()`.
#'
#' @keywords internal
#' @noRd
extract_ml_sampling_params <- function(data, is_polytomous) {
  if (is_polytomous) {
    thresh_mat <- extract_item_thresholds(as.matrix(data))  # from utils-simulation.R
    params_list <- lapply(seq_len(nrow(thresh_mat)), function(i) {
      taus <- as.numeric(thresh_mat[i, !is.na(thresh_mat[i, ])])
      -cumsum(taus)
    })
    names(params_list) <- rownames(thresh_mat)
    params_list
  } else {
    rasch_fit <- psychotools::raschmodel(data, hessian = FALSE)
    # `coef(raschmodel)` returns K-1 parameters under the sum-to-zero
    # constraint, but the sampler needs all K item parameters. Use
    # `itempar()` which returns the full K-vector with proper names.
    -as.numeric(psychotools::itempar(rasch_fit))
  }
}

#' Compute the Martin-Lof likelihood-ratio statistic
#'
#' @keywords internal
#' @noRd
compute_ml_statistic <- function(data, partition_list, is_polytomous) {
  N <- nrow(data)

  # Full-scale CML log-likelihood
  full_fit <- if (is_polytomous) {
    psychotools::pcmodel(data, hessian = FALSE)
  } else {
    psychotools::raschmodel(data, hessian = FALSE)
  }
  ll_full <- as.numeric(stats::logLik(full_fit))

  # Sum of subscale CML log-likelihoods
  ll_subs <- 0
  for (idx in partition_list) {
    sub_data <- data[, idx, drop = FALSE]
    sub_fit <- if (is_polytomous) {
      psychotools::pcmodel(sub_data, hessian = FALSE)
    } else {
      psychotools::raschmodel(sub_data, hessian = FALSE)
    }
    ll_subs <- ll_subs + as.numeric(stats::logLik(sub_fit))
  }

  # Score and joint subscore frequencies
  total_scores <- rowSums(data)
  n_t  <- as.numeric(table(total_scores))

  subscores <- vapply(partition_list,
                      function(idx) rowSums(data[, idx, drop = FALSE]),
                      numeric(N))
  if (is.matrix(subscores)) {
    subscore_keys <- apply(subscores, 1L, function(x) paste(x, collapse = ","))
  } else {
    subscore_keys <- as.character(subscores)
  }
  n_combo <- as.numeric(table(subscore_keys))

  # 0 log 0 = 0 (and table() never returns 0, so safe to use logs directly)
  log_term_combo <- sum(n_combo * log(n_combo / N))
  log_term_total <- sum(n_t    * log(n_t    / N))

  2 * (log_term_combo - log_term_total - ll_full + ll_subs)
}

#' Per-subscale WLE thetas and pairwise correlations
#'
#' For each subscale, fits a CML Rasch / PCM via `psychotools::raschmodel()`
#' / `psychotools::pcmodel()` and extracts Warm's Weighted Likelihood
#' Estimate of theta per person via the patched
#' iarm helper that handles boundary scores cleanly. Returns the WLE
#' matrix and pairwise Pearson correlations with 95% CIs. Robust to
#' subscale fit failures (returns NA WLEs for that subscale).
#'
#' @keywords internal
#' @noRd
compute_subscale_wle_and_correlations <- function(data, partition_list,
                                                  is_polytomous) {
  D <- length(partition_list)
  N <- nrow(data)

  wle_per_subscale <- lapply(seq_along(partition_list), function(d) {
    idx <- partition_list[[d]]
    sub_data <- data[, idx, drop = FALSE]
    tryCatch({
      sub_fit <- if (is_polytomous) {
        psychotools::pcmodel(sub_data, hessian = FALSE)
      } else {
        psychotools::raschmodel(sub_data, hessian = FALSE)
      }
      pe <- iarm_person_estimates(sub_fit, allperson = TRUE)
      as.numeric(pe[, "WLE"])
    }, error = function(e) {
      warning(sprintf("WLE estimation failed for subscale %d: %s",
                      d, conditionMessage(e)), call. = FALSE)
      rep(NA_real_, N)
    })
  })

  wle_scores <- as.data.frame(do.call(cbind, wle_per_subscale))
  names(wle_scores) <- paste0("subscale_", seq_len(D), "_wle")
  rownames(wle_scores) <- NULL

  # Pairwise correlations
  pairs_idx <- if (D == 2L) {
    matrix(c(1L, 2L), ncol = 2L)
  } else {
    t(utils::combn(D, 2L))
  }

  rows <- lapply(seq_len(nrow(pairs_idx)), function(i) {
    a <- pairs_idx[i, 1L]
    b <- pairs_idx[i, 2L]
    x <- wle_scores[[a]]
    y <- wle_scores[[b]]
    finite <- is.finite(x) & is.finite(y)
    if (sum(finite) < 3L) {
      return(data.frame(
        subscale_a = a, subscale_b = b,
        r          = NA_real_,
        ci_lower   = NA_real_,
        ci_upper   = NA_real_,
        p_value    = NA_real_,
        n          = sum(finite),
        stringsAsFactors = FALSE,
        row.names = NULL
      ))
    }
    test <- tryCatch(
      stats::cor.test(x[finite], y[finite], method = "pearson"),
      error = function(e) NULL
    )
    if (is.null(test)) {
      data.frame(
        subscale_a = a, subscale_b = b,
        r          = NA_real_,
        ci_lower   = NA_real_,
        ci_upper   = NA_real_,
        p_value    = NA_real_,
        n          = sum(finite),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    } else {
      data.frame(
        subscale_a = a, subscale_b = b,
        r          = round(unname(test$estimate), 3),
        ci_lower   = round(test$conf.int[1L], 3),
        ci_upper   = round(test$conf.int[2L], 3),
        p_value    = signif(test$p.value, 4),
        n          = sum(finite),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  })
  wle_correlation <- do.call(rbind, rows)
  rownames(wle_correlation) <- NULL

  list(wle_scores = wle_scores, wle_correlation = wle_correlation)
}

#' Sample a dichotomous response vector with target total score
#'
#' Uses the fast algorithm of Christensen & Kreiner (2007, p. 23): pick `s`
#' items without replacement, weighted by `exp(item_params)`.
#'
#' @keywords internal
#' @noRd
sample_dichotomous_at_score <- function(s, item_params) {
  n_items <- length(item_params)
  if (s == 0L)       return(integer(n_items))
  if (s == n_items)  return(rep(1L, n_items))
  weights <- exp(item_params)
  picked <- sample.int(n_items, size = s, prob = weights, replace = FALSE)
  out <- integer(n_items)
  out[picked] <- 1L
  out
}

#' Sample a polytomous response vector with target total score
#'
#' Iteratively samples one item's response at a time. For each item j in the
#' remaining pool, the conditional probability of \eqn{x_j} given the
#' remaining target score \eqn{t'} is proportional to
#' \eqn{\exp(\epsilon_{j,x_j}) \cdot \gamma^{(j)}_{t' - x_j}} where
#' \eqn{\gamma^{(j)}} is the elementary symmetric function over the
#' remaining items excluding item j (Andersen, 1995, eq. 15.22, 15.26).
#' \eqn{\gamma^{(j)}} is recomputed via
#' `psychotools::elementary_symmetric_functions()` for each step.
#'
#' @keywords internal
#' @noRd
sample_polytomous_at_score <- function(t, params_list) {
  n_items   <- length(params_list)
  m_i       <- vapply(params_list, length, integer(1L))   # max category per item
  M_total   <- sum(m_i)

  if (t == 0L)        return(integer(n_items))
  if (t == M_total)   return(m_i)

  x <- integer(n_items)
  remaining_pos    <- seq_len(n_items)        # original positions still in pool
  remaining_params <- params_list
  remaining_max    <- m_i
  remaining_score  <- t

  while (length(remaining_pos) > 0L) {
    if (length(remaining_pos) == 1L) {
      x[remaining_pos] <- remaining_score
      break
    }
    j <- 1L
    item_pos <- remaining_pos[j]
    other_params <- remaining_params[-j]
    M_other      <- sum(remaining_max[-j])

    # gamma^{(j)}: ESF over the items remaining minus item j
    gamma_other <- if (length(other_params) == 0L) {
      1
    } else {
      esf <- psychotools::elementary_symmetric_functions(other_params, order = 0L)
      if (is.list(esf)) esf[[1L]] else esf
    }

    # Conditional log-probabilities for x_j in 0:m_j
    eps_j <- c(0, remaining_params[[j]])  # epsilon_{j,0}=0 prepended
    max_xj <- remaining_max[j]
    log_probs <- rep(-Inf, max_xj + 1L)
    for (xj in 0L:max_xj) {
      target_other <- remaining_score - xj
      if (target_other >= 0L && target_other <= M_other) {
        log_probs[xj + 1L] <- eps_j[xj + 1L] + log(gamma_other[target_other + 1L])
      }
    }
    log_probs <- log_probs - max(log_probs)
    probs <- exp(log_probs)
    if (sum(probs) <= 0 || !is.finite(sum(probs))) {
      # Numerical fallback: should be rare. Pick the most likely value.
      sampled_xj <- which.max(log_probs) - 1L
    } else {
      probs <- probs / sum(probs)
      sampled_xj <- sample(0L:max_xj, size = 1L, prob = probs)
    }
    x[item_pos] <- sampled_xj
    remaining_score  <- remaining_score - sampled_xj
    remaining_pos    <- remaining_pos[-j]
    remaining_params <- remaining_params[-j]
    remaining_max    <- remaining_max[-j]
  }
  x
}

#' Run a single Monte Carlo iteration of the Martin-Lof test
#'
#' @keywords internal
#' @noRd
run_single_ml_iteration <- function(seed, sim_data) {
  set.seed(seed)

  N <- sim_data$N
  scores <- sample(sim_data$score_values, size = N, replace = TRUE,
                   prob = sim_data$score_probs)

  sim_responses <- matrix(0L, nrow = N, ncol = length(sim_data$item_names))
  if (sim_data$is_polytomous) {
    for (i in seq_len(N)) {
      sim_responses[i, ] <- sample_polytomous_at_score(scores[i],
                                                       sim_data$sampling_params)
    }
  } else {
    for (i in seq_len(N)) {
      sim_responses[i, ] <- sample_dichotomous_at_score(scores[i],
                                                        sim_data$sampling_params)
    }
  }

  sim_df <- as.data.frame(sim_responses)
  colnames(sim_df) <- sim_data$item_names

  tryCatch({
    compute_ml_statistic(sim_df, sim_data$partition_list, sim_data$is_polytomous)
  }, error = function(e) NA_real_)
}

#' Parallel runner for ML iterations (mirai)
#'
#' @keywords internal
#' @noRd
run_ml_sim_parallel <- function(iterations, sim_seeds, sim_data_list,
                                n_cores, verbose = FALSE) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  tasks <- lapply(seq_len(iterations), function(i) {
    mirai::mirai(
      { run_single_ml_iteration(seed, data_list) },
      seed                          = sim_seeds[i],
      data_list                     = sim_data_list,
      run_single_ml_iteration       = run_single_ml_iteration,
      compute_ml_statistic          = compute_ml_statistic,
      sample_dichotomous_at_score   = sample_dichotomous_at_score,
      sample_polytomous_at_score    = sample_polytomous_at_score
    )
  })

  results <- vector("list", iterations)
  for (i in seq_len(iterations)) {
    res <- mirai::call_mirai(tasks[[i]])$data
    results[[i]] <- if (inherits(res, "errorValue")) NA_real_ else res
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) { close(pb); message("") }
  results
}

#' Sequential runner for ML iterations
#'
#' @keywords internal
#' @noRd
run_ml_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                  verbose = FALSE) {
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  results <- vector("list", iterations)
  for (i in seq_len(iterations)) {
    results[[i]] <- run_single_ml_iteration(sim_seeds[i], sim_data_list)
    if (verbose) utils::setTxtProgressBar(pb, i)
  }
  if (verbose) { close(pb); message("") }
  results
}

# ===========================================================================
# RMdimMartinLofResiduals: standardised residuals from the joint subscore table
# ===========================================================================

#' Standardised Residuals from the Joint Subscore Distribution
#'
#' Diagnostic accompanying \code{\link{RMdimMartinLof}}: per-cell standardised
#' residuals from the joint distribution of subscores under unidimensionality
#' (Christensen, Bjorner, Kreiner, & Petersen, 2002, eq. 13). Useful for
#' identifying *where* a partition deviates from the unidimensional null
#' rather than just whether it does (which `RMdimMartinLof()` answers).
#'
#' For each cell of the joint subscore table (indexed by
#' \eqn{(t_1, \ldots, t_D)}), the conditional probability under H0 given the
#' total score \eqn{t = \sum_d t_d} is
#' \deqn{p(t_1, \ldots, t_D \mid t) = \prod_d \gamma^{(d)}_{t_d} / \gamma_t,}
#' the expected count is \eqn{e = n_t \cdot p}, and the residual is
#' \eqn{(o - e) / \sqrt{n_t \cdot p \cdot (1 - p)}}. CML estimates from the
#' unidimensional model are used for the \eqn{\gamma}-functions.
#'
#' Reading the table (D = 2): under the unidimensional null, residuals
#' should be patternless and roughly N(0, 1). Multidimensionality with
#' positively correlated dimensions typically shows up as **positive**
#' residuals at the corners of each antidiagonal (high `t1` + low `t2`,
#' low `t1` + high `t2`) and **negative** residuals near the antidiagonal
#' centre (matched subscores). Negatively correlated dimensions show
#' positive residuals at the table corners (high/low and low/high) and
#' negative residuals at high/high and low/low. See Christensen et al.
#' (2002, section7) for a worked example.
#'
#' Cells where the total score has no observed cases (`n_t = 0`) are
#' uninformative and are dropped from the output.
#'
#' @param data A data.frame or matrix of item responses (0-based,
#'   non-negative integers). Rows with any `NA` are dropped.
#' @param partition Same format as in \code{\link{RMdimMartinLof}}: a list of
#'   item-name/index vectors, or a length-`ncol(data)` vector of group
#'   labels. Each subscale must contain at least 2 items.
#' @param output Character. `"kable"` (default) for a 2-D pipe-format
#'   residual table when D = 2 (long-format kable when D > 2), `"dataframe"`
#'   for the underlying long-format data.frame, or `"ggplot"` for a
#'   diverging-fill heatmap (with `facet_wrap` over `t3` for D = 3, error
#'   for D > 3).
#' @param flag_threshold Numeric. Cells with `|residual| > flag_threshold`
#'   are flagged: marked as `**bold**` in the kable, shown in the `flagged`
#'   column of the dataframe. Default `2`.
#' @param color_by Character. For `output = "ggplot"`, what the tile fill
#'   colour encodes. `"residual"` (default) -- diverging red-white-blue scale
#'   centred at 0, the most directly diagnostic. `"n"` -- sequential blue
#'   scale on the observed cell count, useful for spotting whether
#'   large-magnitude residuals are driven by sparse cells. Either way, the
#'   numeric residual is printed inside each cell.
#' @param color_limits Numeric length-2 vector or `NULL`. Caps the colour
#'   scale for `output = "ggplot"`. Default `c(-5, 5)` when
#'   `color_by = "residual"` (sparse-cell residuals from
#'   `(o-e)/sqrt(n*p*(1-p))` can be enormous when expected counts are tiny;
#'   capping the scale stops outliers from compressing the rest of the
#'   plot). The unclipped residual values still appear as cell labels.
#'   `NULL` for `color_by = "n"`, which uses the natural data range.
#' @param min_expected Numeric or `NULL`. If set, cells with `expected <
#'   min_expected` have their residual set to `NA` (they appear grey in the
#'   heatmap and are not flagged). Default `NULL` (no filtering). Setting
#'   `min_expected = 1` (or `5`) is the analogue of the Cochran rule for
#'   sparse-cell chi-square contributions and removes residuals whose
#'   asymptotic standard normal approximation is unreliable.
#'
#' @return
#' * `output = "kable"`: a `knitr_kable` object. For D = 2, a wide table
#'   with rows = `t1`, columns = `t2`, cells = standardised residual
#'   (`**bold**` if flagged, em-dash if NA). For D > 2, long-format with
#'   one row per cell.
#' * `output = "dataframe"`: a long-format data.frame with columns `t1`,
#'   ..., `tD`, `total`, `observed`, `expected`, `residual`, `flagged`.
#' * `output = "ggplot"`: a `geom_tile()` heatmap (D = 2 or 3 only).
#'
#' @references
#' Christensen, K. B., Bjorner, J. B., Kreiner, S., & Petersen, J. H. (2002).
#' Testing unidimensionality in polytomous Rasch models. *Psychometrika,
#' 67*(4), 563-574. \doi{10.1007/BF02295132}
#'
#' @seealso \code{\link{RMdimMartinLof}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- as.data.frame(matrix(sample(0:1, 400 * 8, replace = TRUE),
#'                             nrow = 400, ncol = 8))
#' colnames(dat) <- paste0("I", 1:8)
#'
#' # Wide kable table for D = 2
#' RMdimMartinLofResiduals(dat,
#'                      partition = list(c("I1","I2","I3","I4"),
#'                                       c("I5","I6","I7","I8")))
#'
#' # Heatmap
#' RMdimMartinLofResiduals(dat,
#'                      partition = c(1,1,1,1,2,2,2,2),
#'                      output = "ggplot")
#'
#' # Underlying data.frame for custom analysis
#' df <- RMdimMartinLofResiduals(dat,
#'                            partition = c(1,1,1,1,2,2,2,2),
#'                            output = "dataframe")
#' df[df$flagged, ]
#' }
RMdimMartinLofResiduals <- function(data,
                                 partition,
                                 output         = c("kable", "dataframe", "ggplot"),
                                 flag_threshold = 2,
                                 color_by       = c("residual", "n"),
                                 color_limits   = NULL,
                                 min_expected   = NULL) {

  output   <- match.arg(output)
  color_by <- match.arg(color_by)

  if (!is.null(min_expected) &&
      (!is.numeric(min_expected) || length(min_expected) != 1L ||
       !is.finite(min_expected) || min_expected < 0)) {
    stop("`min_expected` must be NULL or a non-negative numeric scalar.",
         call. = FALSE)
  }
  if (!is.null(color_limits) &&
      (!is.numeric(color_limits) || length(color_limits) != 2L ||
       any(!is.finite(color_limits)) || color_limits[1L] >= color_limits[2L])) {
    stop("`color_limits` must be NULL or a length-2 numeric vector with ",
         "lower < upper.", call. = FALSE)
  }

  validate_response_data(data)

  data <- stats::na.omit(as.data.frame(data))
  if (nrow(data) < 30L) {
    stop("Need at least 30 complete cases.", call. = FALSE)
  }

  partition_list <- normalize_ml_partition(partition, data)
  D <- length(partition_list)
  if (D < 2L) {
    stop("Partition must specify >= 2 subscales.", call. = FALSE)
  }
  if (any(lengths(partition_list) < 2L)) {
    stop("Each subscale must contain at least 2 items.", call. = FALSE)
  }
  if (output == "ggplot" && D > 3L) {
    stop("output = \"ggplot\" supports D = 2 or 3 subscales only. ",
         "Use output = \"dataframe\" for higher D.", call. = FALSE)
  }

  # Restrict to items used in the partition; reindex
  used_idx <- sort(unique(unlist(partition_list)))
  if (length(used_idx) < ncol(data)) {
    data <- data[, used_idx, drop = FALSE]
    remap <- stats::setNames(seq_along(used_idx), used_idx)
    partition_list <- lapply(partition_list,
                             function(idx) unname(remap[as.character(idx)]))
  }

  data_mat      <- as.matrix(data)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L
  N             <- nrow(data)

  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # CML item parameters (vector for RM, list for PCM)
  full_params <- extract_ml_sampling_params(data, is_polytomous)

  # gamma_full over all items in the partition
  gamma_full <- compute_esf_gamma(full_params)

  # gamma per subscale
  gamma_per_sub <- lapply(partition_list, function(idx) {
    sub_params <- if (is_polytomous) full_params[idx] else full_params[idx]
    compute_esf_gamma(sub_params)
  })
  M_d <- vapply(gamma_per_sub, function(g) length(g) - 1L, integer(1L))

  # Observed subscores
  subscores <- vapply(partition_list,
                      function(idx) rowSums(data[, idx, drop = FALSE]),
                      numeric(N))
  if (!is.matrix(subscores)) subscores <- matrix(subscores, nrow = N, ncol = D)

  total_scores  <- rowSums(subscores)
  n_t_table     <- table(total_scores)
  observable_t  <- as.integer(names(n_t_table))

  # All possible cells (cartesian product of 0:M_d)
  cell_indices <- expand.grid(lapply(M_d, function(M) 0:M),
                              KEEP.OUT.ATTRS = FALSE)
  colnames(cell_indices) <- paste0("t", seq_len(D))
  cell_totals <- rowSums(cell_indices)

  # Drop cells with totals not observed (n_t = 0)
  keep_cell <- cell_totals %in% observable_t
  cell_indices <- cell_indices[keep_cell, , drop = FALSE]
  cell_totals  <- cell_totals[keep_cell]

  # Observed counts via key-based lookup
  obs_keys  <- apply(subscores,    1L, function(x) paste(x, collapse = "."))
  cell_keys <- apply(cell_indices, 1L, function(x) paste(x, collapse = "."))
  obs_table <- table(obs_keys)
  observed  <- as.numeric(obs_table[cell_keys])
  observed[is.na(observed)] <- 0

  # n_t per cell
  n_t_per_cell <- as.numeric(n_t_table[as.character(cell_totals)])

  # Conditional probability p(t1,...,tD | t) under H0
  log_p <- numeric(nrow(cell_indices))
  for (d in seq_len(D)) {
    log_p <- log_p + log(gamma_per_sub[[d]][cell_indices[[d]] + 1L])
  }
  log_p <- log_p - log(gamma_full[cell_totals + 1L])
  p <- exp(log_p)
  p[!is.finite(p)] <- 0

  # Expected and standardised residuals
  expected <- n_t_per_cell * p
  se_denom <- sqrt(n_t_per_cell * p * (1 - p))
  residual <- (observed - expected) / se_denom
  residual[!is.finite(residual)] <- NA_real_

  # Optional Cochran-style filter on sparse cells (asymptotic standard
  # normal approximation breaks down when expected is tiny).
  if (!is.null(min_expected)) {
    residual[expected < min_expected] <- NA_real_
  }

  # Result data.frame
  result_df <- data.frame(
    cell_indices,
    total      = as.integer(cell_totals),
    observed   = as.integer(observed),
    expected   = round(expected, 3),
    residual   = round(residual, 3),
    flagged    = !is.na(residual) & abs(residual) > flag_threshold,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )

  # Sort by total then by t1, t2, ...
  ord_cols <- c("total", paste0("t", seq_len(D)))
  result_df <- result_df[do.call(order, result_df[ord_cols]), , drop = FALSE]
  rownames(result_df) <- NULL

  if (output == "dataframe") {
    return(result_df)
  }

  if (output == "kable") {
    if (D == 2L) {
      t1_vals <- 0:M_d[1L]
      t2_vals <- 0:M_d[2L]
      mat <- matrix("", nrow = length(t1_vals), ncol = length(t2_vals),
                    dimnames = list(as.character(t1_vals),
                                    as.character(t2_vals)))
      for (k in seq_len(nrow(result_df))) {
        r1 <- result_df$t1[k] + 1L
        r2 <- result_df$t2[k] + 1L
        v  <- result_df$residual[k]
        if (is.na(v)) {
          mat[r1, r2] <- "--"
        } else {
          formatted <- sprintf("%.2f", v)
          if (isTRUE(result_df$flagged[k])) {
            formatted <- paste0("**", formatted, "**")
          }
          mat[r1, r2] <- formatted
        }
      }
      kable_df <- as.data.frame(mat, stringsAsFactors = FALSE)
      kable_df <- cbind(`t1\\t2` = rownames(mat), kable_df)
      caption <- paste0(
        "Standardised residuals (Christensen et al. 2002, eq. 13). ",
        "Rows = subscale 1 score, columns = subscale 2 score. ",
        "**Bold** = |residual| > ", flag_threshold,
        ". -- = uncomputable. n = ", N, " complete cases."
      )
      return(knitr::kable(kable_df, format = "pipe",
                          row.names = FALSE, caption = caption))
    } else {
      caption <- paste0(
        "Standardised residuals (Christensen et al. 2002, eq. 13) for D = ",
        D, " subscales. n = ", N, " complete cases. ",
        "|residual| > ", flag_threshold, " indicates potential ",
        "dimensionality issue."
      )
      cols <- c(paste0("t", seq_len(D)), "total", "observed", "expected",
                "residual", "flagged")
      return(knitr::kable(result_df[, cols], format = "pipe",
                          caption = caption))
    }
  }

  # ggplot
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\".",
         call. = FALSE)
  }

  plot_df <- result_df
  plot_df$t1 <- as.integer(plot_df$t1)
  plot_df$t2 <- as.integer(plot_df$t2)

  # Stash the chosen fill values in a dedicated column. Avoids any
  # `.data$col` / quosure / aes-evaluation surprises that can swap branches
  # in some ggplot2 setups -- the column we map onto fill is unambiguous.
  plot_df[["fill_value"]] <- if (identical(color_by, "residual")) {
    plot_df[["residual"]]
  } else {
    as.numeric(plot_df[["observed"]])
  }

  # Resolve colour-scale limits and clip fill_value BEFORE building the
  # plot (mutating plot_df after ggplot() captures it has no effect, since
  # ggplot has its own copy of the data).
  if (identical(color_by, "residual")) {
    fill_limits <- if (is.null(color_limits)) c(-5, 5) else color_limits
    plot_df$fill_value <- pmin(pmax(plot_df$fill_value, fill_limits[1L]),
                               fill_limits[2L])
  } else {
    fill_limits <- color_limits  # NULL -> use natural range below
    if (!is.null(fill_limits)) {
      plot_df$fill_value <- pmin(pmax(plot_df$fill_value, fill_limits[1L]),
                                 fill_limits[2L])
    }
  }

  base <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_tile(
      mapping = ggplot2::aes(x = .data[["t2"]], y = .data[["t1"]],
                             fill = .data[["fill_value"]]),
      colour  = "grey90"
    )

  if (identical(color_by, "residual")) {
    base <- base +
      ggplot2::scale_fill_gradient2(
        low      = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, na.value = "grey95",
        limits   = fill_limits,
        name     = "Std.\nresidual"
      )
  } else {
    base <- base +
      ggplot2::scale_fill_gradient(
        low      = "#F7FBFF", high = "#08306B",
        na.value = "grey95",
        limits   = fill_limits,
        name     = "Observed\ncount"
      )
  }

  # Always print residual values inside cells, regardless of fill choice.
  # `inherit.aes = FALSE` so the text layer's mapping is independent of
  # the tile's fill aesthetic.
  base <- base +
    ggplot2::geom_text(
      data    = plot_df[!is.na(plot_df$residual), ],
      mapping = ggplot2::aes(x = .data[["t2"]], y = .data[["t1"]],
                             label = sprintf("%.1f", .data[["residual"]])),
      size    = 3,
      colour  = "grey20",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, max(plot_df$t2), by = 1L), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0L, max(plot_df$t1), by = 1L), expand = c(0, 0)
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(size = 9, hjust = 0)
    ) +
    er2_axis_margins()

  if (D == 2L) {
    return(
      base +
        ggplot2::labs(
          x = expression("Subscale 2 score (" * t[2] * ")"),
          y = expression("Subscale 1 score (" * t[1] * ")"),
          caption = paste0(
            "Standardised residuals (Christensen et al. 2002, eq. 13). ",
            "n = ", N, " complete cases. Cells with |residual| > ",
            flag_threshold, " indicate potential dimensionality issues."
          )
        )
    )
  }

  # D == 3: facet by t3
  base +
    ggplot2::facet_wrap(~ .data$t3, labeller = ggplot2::label_both) +
    ggplot2::labs(
      x = expression("Subscale 2 score (" * t[2] * ")"),
      y = expression("Subscale 1 score (" * t[1] * ")"),
      caption = paste0(
        "Standardised residuals (Christensen et al. 2002, eq. 13). ",
        "n = ", N, " complete cases. Faceted by subscale 3 score (",
        expression(t[3]), "). |residual| > ", flag_threshold,
        " flagged."
      )
    )
}

# ---------------------------------------------------------------------------
# Internal: gamma-function computation that handles both vector (RM) and list
# (PCM) parameter formats and returns a numeric vector gamma_0, ..., gamma_M.
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
compute_esf_gamma <- function(params) {
  if (is.list(params) && length(params) == 0L)   return(1)
  if (is.numeric(params) && length(params) == 0L) return(1)
  esf <- psychotools::elementary_symmetric_functions(params, order = 0L)
  if (is.list(esf)) esf[[1L]] else esf
}

