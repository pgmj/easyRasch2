# Internal helpers shared by RMdifLR() and RMdifAnchor().

#' Check the DIF variable and drop rows that cannot be used
#'
#' Checks that `dif_var` has one value per row of `data`, then drops rows
#' with a missing DIF group and respondents with no responses, each with a
#' message, keeping `data` and the group factor aligned.
#'
#' @param data Response data.frame (already through
#'   `validate_response_data()`).
#' @param dif_var Vector of group memberships, one per row of `data`.
#' @return A list with `data`, `dif_factor` (factor with unused levels
#'   dropped), `n_total` (rows before dropping) and `has_na` (whether the
#'   kept responses have missing values).
#' @keywords internal
#' @noRd
.prepare_dif_data <- function(data, dif_var) {
  data <- as.data.frame(data)
  n_total <- nrow(data)
  if (length(dif_var) != nrow(data)) {
    stop(
      "`dif_var` must have the same length as nrow(data) (", nrow(data),
      "). Got ", length(dif_var), ".",
      call. = FALSE
    )
  }
  dif_factor <- droplevels(as.factor(dif_var))
  na_mask <- is.na(dif_factor)
  if (any(na_mask)) {
    message(sum(na_mask), " row(s) with NA in `dif_var` dropped.")
    data <- data[!na_mask, , drop = FALSE]
    dif_factor <- droplevels(dif_factor[!na_mask])
  }
  # Respondents with no responses contribute nothing to CML. Drop them
  # jointly with `dif_factor`.
  answered <- rowSums(!is.na(as.matrix(data))) > 0L
  data <- .drop_empty_respondents(data)
  dif_factor <- droplevels(dif_factor[answered])
  list(data = data, dif_factor = dif_factor, n_total = n_total,
       has_na = anyNA(data))
}

#' Stop when categories cannot be compared across DIF groups
#'
#' Runs the full-data category checks, then `.check_group_categories()`.
#'
#' @param data Response data.frame.
#' @param dif_factor Factor of DIF groups, without `NA`.
#' @param fn Name of the calling function, for the message.
#' @return Invisibly `TRUE`; stops with an `easyRasch2_category_error`.
#' @keywords internal
#' @noRd
.check_dif_categories <- function(data, dif_factor, fn) {
  .check_response_categories(data)
  ec_full <- .extreme_only_categories(data)
  if (nrow(ec_full) > 0L) {
    stop(.category_error(
      paste0(
        .extreme_only_message(ec_full),
        " ", fn, "() stops here, since these thresholds cannot be ",
        "compared across groups."
      ),
      unique(ec_full$Item)
    ))
  }
  .check_group_categories(data, dif_factor)
}

#' CML fit that names the sample when it fails
#'
#' @param d Response data.frame, already checked.
#' @param where Where the fit is, for the message ("group 'A'").
#' @return A psychotools fit from `.cml_fit()`.
#' @keywords internal
#' @noRd
.cml_fit_or_stop <- function(d, where) {
  tryCatch(
    # pcmodel() for every fit (a dichotomous item is a 2-category PCM; the
    # raschmodel() Hessian fails on some incomplete data). check = FALSE:
    # the data were checked before.
    .cml_fit(d, check = FALSE),
    error = function(e) {
      stop(
        "CML estimation failed in ", where, " (", conditionMessage(e),
        "). Small or sparse groups are the usual cause; inspect the ",
        "counts with RMplotTile(data, group = dif_var).",
        call. = FALSE
      )
    }
  )
}

#' Stop unless every category an item uses occurs in every DIF group
#'
#' psychotools fits an item with an unused category in one group with fewer
#' thresholds there (or, for an empty lowest or interior category, shifts or
#' drops them). Its parameters would then differ in number or meaning across
#' groups, so neither the LR test nor the anchored comparisons are defined. The
#' full data have already passed `.check_response_categories()`, so an item's
#' categories there run 0 to its highest observed one. A category used in a
#' group only by respondents with a zero or perfect score there also stops:
#' CML leaves those respondents out, so the category is unused in that
#' group's fit and its threshold diverges.
#'
#' @param data Response data.frame (items from 0).
#' @param dif_factor Factor of the same length, without `NA`.
#' @return Invisibly `TRUE`; stops with an `easyRasch2_category_error`
#'   naming the items, the categories and the groups.
#' @keywords internal
#' @noRd
.check_group_categories <- function(data, dif_factor) {
  data <- as.matrix(data)
  items <- colnames(data)
  problems <- character(0)
  bad <- character(0)
  unused <- matrix(FALSE, ncol(data), nlevels(dif_factor))
  full_max <- apply(data, 2L, max, na.rm = TRUE)
  for (j in seq_len(ncol(data))) {
    x <- data[, j]
    cats <- sort(unique(x[!is.na(x)]))
    for (gi in seq_along(levels(dif_factor))) {
      g <- levels(dif_factor)[gi]
      xg <- x[dif_factor == g]
      miss <- setdiff(cats, xg[!is.na(xg)])
      if (length(miss) > 0L) {
        unused[j, gi] <- TRUE
        bad <- union(bad, items[j])
        problems <- c(problems, paste0(
          items[j], " (",
          if (length(miss) == length(cats)) {
            "no responses"
          } else {
            paste0(
              if (length(miss) == 1L) "category " else "categories ",
              paste(miss, collapse = ", "), " unused"
            )
          },
          " in group '", g, "')"
        ))
      }
    }
  }
  # Extreme-score-only categories per group, for items whose categories are
  # all used there (else the group's item maxima would differ).
  for (gi in seq_along(levels(dif_factor))) {
    g <- levels(dif_factor)[gi]
    # Perfect is judged against the full-data maxima, so an item whose top
    # category is unused in the group does not make others look extreme-only.
    ec <- .extreme_only_categories(
      data[dif_factor == g, , drop = FALSE],
      item_max = full_max
    )
    ec <- ec[!unused[match(ec$Item, items), gi], , drop = FALSE]
    for (it in unique(ec$Item)) {
      k <- ec$Category[ec$Item == it]
      bad <- union(bad, it)
      problems <- c(problems, paste0(
        it, " (", if (length(k) == 1L) "category " else "categories ",
        paste(k, collapse = ", "), " used only by respondents with a zero ",
        "or perfect score in group '", g, "')"
      ))
    }
  }
  if (length(problems) > 0L) {
    msg <- paste0(
      "Item parameters cannot be compared across DIF groups for: ",
      paste(problems, collapse = ", "),
      ". Every response category of an item must be used in every group, ",
      "by respondents without a zero or perfect score (conditional maximum ",
      "likelihood leaves those out). ",
      "Merge a sparse category with an adjacent one, combine groups, or ",
      "leave the item out. RMplotTile(data, group = dif_var) shows the ",
      "counts."
    )
    stop(.category_error(msg, bad))
  }
  invisible(TRUE)
}


#' Evaluate an expression with a fixed seed, restoring the caller's RNG
#'
#' RMdifAnchor() draws random tie-break orderings for the forward anchor and
#' Monte-Carlo draws for the single-step adjustment. A fixed seed makes the
#' result reproducible without changing the caller's random stream. The
#' generator is pinned to R's defaults, as in the other seeded functions (see
#' `?"easyRasch2-reproducibility"`), so the same seed gives the same result
#' under any `RNGkind()`, and the caller's generator is restored on exit.
#'
#' @keywords internal
#' @noRd
.with_local_seed <- function(seed, expr) {
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had_seed) get(".Random.seed", envir = genv) else NULL
  old_kind <- RNGkind()
  on.exit({
    # RNGkind() first: it reseeds, which the saved state then overwrites.
    suppressWarnings(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]))
    if (had_seed) {
      assign(".Random.seed", old, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  expr
}
