# Prototype: anchored Wald tests for item-level DIF between two groups,
# with multiplicity correction. Development only, not package code.
#
# Each item's DIF estimate is its location difference (group 2 minus group 1,
# in logits) after the two groups' CML scales are linked through an anchor.
# The Wald z is the estimate over its SE, and the p-values are adjusted across
# the tested items.
#
# Anchors:
#   "centred"   All items. Locations centred within each group, as in
#               RMdifLR() today, so every item is tested against the mean
#               item. All k items get a test. Implemented here with
#               psychotools::itempar() and multcomp::glht().
#   "all-other" Each item anchored on all other items (psychotools).
#   "constant"  Constant anchor class, Gini selection (psychotools default).
#   "forward"   Iterative forward anchor class, MTT selection (psychotools
#               default).
# The constant and forward classes leave the first anchor item untested
# (k - 1 results, see ?psychotools::anchortest). Its row gets NA.
#
# Adjustments: "none", "holm" (stats::p.adjust) and "single-step" (multcomp,
# the max-|z| over the joint normal of the tested estimates, which uses their
# correlation).
#
# Requires psychotools and multcomp. Partial credit models (pcmodel()) give
# item locations as mean thresholds, dichotomous items work the same way.

dif_wald_fit <- function(data, group) {
  group <- factor(group)
  if (nlevels(group) != 2L) stop("Prototype handles two groups only.")
  lv <- levels(group)
  data <- as.matrix(data)
  list(
    m1 = psychotools::pcmodel(data[group == lv[1L], , drop = FALSE]),
    m2 = psychotools::pcmodel(data[group == lv[2L], , drop = FALSE]),
    levels = lv,
    items = colnames(data)
  )
}

# Returns one row per item and anchor: Item, Anchor, Diff (group 2 minus
# group 1), SE, z, p, padj_holm, padj_ss, is_anchor.
dif_wald <- function(fit, anchors = c("centred", "all-other", "constant", "forward")) {
  items <- fit$items
  k <- length(items)
  out <- lapply(anchors, function(an) {
    if (an == "centred") {
      i1 <- psychotools::itempar(fit$m1, vcov = TRUE)
      i2 <- psychotools::itempar(fit$m2, vcov = TRUE)
      est <- stats::setNames(as.numeric(stats::coef(i2) - stats::coef(i1)), items)
      V <- stats::vcov(i1) + stats::vcov(i2)
      gl <- multcomp::glht(multcomp::parm(est, V), linfct = diag(k))
      anchor_items <- items
    } else {
      args <- list(fit$m1, fit$m2, adjust = "none")
      args$class <- an
      if (an == "constant") args$select <- "Gini"
      if (an == "forward") args$select <- "MTT"
      a <- do.call(psychotools::anchortest, args)
      gl <- a$final_tests
      anchor_items <- if (an == "all-other") character() else items[a$anchor_items]
    }
    raw <- summary(gl, test = multcomp::adjusted("none"))$test
    ss <- summary(gl, test = multcomp::adjusted("single-step"))$test
    tested <- if (an == "centred") items else names(raw$coefficients)
    # anchortest() reports object minus object2 (group 1 minus group 2).
    sgn <- if (an == "centred") 1 else -1
    res <- data.frame(
      Item = items, Anchor = an,
      Diff = NA_real_, SE = NA_real_, z = NA_real_, p = NA_real_,
      padj_holm = NA_real_, padj_ss = NA_real_,
      is_anchor = items %in% anchor_items,
      stringsAsFactors = FALSE
    )
    idx <- match(tested, items)
    res$Diff[idx] <- sgn * as.numeric(raw$coefficients)
    res$SE[idx] <- as.numeric(raw$sigma)
    res$z[idx] <- sgn * as.numeric(raw$tstat)
    res$p[idx] <- as.numeric(raw$pvalues)
    res$padj_holm[idx] <- stats::p.adjust(res$p[idx], "holm")
    res$padj_ss[idx] <- as.numeric(ss$pvalues)
    res
  })
  do.call(rbind, out)
}
