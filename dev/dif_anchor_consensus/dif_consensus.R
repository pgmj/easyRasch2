# Prototype: consensus forward anchor over random tie-breaks (RMdifAnchor
# development). Development only, not package code.
#
# psychotools::anchor() breaks ties in the item ranking at random, so the
# forward (MTT) anchor can change with the seed on the same data. Here the
# ranking is redone under S seeds (auxiliary tests computed once, see
# dif_thr_anchor_ties() in ../dif_anchor_threshold/dif_thresh.R), and the
# anchor is built from the inclusion frequencies:
#
#   freq_i   share of seeds whose anchor contains item i
#   first_i  share of seeds where item i is the first (untested) anchor item
#   anchor   items with freq_i >= cut (default .5). If none reach the cut,
#            the items with the highest frequency.
#   order    the reference (first, untested) item is the anchor item with the
#            highest first_i, ties broken by freq_i, then by item number.
#            The other anchor items follow by decreasing freq_i.
#
# The result is deterministic in the limit S -> Inf. For finite S it can
# still change with the seed set when some freq_i is close to the cut.

# Consensus from a list of per-seed anchors (each a vector of item indices,
# first element = reference item). k = number of items.
dif_consensus_from <- function(anchors, k, cut = 0.5) {
  S <- length(anchors)
  freq <- tabulate(unlist(anchors), k) / S
  first <- tabulate(vapply(anchors, function(a) as.integer(a[1L]), integer(1)), k) / S
  A <- which(freq >= cut)
  if (!length(A)) A <- which(freq == max(freq))
  ref <- A[order(-first[A], -freq[A], A)][1L]
  rest <- setdiff(A, ref)
  rest <- rest[order(-freq[rest], rest)]
  list(anchor = c(ref, rest), freq = freq, first = first)
}

dif_thr_anchor_consensus <- function(fit, seeds, cut = 0.5) {
  As <- dif_thr_anchor_ties(fit, "forward", seeds)
  dif_consensus_from(As, length(fit$items), cut)
}
