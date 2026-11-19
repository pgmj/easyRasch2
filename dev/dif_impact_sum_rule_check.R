# Sum-score rule: the Study 3 simulations against RMdifImpact().
#
# Study 3 (dif_impact_ext_sim_3b.R) classified a focal person at or above a
# sum-score cut under split when the equated raw score was at or above the
# cut. RMdifImpact() uses the equated score at or above cut - 0.5 (the
# nearest whole score, and the point that defines the true status), since
# the first rule moves the focal cut by a whole raw score for any DIF,
# however small, that makes the items easier for the focal group
# (dev/dif-impact-design.md, "Build notes").
#
# For each of the simulated settings (group sizes do not matter here), this
# script gives the lowest focal raw score classified at or above the cut
# at the true parameters, under keep (the cut itself), the simulated split
# rule and the package's split rule. Deterministic, a few seconds.
#
# Run from the package root:
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript dev/dif_impact_sum_rule_check.R
# Writes dev/dif_impact_ext_sim/sum_rule_slim.rds (+ dated copy), read by
# dev/dif_impact.qmd.

suppressMessages(devtools::load_all(".", quiet = TRUE))  # .dimp_equate()
source("dev/dif_impact_ext_sim_3b.R")  # grid_scenarios() etc.; run block does not fire

cells <- grid_scenarios()
key <- unique(cells[c("type", "k", "scenario", "position", "d", "cut")])
rownames(key) <- NULL

rules <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
  cell <- as.list(key[i, ])
  sc <- make_scale(cell$type, cell$k)
  si <- scenario_items(sc$loc, cell$scenario, cell$position, cell$cut)
  thr <- with_shifts(list(ref = sc$thr, focal = sc$thr), si$j,
                     cell$d * si$sign)
  cut_raw <- cut_raw_for(cell)
  eq <- .dimp_equate(thr$ref, thr$focal)
  r <- seq_along(eq) - 1L
  lowest <- function(x) r[which(x)[1L]]
  data.frame(key[i, ], cut_raw = cut_raw,
             eq_below = eq[r == cut_raw - 1L], eq_at = eq[r == cut_raw],
             keep = cut_raw,
             split_sim = lowest(eq >= cut_raw - 1e-9),  # as in the 3b script
             split_pkg = lowest(eq >= cut_raw - 0.5))
}))

out <- list(
  rules = rules,
  note = paste("Lowest focal raw score at or above the sum-score cut at the",
               "true parameters: keep, the Study 3 split rule (eq >= cut)",
               "and the RMdifImpact() split rule (eq >= cut - 0.5)."),
  created = Sys.Date()
)
dir <- file.path("dev", "dif_impact_ext_sim")
saveRDS(out, file.path(dir, "sum_rule_slim.rds"))
saveRDS(out, file.path(dir, paste0("sum_rule_slim_", Sys.Date(), ".rds")))

print(rules, row.names = FALSE, digits = 3)
cat("settings:", nrow(rules),
    "| split rule differs from keep, Study 3 rule:",
    sum(rules$split_sim != rules$keep),
    "| package rule:", sum(rules$split_pkg != rules$keep), "\n")
