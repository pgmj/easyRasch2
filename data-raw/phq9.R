# Prepare the bundled phq9 dataset.
#
# Source -----------------------------------------------------------------
# Data files were manually downloaded from the September 2024 release of the
# U.S. National Health and Nutrition Examination Survey (NHANES) data
# website:
#   https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Questionnaire&CycleBeginYear=2024
# NHANES microdata are produced by the U.S. CDC/NCHS and released as
# public-use data files in the public domain (see
# https://www.cdc.gov/nchs/policy/data-release-policy.html). No re-
# distribution licence applies.
#
# The full NHANES sample is large; we drew a random subsample of 600
# respondents with complete responses to all nine PHQ-9 items, deliberately
# retaining enough respondents with a sum-score of zero so that the dataset
# can illustrate floor behaviour in a Rasch analysis.
#
# To regenerate the bundled .rda from data-raw/n600.csv, run:
#   Rscript data-raw/phq9.R
# from the package root.

raw <- read.csv("data-raw/n600.csv", stringsAsFactors = FALSE)

# Sanity checks
stopifnot(
  nrow(raw) == 600L,
  all(c(paste0("q", 1:9), "gender", "age", "edu") %in% names(raw)),
  all(vapply(raw[paste0("q", 1:9)], function(x) all(x %in% 0:3), logical(1)))
)

# Item columns kept as integer (PHQ-9 0..3 ordinal)
for (i in 1:9) raw[[paste0("q", i)]] <- as.integer(raw[[paste0("q", i)]])

# Demographics
raw$gender <- factor(raw$gender, levels = c("Female", "Male"))
raw$edu    <- factor(raw$edu,
                     levels = c("Elementary School",
                                "High school",
                                "University"))
raw$age    <- as.integer(raw$age)

phq9 <- raw[, c(paste0("q", 1:9), "gender", "age", "edu")]
rownames(phq9) <- NULL

# Save with maximum compression. usethis::use_data() does this for us,
# but we avoid the build-time dependency by writing directly.
save(phq9, file = "data/phq9.rda", compress = "bzip2", version = 3)

# Optional: confirm round-trip
loaded <- new.env()
load("data/phq9.rda", envir = loaded)
stopifnot(identical(loaded$phq9, phq9))
message("data/phq9.rda written (", nrow(phq9), " rows, ",
        ncol(phq9), " columns).")
