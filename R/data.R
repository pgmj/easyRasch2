#' PHQ-9 Depression Screener (NHANES Subsample)
#'
#' A processed subsample of the Patient Health Questionnaire 9-item (PHQ-9)
#' depression screener from the U.S. National Health and Nutrition
#' Examination Survey (NHANES), September 2024 release. Six hundred
#' respondents were drawn at random from the cycle's PHQ-9 module subject
#' to having complete responses on all nine items, while retaining a
#' realistic share of respondents with a sum-score of zero (n = 8) so
#' that floor behaviour can be illustrated in a Rasch analysis.
#'
#' @format A data frame with 600 rows and 12 variables:
#' \describe{
#'   \item{q1}{Little interest or pleasure in doing things. Integer 0--3.}
#'   \item{q2}{Feeling down, depressed, or hopeless. Integer 0--3.}
#'   \item{q3}{Trouble falling/staying asleep, or sleeping too much.
#'     Integer 0--3.}
#'   \item{q4}{Feeling tired or having little energy. Integer 0--3.}
#'   \item{q5}{Poor appetite or overeating. Integer 0--3.}
#'   \item{q6}{Feeling bad about yourself --- or that you are a failure
#'     or have let yourself or your family down. Integer 0--3.}
#'   \item{q7}{Trouble concentrating on things, such as reading the
#'     newspaper or watching television. Integer 0--3.}
#'   \item{q8}{Moving or speaking so slowly that other people could have
#'     noticed --- or the opposite, being so fidgety or restless that
#'     you have been moving around a lot more than usual. Integer 0--3.}
#'   \item{q9}{Thoughts that you would be better off dead, or of hurting
#'     yourself in some way. Integer 0--3.}
#'   \item{gender}{Self-reported gender, factor with levels
#'     `"Female"` and `"Male"` (31 respondents with missing values).}
#'   \item{age}{Age in years (integer, range 15--85).}
#'   \item{edu}{Highest educational attainment, factor with levels
#'     `"Elementary School"`, `"High school"`, `"University"`.}
#' }
#'
#' Each PHQ-9 item uses a four-point ordinal response scale, scored
#' 0 (`"Not at all"`), 1 (`"Several days"`), 2 (`"More than half the
#' days"`) and 3 (`"Nearly every day"`).
#'
#' @details
#' The dataset is a *processed subsample* intended for teaching and for
#' the package's worked example; it should not be treated as a canonical
#' NHANES microdata file. Users wishing to validate against NCHS-published
#' figures should download the original public-use microdata directly
#' from the NHANES website (see *Source*).
#'
#' @source
#' U.S. Centers for Disease Control and Prevention, National Center for
#' Health Statistics. *National Health and Nutrition Examination Survey*,
#' September 2024 release.
#' \url{https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Questionnaire&CycleBeginYear=2024}.
#' NHANES data are released to the public domain by the U.S. federal
#' government (\url{https://www.cdc.gov/nchs/policy/data-release-policy.html}).
#'
#' @references
#' Kroenke, K., Spitzer, R. L., & Williams, J. B. W. (2001). The PHQ-9:
#'   Validity of a brief depression severity measure. *Journal of General
#'   Internal Medicine*, 16(9), 606--613.
#'   \doi{10.1046/j.1525-1497.2001.016009606.x}
#'
#' @examples
#' data(phq9)
#' str(phq9)
#' summary(rowSums(phq9[, 1:9]))
"phq9"
