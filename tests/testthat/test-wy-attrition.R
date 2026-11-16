# .wy_attrition_message() (review finding DA6)

test_that("the attrition message never reports a floor above 1", {
  expect_silent(easyRasch2:::.wy_attrition_message(100, 100))
  expect_silent(easyRasch2:::.wy_attrition_message(NULL, 100))
  expect_message(easyRasch2:::.wy_attrition_message(50, 100), "1/\\(50\\+1\\)")
  m <- tryCatch(easyRasch2:::.wy_attrition_message(1, 100, floor_mult = 2L),
                message = conditionMessage)
  expect_match(m, "= 1\\.")
  expect_message(easyRasch2:::.wy_attrition_message(0, 100, floor_mult = 2L),
                 "corrected p-values are missing")
})
