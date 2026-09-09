test_that("required columns are accepted", {
  x <- data.frame(id = 1, response = "PR")
  expect_invisible(biomed_check_columns(x, c("id", "response")))
})

test_that("missing columns are reported", {
  x <- data.frame(id = 1)
  expect_error(
    biomed_check_columns(x, c("response", "pfs_months")),
    "Required columns are missing"
  )
})

test_that("invalid inputs fail early", {
  expect_error(biomed_check_columns(matrix(1), "id"), "data frame")
  expect_error(biomed_check_columns(data.frame(id = 1), NA_character_))
})
