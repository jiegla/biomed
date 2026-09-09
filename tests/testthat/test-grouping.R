test_that("makegroup creates binary and tertile groups", {
  d <- data.frame(marker = 1:9)
  binary <- makegroup(d, "marker")
  expect_equal(levels(binary$marker_binary), c("Low", "High"))
  expect_equal(as.character(binary$marker_binary[5]), "Low")
  expect_equal(attr(binary, "cutoff_table")$cutoff, 5)

  ternary <- makegroup(d, "marker", num_group = 3)
  expect_equal(levels(ternary$marker_ternary), c("Low", "Middle", "High"))
  expect_equal(sum(table(ternary$marker_ternary)), 9)
})

test_that("makegroup does not silently encode non-numeric factors", {
  d <- data.frame(marker = factor(c("low", "high")))
  expect_error(makegroup(d, "marker"), "must be numeric")
})

test_that("makegroup records ROC direction and positive class safely", {
  skip_if_not_installed("pROC")
  d <- data.frame(
    response = factor(rep(c("No", "Yes"), each = 10)),
    marker = c(1:10, 11:20)
  )
  out <- makegroup(
    d, "marker", method = "roc_youden", response = "response",
    positive_class = "Yes", direction = "<"
  )
  expect_s3_class(out$marker_binary, "factor")
  expect_gt(attr(out, "cutoff_table")$auc, 0.9)
})
