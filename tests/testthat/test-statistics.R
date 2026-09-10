test_that("batch_ANOVA returns group means and adjusted p-values", {
  d <- data.frame(
    arm = factor(rep(c("A", "B"), each = 6)),
    signal = c(1:6, 10:15),
    noise = rep(c(1, 2), 6)
  )
  out <- batch_ANOVA(d, target = "arm", feature = c("signal", "noise"))
  expect_equal(sort(out$sig_names), c("noise", "signal"))
  expect_true(all(c("p.adj", "mean_A", "mean_B") %in% names(out)))
  expect_lt(out$p.value[out$sig_names == "signal"], 0.001)
})

test_that("categorical tests switch to Fisher for sparse 2 by 2 tables", {
  d <- data.frame(
    feature = factor(c(rep("A", 9), "B")),
    response = factor(c(rep("No", 8), "Yes", "Yes"))
  )
  out <- chi_square_batch_df(d, "feature", response = "response")
  expect_equal(out$Test_Type, "Fisher's Exact Test")
  expect_true(is.finite(out$P_Value))
  expect_true("P_Adjust" %in% names(out))
})

test_that("cophx_batch retains every contrast for a multilevel factor", {
  d <- data.frame(
    time = 1:18,
    status = rep(c(1, 0, 1), 6),
    subtype = factor(rep(c("A", "B", "C"), 6)),
    marker = seq(0.1, 1.8, by = 0.1)
  )
  out <- suppressWarnings(cophx_batch(
    c("subtype", "marker"), d, time = "time", status = "status"
  ))
  expect_equal(sum(out$Marker == "subtype"), 2)
  expect_equal(sum(out$Marker == "marker"), 1)
  expect_true(all(c("HR", "CI_lower", "CI_upper", "p_adjust") %in% names(out)))
})

test_that("roc_batch uses the requested positive class", {
  skip_if_not_installed("pROC")
  d <- data.frame(
    response = factor(rep(c("No", "Yes"), each = 12)),
    marker = c(1:12, 20:31)
  )
  expect_warning(out <- roc_batch(
    "marker", d, "response", positive_class = "Yes", direction = "<"
  ), "AUC == 1")
  expect_equal(out$Positive, 12)
  expect_gt(out$AUC, 0.99)
  expect_true(is.na(out$Error))
})
