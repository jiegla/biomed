test_that("donutPie returns named ggplot objects", {
  d <- data.frame(response = factor(c("CR", "PR", "PR", NA)))
  plots <- donutPie(d, "response", removeNA = FALSE)
  expect_named(plots, "response")
  expect_s3_class(plots$response, "ggplot")
})

test_that("draw_stack_barplot returns statistics without writing by default", {
  d <- data.frame(
    arm = rep(c("A", "B"), each = 6),
    response = c("CR", "CR", "PR", "PR", "SD", "PD",
                 "CR", "PR", "PR", "SD", "SD", "PD")
  )
  result <- draw_stack_barplot(
    d, "arm", "response", return_plot = TRUE, save = FALSE
  )
  expect_s3_class(result$plot$response, "ggplot")
  expect_equal(result$statistics$N, 12)
  expect_true(is.finite(result$statistics$p_value))
})

test_that("plot_roc_batch returns plots and writes nothing when path is null", {
  skip_if_not_installed("pROC")
  d <- data.frame(
    response = factor(rep(c("No", "Yes"), each = 10)),
    marker = c(1:10, 11:20)
  )
  plots <- plot_roc_batch(
    d, "response", "marker", positive_class = "Yes", direction = "<"
  )
  expect_named(plots, "marker")
  expect_s3_class(plots$marker, "ggplot")
})

test_that("vennjgl returns exact intersection members without writing", {
  skip_if_not_installed("VennDiagram")
  result <- vennjgl(
    list(A = c("x", "y"), B = c("y", "z")), showFigure = FALSE,
    save = FALSE, write_xlsx = FALSE
  )
  expect_true(all(c("set", "count", "values") %in% names(result$partition)))
  expect_length(result$image_files, 0)
  expect_null(result$excel_file)
})
