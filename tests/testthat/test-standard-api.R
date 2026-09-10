test_that("standard analysis interfaces preserve the legacy calculations", {
  d <- data.frame(group = rep(c("A", "B"), each = 6),
                  marker = c(1:6, 5:10), status = rep(0:1, 6), time = c(5, 2, 8, 4, 1, 7, 3, 9, 6, 12, 10, 11))
  expect_equal(batch_anova(d, "marker")$p_value, batch_ANOVA(d, feature = "marker")$p.value)
  expect_equal(batch_cox(d, "marker", "time", "status")$hr,
               cophx_batch("marker", d, "time", "status")$HR)
  expect_equal(batch_chi_square(d, "group")$p_value,
               chi_square_batch_df(d, "group")$P_Value)
  expect_equal(make_group(d, "marker"), makegroup(d, "marker"))
})

test_that("ANOVA preprocessing really removes invalid features", {
  d <- data.frame(group = rep(c("A", "B"), each = 4), good = c(1:4, 3:6),
                  constant = 1, missing = c(NA, 2:8), infinite = c(Inf, 2:8))
  expect_warning(out <- batch_anova(d, c("good", "constant", "missing", "infinite"),
                                     preprocess = TRUE), "skipped")
  expect_equal(out$variable, "good")
})

test_that("sparse larger tables use Fisher and match base R", {
  d <- data.frame(x = rep(letters[1:3], each = 3), y = rep(c("a", "a", "b"), 3))
  out <- batch_chi_square(d, "x", "y")
  expect_equal(out$method, "Fisher's Exact Test")
  expect_equal(out$p_value, stats::fisher.test(table(d$x, d$y))$p.value)
})

test_that("tied grouping values and invalid survival input are handled", {
  d <- data.frame(x = c(1, 1, 1, 1, NA), time = 1:5, status = c(0, 1, 2, 0, 1))
  out <- make_group(d, "x", n_groups = 3)
  expect_equal(as.character(out$x_ternary), c(rep("Low", 4), NA))
  expect_error(batch_cox(d, "x", "time", "status"), "coded 0")
  expect_error(batch_cox(d, "x", "time", "status", conf_level = NA))
  expect_error(make_group(d, "x", n_groups = NA))
})

test_that("numeric factors keep their numeric values in ROC", {
  skip_if_not_installed("pROC")
  d <- data.frame(y = rep(c("No", "Yes"), each = 5),
                  x = factor(c(1, 2, 10, 3, 5, 4, 8, 11, 12, 7)))
  out <- batch_roc(d, "x", "y", positive_class = "Yes", direction = "<")
  reference <- pROC::roc(d$y, as.numeric(as.character(d$x)),
                         levels = c("No", "Yes"), direction = "<", quiet = TRUE)
  expect_equal(out$auc, as.numeric(pROC::auc(reference)))
  expect_equal(out$positive_class, "Yes")
  expect_true(is.na(out$error))
  curve <- plot_roc(d, "x", "y", positive_class = "Yes", direction = "<")$x
  expect_equal(curve$data$false_positive_rate, 1 - reference$specificities)
  expect_equal(curve$data$true_positive_rate, reference$sensitivities)
  for (method in c("fixed_sensitivity", "fixed_specificity")) {
    grouped <- make_group(d, "x", method = method, response = "y",
                           positive_class = "Yes", direction = "<", fixed_value = 0.7)
    cut <- attr(grouped, "cutoff_table")
    expect_true(is.finite(cut$cutoff))
    metric <- if (method == "fixed_sensitivity") cut$sensitivity else cut$specificity
    expect_gte(metric, 0.7)
  }
})

test_that("plots build and outputs are readable", {
  skip_if_not_installed("openxlsx")
  d <- data.frame(group = rep(c("A", "B"), each = 4),
                  response = rep(c("Yes", "No"), 4), marker = c(1:4, 3:6))
  path <- tempfile()
  on.exit(unlink(path, recursive = TRUE))
  stats <- batch_anova(d, "marker", output_dir = path)
  expect_equal(openxlsx::read.xlsx(file.path(path, "anova.xlsx"))$p_value, stats$p_value)
  plots <- plot_donut(d, "response", output_dir = path)
  expect_silent(ggplot2::ggplot_build(plots$response))
  expect_true(all(file.exists(file.path(path, paste0("response_donut.", c("pdf", "png"))))))
  stacked <- plot_stacked_bar(d, "response", "group", output_dir = path, dpi = 72)
  expect_silent(ggplot2::ggplot_build(stacked$plots$response))
  expect_true(file.exists(file.path(path, "group_statistics.xlsx")))
})

test_that("Venn works for two through six sets and retains exact members", {
  skip_if_not_installed("VennDiagram")
  for (n in 2:6) {
    sets <- stats::setNames(lapply(seq_len(n), function(i) c("shared", paste0("unique", i))),
                            LETTERS[seq_len(n)])
    out <- plot_venn(sets, show_plot = FALSE)
    expect_equal(sum(out$partition$count), length(unique(unlist(sets))))
    expect_true(any(out$partition$values == "shared"))
    expect_length(out$image_files, 0)
  }
})

test_that("unequal Venn partitions preserve all members and export correctly", {
  skip_if_not_installed("VennDiagram")
  skip_if_not_installed("openxlsx")
  path <- tempfile()
  on.exit(unlink(path, recursive = TRUE))
  out <- plot_venn(list(A = letters[1:4], B = letters[3:6]), show_plot = FALSE,
                   output_dir = path, file_type = "pdf")
  expect_equal(sum(out$partition$count), 6L)
  expect_equal(out$partition$values[out$partition$A & out$partition$B], "c, d")
  expect_true(all(file.exists(c(out$image_files, out$excel_file))))
  expect_equal(sum(openxlsx::read.xlsx(out$excel_file)$count), 6L)
})
