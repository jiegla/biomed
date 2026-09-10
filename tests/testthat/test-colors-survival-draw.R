test_that("the built-in palette preserves the exact supplied RGBA colors", {
  expected <- c("#E64B35CC", "#4DBBD5CC", "#00A087CC", "#3C5488CC", "#F39B7FCC",
                "#8491B4CC", "#91D1C2CC", "#DC0000CC", "#7E6148CC", "#5050FFCC",
                "#CE3D32CC", "#749B58CC", "#F0E685CC", "#466983CC")
  expect_identical(biomed_colors(), expected)
  expect_identical(biomed_colors(3), expected[1:3])
  expect_identical(biomed:::.biomed_palette(14), expected)
  expect_length(biomed_colors(0), 0)
  expect_length(biomed_colors(20), 20)
  expect_true(all(grDevices::col2rgb(biomed_colors(20), alpha = TRUE)[4, ] == 204))
  expect_error(biomed_colors(-1))
  expect_error(biomed_colors(1.5))
  d <- data.frame(x = 1:8, y = c(2, 1, 4, 3, 6, 5, 8, 7), g = rep(c("A", "B"), 4))
  p <- plot_correlation(d, "x", "y", group = "g", add.regress = FALSE)
  expect_equal(unname(p$scales$get_scales("colour")$palette(2)), expected[1:2])
})

test_that("registered grid drawing lets ggsave export whole survival objects", {
  skip_if_not_installed("survminer")
  set.seed(21)
  d <- data.frame(time = stats::rexp(60), status = rep(c(1, 1, 0), 20),
                   group = rep(c("A", "B"), 30))
  p <- surv_fig_hr("group", d, time = "time", status = "status", output_dir = NULL)
  expect_s3_class(p$table, "ggplot")
  expect_identical(utils::getS3method("grid.draw", "ggsurvplot"),
                   get("grid.draw.ggsurvplot", envir = asNamespace("biomed")))
  expect_equal(unname(p$plot$scales$get_scales("colour")$palette(2)), biomed_colors(2))
  path <- tempfile()
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  for (ext in c("pdf", "png")) {
    f <- file.path(path, paste0("survival.", ext))
    expect_no_error(ggplot2::ggsave(f, plot = p, width = 8, height = 7, dpi = 72))
    expect_gt(file.info(f)$size, 1000)
  }
})
