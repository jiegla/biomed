survival_fixture <- function() {
  set.seed(412)
  n <- 120
  d <- data.frame(marker = stats::rnorm(n), group = factor(rep(c("A", "B", "C"), each = 40)))
  d$time <- stats::rexp(n, rate = exp(.4 * d$marker + .5 * as.numeric(d$group)))
  d$status <- rep(c(1, 1, 0), 40)
  d
}

test_that("file names are portable and missing values stay missing", {
  expect_equal(sanitize_filename(c("a/b\\c", "a b\tc", "a:*?\"<>|", "..", "CON", NA)),
               c("a_b_c", "a_b_c", "a", "output", "_CON", NA))
  expect_equal(nchar(sanitize_filename(strrep("x", 250))), 200L)
  expect_error(sanitize_filename("x", max_length = 0))
})

test_that("status conversion is explicit and never silently loses unknown labels", {
  expect_equal(biomed:::.biomed_status01(c("alive", "dead", NA), "labels"), c(0, 1, NA))
  expect_equal(biomed:::.biomed_status01(c(1, 1), "12"), c(0, 0))
  expect_equal(biomed:::.biomed_status01(c(1, 1), "01"), c(1, 1))
  expect_equal(biomed:::.biomed_status01(factor(c("1", "2")), "auto"), c(0, 1))
  expect_error(biomed:::.biomed_status01(c("alive", "unknown")), "Unrecognized")
  expect_error(biomed:::.biomed_status01(c(0, 2), "01"), "encoding")
  expect_equal(biomed:::.biomed_status01(c(TRUE, FALSE, NA), "01"), c(1, 0, NA))
})

test_that("batch survival retains all contrasts and reports failures", {
  d <- survival_fixture()
  out <- batch_survival(d, c("marker", "group", "missing"))
  expect_equal(nrow(out$results), 3L)
  expect_equal(out$failed$variable, "missing")
  fit <- survival::coxph(survival::Surv(time, status) ~ marker, data = d)
  expect_equal(out$results$hr[out$results$variable == "marker"], unname(exp(stats::coef(fit))))
  expect_equal(out$results$p_adjusted, stats::p.adjust(out$results$p_value, "BH"))
  legacy <- batch_surv(d, c("marker", "group"), verbose = FALSE)
  expect_equal(legacy$HR, out$results$hr)
  expect_s3_class(attr(legacy, "failed"), "data.frame")
  d$status <- d$status + 1
  expect_equal(batch_survival(d, "marker", status_encoding = "12")$results$hr,
               out$results$hr[out$results$variable == "marker"])
  expect_error(batch_survival(d, "marker"), "encoding")
  expect_error(batch_survival(d, "marker", min_prop = .9), "minprop")
  d$time[1] <- -1
  expect_error(batch_survival(d, "marker", status_encoding = "12"), "nonnegative")
})

test_that("cutoffs respect group sizes and preserve existing binary columns", {
  skip_if_not_installed("survminer")
  d <- survival_fixture()
  d$marker_binary <- stats::rnorm(nrow(d))
  cut <- find_survival_cutoff(d, "marker", min_prop = .2, fallback = "median")
  expect_true(is.finite(cut$cutoff))
  expect_equal(cut$data$marker_binary, d$marker_binary)
  expect_false(identical(cut$group_column, "marker_binary"))
  expect_true(all(table(cut$data[[cut$group_column]]) >= ceiling(.2 * nrow(d))))
  expect_equal(as.character(cut$data[[cut$group_column]]), ifelse(d$marker <= cut$cutoff, "Low", "High"))
  out <- batch_survival(d, c("marker", "marker_binary"), best_cutoff = TRUE,
                         fallback = "median", min_event_group = 0)
  expect_setequal(out$results$variable, c("marker", "marker_binary"))
  one <- batch_survival(d, "marker_binary", best_cutoff = TRUE,
                         fallback = "median", min_event_group = 0)
  expect_equal(out$results$hr[out$results$variable == "marker_binary"], one$results$hr)
  d$marker <- 1
  fail <- find_survival_cutoff(d, "marker")
  expect_true(is.na(fail$cutoff))
  expect_match(fail$reason, "unique")
})

test_that("correlation filters original zeros and missing groups before scaling", {
  d <- data.frame(x = c(0, 1, 2, 3, 8, 10), y = c(7, 3, 9, 2, 8, 5),
                   g = c("a", "a", "b", NA, "a", "b"))
  p <- plot_correlation(d, "x", "y", group = "g", remove_x_zero = TRUE,
                         remove_group_na = TRUE, colors = c(a = "red", b = "blue"))
  used <- d[d$x != 0 & !is.na(d$g), ]
  expect_equal(attr(p, "correlation")$n, nrow(used))
  expect_equal(attr(p, "correlation")$estimate,
               unname(stats::cor.test(used$x, used$y, method = "spearman", exact = FALSE)$estimate))
  expect_equal(p$data$x, as.numeric(scale(used$x)))
  expect_equal(p$data$y, as.numeric(scale(used$y)))
  collision <- data.frame(categorys = 1:8, y = c(2, 4, 1, 3, 8, 7, 5, 6), g = rep(c("a", "b"), 4))
  pp <- plot_correlation(collision, "categorys", "y", group = "g", scale = FALSE,
                          colors = c(a = "red", b = "blue"), add.regress = FALSE)
  expect_equal(pp$data$categorys, collision$categorys)
  expect_silent(ggplot2::ggplot_build(pp))
  for (method in c("pearson", "spearman", "kendall")) {
    p <- plot_correlation(d, 1, 2, method = method, scale = FALSE, add.regress = FALSE)
    expect_equal(attr(p, "correlation")$p_value,
                 stats::cor.test(d$x, d$y, method = method, exact = FALSE)$p.value)
    expect_silent(ggplot2::ggplot_build(p))
  }
  expect_warning(plot_correlation(data.frame(x = 1, y = 1), "x", "y"), "Insufficient")
  expect_error(plot_correlation(d, "missing", "y"), "not found")
})

test_that("pairwise HRs retain their direction and reference", {
  d <- survival_fixture()
  d <- d[d$group %in% c("A", "B"), ]
  d$group <- droplevels(d$group)
  out <- pairwise_survival(d, "group", reference = "A")
  fit <- survival::coxph(survival::Surv(time, status) ~ group, data = d)
  expect_equal(out$hr, unname(exp(stats::coef(fit))))
  expect_equal(out$reference, "A")
  reverse <- pairwise_survival(d, "group", reference = "B")
  expect_equal(reverse$hr, 1 / out$hr)
  expect_equal(reverse$conf_low, 1 / out$conf_high)
  expect_equal(reverse$logrank_p, out$logrank_p)
  expect_equal(nrow(pairwise_survival(d[d$group == "A", ], "group")), 0L)
  names(d) <- c("marker", "treatment group", "follow up", "event status")
  expect_equal(pairwise_survival(d, "treatment group", "follow up", "event status")$hr, out$hr)
})

test_that("survival plots retain risk tables, medians and named comparisons", {
  skip_if_not_installed("survminer")
  d <- survival_fixture()
  for (ng in 1:3) {
    sub <- d[d$group %in% LETTERS[seq_len(ng)], ]
    p <- plot_survival(sub, "group")
    expect_s3_class(p, "ggsurvplot")
    expect_s3_class(p$table, "ggplot")
    expect_equal(nrow(attr(p, "statistics")$medians), ng)
    expect_equal(nrow(attr(p, "statistics")$comparisons), choose(ng, 2))
    expect_equal(attr(p, "statistics")$n, nrow(sub))
  }
  two <- d[d$group %in% c("A", "B"), ]
  p <- plot_survival(two, "group", reference = "A")
  expect_equal(attr(p, "statistics")$comparisons$hr, pairwise_survival(two, "group")$hr)
  old <- surv_fig_hr("group", two, time = "time", status = "status", output_dir = NULL)
  expect_s3_class(old, "ggsurvplot")
})

test_that("new exports write readable results and composite survival plots", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("openxlsx")
  path <- tempfile()
  on.exit(unlink(path, recursive = TRUE))
  d <- survival_fixture()
  p <- plot_survival(d, "group", output_dir = path, formats = "pdf")
  expect_true(file.exists(file.path(path, "group.pdf")))
  expect_true(file.info(file.path(path, "group.pdf"))$size > 1000)
  expect_true(file.exists(file.path(path, "group_statistics.txt")))
  out <- batch_survival(d, "marker", output_dir = path)
  expect_equal(openxlsx::read.xlsx(file.path(path, "survival.xlsx"))$hr, out$results$hr)
  p <- plot_correlation(d, "marker", "time", output_dir = path, fig.format = "pdf", add.regress = FALSE)
  expect_true(file.exists(file.path(path, "1-time-marker-correlation.pdf")))
  expect_length(list.files(path, pattern = "RData$"), 0)
})
