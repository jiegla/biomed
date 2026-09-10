#' Plot Kaplan-Meier curves, risk tables and directed hazard ratios
#' @param data Data frame.
#' @param group,time,status Column names.
#' @param max_time Optional upper display limit; does not truncate model fitting.
#' @param reference Optional reference group; otherwise the first factor level.
#' @param status_encoding Outcome coding: "01", "12", "labels", or "auto".
#' @param show_hr Show directed two-group HR and confidence interval.
#' @param colors Optional group colors, optionally named.
#' @param x_label Time-axis units.
#' @param risk_table Include the number-at-risk table.
#' @param sort_by_median Order display groups by decreasing median survival.
#'   An explicit reference is kept first; comparisons always name their direction.
#' @param show_plot Print the plot.
#' @param output_dir Optional output directory; NULL writes no files.
#' @param filename Output base name.
#' @param formats Output formats: "pdf", "png", or both.
#' @param width,height Export size in inches, including the risk table.
#' @return A survminer ggsurvplot object. Attribute "statistics" contains medians,
#'   pairwise comparisons, global log-rank P value, sample and event counts.
#' @export
plot_survival <- function(data, group, time = "time", status = "status",
                          max_time = NULL, reference = NULL, status_encoding = "01",
                          show_hr = TRUE, colors = NULL, x_label = "Months",
                          risk_table = TRUE, sort_by_median = FALSE, show_plot = FALSE,
                          output_dir = NULL, filename = group, formats = c("pdf", "png"),
                          width = 8, height = 7) {
  .biomed_require("survminer")
  d <- .biomed_survival_data(data, group, time, status, status_encoding, reference)
  ng <- nlevels(d$.group)
  if (is.null(max_time)) max_time <- max(d$.time)
  if (length(max_time) != 1L || !is.numeric(max_time) || !is.finite(max_time) || max_time <= 0) {
    stop("max_time must be a positive finite number.")
  }
  fit <- if (ng > 1L) {
    survival::survfit(survival::Surv(.time, .status) ~ .group, data = d)
  } else {
    survival::survfit(survival::Surv(.time, .status) ~ 1, data = d)
  }
  medians <- survminer::surv_median(fit)
  if (ng > 1L && sort_by_median) {
    lv <- levels(d$.group)[order(-medians$median, na.last = TRUE)]
    if (!is.null(reference)) lv <- c(reference, setdiff(lv, reference))
    d$.group <- factor(d$.group, levels = lv)
    fit <- survival::survfit(survival::Surv(.time, .status) ~ .group, data = d)
    medians <- survminer::surv_median(fit)
  }
  medians$strata <- levels(d$.group)
  comparisons <- .biomed_pairwise_survival(d)
  logrank_p <- if (ng > 1L && sum(d$.status) > 0L) {
    lr <- survival::survdiff(survival::Surv(.time, .status) ~ .group, data = d)
    stats::pchisq(lr$chisq, ng - 1L, lower.tail = FALSE)
  } else NA_real_
  color_values <- unname(.biomed_colour_values(levels(d$.group), colors))
  # The formula stored by survfit remains evaluable by survminer.
  gp <- survminer::ggsurvplot(
    fit, data = d, pval = if (is.finite(logrank_p)) paste0("Log-rank P = ", format.pval(logrank_p, digits = 3)) else FALSE,
    conf.int = FALSE, risk.table = risk_table, risk.table.y.text = FALSE,
    risk.table.title = "At risk", risk.table.col = "strata",
    surv.median.line = if (any(is.finite(medians$median))) "hv" else "none",
    ggtheme = ggplot2::theme_classic(), tables.theme = survminer::theme_cleantable(),
    break.time.by = max_time / 5, xlim = c(0, max_time),
    xlab = paste0(time, " (", x_label, ")"), palette = color_values,
    legend.labs = levels(d$.group), legend.title = group, title = group,
    risk.table.height = 0.22, surv.scale = "percent")
  fmt <- function(x) ifelse(is.finite(x), sprintf("%.1f", x), "NR")
  labels <- paste0(medians$strata, ": ", fmt(medians$median), " (95% CI ",
                   fmt(medians$lower), "-", fmt(medians$upper), ")")
  if (show_hr && ng == 2L) {
    z <- comparisons[1L, ]
    hr_label <- sprintf("%s vs %s: HR %.2f (95%% CI %.2f-%.2f)",
                        z$group, z$reference, z$hr, z$conf_low, z$conf_high)
    labels <- c(hr_label, labels)
  }
  gp$plot <- gp$plot + ggplot2::annotate("text", x = max_time * .98, y = .85,
                                        label = paste(labels, collapse = "\n"),
                                        hjust = 1, vjust = 1, size = 3)
  if (ng == 1L) {
    gp$plot <- gp$plot + ggplot2::theme(legend.position = "none")
    if (risk_table) gp$table <- gp$table + ggplot2::theme(legend.position = "none")
  }
  attr(gp, "statistics") <- list(medians = medians, comparisons = comparisons,
                                  logrank_p = logrank_p, n = nrow(d), events = sum(d$.status),
                                  excluded = nrow(data) - nrow(d))
  attr(gp, "survival_fit") <- fit
  .biomed_save_survival(gp, output_dir, filename, formats, width, height)
  if (show_plot) print(gp)
  gp
}
