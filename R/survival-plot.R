.biomed_survival_data <- function(data, group, time, status, status_encoding, reference = NULL) {
  .biomed_required_columns(data, c(group, time, status))
  d <- data.frame(.time = .biomed_as_numeric(data[[time]], time),
                  .status = .status01_jgl(data[[status]], status_encoding),
                  .group = data[[group]])
  if (any(d$.time < 0, na.rm = TRUE)) stop("Time must be nonnegative.")
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (!nrow(d)) stop("No complete survival observations.")
  d$.group <- droplevels(factor(d$.group))
  if (!is.null(reference)) {
    if (length(reference) != 1L || !reference %in% levels(d$.group)) stop("reference must name an observed group.")
    d$.group <- stats::relevel(d$.group, reference)
  }
  d
}

.biomed_pairwise_survival <- function(d, p_adjust = "BH") {
  lv <- levels(d$.group)
  empty <- data.frame(group = character(), reference = character(), n = integer(),
                      events = integer(), hr = double(), conf_low = double(),
                      conf_high = double(), p_value = double(), logrank_p = double(),
                      warning = character(), p_adjusted = double())
  if (length(lv) < 2L) return(empty)
  pairs <- utils::combn(lv, 2L, simplify = FALSE)
  result <- lapply(pairs, function(pair) {
    sub <- d[d$.group %in% pair, , drop = FALSE]
    sub$.group <- factor(sub$.group, levels = pair)
    warnings <- character()
    capture <- function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
    fit <- withCallingHandlers(tryCatch(
      survival::coxph(survival::Surv(.time, .status) ~ .group, data = sub),
      error = function(e) { warnings <<- c(warnings, conditionMessage(e)); NULL }),
      warning = capture)
    vals <- rep(NA_real_, 4L)
    if (!is.null(fit)) {
      s <- summary(fit)
      vals <- c(s$conf.int[1L, "exp(coef)"], s$conf.int[1L, "lower .95"],
                s$conf.int[1L, "upper .95"], s$coefficients[1L, "Pr(>|z|)"])
    }
    lp <- withCallingHandlers(tryCatch({
      lr <- survival::survdiff(survival::Surv(.time, .status) ~ .group, data = sub)
      stats::pchisq(lr$chisq, 1L, lower.tail = FALSE)
    }, error = function(e) NA_real_), warning = capture)
    data.frame(group = pair[2L], reference = pair[1L], n = nrow(sub),
               events = sum(sub$.status), hr = vals[1L], conf_low = vals[2L],
               conf_high = vals[3L], p_value = vals[4L], logrank_p = lp,
               warning = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_)
  })
  result <- do.call(rbind, result)
  result$p_adjusted <- stats::p.adjust(result$logrank_p, method = p_adjust)
  rownames(result) <- NULL
  result
}

#' Pairwise survival comparisons with explicit HR direction
#' @param data Data frame.
#' @param group,time,status Column names for groups, follow-up time and outcome.
#' @param reference Optional group to place first (the comparison reference).
#' @param status_encoding One of "01", "12", "labels", or "auto".
#' @param p_adjust Adjustment method for pairwise log-rank P values.
#' @return Data frame with directed HRs, Cox Wald P values, log-rank P values,
#'   adjusted log-rank P values and fit warnings. HR is group versus reference.
#' @export
pairwise_survival <- function(data, group, time = "time", status = "status",
                              reference = NULL, status_encoding = "01", p_adjust = "BH") {
  d <- .biomed_survival_data(data, group, time, status, status_encoding, reference)
  .biomed_pairwise_survival(d, p_adjust)
}

.biomed_save_survival <- function(plot, output_dir, filename, formats, width, height) {
  if (is.null(output_dir)) return(invisible(NULL))
  formats <- match.arg(formats, c("pdf", "png"), several.ok = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (ext in formats) {
    file <- file.path(output_dir, paste0(sanitize_filename(filename), ".", ext))
    if (ext == "pdf") grDevices::pdf(file, width = width, height = height)
    else grDevices::png(file, width = width, height = height, units = "in", res = 150)
    tryCatch(print(plot), finally = grDevices::dev.off())
  }
  stats <- attr(plot, "statistics")
  writeLines(utils::capture.output(print(stats)),
             file.path(output_dir, paste0(sanitize_filename(filename), "_statistics.txt")))
  invisible(NULL)
}

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
  fit <- if (ng > 1L) survival::survfit(survival::Surv(.time, .status) ~ .group, data = d)
         else survival::survfit(survival::Surv(.time, .status) ~ 1, data = d)
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

#' Legacy survival plot interface
#' @param group_var,df,max_time,time,status,hr,cols,x.label Legacy arguments.
#' @param output_dir Directory for the legacy text log; NULL disables writing.
#' @param ... Additional arguments to plot_survival.
#' @return A ggsurvplot object with a statistics attribute.
#' @export
surv_fig_hr <- function(group_var, df, max_time = 60, time = "OS",
                        status = "Survival_status", hr = "yes", cols = NULL,
                        x.label = "Months", output_dir = "Surv_Output", ...) {
  if (is.null(cols) || !length(cols) || all(is.na(cols))) cols <- NULL
  else cols <- rep(cols, length.out = length(unique(stats::na.omit(df[[group_var]]))))
  gp <- plot_survival(df, group_var, time, status, max_time = max_time,
                      show_hr = identical(hr, "yes"), colors = cols, x_label = x.label,
                      status_encoding = "auto", sort_by_median = TRUE, ...)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    file <- file.path(output_dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"),
                                         sanitize_filename(group_var), "_Output.txt"))
    writeLines(utils::capture.output(print(attr(gp, "statistics"))), file)
  }
  gp
}
