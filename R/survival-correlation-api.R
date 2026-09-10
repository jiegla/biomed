#' Find a survival cutoff with an optional median fallback
#' @param data Data frame.
#' @param variable Numeric predictor column.
#' @param time,status Follow-up and outcome column names.
#' @param min_prop Minimum fraction in each cutoff group, in (0, 0.5].
#' @param fallback One of "none" or "median"; median chooses the closest valid observed cutoff.
#' @param status_encoding One of "01", "12", "labels", or "auto".
#' @param verbose Print the selected cutoff.
#' @return List containing data, cutoff, method, group_column and reason.
#' @export
find_survival_cutoff <- function(data, variable, time = "time", status = "status",
                                 min_prop = .1, fallback = c("none", "median"),
                                 status_encoding = "01", verbose = FALSE) {
  out <- best_cutoff_jgl(data, variable, time, status, print_result = verbose,
                         minprop = min_prop, fallback_cutoff = match.arg(fallback),
                         status_encoding = status_encoding)
  names(out)[names(out) == "pdata"] <- "data"
  names(out)[names(out) == "best_cutoff"] <- "cutoff"
  names(out)[names(out) == "cutoff_method"] <- "method"
  out
}

#' Batch survival analysis with screening and failure reports
#' @param data Data frame.
#' @param variables Predictor column names.
#' @param time,status Follow-up and outcome column names.
#' @param best_cutoff Optimize a binary cutoff before fitting each model.
#' @param min_sample,min_event,min_censored Overall and per-predictor thresholds.
#' @param min_prop Minimum cutoff group fraction, in (0, 0.5].
#' @param min_event_group,min_censored_group Minimum group counts in cutoff mode.
#' @param fallback Cutoff fallback: "none" or "median".
#' @param p_adjust Apply BH adjustment to Cox Wald P values over all returned contrasts.
#' @param status_encoding One of "01", "12", "labels", or "auto".
#' @param verbose Show progress.
#' @param output_dir Optional directory for survival.xlsx.
#' @return List with results, failed and cutoffs tables, using snake_case columns.
#' @export
batch_survival <- function(data, variables, time = "time", status = "status",
                            best_cutoff = FALSE, min_sample = 20, min_event = 5,
                            min_censored = 5, min_prop = .1, min_event_group = 1,
                            min_censored_group = 0, fallback = c("none", "median"),
                            p_adjust = TRUE, status_encoding = "01", verbose = FALSE,
                            output_dir = NULL) {
  out <- batch_surv_jgl(data, variables, time, status, best_cutoff = best_cutoff,
                        min_sample = min_sample, min_event = min_event,
                        min_censored = min_censored, minprop = min_prop,
                        min_event_group = min_event_group,
                        min_censored_group = min_censored_group,
                        fallback_cutoff = match.arg(fallback), p_adjust = p_adjust,
                        return_failed = TRUE, verbose = verbose,
                        status_encoding = status_encoding)
  mapping <- c(ID = "variable", HR = "hr", lower95 = "conf_low", upper95 = "conf_high",
                P = "p_value", N = "n", Event = "events", Censored = "censored",
                FDR = "p_adjusted", high = "n_high", low = "n_low",
                high_m = "median_high", low_m = "median_low")
  result <- list(results = .biomed_standard_table(out$result, mapping),
                  failed = .biomed_standard_table(out$failed, mapping),
                  cutoffs = if (is.null(out$cutoff)) data.frame() else .biomed_standard_table(out$cutoff, mapping))
  attr(result$results, "failed") <- attr(result$results, "cutoff") <- NULL
  .biomed_save_table(result, output_dir, "survival.xlsx")
}

#' Plot correlation with consistent filtering and statistics
#' @param data Data frame or matrix.
#' @param x,y Column names or positions.
#' @param group Optional grouping column for colors and regression lines.
#' @param method Correlation method: spearman, pearson or kendall.
#' @param scale Standardize both columns after filtering.
#' @param colors Optional group colors.
#' @param remove_group_na Drop missing groups before calculation.
#' @param remove_x_zero,remove_y_zero Remove original zero values before scaling.
#' @param show_plot Print the plot.
#' @param output_dir Optional output directory.
#' @param save_data Save analyzed plot data as RData in addition to the figure.
#' @param ... Additional styling and export arguments to get_cor_jgl.
#' @return A ggplot, invisibly, with a correlation attribute containing n,
#'   estimate and P value. Insufficient or constant data warns and returns NULL.
#' @export
plot_correlation <- function(data, x, y, group = NULL,
                              method = c("spearman", "pearson", "kendall"),
                              scale = TRUE, colors = NULL, remove_group_na = FALSE,
                              remove_x_zero = FALSE, remove_y_zero = FALSE,
                              show_plot = FALSE, output_dir = NULL, save_data = FALSE, ...) {
  get_cor_jgl(df = data, var1 = x, var2 = y, subtype = group,
              method = match.arg(method), scale = scale, color_subtype = colors,
              na.subtype.rm = remove_group_na, remove.x.zero = remove_x_zero,
              remove.y.zero = remove_y_zero, show_plot = show_plot,
              save_plot = !is.null(output_dir), path = output_dir, save_data = save_data, ...)
}
