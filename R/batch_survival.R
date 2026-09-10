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
  out <- batch_surv(data, variables, time, status, best_cutoff = best_cutoff,
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
