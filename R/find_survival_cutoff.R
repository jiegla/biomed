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
  out <- best_cutoff(data, variable, time, status, print_result = verbose,
                         minprop = min_prop, fallback_cutoff = match.arg(fallback),
                         status_encoding = status_encoding)
  names(out)[names(out) == "pdata"] <- "data"
  names(out)[names(out) == "best_cutoff"] <- "cutoff"
  names(out)[names(out) == "cutoff_method"] <- "method"
  out
}
