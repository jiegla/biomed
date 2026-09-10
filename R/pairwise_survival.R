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
