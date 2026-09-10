#' @rdname biomed_standard_api
#' @export
batch_cox <- function(data, variables, time = "PFS", status = "PFS_status",
                      conf_level = 0.95, output_dir = NULL) {
  result <- cox_batch(variables, data, time, status, conf_level)
  result <- .biomed_standard_table(result, c(
    Marker = "variable", Term = "term", Level = "level", Reference = "reference",
    N = "n", Events = "events", HR = "hr", CI_lower = "conf_low",
    CI_upper = "conf_high", Error = "error"
  ))
  .biomed_save_table(result, output_dir, "cox.xlsx")
}
