#' @rdname biomed_standard_api
#' @export
batch_roc <- function(data, variables, response, method = "best",
                      positive_class = NULL, direction = "auto",
                      conf_level = 0.95, output_dir = NULL) {
  result <- roc_batch(variables, data, response, method, positive_class,
                      direction, conf_level)
  names(result) <- tolower(names(result))
  result <- .biomed_standard_table(result, c(
    markers = "variable", auc_ci_lower = "conf_low", auc_ci_upper = "conf_high"
  ))
  .biomed_save_table(result, output_dir, "roc.xlsx")
}
