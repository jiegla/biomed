#' Check required columns
#'
#' Verifies that a data frame contains the fields required by an analysis.
#' The error message lists every missing field so problems can be fixed in one
#' pass.
#'
#' @param data A data frame or data-frame-like object.
#' @param required A character vector of required column names.
#'
#' @return The input `data`, invisibly.
#' @export
#'
#' @examples
#' cohort <- data.frame(patient_id = c("P001", "P002"), response = c(1, 0))
#' biomed_check_columns(cohort, c("patient_id", "response"))
biomed_check_columns <- function(data, required) {
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  if (!is.character(required) || anyNA(required) || any(!nzchar(required))) {
    cli::cli_abort(
      "{.arg required} must be a character vector of non-empty names."
    )
  }

  missing <- setdiff(unique(required), names(data))

  if (length(missing) > 0L) {
    cli::cli_abort(c(
      "Required columns are missing.",
      "x" = "Missing: {missing}."
    ))
  }

  invisible(data)
}
