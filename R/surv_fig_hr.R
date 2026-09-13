#' Legacy survival plot interface
#' @param group_var,df,max_time,time,status,hr,cols,x.label Legacy arguments.
#' @param output_dir Directory for the legacy text log; NULL disables writing.
#' @param ... Additional arguments to plot_survival.
#' @return A ggsurvplot object with a statistics attribute.
#' @export
surv_fig_hr <- function(group_var, df, max_time = 60, time = "OS",
                        status = "Survival_status", hr = "yes", cols = NULL,
                        x.label = "Months", output_dir = "Surv_Output", ...) {
  column_args <- list(group_var = group_var, time = time, status = status)
  valid_column <- vapply(column_args, function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  }, logical(1))
  if (any(!valid_column)) {
    bad <- names(column_args)[!valid_column]
    stop("Column-name arguments must be nonempty. Invalid argument(s): ",
         paste(bad, collapse = ", "), ".")
  }
  .biomed_required_columns(df, unname(unlist(column_args, use.names = FALSE)))
  if (is.null(cols) || !length(cols) || all(is.na(cols))) {
    cols <- NULL
  } else {
    cols <- rep(cols, length.out = length(unique(stats::na.omit(df[[group_var]]))))
  }
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

#' Draw a complete survival plot on the current graphics page
#'
#' Registered automatically for grid drawing and ggplot2::ggsave. Uses the
#' survminer print method through standard S3 dispatch, retaining risk tables.
#' @param x A survminer ggsurvplot object.
#' @param recording Compatibility argument for grid.draw.
#' @rdname surv_fig_hr
#' @importFrom grid grid.draw
#' @export
grid.draw.ggsurvplot <- function(x, recording = TRUE) {
  .biomed_require("survminer")
  print(x, newpage = FALSE)
  invisible(x)
}
