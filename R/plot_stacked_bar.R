#' @rdname biomed_standard_api
#' @export
plot_stacked_bar <- function(data, variables, group, colors = NULL,
                             remove_na = TRUE, output_dir = NULL, ...) {
  result <- draw_stack_barplot(
    data, group, variables, color = colors, remove_na = remove_na,
    out_dir = if (is.null(output_dir)) "StackBar" else output_dir,
    save = !is.null(output_dir), return_plot = TRUE, ...
  )
  names(result)[names(result) == "plot"] <- "plots"
  result$statistics <- .biomed_standard_table(result$statistics, c(N = "n"))
  result
}
