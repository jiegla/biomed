#' @rdname biomed_standard_api
#' @export
plot_venn <- function(sets, filename = "venn", colors = NULL,
                      show_plot = TRUE, output_dir = NULL, ...) {
  result <- venn_plot(
    sets, venn_name = filename, col = colors, fill_col = colors, cat.col = colors,
    showFigure = show_plot, out_dir = if (is.null(output_dir)) "." else output_dir,
    save = !is.null(output_dir), write_xlsx = !is.null(output_dir), ...
  )
  invisible(result)
}
