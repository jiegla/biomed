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
#' @param ... Additional styling and export arguments to get_cor.
#' @return A ggplot, invisibly, with a correlation attribute containing n,
#'   estimate and P value. Insufficient or constant data warns and returns NULL.
#' @export
plot_correlation <- function(data, x, y, group = NULL,
                              method = c("spearman", "pearson", "kendall"),
                              scale = TRUE, colors = NULL, remove_group_na = FALSE,
                              remove_x_zero = FALSE, remove_y_zero = FALSE,
                              show_plot = FALSE, output_dir = NULL, save_data = FALSE, ...) {
  get_cor(df = data, var1 = x, var2 = y, subtype = group,
              method = match.arg(method), scale = scale, color_subtype = colors,
              na.subtype.rm = remove_group_na, remove.x.zero = remove_x_zero,
              remove.y.zero = remove_y_zero, show_plot = show_plot,
              save_plot = !is.null(output_dir), path = output_dir, save_data = save_data, ...)
}
