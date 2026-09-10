\name{plot_correlation}
\alias{plot_correlation}
\title{Correlation plots with consistent filtering}
\usage{
plot_correlation(data, x, y, group = NULL,
                              method = c("spearman", "pearson", "kendall"),
                              scale = TRUE, colors = NULL, remove_group_na = FALSE,
                              remove_x_zero = FALSE, remove_y_zero = FALSE,
                              show_plot = FALSE, output_dir = NULL,
                              save_data = FALSE, ...)
}
\arguments{
\item{data}{Data frame; plot_correlation also accepts a matrix.}
\item{x}{Character file names for sanitize_filename; an x-axis column name or position for plot_correlation.}
\item{y}{Y-axis column name or position.}
\item{group}{Group column; optional for correlation.}
\item{method}{Correlation method: spearman, pearson or kendall.}
\item{scale}{Standardize x and y after filtering.}
\item{colors}{Optional group colors, optionally named.}
\item{remove_group_na}{Drop missing group values before calculation and plotting.}
\item{remove_x_zero}{Remove original x values equal to zero before scaling.}
\item{remove_y_zero}{Remove original y values equal to zero before scaling.}
\item{show_plot}{Print the plot.}
\item{output_dir}{Output directory; NULL writes no files. The legacy surv_fig_hr default writes a text log in Surv_Output.}
\item{save_data}{Also export analyzed plot data as RData. Legacy get_cor_jgl defaults to TRUE when saving a plot; plot_correlation defaults to FALSE.}
\item{...}{Additional plot styling and export arguments passed to the underlying plotting function. Do not repeat arguments already mapped by the wrapper.}
}
\description{Compute correlation and plot the same filtered observations. Supports Pearson, Spearman and Kendall tests, group colors, labels, regression lines, density contours and optional export.}
\value{A ggplot invisibly, with a correlation attribute containing x, y, method, n, estimate and p_value. Insufficient or constant data warns and returns NULL.}
\details{Correlation is computed across all retained observations. In grouped plots, regression lines are fitted separately within each group. Numeric factor labels are converted to their numeric values. Non-numeric values and infinities cause an error.}
