\name{get_cor}
\alias{get_cor}
\title{Legacy correlation plot}
\usage{
get_cor(
    df = NULL,
    var1,
    var2,
    is.matrix = FALSE,
    scale = TRUE,
    subtype = NULL,
    na.subtype.rm = FALSE,
    color_subtype = NULL,
    palette = "biomed",
    index = NULL,
    method = c("spearman", "pearson", "kendall"),
    show_cor_result = FALSE,
    col_line = NULL,
    id = NULL,
    show_label = FALSE,
    point_size = 4,
    title = NULL,
    alpha = 0.5,
    title_size = 1.5,
    text_size = 10,
    axis_angle = 0,
    hjust = 0,
    show_plot = TRUE,
    save_plot = FALSE,
    path = NULL,
    fig.format = "png",
    fig.width = 7,
    fig.height = 7.3,
    add.hdr.line = FALSE,
    remove.x.zero = FALSE,
    remove.y.zero = FALSE,
    point.fill = "lightblue",
    point.color = "black",
    point.shape = 21,
    stroke = 0.5,
    add.regress = TRUE,
    save_data = TRUE
)
}
\arguments{
\item{df}{Data frame or matrix for correlation; data frame for survival.}
\item{var1}{X-axis column name or numeric position.}
\item{var2}{Y-axis column name or numeric position.}
\item{is.matrix}{Compatibility argument; matrix input is converted to a data frame.}
\item{scale}{Standardize x and y after filtering.}
\item{subtype}{Optional grouping column for colors and group-specific regression lines.}
\item{na.subtype.rm}{Drop missing group values before correlation calculation and plotting.}
\item{color_subtype}{Optional group colors, optionally named.}
\item{palette}{Built-in biomed (default), or a ggsci palette: jama, npg, aaas, nejm, lancet, jco or d3; used when group colors are not supplied.}
\item{index}{Optional output filename prefix; defaults to 1.}
\item{method}{Correlation method: spearman, pearson or kendall.}
\item{show_cor_result}{Print the cor.test result.}
\item{col_line}{Regression-line color; default depends on overall correlation sign.}
\item{id}{Optional text-label column.}
\item{show_label}{Show text labels using id.}
\item{point_size}{Point size.}
\item{title}{Plot title; defaults to the compared column names.}
\item{alpha}{Point opacity.}
\item{title_size}{Title size relative to text_size.}
\item{text_size}{Base text size.}
\item{axis_angle}{X-axis text angle.}
\item{hjust}{Horizontal justification for title and x-axis text.}
\item{show_plot}{Print the plot.}
\item{save_plot}{Save the plot when path is supplied.}
\item{path}{Optional plot output directory.}
\item{fig.format}{ggsave output format, e.g. png or pdf.}
\item{fig.width}{Figure width in inches.}
\item{fig.height}{Figure height in inches.}
\item{add.hdr.line}{Add highest-density-region contours using optional ggdensity.}
\item{remove.x.zero}{Remove original x values equal to zero before scaling.}
\item{remove.y.zero}{Remove original y values equal to zero before scaling.}
\item{point.fill}{Ungrouped point fill.}
\item{point.color}{Ungrouped point outline color.}
\item{point.shape}{Ungrouped point shape.}
\item{stroke}{Ungrouped point outline width.}
\item{add.regress}{Add linear regression line(s) with confidence bands, regardless of correlation method; grouped plots fit a line within each group.}
\item{save_data}{Also export analyzed plot data as RData. Legacy get_cor defaults to TRUE when saving a plot; plot_correlation defaults to FALSE.}
}
\description{Compute correlation and plot the same filtered observations. Original zero values and missing groups are filtered before scaling and testing. Prefer plot_correlation for standardized names.}
\value{A ggplot invisibly, with a correlation attribute containing x, y, method, n, estimate and p_value. Insufficient or constant data warns and returns NULL. NULL df returns NULL.}
\details{Correlation is computed across all retained observations. In grouped plots, regression lines are fitted separately within each group. Numeric factor labels are converted to their numeric values. Non-numeric values and infinities cause an error.}
\examples{
\dontrun{
biomarker <- data.frame(
  PD_L1 = c(1, 5, 10, 20, 30, 50, 70, 90),
  TMB = c(2.1, 2.4, 3.2, 4.0, 5.1, 6.2, 7.0, 8.4),
  arm = rep(c("A", "B"), 4)
)
get_cor(biomarker, var1 = "PD_L1", var2 = "TMB", subtype = "arm",
  scale = FALSE, show_plot = FALSE, save_plot = FALSE)
}
}
