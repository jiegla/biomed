\name{draw_stack_barplot}
\alias{draw_stack_barplot}
\title{Draw proportional stacked bar charts}
\usage{
draw_stack_barplot(dataframe, x_var, y_vars, out_dir = "StackBar",
  color = NULL, test_method = c("auto", "chi_square", "fisher"),
  show_label = TRUE, show_n = TRUE, show_p = TRUE,
  percent_accuracy = 0.1, palette_name = NULL, width = 4.5,
  height = 5.5, dpi = 600, horizontal = FALSE, remove_na = TRUE,
  return_plot = FALSE, theme_size = 16, save = FALSE,
  statistics_format = c("xlsx", "csv"))
}
\arguments{
  \item{dataframe}{A data frame.}
  \item{x_var}{Grouping column.}
  \item{y_vars}{Categorical outcome columns.}
  \item{out_dir}{Output directory.}
  \item{color}{Optional colour vector.}
  \item{test_method}{Use automatic, chi-square, or Fisher testing.}
  \item{show_label, show_n, show_p}{Control plot annotations.}
  \item{percent_accuracy}{Accuracy passed to \code{scales::percent()}.}
  \item{palette_name}{Deprecated compatibility argument.}
  \item{width, height}{Figure dimensions in inches.}
  \item{dpi}{PNG resolution.}
  \item{horizontal}{Flip coordinates.}
  \item{remove_na}{Remove incomplete observations.}
  \item{return_plot}{Return plots alongside the statistics.}
  \item{theme_size}{Base font size.}
  \item{save}{Save PDF and PNG figures and a statistics table.}
  \item{statistics_format}{Save statistics as XLSX or CSV.}
}
\value{A statistics data frame, optionally bundled with named ggplot objects.}
\description{Draws proportional stacked bars and tests categorical association.}
\examples{
d <- data.frame(arm = rep(c("A", "B"), each = 4),
  response = c("CR", "PR", "SD", "PD", "CR", "CR", "PR", "SD"))
draw_stack_barplot(d, "arm", "response", return_plot = TRUE)
}
