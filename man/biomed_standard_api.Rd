\name{biomed_standard_api}
\alias{biomed_standard_api}
\alias{make_group}
\alias{batch_anova}
\alias{batch_chi_square}
\alias{batch_cox}
\alias{batch_roc}
\alias{plot_roc}
\alias{plot_donut}
\alias{plot_stacked_bar}
\alias{plot_venn}
\title{Standard clinical analysis interfaces}
\usage{
make_group(data, variables, n_groups = 2, method = "median", ...)
batch_anova(data, variables = NULL, group = "group",
  preprocess = FALSE, output_dir = NULL)
batch_chi_square(data, variables, response = "status",
  test_method = "auto", correct = FALSE, output_dir = NULL)
batch_cox(data, variables, time = "PFS", status = "PFS_status",
  conf_level = 0.95, output_dir = NULL)
batch_roc(data, variables, response, method = "best",
  positive_class = NULL, direction = "auto", conf_level = 0.95,
  output_dir = NULL)
plot_roc(data, variables, response, output_dir = NULL,
  positive_class = NULL, direction = "auto", ...)
plot_donut(data, variables, colors = NULL, remove_na = TRUE,
  output_dir = NULL, ...)
plot_stacked_bar(data, variables, group, colors = NULL,
  remove_na = TRUE, output_dir = NULL, ...)
plot_venn(sets, filename = "venn", colors = NULL, show_plot = TRUE,
  output_dir = NULL, ...)
}
\arguments{
\item{data}{A data frame.}
\item{variables}{Columns to analyse.}
\item{group}{Grouping column.}
\item{response}{Binary or categorical outcome column.}
\item{time, status}{Survival time and event columns; status must be 0 (censored) or 1 (event).}
\item{method}{Grouping or ROC threshold method.}
\item{n_groups}{Number of groups, two or three.}
\item{preprocess}{Filter incomplete, non-numeric, infinite and constant features.}
\item{positive_class}{Positive outcome class; NULL selects the second observed factor level or sorted label.}
\item{direction}{ROC direction.}
\item{conf_level}{Confidence level.}
\item{output_dir}{Optional output directory; NULL writes no files.}
\item{test_method}{Automatic, chi-square or Fisher testing.}
\item{correct}{Apply Yates' correction.}
\item{colors}{Colour vector, or for donuts a list by variable.}
\item{remove_na}{Remove missing values.}
\item{sets}{Named list of sets.}
\item{filename}{Venn filename stem.}
\item{show_plot}{Draw the Venn diagram.}
\item{...}{Additional options passed to the corresponding original function.}
}
\value{Analysis functions return data frames with snake_case columns. Grouping
returns augmented data with a cutoff_table attribute. Donut and ROC plots return
named lists; stacked bars return plots and statistics. Venn returns the drawing,
partition table and output paths.}
\description{Consistent interfaces for the nine original analysis utilities.
Original entry points remain available.}
\details{Batch tables can be exported as XLSX. Figures use PDF/PNG by default;
Venn also supports TIFF. Optimized cutoffs and same-data ROC estimates are
exploratory and require independent validation. Additional options retain the
original names when passed through \code{...}. See README for migration details.}
\examples{
d <- data.frame(group = rep(c("A", "B"), each = 5), marker = c(1:5, 3:7))
batch_anova(d, "marker")
make_group(d, "marker")
}
