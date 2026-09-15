\name{analysis_utilities}
\alias{makegroup}
\alias{anova_batch}
\alias{chi_square_batch_df}
\alias{cox_batch}
\alias{roc_batch}
\alias{plot_roc_batch}
\title{Clinical and biomarker analysis utilities}
\usage{
makegroup(pdata, variables, num_group = 2,
  method = c("median", "survival_best", "roc_youden", "roc_closest",
    "fixed_sensitivity", "fixed_specificity"), time = "time",
  status = "status", response = "group", fixed_value = 0.8,
  positive_class = NULL, direction = "auto",
  add_method_suffix = FALSE, verbose = FALSE)

anova_batch(data, target = "group", feature = NULL,
  feature_manipulation = FALSE)

chi_square_batch_df(pdata, variables, response = "status",
  test_method = c("auto", "chi_square", "fisher"), correct = FALSE)

cox_batch(var, pdata, time = "PFS", status = "PFS_status",
  conf_level = 0.95)

roc_batch(vars, pdata, outcome, method = c("best", "median"),
  positive_class = "1", direction = "auto", conf_level = 0.95)

plot_roc_batch(data, response, variables, fig.path = NULL, main = NULL,
  alpha = 1, smooth = FALSE, positive_class = NULL, direction = "auto",
  file_type = c("png", "pdf"), width = 5, height = 5, dpi = 300)
}
\arguments{
  \item{pdata, data}{A data frame.}
  \item{variables, vars, var, feature}{Columns to analyse.}
  \item{num_group}{Create two or three groups.}
  \item{method}{Cutoff or threshold method.}
  \item{time, status}{Survival time and event columns.}
  \item{response, outcome, target}{Outcome or grouping column.}
  \item{fixed_value}{Target sensitivity or specificity.}
  \item{positive_class}{Outcome value treated as the positive event.}
  \item{direction}{ROC direction passed to \code{pROC::roc()}.}
  \item{add_method_suffix}{Include the method in new column names.}
  \item{verbose}{Print cutoff summaries.}
  \item{feature_manipulation}{Filter incomplete, non-numeric, infinite and constant features.}
  \item{test_method}{Association-test selection.}
  \item{correct}{Use Yates' correction for 2-by-2 chi-square tests.}
  \item{conf_level}{Confidence level.}
  \item{fig.path}{Optional figure output directory.}
  \item{main}{Optional plot title prefix.}
  \item{alpha}{ROC line opacity.}
  \item{smooth}{Smooth ROC curves.}
  \item{file_type}{Figure file formats.}
  \item{width, height}{Figure dimensions in inches.}
  \item{dpi}{PNG resolution.}
}
\value{
\code{makegroup()} returns the input with grouping columns and a cutoff-table
attribute. The batch analysis functions return data frames. \code{plot_roc_batch()}
returns an invisible named list of ggplot objects.
}
\description{
Utilities for reproducible grouping, one-way ANOVA, categorical association
tests, univariable Cox regression, and ROC analysis.
}
\examples{
clinical <- data.frame(
  group = factor(rep(c("A", "B"), each = 10)),
  response = factor(rep(c(0, 1), 10)),
  time = seq_len(20),
  status = rep(c(1, 0, 1, 1, 0), 4),
  marker1 = c(1:10, 5:14),
  marker2 = c(10:1, 8:17)
)

# Create Low/High groups from a continuous biomarker
makegroup(clinical, "marker1", method = "median")

# Batch ANOVA and categorical association
anova_batch(clinical, target = "group", feature = c("marker1", "marker2"))
chi_square_batch_df(clinical, "group", response = "response")

# Univariable Cox regression
cox_batch(c("marker1", "marker2"), clinical, time = "time", status = "status")

\dontrun{
# ROC summaries and plots require optional ROC/plotting dependencies
roc_batch(c("marker1", "marker2"), clinical, outcome = "response")
plot_roc_batch(clinical, response = "response",
  variables = c("marker1", "marker2"), fig.path = "ROC")
}
}
