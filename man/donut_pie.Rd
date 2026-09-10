\name{donut_pie}
\alias{donut_pie}
\title{Draw donut charts for categorical variables}
\usage{
donut_pie(df, vars, blank = 0.8, linewidth = 2.5, showtext = TRUE,
  text_type = c("both", "percent", "count"), show_center_n = TRUE,
  removeNA = TRUE, na_label = "NA", color_list = NULL, legend = TRUE)
}
\arguments{
  \item{df}{A data frame.}
  \item{vars}{Categorical columns to plot.}
  \item{blank}{Inner-radius control in the interval \code{[0, 2.5)}.}
  \item{linewidth}{Width of borders between slices.}
  \item{showtext}{Show labels in slices.}
  \item{text_type}{Label counts, percentages, or both.}
  \item{show_center_n}{Show analysed sample count in the centre.}
  \item{removeNA}{Remove missing values.}
  \item{na_label}{Displayed missing-value label.}
  \item{color_list}{A shared colour vector or named list by variable.}
  \item{legend}{Display a legend.}
}
\value{A named list of ggplot objects.}
\description{Creates one donut chart for each requested categorical variable.}
\examples{
d <- data.frame(response = factor(c("CR", "PR", "SD", "PR")))
donut_pie(d, "response")
}
