\name{best_cutoff}
\alias{best_cutoff}
\title{Legacy survival cutoff selection}
\usage{
best_cutoff(pdata, variable, time = "time",
                            status = "status", print_result = TRUE,
                            minprop = 0.1,
                            fallback_cutoff = c("none", "median"),
                            status_encoding = "auto")
}
\arguments{
\item{pdata}{Data frame.}
\item{variable}{Predictor column name for cutoff selection; a vector of predictor names for batch_surv.}
\item{time}{Follow-up time column, containing nonnegative numeric values.}
\item{status}{Outcome column.}
\item{print_result}{Print the selected cutoff.}
\item{minprop}{Minimum fraction per cutoff group, in (0, 0.5].}
\item{fallback_cutoff}{One of "none" or "median". Median fallback selects the valid observed cutoff closest to the median only when optimized selection fails.}
\item{status_encoding}{One of "01" (0 censored, 1 event), "12" (1 censored, 2 event), "labels", or "auto". Auto tries 0/1 first: an all-1 cohort is treated as all events. Select "12" explicitly for an all-censored 1/2-coded cohort. Labels include alive/dead, censored/event, false/true, no/yes and progression/recurrence terms. Unknown nonmissing labels cause an error.}
}
\description{Find a maximally selected survival cutoff, with optional fallback to a valid observed value nearest the median. Prefer find_survival_cutoff for a standardized interface.}
\value{A list with pdata, best_cutoff, cutoff_method, reason and, on success, group_column. Clean complete observations are returned. Low includes values equal to the cutoff. Existing columns are preserved by choosing a unique group column name.}
\details{Data-driven cutoffs are exploratory. Cox P values and confidence intervals computed after selecting a cutoff are not adjusted for cutoff selection; BH correction across contrasts does not remove that selection bias. Validate selected cutoffs in independent data. Batch fit warnings are retained in the warning column.}
\examples{
\dontrun{
surv_data <- data.frame(
  time = c(6, 8, 12, 15, 18, 21, 24, 30, 36, 42),
  status = c(1, 1, 0, 1, 0, 1, 0, 1, 0, 0),
  biomarker = c(2.1, 3.0, 1.7, 5.2, 4.1, 6.3, 2.8, 7.1, 3.5, 5.9)
)
best_cutoff(surv_data, "biomarker", time = "time", status = "status")
}
}
