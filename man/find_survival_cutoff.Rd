\name{find_survival_cutoff}
\alias{find_survival_cutoff}
\title{Find a survival cutoff}
\usage{
find_survival_cutoff(data, variable, time = "time", status = "status",
                                 min_prop = .1, fallback = c("none", "median"),
                                 status_encoding = "01", verbose = FALSE)
}
\arguments{
\item{data}{Data frame; plot_correlation also accepts a matrix.}
\item{variable}{Predictor column name for cutoff selection; a vector of predictor names for batch_surv.}
\item{time}{Follow-up time column, containing nonnegative numeric values.}
\item{status}{Outcome column.}
\item{min_prop}{Minimum fraction per cutoff group, in (0, 0.5].}
\item{fallback}{One of "none" or "median"; fallback is used only if optimized selection fails.}
\item{status_encoding}{One of "01" (0 censored, 1 event), "12" (1 censored, 2 event), "labels", or "auto". Auto tries 0/1 first: an all-1 cohort is treated as all events. Select "12" explicitly for an all-censored 1/2-coded cohort. Labels include alive/dead, censored/event, false/true, no/yes and progression/recurrence terms. Unknown nonmissing labels cause an error.}
\item{verbose}{Print cutoff details or batch progress.}
}
\description{Find a maximally selected survival cutoff with minimum group sizes and an optional median fallback.}
\value{A list with data, cutoff, method, reason and, on success, group_column. The added factor has levels Low and High, with Low including the cutoff. Failure returns an NA cutoff and a reason.}
\details{Data-driven cutoffs are exploratory. Cox P values and confidence intervals computed after selecting a cutoff are not adjusted for cutoff selection; BH correction across contrasts does not remove that selection bias. Validate selected cutoffs in independent data. Batch fit warnings are retained in the warning column.}
