\name{pairwise_survival}
\alias{pairwise_survival}
\title{Directed pairwise survival comparisons}
\usage{
pairwise_survival(data, group, time = "time", status = "status",
                              reference = NULL, status_encoding = "01", p_adjust = "BH")
}
\arguments{
\item{data}{Data frame; plot_correlation also accepts a matrix.}
\item{group}{Group column; optional for correlation.}
\item{time}{Follow-up time column, containing nonnegative numeric values.}
\item{status}{Outcome column.}
\item{reference}{Optional observed group placed first; HR is always group versus reference as named in the result.}
\item{status_encoding}{One of "01" (0 censored, 1 event), "12" (1 censored, 2 event), "labels", or "auto". Auto tries 0/1 first: an all-1 cohort is treated as all events. Select "12" explicitly for an all-censored 1/2-coded cohort. Labels include alive/dead, censored/event, false/true, no/yes and progression/recurrence terms. Unknown nonmissing labels cause an error.}
\item{p_adjust}{For batch functions, logical: apply BH correction to all returned Cox Wald P values. For pairwise_survival, a p.adjust method applied to pairwise log-rank P values.}
}
\description{Compare each pair of groups with explicitly named HR direction and separate Cox Wald and log-rank P values.}
\value{A data frame with group, reference, n, events, hr, conf_low, conf_high, p_value (Cox Wald), logrank_p, warning and p_adjusted (adjusted log-rank). A single group yields a zero-row table.}
