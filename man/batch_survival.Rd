\name{batch_survival}
\alias{batch_survival}
\title{Batch survival analysis with failure reports}
\usage{
batch_survival(data, variables, time = "time", status = "status",
                            best_cutoff = FALSE, min_sample = 20, min_event = 5,
                            min_censored = 5, min_prop = .1, min_event_group = 1,
                            min_censored_group = 0, fallback = c("none", "median"),
                            p_adjust = TRUE, status_encoding = "01", verbose = FALSE,
                            output_dir = NULL)
}
\arguments{
\item{data}{Data frame; plot_correlation also accepts a matrix.}
\item{variables}{Character vector of predictor column names.}
\item{time}{Follow-up time column, containing nonnegative numeric values.}
\item{status}{Outcome column.}
\item{best_cutoff}{Optimize and dichotomize each numeric predictor before fitting Cox models.}
\item{min_sample}{Minimum complete sample count, overall and per predictor.}
\item{min_event}{Minimum event count, overall and per predictor.}
\item{min_censored}{Minimum censored count, overall and per predictor.}
\item{min_prop}{Minimum fraction per cutoff group, in (0, 0.5].}
\item{min_event_group}{Minimum events per cutoff group.}
\item{min_censored_group}{Minimum censored observations per cutoff group.}
\item{fallback}{One of "none" or "median"; fallback is used only if optimized selection fails.}
\item{p_adjust}{For batch functions, logical: apply BH correction to all returned Cox Wald P values. For pairwise_survival, a p.adjust method applied to pairwise log-rank P values.}
\item{status_encoding}{One of "01" (0 censored, 1 event), "12" (1 censored, 2 event), "labels", or "auto". Auto tries 0/1 first: an all-1 cohort is treated as all events. Select "12" explicitly for an all-censored 1/2-coded cohort. Labels include alive/dead, censored/event, false/true, no/yes and progression/recurrence terms. Unknown nonmissing labels cause an error.}
\item{verbose}{Print cutoff details or batch progress.}
\item{output_dir}{Output directory; NULL writes no files. The legacy surv_fig_hr default writes a text log in Surv_Output.}
}
\description{Fit univariable Cox models for multiple predictors with eligibility thresholds, optional cutoff selection, multiple-testing adjustment and detailed failure reporting.}
\value{A list with results, failed and cutoffs tables. Results use variable, hr, conf_low, conf_high, p_value, n, events, censored, warning and optional p_adjusted/cutoff columns. A factor produces one row per nonreference contrast. Optional output is survival.xlsx.}
\details{Data-driven cutoffs are exploratory. Cox P values and confidence intervals computed after selecting a cutoff are not adjusted for cutoff selection; BH correction across contrasts does not remove that selection bias. Validate selected cutoffs in independent data. Batch fit warnings are retained in the warning column.}
