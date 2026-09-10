\name{batch_surv}
\alias{batch_surv}
\title{Legacy batch survival analysis}
\usage{
batch_surv(pdata, variable, time = "time", status = "status",
                           best_cutoff = FALSE,
                           min_sample = 20,
                           min_event = 5,
                           min_censored = 5,
                           minprop = 0.1,
                           min_event_group = 1,
                           min_censored_group = 0,
                           fallback_cutoff = c("none", "median"),
                           p_adjust = TRUE,
                           return_failed = FALSE,
                           verbose = TRUE, status_encoding = "auto")
}
\arguments{
\item{pdata}{Data frame.}
\item{variable}{Predictor column name for cutoff selection; a vector of predictor names for batch_surv.}
\item{time}{Follow-up time column, containing nonnegative numeric values.}
\item{status}{Outcome column.}
\item{best_cutoff}{Optimize and dichotomize each numeric predictor before fitting Cox models.}
\item{min_sample}{Minimum complete sample count, overall and per predictor.}
\item{min_event}{Minimum event count, overall and per predictor.}
\item{min_censored}{Minimum censored count, overall and per predictor.}
\item{minprop}{Minimum fraction per cutoff group, in (0, 0.5].}
\item{min_event_group}{Minimum events per cutoff group.}
\item{min_censored_group}{Minimum censored observations per cutoff group.}
\item{fallback_cutoff}{One of "none" or "median". Median fallback selects the valid observed cutoff closest to the median only when optimized selection fails.}
\item{p_adjust}{For batch functions, logical: apply BH correction to all returned Cox Wald P values. For pairwise_survival, a p.adjust method applied to pairwise log-rank P values.}
\item{return_failed}{Return a list of result, failed and (when applicable) cutoff tables; otherwise return a tibble with failed and cutoff attributes.}
\item{verbose}{Print cutoff details or batch progress.}
\item{status_encoding}{One of "01" (0 censored, 1 event), "12" (1 censored, 2 event), "labels", or "auto". Auto tries 0/1 first: an all-1 cohort is treated as all events. Select "12" explicitly for an all-censored 1/2-coded cohort. Labels include alive/dead, censored/event, false/true, no/yes and progression/recurrence terms. Unknown nonmissing labels cause an error.}
}
\description{Fit univariable Cox models with eligibility thresholds, optional cutoff selection, multiple-testing adjustment and detailed failure reporting. Prefer batch_survival for standardized names.}
\value{A tibble with ID, cox_variable, term, HR, lower95, upper95, P, N, Event, Censored, warning and optional FDR/cutoff summaries. Multiple factor contrasts are retained. Failures are available as attributes or list elements with return_failed = TRUE.}
\details{Data-driven cutoffs are exploratory. Cox P values and confidence intervals computed after selecting a cutoff are not adjusted for cutoff selection; BH correction across contrasts does not remove that selection bias. Validate selected cutoffs in independent data. Batch fit warnings are retained in the warning column.}
