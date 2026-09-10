\name{plot_survival}
\alias{plot_survival}
\title{Kaplan-Meier curves with directed hazard ratios}
\usage{
plot_survival(data, group, time = "time", status = "status",
                          max_time = NULL, reference = NULL, status_encoding = "01",
                          show_hr = TRUE, colors = NULL, x_label = "Months",
                          risk_table = TRUE, sort_by_median = FALSE, show_plot = FALSE,
                          output_dir = NULL, filename = group, formats = c("pdf", "png"),
                          width = 8, height = 7)
}
\arguments{
\item{data}{Data frame; plot_correlation also accepts a matrix.}
\item{group}{Group column; optional for correlation.}
\item{time}{Follow-up time column, containing nonnegative numeric values.}
\item{status}{Outcome column.}
\item{max_time}{Positive upper display limit; does not truncate observations used for fitting. NULL in plot_survival uses the maximum observed time.}
\item{reference}{Optional observed group placed first; HR is always group versus reference as named in the result.}
\item{status_encoding}{One of "01" (0 censored, 1 event), "12" (1 censored, 2 event), "labels", or "auto". Auto tries 0/1 first: an all-1 cohort is treated as all events. Select "12" explicitly for an all-censored 1/2-coded cohort. Labels include alive/dead, censored/event, false/true, no/yes and progression/recurrence terms. Unknown nonmissing labels cause an error.}
\item{show_hr}{Show directed two-group HR and 95 percent confidence interval.}
\item{colors}{Optional group colors, optionally named.}
\item{x_label}{Time-axis unit label.}
\item{risk_table}{Include the number-at-risk table.}
\item{sort_by_median}{Order display groups by decreasing median survival, with unreached medians last. An explicit reference is kept first. Result comparison labels identify HR direction.}
\item{show_plot}{Print the plot.}
\item{output_dir}{Output directory; NULL writes no files. The legacy surv_fig_hr default writes a text log in Surv_Output.}
\item{filename}{Base output filename.}
\item{formats}{One or both of pdf and png; survival export includes the risk table and a statistics text file.}
\item{width}{Export width in inches.}
\item{height}{Export height in inches.}
}
\description{Plot survival curves, risk tables, median survival, global log-rank P values and explicitly directed two-group HRs. Pairwise statistics are retained for more than two groups.}
\value{A survminer ggsurvplot object. Attribute statistics contains medians, pairwise comparisons, global logrank_p, n, events and excluded counts. Attribute survival_fit stores the fitted curve. NR labels indicate a median or bound was not reached.}
