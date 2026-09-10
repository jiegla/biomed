\name{surv_fig_hr}
\alias{surv_fig_hr}
\alias{grid.draw.ggsurvplot}
\title{Legacy Kaplan-Meier plot interface}
\usage{
surv_fig_hr(group_var, df, max_time = 60, time = "OS",
                        status = "Survival_status", hr = "yes", cols = NULL,
                        x.label = "Months", output_dir = "Surv_Output", ...)

\method{grid.draw}{ggsurvplot}(x, recording = TRUE)
}
\arguments{
\item{x}{A survminer ggsurvplot object.}
\item{recording}{Compatibility argument for grid.draw.}
\item{group_var}{Grouping column name.}
\item{df}{Data frame or matrix for correlation; data frame for survival.}
\item{max_time}{Positive upper display limit; does not truncate observations used for fitting. NULL in plot_survival uses the maximum observed time.}
\item{time}{Follow-up time column, containing nonnegative numeric values.}
\item{status}{Outcome column.}
\item{hr}{Legacy switch: "yes" shows the directed two-group HR.}
\item{cols}{Legacy group palette; recycled to group count.}
\item{x.label}{Legacy time-axis unit label.}
\item{output_dir}{Output directory; NULL writes no files. The legacy surv_fig_hr default writes a text log in Surv_Output.}
\item{...}{Additional plot styling and export arguments passed to the underlying plotting function. Do not repeat arguments already mapped by the wrapper.}
}
\description{Compatibility interface for plot_survival, retaining the original text log and median-based display ordering. HRs are always labeled with their actual comparison direction.}
\value{A survminer ggsurvplot object with statistics and survival_fit attributes, as in plot_survival.}

\details{The registered grid.draw method uses standard print dispatch with newpage = FALSE, allowing ggplot2::ggsave to save the whole survival plot including its risk table. No global helper definition is needed.}
