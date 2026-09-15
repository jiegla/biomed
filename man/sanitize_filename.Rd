\name{sanitize_filename}
\alias{sanitize_filename}
\title{Sanitize portable file names}
\usage{
sanitize_filename(x, max_length = 200L)
}
\arguments{
\item{x}{Character file names for sanitize_filename; an x-axis column name or position for plot_correlation.}
\item{max_length}{Maximum file-name length in characters, an integer from 1 to 200.}
}
\description{Remove invalid file-name characters without depending on attached packages. Missing values are preserved; empty and Windows reserved names are replaced. Names are not guaranteed unique; the limit counts characters, not bytes.}
\value{A character vector of safe file names, preserving missing values.}
\examples{
sanitize_filename(c("Figure 1: PFS/OS", "PD-L1 > 50%.pdf"))
}
