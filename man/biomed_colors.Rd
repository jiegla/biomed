\name{biomed_colors}
\alias{biomed_colors}
\title{Built-in biomedical research colors}
\usage{biomed_colors(n = 14L)}
\arguments{\item{n}{Nonnegative integer number of colors; defaults to 14.}}
\description{Return the supplied 14 colors in their original order, preserving the CC alpha channel (80 percent opacity). More than 14 colors are interpolated with alpha preserved.}
\value{A hexadecimal RGBA character vector. Zero returns character(0).}
\examples{
biomed_colors()
biomed_colors(3)
}
