\name{vennjgl}
\alias{vennjgl}
\title{Draw a Venn diagram and extract partitions}
\usage{
vennjgl(list, venn_name = "test", col = NULL,
  palette = c("npg", "aaas", "nejm", "lancet", "jama"),
  fill_col = NULL, cat.col = NULL, showFigure = TRUE,
  file_type = c("tif", "png", "pdf"), width = 6, height = 5,
  units = "in", res = 300, alpha = 0.3, cex = 1.5, cat.cex = 2,
  lwd = 2, lty = 1, main = NULL, main.cex = 2, sub = NULL,
  sub.cex = 1, rotation.degree = 0, reverse = TRUE, out_dir = ".",
  save = TRUE, write_xlsx = TRUE)
}
\arguments{
  \item{list}{A named list of two to six vectors.}
  \item{venn_name}{Output filename stem.}
  \item{col, fill_col, cat.col}{Outline, fill, and label colours.}
  \item{palette}{Fallback scientific colour palette.}
  \item{showFigure}{Draw on the current device.}
  \item{file_type}{Any of \code{tif}, \code{png}, and \code{pdf}.}
  \item{width, height, units, res}{Graphics-device settings.}
  \item{alpha, cex, cat.cex, lwd, lty, main, main.cex, sub, sub.cex,
    rotation.degree, reverse}{Arguments passed to \code{VennDiagram::venn.diagram()}.}
  \item{out_dir}{Output directory.}
  \item{save}{Save figure files.}
  \item{write_xlsx}{Save the partition table as XLSX.}
}
\value{Invisibly, a list containing the grob, partition table, and output paths.}
\description{Creates Venn diagrams and an auditable table of set partitions.}
\examples{
if (requireNamespace("VennDiagram", quietly = TRUE)) {
  vennjgl(list(A = letters[1:4], B = letters[3:6]),
    showFigure = FALSE, save = FALSE, write_xlsx = FALSE)
}
}
