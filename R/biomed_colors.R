#' Built-in biomedical research colors
#'
#' The 14 supplied colors are returned in their original order, including
#' the CC alpha channel (80 percent opacity). Larger palettes are interpolated.
#' @param n Number of colors. Defaults to all 14; zero returns character(0).
#' @return A character vector of hexadecimal RGBA colors.
#' @export
#' @examples
#' biomed_colors()
#' biomed_colors(3)
biomed_colors <- function(n = 14L) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 0 || n != floor(n)) stop("n must be a nonnegative integer.")
  colors <- c(
    "#E64B35CC", "#4DBBD5CC", "#00A087CC", "#3C5488CC", "#F39B7FCC",
    "#8491B4CC", "#91D1C2CC", "#DC0000CC", "#7E6148CC", "#5050FFCC",
    "#CE3D32CC", "#749B58CC", "#F0E685CC", "#466983CC"
  )
  if (n <= length(colors)) return(colors[seq_len(n)])
  grDevices::colorRampPalette(colors, alpha = TRUE)(n)
}
