#' Sanitize a file name without requiring attached packages
#' @param x Character vector of file names (not full paths).
#' @param max_length Maximum number of characters, between 1 and 200.
#' @return A character vector; missing inputs stay missing. Empty and reserved
#'   names become safe names. This function does not guarantee uniqueness.
#' @export
sanitize_filename <- function(x, max_length = 200L) {
  if (!is.character(x)) stop("x must be a character vector.")
  if (length(max_length) != 1L || !is.numeric(max_length) ||
      !is.finite(max_length) || max_length < 1 || max_length > 200 ||
      max_length != floor(max_length)) stop("max_length must be an integer from 1 to 200.")
  # Preserve the original space-to-underscore behavior.
  out <- gsub("[/\\\\]", "_", x)
  out <- gsub('[:*?"<>|]', "", out)
  out <- gsub("[[:space:]]+", "_", out)
  out <- gsub("[^[:alnum:]_.-]", "", out)
  out <- substr(out, 1L, max_length)
  out <- sub("[. ]+$", "", out)
  empty <- !is.na(out) & !nzchar(out)
  out[empty] <- substr("output", 1L, max_length)
  reserved <- !is.na(out) & grepl("^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])($|\\.)", out, ignore.case = TRUE)
  out[reserved] <- substr(paste0("_", out[reserved]), 1L, max_length)
  out
}
