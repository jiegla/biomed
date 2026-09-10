#' @rdname biomed_standard_api
#' @export
plot_roc <- function(data, variables, response, output_dir = NULL,
                     positive_class = NULL, direction = "auto", ...) {
  plot_roc_batch(data, response, variables, fig.path = output_dir,
                 positive_class = positive_class, direction = direction, ...)
}
