#' @rdname biomed_standard_api
#' @export
plot_donut <- function(data, variables, colors = NULL, remove_na = TRUE,
                       output_dir = NULL, ...) {
  plots <- donut_pie(data, variables, color_list = colors, removeNA = remove_na, ...)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stems <- make.unique(.biomed_safe_name(names(plots)))
    for (i in seq_along(plots)) {
      for (type in c("pdf", "png")) {
        ggplot2::ggsave(file.path(output_dir, paste0(stems[i], "_donut.", type)),
                        plots[[i]], width = 5, height = 5, dpi = 300, bg = "white")
      }
    }
  }
  plots
}
