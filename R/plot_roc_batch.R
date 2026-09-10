#' Plot ROC curves in batch
#'
#' @param data A data frame.
#' @param response Binary outcome column.
#' @param variables Continuous predictor columns.
#' @param fig.path Optional output directory. If `NULL`, files are not saved.
#' @param main Optional title prefix.
#' @param alpha Line opacity.
#' @param smooth Smooth ROC curves before plotting.
#' @param positive_class Outcome value treated as the event.
#' @param direction Passed to [pROC::roc()].
#' @param file_type Any of `"png"` and `"pdf"`.
#' @param width,height Figure dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return An invisible named list of ggplot objects.
#' @export
plot_roc_batch <- function(
    data,
    response,
    variables,
    fig.path = NULL,
    main = NULL,
    alpha = 1,
    smooth = FALSE,
    positive_class = NULL,
    direction = "auto",
    file_type = c("png", "pdf"),
    width = 5,
    height = 5,
    dpi = 300) {
  .biomed_require("pROC")
  data <- as.data.frame(data)
  .biomed_required_columns(data, c(response, variables))
  .biomed_validate_columns(variables)
  file_type <- match.arg(file_type, c("png", "pdf"), several.ok = TRUE)
  if (!is.null(fig.path)) dir.create(fig.path, recursive = TRUE, showWarnings = FALSE)

  plots <- lapply(variables, function(var) {
    prepared <- .biomed_roc_data(data, response, var, positive_class)
    roc_obj <- pROC::roc(
      prepared$data$.response, prepared$data$.predictor,
      levels = c(prepared$negative, prepared$positive),
      direction = direction, quiet = TRUE
    )
    if (smooth) roc_obj <- pROC::smooth(roc_obj)
    ci <- as.numeric(pROC::ci.auc(roc_obj))
    title <- paste0(
      if (is.null(main)) paste0(var, ": ") else paste0(main, " - ", var, ": "),
      "AUC = ", sprintf("%.3f", as.numeric(pROC::auc(roc_obj))),
      " (95% CI ", sprintf("%.3f", ci[1L]), "-", sprintf("%.3f", ci[3L]), ")"
    )
    curve_data <- data.frame(
      false_positive_rate = 1 - roc_obj$specificities,
      true_positive_rate = roc_obj$sensitivities
    )
    p <- ggplot2::ggplot(curve_data, ggplot2::aes(
      x = .data[["false_positive_rate"]], y = .data[["true_positive_rate"]]
    )) +
      ggplot2::geom_path(alpha = alpha) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
      ggplot2::labs(title = title, x = "False positive rate", y = "True positive rate") +
      ggplot2::coord_equal() +
      ggplot2::theme_bw()

    if (!is.null(fig.path)) {
      stem <- file.path(fig.path, paste0(.biomed_safe_name(var), "_ROC"))
      for (type in file_type) {
        args <- list(
          filename = paste0(stem, ".", type), plot = p,
          width = width, height = height, bg = "white"
        )
        if (type == "png") args$dpi <- dpi
        do.call(ggplot2::ggsave, args)
      }
    }
    p
  })
  names(plots) <- variables
  invisible(plots)
}
