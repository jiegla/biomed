.biomed_roc_data <- function(data, response, predictor, positive_class) {
  .biomed_required_columns(data, c(response, predictor))
  dat <- data.frame(
    .response = data[[response]],
    .predictor = .biomed_as_numeric(data[[predictor]], predictor)
  )
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  y <- as.character(dat$.response)
  observed <- if (is.factor(dat$.response)) levels(droplevels(dat$.response)) else sort(unique(y))
  if (length(observed) != 2L) {
    cli::cli_abort("{.arg response} must contain exactly two observed classes.")
  }
  positive <- if (is.null(positive_class)) observed[2L] else as.character(positive_class)
  if (length(positive) != 1L || is.na(positive)) {
    cli::cli_abort("{.arg positive_class} must identify exactly one observed class.")
  }
  if (!positive %in% observed) {
    cli::cli_abort("Positive class {.val {positive}} was not found in {.arg response}.")
  }
  negative <- setdiff(observed, positive)
  dat$.response <- factor(y, levels = c(negative, positive))
  list(data = dat, negative = negative, positive = positive)
}

.biomed_threshold_metrics <- function(roc_obj, dat, threshold) {
  positive <- levels(dat$.response)[2L]
  predicted <- if (roc_obj$direction == "<") {
    dat$.predictor >= threshold
  } else {
    dat$.predictor <= threshold
  }
  actual <- dat$.response == positive
  tp <- sum(predicted & actual)
  tn <- sum(!predicted & !actual)
  fp <- sum(predicted & !actual)
  fn <- sum(!predicted & actual)
  safe_div <- function(x, y) if (y == 0) NA_real_ else x / y
  c(
    Specificity = safe_div(tn, tn + fp),
    Sensitivity = safe_div(tp, tp + fn),
    Accuracy = safe_div(tp + tn, length(actual)),
    Precision = safe_div(tp, tp + fp),
    Recall = safe_div(tp, tp + fn)
  )
}

#' Batch ROC analysis
#'
#' Computes AUC, confidence intervals, a selected threshold and classification
#' metrics for multiple continuous predictors.
#'
#' @param vars Predictor columns.
#' @param pdata A data frame.
#' @param outcome Binary outcome column.
#' @param method Use the Youden-optimal (`"best"`) or median threshold.
#' @param positive_class Outcome value treated as the event.
#' @param direction Passed to [pROC::roc()].
#' @param conf_level Confidence level for AUC intervals.
#'
#' @return A data frame with one row per predictor.
#' @export
roc_batch <- function(
    vars,
    pdata,
    outcome,
    method = c("best", "median"),
    positive_class = "1",
    direction = "auto",
    conf_level = 0.95) {
  .biomed_require("pROC")
  method <- match.arg(method)
  pdata <- as.data.frame(pdata)
  .biomed_required_columns(pdata, c(outcome, vars))
  .biomed_validate_columns(vars)
  .biomed_probability(conf_level, "conf_level")

  rows <- lapply(vars, function(var) {
    tryCatch({
      prepared <- .biomed_roc_data(pdata, outcome, var, positive_class)
      dat <- prepared$data
      roc_obj <- pROC::roc(
        dat$.response, dat$.predictor,
        levels = c(prepared$negative, prepared$positive),
        direction = direction, quiet = TRUE
      )
      threshold <- if (method == "best") {
        coords <- pROC::coords(
          roc_obj, "best", best.method = "youden", ret = "threshold",
          transpose = FALSE
        )
        values <- as.numeric(as.data.frame(coords)$threshold)
        values <- values[is.finite(values)]
        if (!length(values)) NA_real_ else values[1L]
      } else {
        stats::median(dat$.predictor)
      }
      metrics <- if (is.finite(threshold)) {
        .biomed_threshold_metrics(roc_obj, dat, threshold)
      } else {
        stats::setNames(rep(NA_real_, 5L),
                        c("Specificity", "Sensitivity", "Accuracy", "Precision", "Recall"))
      }
      ci <- as.numeric(pROC::ci.auc(roc_obj, conf.level = conf_level))
      data.frame(
        Markers = var, AUC = as.numeric(pROC::auc(roc_obj)),
        AUC_CI_lower = ci[1L], AUC_CI_upper = ci[3L], Threshold = threshold,
        Specificity = unname(metrics["Specificity"]),
        Sensitivity = unname(metrics["Sensitivity"]),
        Accuracy = unname(metrics["Accuracy"]),
        Precision = unname(metrics["Precision"]),
        Recall = unname(metrics["Recall"]), N = nrow(dat),
        Positive = sum(dat$.response == prepared$positive),
        Negative = sum(dat$.response == prepared$negative),
        Positive_Class = prepared$positive, Negative_Class = prepared$negative,
        Direction = roc_obj$direction, Error = NA_character_,
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }, error = function(e) {
      data.frame(
        Markers = var, AUC = NA_real_, AUC_CI_lower = NA_real_,
        AUC_CI_upper = NA_real_, Threshold = NA_real_, Specificity = NA_real_,
        Sensitivity = NA_real_, Accuracy = NA_real_, Precision = NA_real_,
        Recall = NA_real_, N = NA_integer_, Positive = NA_integer_,
        Negative = NA_integer_, Positive_Class = NA_character_,
        Negative_Class = NA_character_, Direction = NA_character_,
        Error = conditionMessage(e), stringsAsFactors = FALSE
      )
    })
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

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
