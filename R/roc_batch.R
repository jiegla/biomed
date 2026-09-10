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
