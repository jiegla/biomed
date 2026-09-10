#' Group continuous variables using reproducible cutoffs
#'
#' Creates binary or three-level factor columns from continuous variables.
#' Cutoffs may be based on the median, maximally selected survival statistics,
#' or ROC criteria. A summary is stored in the `cutoff_table` attribute.
#'
#' @param pdata A data frame.
#' @param variables Continuous columns to group.
#' @param num_group Either `2` or `3`.
#' @param method Cutoff method. Three groups currently require `"median"` and
#'   use the 33rd and 66th percentiles.
#' @param time,status Survival time and event columns.
#' @param response Binary response column for ROC methods.
#' @param fixed_value Target sensitivity or specificity in `(0, 1)`.
#' @param positive_class Event class for ROC analysis. If `NULL`, the second
#'   factor level is used.
#' @param direction Passed to [pROC::roc()].
#' @param add_method_suffix Include the method in generated column names.
#' @param verbose Print cutoff summaries.
#'
#' @return `pdata` with new factor columns and a `cutoff_table` attribute.
#' @export
makegroup <- function(
    pdata,
    variables,
    num_group = 2,
    method = c(
      "median", "survival_best", "roc_youden", "roc_closest",
      "fixed_sensitivity", "fixed_specificity"
    ),
    time = "time",
    status = "status",
    response = "group",
    fixed_value = 0.8,
    positive_class = NULL,
    direction = "auto",
    add_method_suffix = FALSE,
    verbose = FALSE) {
  method <- match.arg(method)
  pdata <- as.data.frame(pdata)
  .biomed_validate_columns(variables)

  if (!is.character(variables) || length(variables) == 0L) {
    cli::cli_abort("{.arg variables} must be a non-empty character vector.")
  }
  .biomed_required_columns(pdata, variables)

  if (length(num_group) != 1L || is.na(num_group) || !num_group %in% c(2, 3)) {
    cli::cli_abort("{.arg num_group} must be 2 or 3.")
  }
  if (num_group == 3 && method != "median") {
    cli::cli_abort("Three groups are supported only with {.val median}.")
  }

  roc_methods <- c(
    "roc_youden", "roc_closest", "fixed_sensitivity", "fixed_specificity"
  )

  if (method == "survival_best") {
    .biomed_required_columns(pdata, c(time, status))
    .biomed_require("survminer", "It is used for maximally selected cutoffs.")
  }
  if (method %in% roc_methods) {
    .biomed_required_columns(pdata, response)
    .biomed_require("pROC", "It is used for ROC-based cutoffs.")
  }
  if (method %in% c("fixed_sensitivity", "fixed_specificity") &&
      (!is.numeric(fixed_value) || length(fixed_value) != 1L ||
       is.na(fixed_value) || fixed_value <= 0 || fixed_value >= 1)) {
    cli::cli_abort("{.arg fixed_value} must be one number between 0 and 1.")
  }

  survival_cutoff <- function(feature) {
    dat <- data.frame(
      .time = .biomed_as_numeric(pdata[[time]], time),
      .status = .biomed_as_numeric(pdata[[status]], status),
      .feature = pdata[[feature]]
    )
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    if (nrow(dat) < 5L || length(unique(dat$.feature)) < 2L) {
      cli::cli_abort("Too few informative observations for {.val {feature}}.")
    }
    fit <- survminer::surv_cutpoint(
      dat, time = ".time", event = ".status", variables = ".feature",
      progressbar = FALSE
    )
    as.numeric(fit$cutpoint[1L, "cutpoint"])
  }

  roc_cutoff <- function(feature) {
    dat <- data.frame(
      .response = pdata[[response]],
      .feature = pdata[[feature]]
    )
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    y <- as.character(dat$.response)
    levels_y <- if (is.factor(dat$.response)) levels(droplevels(dat$.response)) else sort(unique(y))
    if (length(levels_y) != 2L) {
      cli::cli_abort("{.arg response} must contain exactly two classes.")
    }
    positive <- if (is.null(positive_class)) levels_y[2L] else as.character(positive_class)
    if (length(positive) != 1L || is.na(positive)) {
      cli::cli_abort("{.arg positive_class} must identify exactly one class.")
    }
    if (!positive %in% levels_y) {
      cli::cli_abort("{.arg positive_class} was not found in {.arg response}.")
    }
    negative <- setdiff(levels_y, positive)
    roc_obj <- pROC::roc(
      factor(y, levels = c(negative, positive)), dat$.feature,
      levels = c(negative, positive), direction = direction, quiet = TRUE
    )

    if (method == "roc_youden") {
      coord <- pROC::coords(
        roc_obj, "best", best.method = "youden",
        ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE
      )
    } else if (method == "roc_closest") {
      coord <- pROC::coords(
        roc_obj, "best", best.method = "closest.topleft",
        ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE
      )
    } else {
      target <- if (method == "fixed_sensitivity") "sensitivity" else "specificity"
      coord <- pROC::coords(
        roc_obj, "all",
        ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE
      )
      coord <- as.data.frame(coord)
      coord <- coord[is.finite(coord$threshold) & coord[[target]] >= fixed_value, , drop = FALSE]
      other <- if (target == "sensitivity") "specificity" else "sensitivity"
      coord <- coord[order(-coord[[other]], -coord[[target]], coord$threshold), , drop = FALSE]
      coord <- utils::head(coord, 1L)
    }
    coord <- as.data.frame(coord)
    thresholds <- as.numeric(coord$threshold)
    thresholds <- thresholds[is.finite(thresholds)]
    if (length(thresholds) == 0L) {
      cli::cli_abort("No finite ROC threshold was found for {.val {feature}}.")
    }
    list(
      cutoff = thresholds[1L],
      auc = as.numeric(pROC::auc(roc_obj)),
      sensitivity = as.numeric(pROC::coords(roc_obj, thresholds[1L], ret = "sensitivity", transpose = FALSE)[[1L]]),
      specificity = as.numeric(pROC::coords(roc_obj, thresholds[1L], ret = "specificity", transpose = FALSE)[[1L]])
    )
  }

  summaries <- list()
  for (variable in variables) {
    pdata[[variable]] <- .biomed_as_numeric(pdata[[variable]], variable)

    if (num_group == 2L) {
      info <- if (method == "median") {
        list(
          cutoff = stats::median(pdata[[variable]], na.rm = TRUE),
          auc = NA_real_, sensitivity = NA_real_, specificity = NA_real_
        )
      } else if (method == "survival_best") {
        list(
          cutoff = survival_cutoff(variable),
          auc = NA_real_, sensitivity = NA_real_, specificity = NA_real_
        )
      } else {
        roc_cutoff(variable)
      }

      if (!is.finite(info$cutoff)) {
        cli::cli_warn("No finite cutoff for {.val {variable}}; it was skipped.")
        next
      }
      group_col <- if (add_method_suffix) {
        paste0(variable, "_", method)
      } else {
        paste0(variable, "_binary")
      }
      pdata[[group_col]] <- factor(
        ifelse(is.na(pdata[[variable]]), NA_character_,
               ifelse(pdata[[variable]] <= info$cutoff, "Low", "High")),
        levels = c("Low", "High")
      )
      counts <- table(pdata[[group_col]])
      summaries[[variable]] <- data.frame(
        feature = variable, group_col = group_col, method = method,
        cutoff = info$cutoff, low_cutoff = NA_real_, high_cutoff = NA_real_,
        auc = info$auc, sensitivity = info$sensitivity,
        specificity = info$specificity,
        n_low = unname(counts["Low"]), n_middle = NA_integer_,
        n_high = unname(counts["High"]), stringsAsFactors = FALSE
      )
    } else {
      cutoffs <- stats::quantile(
        pdata[[variable]], c(0.33, 0.66), na.rm = TRUE, names = FALSE
      )
      if (any(!is.finite(cutoffs))) {
        cli::cli_warn("No finite tertile cutoffs for {.val {variable}}; it was skipped.")
        next
      }
      group_col <- paste0(variable, "_ternary")
      pdata[[group_col]] <- factor(
        ifelse(pdata[[variable]] <= cutoffs[1L], "Low",
               ifelse(pdata[[variable]] <= cutoffs[2L], "Middle", "High")),
        levels = c("Low", "Middle", "High")
      )
      counts <- table(pdata[[group_col]])
      summaries[[variable]] <- data.frame(
        feature = variable, group_col = group_col, method = "tertile",
        cutoff = NA_real_, low_cutoff = cutoffs[1L], high_cutoff = cutoffs[2L],
        auc = NA_real_, sensitivity = NA_real_, specificity = NA_real_,
        n_low = unname(counts["Low"]), n_middle = unname(counts["Middle"]),
        n_high = unname(counts["High"]), stringsAsFactors = FALSE
      )
    }

    if (verbose) {
      print(summaries[[variable]])
    }
  }

  cutoff_table <- if (length(summaries)) {
    do.call(rbind, unname(summaries))
  } else {
    data.frame()
  }
  rownames(cutoff_table) <- NULL
  attr(pdata, "cutoff_table") <- cutoff_table
  pdata
}
