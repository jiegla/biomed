#' Legacy survival cutoff selection
#'
#' Find a maximally selected survival cutoff, with optional fallback to a valid observed value nearest the median. Prefer find_survival_cutoff for a standardized interface.
#' @param pdata Data frame.
#' @param variable Predictor column name for cutoff selection; a vector of predictor names for batch_surv.
#' @param time Follow-up time column, containing nonnegative numeric values.
#' @param status Outcome column.
#' @param print_result Print the selected cutoff.
#' @param minprop Minimum fraction per cutoff group, in (0, 0.5].
#' @param fallback_cutoff One of "none" or "median". Median fallback selects the valid observed cutoff closest to the median only when optimized selection fails.
#' @param status_encoding One of "01" (0 censored, 1 event), "12" (1 censored, 2 event), "labels", or "auto". Auto tries 0/1 first: an all-1 cohort is treated as all events. Select "12" explicitly for an all-censored 1/2-coded cohort. Labels include alive/dead, censored/event, false/true, no/yes and progression/recurrence terms. Unknown nonmissing labels cause an error.
#' @return A list with pdata, best_cutoff, cutoff_method, reason and, on success, group_column. Clean complete observations are returned. Low includes values equal to the cutoff. Existing columns are preserved by choosing a unique group column name.
#' @export
best_cutoff <- function(pdata, variable, time = "time",
                            status = "status", print_result = TRUE,
                            minprop = 0.1,
                            fallback_cutoff = c("none", "median"),
                            status_encoding = "auto") {

  fallback_cutoff <- match.arg(fallback_cutoff)
  .biomed_survival_options(minprop)

  if (!is.data.frame(pdata)) {
    stop("pdata must be a data frame.")
  }

  if (!variable %in% colnames(pdata)) {
    stop(paste0("Variable not found in pdata: ", variable))
  }

  if (!time %in% colnames(pdata)) {
    stop(paste0("Time column not found in pdata: ", time))
  }

  if (!status %in% colnames(pdata)) {
    stop(paste0("Status column not found in pdata: ", status))
  }

  pdata <- as.data.frame(pdata)

  time_tmp <- .biomed_unique_name("..time_cut_biomed..", colnames(pdata))
  status_tmp <- .biomed_unique_name("..status_cut_biomed..", c(colnames(pdata), time_tmp))
  x_tmp <- .biomed_unique_name("..x_cut_biomed..", c(colnames(pdata), time_tmp, status_tmp))

  pdata[[time_tmp]] <- .biomed_as_numeric(pdata[[time]], time)
  if (any(pdata[[time_tmp]] < 0, na.rm = TRUE)) stop("Time must be nonnegative.")
  pdata[[status_tmp]] <- .biomed_status01(pdata[[status]], status_encoding)
  pdata[[x_tmp]] <- .biomed_numeric_values(pdata[[variable]])

  pdata <- pdata[
    !is.na(pdata[[time_tmp]]) &
      is.finite(pdata[[time_tmp]]) &
      !is.na(pdata[[status_tmp]]) &
      pdata[[status_tmp]] %in% c(0, 1) &
      !is.na(pdata[[x_tmp]]) &
      is.finite(pdata[[x_tmp]]),
    ,
    drop = FALSE
  ]

  clean_return <- function(df) {
    df[[time_tmp]] <- NULL
    df[[status_tmp]] <- NULL
    df[[x_tmp]] <- NULL
    df
  }

  fail_return <- function(reason) {
    list(
      pdata = clean_return(pdata),
      best_cutoff = NA_real_,
      cutoff_method = NA_character_,
      reason = reason
    )
  }

  if (nrow(pdata) == 0) {
    return(fail_return("no complete data"))
  }

  if (sum(pdata[[status_tmp]] == 1, na.rm = TRUE) == 0) {
    return(fail_return("no event"))
  }

  if (sum(pdata[[status_tmp]] == 0, na.rm = TRUE) == 0) {
    return(fail_return("no censored sample"))
  }

  if (length(unique(stats::na.omit(pdata[[x_tmp]]))) < 3 &&
      fallback_cutoff == "none") {
    return(fail_return("less than 3 unique expression values"))
  }

  valid_cutoffs <- .biomed_valid_cutoffs(pdata[[x_tmp]], minprop = minprop)

  if (length(valid_cutoffs) == 0) {
    return(fail_return("no cutoff satisfies minprop"))
  }

  cutoff_value <- NA_real_
  cutoff_method <- NA_character_
  cutoff_error <- NA_character_

  if (requireNamespace("survminer", quietly = TRUE)) {

    dat_cut <- data.frame(
      time_iobr = pdata[[time_tmp]],
      status_iobr = pdata[[status_tmp]],
      x_iobr = pdata[[x_tmp]]
    )

    iscutoff <- tryCatch(
      {
        survminer::surv_cutpoint(
          dat_cut,
          time = "time_iobr",
          event = "status_iobr",
          variables = "x_iobr",
          minprop = minprop
        )
      },
      error = function(e) {
        cutoff_error <<- conditionMessage(e)
        return(NULL)
      }
    )

    if (!is.null(iscutoff)) {
      cutoff_value <- suppressWarnings(
        as.numeric(iscutoff$cutpoint$cutpoint[1])
      )

      if (length(cutoff_value) != 1 ||
          is.na(cutoff_value) ||
          !is.finite(cutoff_value)) {
        cutoff_value <- NA_real_
      } else {
        n_low <- sum(pdata[[x_tmp]] <= cutoff_value, na.rm = TRUE)
        n_high <- sum(pdata[[x_tmp]] > cutoff_value, na.rm = TRUE)
        min_n <- ceiling(nrow(pdata) * minprop)

        if (n_low < min_n || n_high < min_n) {
          cutoff_value <- NA_real_
          cutoff_error <- "surv_cutpoint returned invalid grouping"
        } else {
          cutoff_method <- "surv_cutpoint"
        }
      }
    }
  } else {
    cutoff_error <- "package survminer is not installed"
  }

  if ((is.na(cutoff_value) || !is.finite(cutoff_value)) &&
      fallback_cutoff == "median") {

    med <- stats::median(pdata[[x_tmp]], na.rm = TRUE)
    cutoff_value <- valid_cutoffs[which.min(abs(valid_cutoffs - med))]
    cutoff_method <- "median_valid"
  }

  if (is.na(cutoff_value) || !is.finite(cutoff_value)) {
    return(fail_return(
      paste0("cutoff failed", ifelse(is.na(cutoff_error), "", paste0(": ", cutoff_error)))
    ))
  }

  variable2 <- .biomed_unique_name(paste0(variable, "_binary"), names(pdata))

  pdata[[variable2]] <- ifelse(
    pdata[[x_tmp]] <= cutoff_value,
    "Low",
    "High"
  )

  pdata[[variable2]] <- factor(
    pdata[[variable2]],
    levels = c("Low", "High")
  )

  if (print_result) {
    message(
      "Best cutoff for ",
      variable,
      ": ",
      round(cutoff_value, 3),
      " [",
      cutoff_method,
      "]"
    )
  }

  list(
    pdata = clean_return(pdata),
    best_cutoff = cutoff_value,
    cutoff_method = cutoff_method,
    group_column = variable2,
    reason = NA_character_
  )
}
