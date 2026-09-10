.num_jgl <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

.status01_jgl <- function(x, encoding = "auto") {
  encoding <- match.arg(encoding, c("auto", "01", "12", "labels"))
  if (is.logical(x) && encoding %in% c("01", "auto", "labels")) return(as.numeric(x))
  raw <- trimws(tolower(as.character(x)))
  numeric_x <- suppressWarnings(as.numeric(raw))
  present <- !is.na(x)
  if (encoding %in% c("01", "12")) {
    allowed <- if (encoding == "01") c(0, 1) else c(1, 2)
    if (any(present & (is.na(numeric_x) | !numeric_x %in% allowed))) {
      stop("Status must use the selected ", encoding, " encoding.")
    }
    return(if (encoding == "12") numeric_x - 1 else numeric_x)
  }
  if (encoding == "auto" && all(!present | !is.na(numeric_x))) {
    observed <- unique(numeric_x[present])
    if (all(observed %in% c(0, 1))) return(numeric_x)
    if (all(observed %in% c(1, 2))) return(numeric_x - 1)
    stop("Unsupported numeric status; use 0/1 or 1/2.")
  }
  events <- c("1", "true", "t", "yes", "y", "event", "dead", "death",
              "deceased", "progression", "progressed", "pd", "relapse",
              "recurred", "recurrence")
  censored <- c("0", "false", "f", "no", "n", "censored", "alive", "living",
                "no event", "noevent", "non-event", "non event", "no progression")
  out <- rep(NA_real_, length(x))
  out[raw %in% events] <- 1
  out[raw %in% censored] <- 0
  if (any(present & is.na(out))) stop("Unrecognized status labels.")
  out
}

.biomed_survival_options <- function(minprop, counts = numeric()) {
  if (length(minprop) != 1L || !is.numeric(minprop) || !is.finite(minprop) ||
      minprop <= 0 || minprop > 0.5) stop("minprop must be in (0, 0.5].")
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts))) stop("Sample/event thresholds must be nonnegative integers.")
}

.valid_cutoffs_jgl <- function(x, minprop = 0.1) {
  x <- x[!is.na(x)]
  n <- length(x)

  if (n == 0) return(numeric(0))

  ux <- sort(unique(x))

  if (length(ux) < 2) return(numeric(0))

  cuts <- ux[-length(ux)]
  min_n <- ceiling(n * minprop)

  cuts[vapply(
    cuts,
    function(cut) {
      n_low <- sum(x <= cut, na.rm = TRUE)
      n_high <- sum(x > cut, na.rm = TRUE)
      n_low >= min_n && n_high >= min_n
    },
    logical(1)
  )]
}

.extract_cox_jgl <- function(fit) {
  s <- summary(fit)

  co <- as.data.frame(s$coefficients, check.names = FALSE)
  ci <- as.data.frame(s$conf.int, check.names = FALSE)

  low_col <- grep("lower", colnames(ci), value = TRUE)[1]
  high_col <- grep("upper", colnames(ci), value = TRUE)[1]
  p_col <- grep("^Pr", colnames(co), value = TRUE)[1]

  data.frame(
    term = rownames(co),
    HR = ci[["exp(coef)"]],
    lower95 = ci[[low_col]],
    upper95 = ci[[high_col]],
    P = co[[p_col]],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.add_fail_jgl <- function(fail_df, ID, stage, reason) {
  rbind(
    fail_df,
    data.frame(
      ID = as.character(ID),
      stage = as.character(stage),
      reason = as.character(reason),
      stringsAsFactors = FALSE
    )
  )
}

.tmp_name_jgl <- function(base, nms) {
  nm <- base
  while (nm %in% nms) {
    nm <- paste0(".", nm)
  }
  nm
}

#' Legacy survival cutoff selection
#'
#' Find a maximally selected survival cutoff, with optional fallback to a valid observed value nearest the median. Prefer find_survival_cutoff for a standardized interface.
#' @param pdata Data frame.
#' @param variable Predictor column name for cutoff selection; a vector of predictor names for batch_surv_jgl.
#' @param time Follow-up time column, containing nonnegative numeric values.
#' @param status Outcome column.
#' @param print_result Print the selected cutoff.
#' @param minprop Minimum fraction per cutoff group, in (0, 0.5].
#' @param fallback_cutoff One of "none" or "median". Median fallback selects the valid observed cutoff closest to the median only when optimized selection fails.
#' @param status_encoding One of "01" (0 censored, 1 event), "12" (1 censored, 2 event), "labels", or "auto". Auto tries 0/1 first: an all-1 cohort is treated as all events. Select "12" explicitly for an all-censored 1/2-coded cohort. Labels include alive/dead, censored/event, false/true, no/yes and progression/recurrence terms. Unknown nonmissing labels cause an error.
#' @return A list with pdata, best_cutoff, cutoff_method, reason and, on success, group_column. Clean complete observations are returned. Low includes values equal to the cutoff. Existing columns are preserved by choosing a unique group column name.
#' @export
best_cutoff_jgl <- function(pdata, variable, time = "time",
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

  time_tmp <- .tmp_name_jgl("..time_cut_jgl..", colnames(pdata))
  status_tmp <- .tmp_name_jgl("..status_cut_jgl..", c(colnames(pdata), time_tmp))
  x_tmp <- .tmp_name_jgl("..x_cut_jgl..", c(colnames(pdata), time_tmp, status_tmp))

  pdata[[time_tmp]] <- .biomed_as_numeric(pdata[[time]], time)
  if (any(pdata[[time_tmp]] < 0, na.rm = TRUE)) stop("Time must be nonnegative.")
  pdata[[status_tmp]] <- .status01_jgl(pdata[[status]], status_encoding)
  pdata[[x_tmp]] <- .num_jgl(pdata[[variable]])

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

  valid_cutoffs <- .valid_cutoffs_jgl(pdata[[x_tmp]], minprop = minprop)

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

  variable2 <- .tmp_name_jgl(paste0(variable, "_binary"), names(pdata))

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

#' Legacy batch survival analysis
#'
#' Fit univariable Cox models with eligibility thresholds, optional cutoff selection, multiple-testing adjustment and detailed failure reporting. Prefer batch_survival for standardized names.
#' @param pdata Data frame.
#' @param variable Predictor column name for cutoff selection; a vector of predictor names for batch_surv_jgl.
#' @param time Follow-up time column, containing nonnegative numeric values.
#' @param status Outcome column.
#' @param best_cutoff Optimize and dichotomize each numeric predictor before fitting Cox models.
#' @param min_sample Minimum complete sample count, overall and per predictor.
#' @param min_event Minimum event count, overall and per predictor.
#' @param min_censored Minimum censored count, overall and per predictor.
#' @param minprop Minimum fraction per cutoff group, in (0, 0.5].
#' @param min_event_group Minimum events per cutoff group.
#' @param min_censored_group Minimum censored observations per cutoff group.
#' @param fallback_cutoff One of "none" or "median". Median fallback selects the valid observed cutoff closest to the median only when optimized selection fails.
#' @param p_adjust For batch functions, logical: apply BH correction to all returned Cox Wald P values. For pairwise_survival, a p.adjust method applied to pairwise log-rank P values.
#' @param return_failed Return a list of result, failed and (when applicable) cutoff tables; otherwise return a tibble with failed and cutoff attributes.
#' @param verbose Print cutoff details or batch progress.
#' @param status_encoding One of "01" (0 censored, 1 event), "12" (1 censored, 2 event), "labels", or "auto". Auto tries 0/1 first: an all-1 cohort is treated as all events. Select "12" explicitly for an all-censored 1/2-coded cohort. Labels include alive/dead, censored/event, false/true, no/yes and progression/recurrence terms. Unknown nonmissing labels cause an error.
#' @return A tibble with ID, cox_variable, term, HR, lower95, upper95, P, N, Event, Censored, warning and optional FDR/cutoff summaries. Multiple factor contrasts are retained. Failures are available as attributes or list elements with return_failed = TRUE.
#' @export
batch_surv_jgl <- function(pdata, variable, time = "time", status = "status",
                           best_cutoff = FALSE,
                           min_sample = 20,
                           min_event = 5,
                           min_censored = 5,
                           minprop = 0.1,
                           min_event_group = 1,
                           min_censored_group = 0,
                           fallback_cutoff = c("none", "median"),
                           p_adjust = TRUE,
                           return_failed = FALSE,
                           verbose = TRUE, status_encoding = "auto") {

  fallback_cutoff <- match.arg(fallback_cutoff)
  .biomed_survival_options(minprop)

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }

  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.")
  }

  pdata <- as.data.frame(pdata)
  variable <- unique(as.character(variable))
  .biomed_validate_columns(variable)
  .biomed_survival_options(minprop, c(min_sample, min_event, min_censored,
                                     min_event_group, min_censored_group))

  fail_df <- data.frame(
    ID = character(),
    stage = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )

  cutoff_df <- NULL

  make_output <- function(result, fail_df, cutoff_df = NULL) {
    result <- tibble::as_tibble(result)
    attr(result, "failed") <- tibble::as_tibble(fail_df)

    if (!is.null(cutoff_df)) {
      attr(result, "cutoff") <- tibble::as_tibble(cutoff_df)
    }

    if (return_failed) {
      out <- list(
        result = result,
        failed = tibble::as_tibble(fail_df)
      )

      if (!is.null(cutoff_df)) {
        out$cutoff <- tibble::as_tibble(cutoff_df)
      }

      return(out)
    }

    result
  }

  if (!time %in% colnames(pdata)) {
    stop(paste0("Time column not found: ", time))
  }

  if (!status %in% colnames(pdata)) {
    stop(paste0("Status column not found: ", status))
  }

  missing_var <- setdiff(variable, colnames(pdata))

  if (length(missing_var) > 0) {
    for (v in missing_var) {
      fail_df <- .add_fail_jgl(
        fail_df,
        ID = v,
        stage = "input",
        reason = "variable not found in pdata"
      )
    }
  }

  variable <- intersect(variable, colnames(pdata))

  if (length(variable) == 0) {
    warning("No valid variables found in pdata.")
    return(make_output(tibble::tibble(), fail_df, cutoff_df))
  }

  row_col <- .tmp_name_jgl("..row_id_jgl..", colnames(pdata))
  time_col <- .tmp_name_jgl("..time_jgl..", c(colnames(pdata), row_col))
  status_col <- .tmp_name_jgl("..status_jgl..", c(colnames(pdata), row_col, time_col))

  pdata[[row_col]] <- seq_len(nrow(pdata))
  pdata[[time_col]] <- .biomed_as_numeric(pdata[[time]], time)
  if (any(pdata[[time_col]] < 0, na.rm = TRUE)) stop("Time must be nonnegative.")
  pdata[[status_col]] <- .status01_jgl(pdata[[status]], status_encoding)

  pdata <- pdata[
    !is.na(pdata[[time_col]]) &
      is.finite(pdata[[time_col]]) &
      !is.na(pdata[[status_col]]) &
      pdata[[status_col]] %in% c(0, 1),
    ,
    drop = FALSE
  ]

  if (nrow(pdata) < min_sample) {
    fail_df <- .add_fail_jgl(
      fail_df,
      ID = "ALL",
      stage = "input",
      reason = "sample size after removing missing time/status is less than min_sample"
    )
    return(make_output(tibble::tibble(), fail_df, cutoff_df))
  }

  if (sum(pdata[[status_col]] == 1, na.rm = TRUE) < min_event) {
    fail_df <- .add_fail_jgl(
      fail_df,
      ID = "ALL",
      stage = "input",
      reason = "total event number is less than min_event"
    )
    return(make_output(tibble::tibble(), fail_df, cutoff_df))
  }

  if (sum(pdata[[status_col]] == 0, na.rm = TRUE) < min_censored) {
    fail_df <- .add_fail_jgl(
      fail_df,
      ID = "ALL",
      stage = "input",
      reason = "total censored number is less than min_censored"
    )
    return(make_output(tibble::tibble(), fail_df, cutoff_df))
  }

  if (best_cutoff) {

    cutoff_df <- data.frame(
      ID = variable,
      cutoff = NA_real_,
      cutoff_method = NA_character_,
      stringsAsFactors = FALSE
    )

    valid_variable <- character()
    binary_names <- stats::setNames(character(length(variable)), variable)
    used_names <- names(pdata)
    for (v in variable) {
      binary_names[[v]] <- .tmp_name_jgl(paste0(v, "_binary"), used_names)
      used_names <- c(used_names, binary_names[[v]])
    }

    if (verbose) {
      pb <- utils::txtProgressBar(
        min = 0,
        max = length(variable),
        style = 3
      )
    }

    for (i in seq_along(variable)) {

      var <- variable[i]
      pdata[[var]] <- .num_jgl(pdata[[var]])

      tmp <- pdata[
        !is.na(pdata[[var]]) &
          is.finite(pdata[[var]]),
        ,
        drop = FALSE
      ]

      fail_reason <- NULL

      if (nrow(tmp) < min_sample) {
        fail_reason <- "valid sample size is less than min_sample"
      } else if (sum(tmp[[status_col]] == 1, na.rm = TRUE) < min_event) {
        fail_reason <- "event number is less than min_event"
      } else if (sum(tmp[[status_col]] == 0, na.rm = TRUE) < min_censored) {
        fail_reason <- "censored number is less than min_censored"
      } else if (length(unique(stats::na.omit(tmp[[var]]))) < 3 &&
                 fallback_cutoff == "none") {
        fail_reason <- "less than 3 unique expression values"
      } else if (length(.valid_cutoffs_jgl(tmp[[var]], minprop = minprop)) == 0) {
        fail_reason <- "no cutoff satisfies minprop"
      }

      if (!is.null(fail_reason)) {
        fail_df <- .add_fail_jgl(
          fail_df,
          ID = var,
          stage = "best_cutoff",
          reason = fail_reason
        )
        if (verbose) utils::setTxtProgressBar(pb, i)
        next
      }

      cutoff_res <- tryCatch(
        {
          best_cutoff_jgl(
            pdata = tmp,
            time = time_col,
            status = status_col,
            variable = var,
            print_result = FALSE,
            minprop = minprop,
            fallback_cutoff = fallback_cutoff
          )
        },
        error = function(e) {
          list(
            pdata = tmp,
            best_cutoff = NA_real_,
            cutoff_method = NA_character_,
            reason = conditionMessage(e)
          )
        }
      )

      if (is.na(cutoff_res$best_cutoff)) {
        fail_df <- .add_fail_jgl(
          fail_df,
          ID = var,
          stage = "best_cutoff",
          reason = cutoff_res$reason
        )
        if (verbose) utils::setTxtProgressBar(pb, i)
        next
      }

      binary_col <- binary_names[[var]]
      pdata[[binary_col]] <- NA_real_

      idx <- match(cutoff_res$pdata[[row_col]], pdata[[row_col]])

      pdata[[binary_col]][idx] <- ifelse(
        cutoff_res$pdata[[cutoff_res$group_column]] == "High",
        1,
        0
      )

      tmp_b <- pdata[
        !is.na(pdata[[binary_col]]),
        ,
        drop = FALSE
      ]

      high_n <- sum(tmp_b[[binary_col]] == 1, na.rm = TRUE)
      low_n <- sum(tmp_b[[binary_col]] == 0, na.rm = TRUE)
      min_n <- ceiling(nrow(tmp_b) * minprop)

      high_event <- sum(
        tmp_b[[binary_col]] == 1 & tmp_b[[status_col]] == 1,
        na.rm = TRUE
      )

      low_event <- sum(
        tmp_b[[binary_col]] == 0 & tmp_b[[status_col]] == 1,
        na.rm = TRUE
      )

      high_censored <- sum(
        tmp_b[[binary_col]] == 1 & tmp_b[[status_col]] == 0,
        na.rm = TRUE
      )

      low_censored <- sum(
        tmp_b[[binary_col]] == 0 & tmp_b[[status_col]] == 0,
        na.rm = TRUE
      )

      if (high_n < min_n || low_n < min_n) {
        fail_df <- .add_fail_jgl(
          fail_df,
          ID = var,
          stage = "best_cutoff",
          reason = "high/low sample size does not satisfy minprop after cutoff"
        )
        if (verbose) utils::setTxtProgressBar(pb, i)
        next
      }

      if (min_event_group > 0 &&
          (high_event < min_event_group || low_event < min_event_group)) {
        fail_df <- .add_fail_jgl(
          fail_df,
          ID = var,
          stage = "best_cutoff",
          reason = "event number in high or low group is too small"
        )
        if (verbose) utils::setTxtProgressBar(pb, i)
        next
      }

      if (min_censored_group > 0 &&
          (high_censored < min_censored_group || low_censored < min_censored_group)) {
        fail_df <- .add_fail_jgl(
          fail_df,
          ID = var,
          stage = "best_cutoff",
          reason = "censored number in high or low group is too small"
        )
        if (verbose) utils::setTxtProgressBar(pb, i)
        next
      }

      hit <- match(var, cutoff_df$ID)
      cutoff_df$cutoff[hit] <- cutoff_res$best_cutoff
      cutoff_df$cutoff_method[hit] <- cutoff_res$cutoff_method

      valid_variable <- c(valid_variable, var)

      if (verbose) utils::setTxtProgressBar(pb, i)
    }

    if (verbose) close(pb)

    variable_for_cox <- unname(binary_names[valid_variable])

    if (length(variable_for_cox) == 0) {
      warning("No variables passed best cutoff filtering.")
      return(make_output(tibble::tibble(), fail_df, cutoff_df))
    }

  } else {
    variable_for_cox <- variable
  }

  result_list <- list()

  for (var in variable_for_cox) {

    if (!var %in% colnames(pdata)) {
      fail_df <- .add_fail_jgl(
        fail_df,
        ID = var,
        stage = "cox",
        reason = "variable not found before Cox"
      )
      next
    }

    tmp_cox <- data.frame(
      time_iobr = pdata[[time_col]],
      status_iobr = pdata[[status_col]],
      x = pdata[[var]],
      stringsAsFactors = FALSE
    )

    tmp_cox <- tmp_cox[
      !is.na(tmp_cox$time_iobr) &
        is.finite(tmp_cox$time_iobr) &
        !is.na(tmp_cox$status_iobr) &
        tmp_cox$status_iobr %in% c(0, 1) &
        !is.na(tmp_cox$x),
      ,
      drop = FALSE
    ]

    if (is.numeric(tmp_cox$x)) tmp_cox <- tmp_cox[is.finite(tmp_cox$x), , drop = FALSE]

    if (is.character(tmp_cox$x)) {
      tmp_cox$x <- factor(tmp_cox$x)
    }

    if (is.factor(tmp_cox$x)) {
      tmp_cox$x <- droplevels(tmp_cox$x)
    }

    original_id <- if (best_cutoff) names(binary_names)[match(var, binary_names)] else var

    if (nrow(tmp_cox) < min_sample) {
      fail_df <- .add_fail_jgl(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = "valid sample size is less than min_sample"
      )
      next
    }

    if (sum(tmp_cox$status_iobr == 1, na.rm = TRUE) < min_event) {
      fail_df <- .add_fail_jgl(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = "event number is less than min_event"
      )
      next
    }

    if (sum(tmp_cox$status_iobr == 0, na.rm = TRUE) < min_censored) {
      fail_df <- .add_fail_jgl(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = "censored number is less than min_censored"
      )
      next
    }

    if (length(unique(tmp_cox$x)) < 2) {
      fail_df <- .add_fail_jgl(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = "only one group/value available"
      )
      next
    }

    warn_msg <- character()
    err_msg <- NULL

    fit <- withCallingHandlers(
      tryCatch(
        {
          survival::coxph(
            survival::Surv(time_iobr, status_iobr) ~ x,
            data = tmp_cox
          )
        },
        error = function(e) {
          err_msg <<- conditionMessage(e)
          return(NULL)
        }
      ),
      warning = function(w) {
        warn_msg <<- c(warn_msg, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    if (is.null(fit)) {
      fail_df <- .add_fail_jgl(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = paste0("coxph failed: ", err_msg)
      )
      next
    }

    result1 <- tryCatch(
      {
        .extract_cox_jgl(fit)
      },
      error = function(e) {
        err_msg <<- conditionMessage(e)
        return(NULL)
      }
    )

    if (is.null(result1) || nrow(result1) == 0) {
      fail_df <- .add_fail_jgl(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = paste0("extract Cox result failed: ", err_msg)
      )
      next
    }

    result1$ID <- original_id
    result1$cox_variable <- var
    result1$N <- nrow(tmp_cox)
    result1$Event <- sum(tmp_cox$status_iobr == 1, na.rm = TRUE)
    result1$Censored <- sum(tmp_cox$status_iobr == 0, na.rm = TRUE)
    result1$warning <- ifelse(
      length(warn_msg) == 0,
      NA_character_,
      paste(unique(warn_msg), collapse = " | ")
    )

    result_list[[length(result_list) + 1]] <- result1
  }

  if (length(result_list) == 0) {
    warning("No variables were available for Cox analysis.")
    return(make_output(tibble::tibble(), fail_df, cutoff_df))
  }

  result <- do.call(rbind, result_list)

  result <- result[
    ,
    c(
      "ID", "cox_variable", "term",
      "HR", "lower95", "upper95", "P",
      "N", "Event", "Censored", "warning"
    ),
    drop = FALSE
  ]

  if (p_adjust) {
    result$FDR <- stats::p.adjust(result$P, method = "BH")
  }

  if (best_cutoff) {

    result$high <- NA_real_
    result$low <- NA_real_
    result$high_m <- NA_real_
    result$low_m <- NA_real_
    result$cutoff <- NA_real_
    result$cutoff_method <- NA_character_

    for (i in seq_len(nrow(result))) {

      original_var <- result$ID[i]
      binary_col <- binary_names[[original_var]]

      if (!binary_col %in% colnames(pdata)) next
      if (!original_var %in% colnames(pdata)) next

      expr <- .num_jgl(pdata[[original_var]])
      bin <- pdata[[binary_col]]

      result$high[i] <- sum(bin == 1, na.rm = TRUE)
      result$low[i] <- sum(bin == 0, na.rm = TRUE)

      result$high_m[i] <- stats::median(expr[bin == 1], na.rm = TRUE)
      result$low_m[i] <- stats::median(expr[bin == 0], na.rm = TRUE)

      hit <- match(original_var, cutoff_df$ID)

      if (!is.na(hit)) {
        result$cutoff[i] <- cutoff_df$cutoff[hit]
        result$cutoff_method[i] <- cutoff_df$cutoff_method[hit]
      }
    }
  }

  result <- result[order(result$P, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL

  make_output(result, fail_df, cutoff_df)
}
