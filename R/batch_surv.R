#' Legacy batch survival analysis
#'
#' Fit univariable Cox models with eligibility thresholds, optional cutoff selection, multiple-testing adjustment and detailed failure reporting. Prefer batch_survival for standardized names.
#' @param pdata Data frame.
#' @param variable Predictor column name for cutoff selection; a vector of predictor names for batch_surv.
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
batch_surv <- function(pdata, variable, time = "time", status = "status",
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
      fail_df <- .biomed_add_failure(
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

  row_col <- .biomed_unique_name("..row_id_biomed..", colnames(pdata))
  time_col <- .biomed_unique_name("..time_biomed..", c(colnames(pdata), row_col))
  status_col <- .biomed_unique_name("..status_biomed..", c(colnames(pdata), row_col, time_col))

  pdata[[row_col]] <- seq_len(nrow(pdata))
  pdata[[time_col]] <- .biomed_as_numeric(pdata[[time]], time)
  if (any(pdata[[time_col]] < 0, na.rm = TRUE)) stop("Time must be nonnegative.")
  pdata[[status_col]] <- .biomed_status01(pdata[[status]], status_encoding)

  pdata <- pdata[
    !is.na(pdata[[time_col]]) &
      is.finite(pdata[[time_col]]) &
      !is.na(pdata[[status_col]]) &
      pdata[[status_col]] %in% c(0, 1),
    ,
    drop = FALSE
  ]

  if (nrow(pdata) < min_sample) {
    fail_df <- .biomed_add_failure(
      fail_df,
      ID = "ALL",
      stage = "input",
      reason = "sample size after removing missing time/status is less than min_sample"
    )
    return(make_output(tibble::tibble(), fail_df, cutoff_df))
  }

  if (sum(pdata[[status_col]] == 1, na.rm = TRUE) < min_event) {
    fail_df <- .biomed_add_failure(
      fail_df,
      ID = "ALL",
      stage = "input",
      reason = "total event number is less than min_event"
    )
    return(make_output(tibble::tibble(), fail_df, cutoff_df))
  }

  if (sum(pdata[[status_col]] == 0, na.rm = TRUE) < min_censored) {
    fail_df <- .biomed_add_failure(
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
      binary_names[[v]] <- .biomed_unique_name(paste0(v, "_binary"), used_names)
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
      pdata[[var]] <- .biomed_numeric_values(pdata[[var]])

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
      } else if (length(.biomed_valid_cutoffs(tmp[[var]], minprop = minprop)) == 0) {
        fail_reason <- "no cutoff satisfies minprop"
      }

      if (!is.null(fail_reason)) {
        fail_df <- .biomed_add_failure(
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
          best_cutoff(
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
        fail_df <- .biomed_add_failure(
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
        fail_df <- .biomed_add_failure(
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
        fail_df <- .biomed_add_failure(
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
        fail_df <- .biomed_add_failure(
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
      fail_df <- .biomed_add_failure(
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
      fail_df <- .biomed_add_failure(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = "valid sample size is less than min_sample"
      )
      next
    }

    if (sum(tmp_cox$status_iobr == 1, na.rm = TRUE) < min_event) {
      fail_df <- .biomed_add_failure(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = "event number is less than min_event"
      )
      next
    }

    if (sum(tmp_cox$status_iobr == 0, na.rm = TRUE) < min_censored) {
      fail_df <- .biomed_add_failure(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = "censored number is less than min_censored"
      )
      next
    }

    if (length(unique(tmp_cox$x)) < 2) {
      fail_df <- .biomed_add_failure(
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
      fail_df <- .biomed_add_failure(
        fail_df,
        ID = original_id,
        stage = "cox",
        reason = paste0("coxph failed: ", err_msg)
      )
      next
    }

    result1 <- tryCatch(
      {
        .biomed_extract_cox(fit)
      },
      error = function(e) {
        err_msg <<- conditionMessage(e)
        return(NULL)
      }
    )

    if (is.null(result1) || nrow(result1) == 0) {
      fail_df <- .biomed_add_failure(
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

      expr <- .biomed_numeric_values(pdata[[original_var]])
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
