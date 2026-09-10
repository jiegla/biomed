#' Batch one-way ANOVA
#'
#' Runs one-way ANOVA for multiple numeric features, reports group means and
#' controls the false-discovery rate with the Benjamini-Hochberg method.
#'
#' @param data A data frame.
#' @param target Grouping column.
#' @param feature Numeric feature columns. If `NULL`, all numeric columns other
#'   than `target` are used.
#' @param feature_manipulation Remove incomplete, non-numeric, infinite and
#'   constant features before testing.
#'
#' @return A data frame ordered by raw p-value.
#' @export
batch_ANOVA <- function(
    data,
    target = "group",
    feature = NULL,
    feature_manipulation = FALSE) {
  if (is.null(data)) {
    return(NULL)
  }
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }
  .biomed_required_columns(data, target)

  if (is.null(feature)) {
    feature <- setdiff(names(data)[vapply(data, is.numeric, logical(1))], target)
  }
  .biomed_validate_columns(feature)
  .biomed_required_columns(data, feature)
  numeric_feature <- feature[vapply(data[feature], is.numeric, logical(1))]
  if (isTRUE(feature_manipulation)) {
    numeric_feature <- numeric_feature[vapply(data[numeric_feature], function(x) {
      length(x) > 1L && all(is.finite(x)) && stats::sd(x) > 0
    }, logical(1))]
  }
  skipped <- setdiff(feature, numeric_feature)
  if (length(skipped)) {
    cli::cli_warn("Non-numeric features were skipped: {skipped}.")
  }
  if (!length(numeric_feature)) {
    cli::cli_abort("No numeric features are available for ANOVA.")
  }

  group <- droplevels(factor(data[[target]]))
  if (nlevels(group) < 2L) {
    cli::cli_abort("ANOVA requires at least two groups.")
  }

  rows <- lapply(numeric_feature, function(feat) {
    dat <- data.frame(.value = data[[feat]], .group = group)
    dat <- dat[stats::complete.cases(dat) & is.finite(dat$.value), , drop = FALSE]
    dat$.group <- droplevels(dat$.group)
    if (nrow(dat) < 3L || nlevels(dat$.group) < 2L) {
      return(data.frame(
        sig_names = feat, statistic = NA_real_, p.value = NA_real_,
        n = nrow(dat), error = "Too few complete observations",
        stringsAsFactors = FALSE
      ))
    }

    fit <- tryCatch(stats::aov(.value ~ .group, data = dat), error = identity)
    if (inherits(fit, "error")) {
      return(data.frame(
        sig_names = feat, statistic = NA_real_, p.value = NA_real_,
        n = nrow(dat), error = conditionMessage(fit), stringsAsFactors = FALSE
      ))
    }
    tab <- summary(fit)[[1L]]
    means <- tapply(dat$.value, dat$.group, mean, na.rm = TRUE)
    mean_names <- paste0("mean_", make.names(names(means), unique = TRUE))
    mean_row <- as.data.frame(as.list(stats::setNames(as.numeric(means), mean_names)))
    base <- data.frame(
      sig_names = feat,
      statistic = unname(tab[["F value"]][1L]),
      p.value = unname(tab[["Pr(>F)"]][1L]),
      n = nrow(dat), error = NA_character_, stringsAsFactors = FALSE
    )
    if (length(means) == 2L) {
      base$mean_diff <- unname(means[1L] - means[2L])
    } else {
      base$mean_diff <- NA_real_
    }
    data.frame(base, mean_row, check.names = FALSE)
  })

  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    for (nm in setdiff(all_names, names(x))) x[[nm]] <- NA
    x[all_names]
  })
  out <- do.call(rbind, rows)
  out$p.adj <- stats::p.adjust(out$p.value, method = "BH")
  out$log10pvalue <- -log10(out$p.value)
  out$stars <- .biomed_significance(out$p.value)
  out[order(out$p.value, na.last = TRUE), , drop = FALSE]
}
