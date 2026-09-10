#' Batch univariable Cox regression
#'
#' Fits one Cox model per predictor and returns every model coefficient. This
#' correctly represents multi-level factors instead of silently retaining only
#' their first contrast.
#'
#' @param var Predictor columns.
#' @param pdata A data frame.
#' @param time,status Survival time and event indicator columns.
#' @param conf_level Confidence level for hazard-ratio intervals.
#'
#' @return A data frame with hazard ratios, confidence intervals and p-values.
#' @export
cophx_batch <- function(
    var,
    pdata,
    time = "PFS",
    status = "PFS_status",
    conf_level = 0.95) {
  pdata <- as.data.frame(pdata)
  .biomed_validate_columns(var)
  .biomed_required_columns(pdata, c(time, status, var))
  event <- .biomed_as_numeric(pdata[[status]], status)
  followup <- .biomed_as_numeric(pdata[[time]], time)
  if (any(!is.na(event) & !event %in% c(0, 1))) {
    cli::cli_abort("Event status must be coded 0 (censored) or 1 (event).")
  }
  if (any(followup < 0, na.rm = TRUE)) {
    cli::cli_abort("Survival time must be non-negative.")
  }
  .biomed_probability(conf_level, "conf_level")
  if (!is.character(var) || !length(var)) {
    cli::cli_abort("{.arg var} must be a non-empty character vector.")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    cli::cli_abort("{.arg conf_level} must be one number between 0 and 1.")
  }

  rows <- lapply(var, function(marker) {
    predictor <- pdata[[marker]]
    if (is.character(predictor) || is.logical(predictor)) {
      predictor <- factor(predictor)
    }
    dat <- data.frame(
      .time = .biomed_as_numeric(pdata[[time]], time),
      .status = .biomed_as_numeric(pdata[[status]], status),
      .marker = predictor
    )
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    if (is.factor(dat$.marker)) dat$.marker <- droplevels(dat$.marker)
    events <- sum(dat$.status != 0)
    reference <- if (is.factor(dat$.marker)) levels(dat$.marker)[1L] else NA_character_
    if (nrow(dat) < 3L || events < 1L || length(unique(dat$.marker)) < 2L) {
      return(data.frame(
        Marker = marker, Term = NA_character_, Level = NA_character_,
        Reference = reference, N = nrow(dat), Events = events,
        HR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
        p_value = NA_real_, Error = "Too few informative observations",
        stringsAsFactors = FALSE
      ))
    }

    fit <- tryCatch(
      survival::coxph(survival::Surv(.time, .status) ~ .marker, data = dat),
      error = identity
    )
    if (inherits(fit, "error")) {
      return(data.frame(
        Marker = marker, Term = NA_character_, Level = NA_character_,
        Reference = reference, N = nrow(dat), Events = events,
        HR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
        p_value = NA_real_, Error = conditionMessage(fit),
        stringsAsFactors = FALSE
      ))
    }

    sm <- summary(fit)
    coef_tab <- as.matrix(sm$coefficients)
    ci_tab <- exp(stats::confint(fit, level = conf_level))
    if (is.null(dim(ci_tab))) ci_tab <- matrix(ci_tab, nrow = 1L)
    terms <- rownames(coef_tab)
    level <- if (is.factor(dat$.marker)) sub("^\\.marker", "", terms) else NA_character_
    data.frame(
      Marker = marker, Term = terms, Level = level, Reference = reference,
      N = nrow(dat), Events = events,
      HR = exp(coef_tab[, "coef"]),
      CI_lower = ci_tab[, 1L], CI_upper = ci_tab[, 2L],
      p_value = coef_tab[, "Pr(>|z|)"], Error = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$p_adjust <- stats::p.adjust(out$p_value, method = "BH")
  out
}
