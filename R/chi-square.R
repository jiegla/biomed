#' Batch association tests for categorical variables
#'
#' Uses Pearson's chi-squared test and automatically switches to Fisher's exact
#' test whenever any expected cell count is below five.
#'
#' @param pdata A data frame.
#' @param variables Categorical columns to test.
#' @param response Categorical outcome column.
#' @param test_method One of `"auto"`, `"chi_square"`, or `"fisher"`.
#' @param correct Apply Yates' correction to 2-by-2 chi-squared tests.
#'
#' @return A data frame with test statistics, raw and adjusted p-values.
#' @export
chi_square_batch_df <- function(
    pdata,
    variables,
    response = "status",
    test_method = c("auto", "chi_square", "fisher"),
    correct = FALSE) {
  test_method <- match.arg(test_method)
  pdata <- as.data.frame(pdata)
  .biomed_validate_columns(variables)
  if (!is.character(variables) || !length(variables)) {
    cli::cli_abort("{.arg variables} must be a non-empty character vector.")
  }
  .biomed_required_columns(pdata, c(response, variables))

  rows <- lapply(variables, function(variable) {
    dat <- pdata[, c(variable, response), drop = FALSE]
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    tab <- table(droplevels(factor(dat[[variable]])),
                 droplevels(factor(dat[[response]])))
    if (any(dim(tab) < 2L)) {
      return(data.frame(
        Variable = variable, Chi_Square_Value = NA_real_, DF = NA_real_,
        P_Value = NA_real_, Test_Type = NA_character_, N = sum(tab),
        Cramers_V = NA_real_, Error = "Fewer than two observed levels",
        stringsAsFactors = FALSE
      ))
    }

    chi <- suppressWarnings(stats::chisq.test(tab, correct = correct))
    use_fisher <- test_method == "fisher" ||
      (test_method == "auto" && any(chi$expected < 5))
    test <- tryCatch(
      if (use_fisher) stats::fisher.test(tab) else chi,
      error = identity
    )
    if (inherits(test, "error")) {
      return(data.frame(
        Variable = variable, Chi_Square_Value = NA_real_, DF = NA_real_,
        P_Value = NA_real_, Test_Type = NA_character_, N = sum(tab),
        Cramers_V = NA_real_, Error = conditionMessage(test),
        stringsAsFactors = FALSE
      ))
    }
    n <- sum(tab)
    cramer <- sqrt(unname(chi$statistic) / (n * min(dim(tab) - 1L)))
    data.frame(
      Variable = variable,
      Chi_Square_Value = if (use_fisher) NA_real_ else unname(chi$statistic),
      DF = if (use_fisher) NA_real_ else unname(chi$parameter),
      P_Value = unname(test$p.value),
      Test_Type = if (use_fisher) "Fisher's Exact Test" else "Chi-squared Test",
      N = n, Cramers_V = cramer, Error = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$P_Adjust <- stats::p.adjust(out$P_Value, method = "BH")
  out[order(out$P_Value, na.last = TRUE), , drop = FALSE]
}
