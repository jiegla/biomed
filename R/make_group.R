#' Standard clinical analysis interfaces
#'
#' These interfaces use snake_case names and consistent data, variables,
#' output_dir and colors arguments. Legacy interfaces use the neutral names listed in the README migration table.
#' @name biomed_standard_api
#' @param data A data frame.
#' @param variables Columns to analyse.
#' @param group Grouping column.
#' @param response Binary or categorical response column.
#' @param time,status Survival time and 0/1 event columns.
#' @param method Grouping or ROC threshold method.
#' @param n_groups Number of groups (2 or 3).
#' @param preprocess Filter non-numeric, incomplete, infinite or constant features.
#' @param positive_class Event class. NULL uses the second observed factor level
#'   or second sorted label, recorded in ROC results.
#' @param direction ROC direction; use a prespecified direction for validation.
#' @param conf_level Confidence level.
#' @param output_dir Optional output directory; NULL writes no files.
#' @param test_method Association test selection.
#' @param correct Apply Yates' correction.
#' @param colors Colour vector, or for donuts a list by variable.
#' @param remove_na Remove missing values.
#' @param sets Named list of sets.
#' @param filename Venn output filename stem.
#' @param show_plot Draw the Venn diagram on the active device.
#' @param ... Additional options of the corresponding original function.
#' @return Analysis functions return data frames with snake_case columns;
#'   grouping returns augmented data with a cutoff_table attribute; plotting
#'   functions return named plot lists or a list with plots and statistics.
#' @details Tables are optionally exported as XLSX; figures use PDF/PNG by
#'   default. Legacy arguments supplied through ... retain their original names.
#'   Old entry points preserve their result column names. Optimized cutoffs and
#'   ROC estimates computed on the same data are exploratory, not validation.
NULL

#' @rdname biomed_standard_api
#' @export
make_group <- function(data, variables, n_groups = 2, method = "median", ...) {
  makegroup(data, variables, num_group = n_groups, method = method, ...)
}
