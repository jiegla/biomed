#' Standard clinical analysis interfaces
#'
#' These interfaces use snake_case names and consistent data, variables,
#' output_dir and colors arguments. Original functions remain available.
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

#' @rdname biomed_standard_api
#' @export
batch_anova <- function(data, variables = NULL, group = "group",
                        preprocess = FALSE, output_dir = NULL) {
  result <- batch_ANOVA(data, target = group, feature = variables,
                        feature_manipulation = preprocess)
  if (is.null(result)) return(NULL)
  result <- .biomed_standard_table(result, c(
    sig_names = "variable", p.value = "p_value", p.adj = "p_adjust",
    log10pvalue = "neg_log10_p", stars = "significance"
  ))
  .biomed_save_table(result, output_dir, "anova.xlsx")
}

#' @rdname biomed_standard_api
#' @export
batch_chi_square <- function(data, variables, response = "status",
                            test_method = "auto", correct = FALSE,
                            output_dir = NULL) {
  result <- chi_square_batch_df(data, variables, response, test_method, correct)
  result <- .biomed_standard_table(result, c(
    Variable = "variable", Chi_Square_Value = "statistic", DF = "df",
    P_Value = "p_value", Test_Type = "method", N = "n",
    Cramers_V = "cramers_v", Error = "error", P_Adjust = "p_adjust"
  ))
  .biomed_save_table(result, output_dir, "chi_square.xlsx")
}

#' @rdname biomed_standard_api
#' @export
batch_cox <- function(data, variables, time = "PFS", status = "PFS_status",
                      conf_level = 0.95, output_dir = NULL) {
  result <- cophx_batch(variables, data, time, status, conf_level)
  result <- .biomed_standard_table(result, c(
    Marker = "variable", Term = "term", Level = "level", Reference = "reference",
    N = "n", Events = "events", HR = "hr", CI_lower = "conf_low",
    CI_upper = "conf_high", Error = "error"
  ))
  .biomed_save_table(result, output_dir, "cox.xlsx")
}

#' @rdname biomed_standard_api
#' @export
batch_roc <- function(data, variables, response, method = "best",
                      positive_class = NULL, direction = "auto",
                      conf_level = 0.95, output_dir = NULL) {
  result <- roc_batch(variables, data, response, method, positive_class,
                      direction, conf_level)
  names(result) <- tolower(names(result))
  result <- .biomed_standard_table(result, c(
    markers = "variable", auc_ci_lower = "conf_low", auc_ci_upper = "conf_high"
  ))
  .biomed_save_table(result, output_dir, "roc.xlsx")
}

#' @rdname biomed_standard_api
#' @export
plot_roc <- function(data, variables, response, output_dir = NULL,
                     positive_class = NULL, direction = "auto", ...) {
  plot_roc_batch(data, response, variables, fig.path = output_dir,
                 positive_class = positive_class, direction = direction, ...)
}

#' @rdname biomed_standard_api
#' @export
plot_donut <- function(data, variables, colors = NULL, remove_na = TRUE,
                       output_dir = NULL, ...) {
  plots <- donutPie(data, variables, color_list = colors, removeNA = remove_na, ...)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stems <- make.unique(.biomed_safe_name(names(plots)))
    for (i in seq_along(plots)) {
      for (type in c("pdf", "png")) {
        ggplot2::ggsave(file.path(output_dir, paste0(stems[i], "_donut.", type)),
                        plots[[i]], width = 5, height = 5, dpi = 300, bg = "white")
      }
    }
  }
  plots
}

#' @rdname biomed_standard_api
#' @export
plot_stacked_bar <- function(data, variables, group, colors = NULL,
                             remove_na = TRUE, output_dir = NULL, ...) {
  result <- draw_stack_barplot(
    data, group, variables, color = colors, remove_na = remove_na,
    out_dir = if (is.null(output_dir)) "StackBar" else output_dir,
    save = !is.null(output_dir), return_plot = TRUE, ...
  )
  names(result)[names(result) == "plot"] <- "plots"
  result
}

#' @rdname biomed_standard_api
#' @export
plot_venn <- function(sets, filename = "venn", colors = NULL,
                      show_plot = TRUE, output_dir = NULL, ...) {
  result <- vennjgl(
    sets, venn_name = filename, col = colors, fill_col = colors, cat.col = colors,
    showFigure = show_plot, out_dir = if (is.null(output_dir)) "." else output_dir,
    save = !is.null(output_dir), write_xlsx = !is.null(output_dir), ...
  )
  invisible(result)
}
