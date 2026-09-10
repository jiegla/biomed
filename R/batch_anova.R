#' @rdname biomed_standard_api
#' @export
batch_anova <- function(data, variables = NULL, group = "group",
                        preprocess = FALSE, output_dir = NULL) {
  result <- anova_batch(data, target = group, feature = variables,
                        feature_manipulation = preprocess)
  if (is.null(result)) return(NULL)
  result <- .biomed_standard_table(result, c(
    sig_names = "variable", p.value = "p_value", p.adj = "p_adjust",
    log10pvalue = "neg_log10_p", stars = "significance"
  ))
  .biomed_save_table(result, output_dir, "anova.xlsx")
}
