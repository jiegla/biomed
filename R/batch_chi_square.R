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
