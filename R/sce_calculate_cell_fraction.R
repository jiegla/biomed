#' Calculate cell-type fractions within samples
#'
#' Count cells in each observed sample and cell-type combination and write an Excel table.
#'
#' @param sce A Seurat object with sample and cell-type metadata.
#' @param output_file Path of the Excel workbook to write.
#' @param sample_col Metadata column identifying biological samples.
#' @param celltype_col Metadata column identifying cell types. For DEG, NULL or 'all' combines all cells.
#' @details The original fixed output column names are retained even when custom metadata columns are supplied. Only observed combinations are reported; missing labels form their own groups, as in the contributed script. The Excel file is written by default.
#' @return A tibble with orig.ident, celltype_major, cell_number and fraction columns.
#' @export
sce_calculate_cell_fraction <- function(sce, output_file = "celltype_fraction_by_sample.xlsx",
                                        sample_col = "orig.ident",
                                        celltype_col = "celltype_major") {
  .biomed_require("dplyr")
  .biomed_require("openxlsx")
  if (!inherits(sce, "Seurat")) stop("sce must be a Seurat object.")
  meta <- sce[[]]
  .biomed_required_columns(meta, c(sample_col, celltype_col))
  if (identical(sample_col, celltype_col)) stop("sample_col and celltype_col must differ.")
  result <- meta |>
    dplyr::group_by(.data[[sample_col]], .data[[celltype_col]]) |>
    dplyr::summarise(cell_number = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::mutate(fraction = .data$cell_number / sum(.data$cell_number)) |>
    dplyr::ungroup()
  # Preserve the contributed table schema, including custom input column names.
  colnames(result)[1:2] <- c("orig.ident", "celltype_major")
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  openxlsx::write.xlsx(result, file = output_file, rowNames = FALSE)
  message("Cell fractions saved to: ", normalizePath(output_file))
  result
}

