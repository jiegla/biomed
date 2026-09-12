#' Score single-cell gene signatures
#'
#' Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.
#' @param sce A Seurat object; for SCENIC, also accepts the path to an RDS object.
#' @param gene_sets Nonempty list of gene vectors. Supply unique names for stable score column names.
#' @param method Scoring backend: UCell, GSVA, or AUCell.
#' @param assay Assay to use. NULL uses the default assay; SCENIC defaults to RNA.
#' @param prefix Prefix for score metadata columns; NULL or an empty string omits it.
#' @param return_type Return the updated Seurat object, a metadata data frame, or a cells-by-signatures matrix.
#' @param layer Expression layer to use; data denotes normalized expression. Split Seurat v5 layers are joined on a local copy.
#' @return An updated Seurat object, data frame, or numeric matrix as selected by return_type. Rows in score outputs follow the original cell order.
#' @export
sce_run_scoring <- function(sce, gene_sets, method = c("UCell", "GSVA", "AUCell"),
                            assay = NULL, prefix = "Pathway",
                            return_type = c("sce", "meta", "matrix"), layer = "data") {
  method <- match.arg(method)
  return_type <- match.arg(return_type)
  if (!is.list(gene_sets) || !length(gene_sets) ||
      !all(vapply(gene_sets, function(x) is.character(x) && length(x) > 0L &&
                  !anyNA(x) && all(nzchar(x)), logical(1)))) {
    stop("`gene_sets` must be a nonempty list of character gene vectors.")
  }
  if (is.null(names(gene_sets))) names(gene_sets) <- paste0("Set", seq_along(gene_sets))
  if (anyNA(names(gene_sets)) || any(!nzchar(names(gene_sets))) || anyDuplicated(names(gene_sets))) {
    stop("Gene-set names must be unique and nonempty.")
  }
  if (!is.null(prefix) && (length(prefix) != 1L || is.na(prefix))) stop("Invalid prefix.")
  expr <- .biomed_sce_expression(sce, assay, layer)
  .biomed_require(method)
  scores <- switch(method,
    UCell = as.matrix(UCell::ScoreSignatures_UCell(expr, features = gene_sets,
      maxRank = min(1500L, nrow(expr)))),
    GSVA = {
      params <- GSVA::gsvaParam(as.matrix(expr), gene_sets, kcdf = "Gaussian")
      t(GSVA::gsva(params, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
    },
    AUCell = {
      ranks <- AUCell::AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
      t(AUCell::getAUC(AUCell::AUCell_calcAUC(gene_sets, ranks)))
    }
  )
  if (!all(colnames(sce) %in% rownames(scores))) stop("Scoring output is missing cells.")
  scores <- scores[colnames(sce), , drop = FALSE]
  if (!is.null(prefix) && nzchar(prefix)) colnames(scores) <- paste0(prefix, "_", colnames(scores))
  if (anyDuplicated(colnames(scores))) stop("Scoring returned duplicate column names.")
  if (return_type == "matrix") return(scores)
  meta <- as.data.frame(scores)
  if (return_type == "meta") return(meta)
  if (!is.null(assay)) SeuratObject::DefaultAssay(sce) <- assay
  SeuratObject::AddMetaData(sce, metadata = meta)
}
