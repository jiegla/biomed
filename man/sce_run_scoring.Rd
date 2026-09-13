\name{sce_run_scoring}
\alias{sce_run_scoring}
\title{Score single-cell gene signatures}
\usage{
sce_run_scoring(sce, gene_sets, method = c("UCell", "GSVA", "AUCell"),
                            assay = NULL, prefix = "Pathway",
                            return_type = c("sce", "meta", "matrix"), layer = "data")
}
\arguments{
\item{sce}{A Seurat object; for SCENIC, also accepts the path to an RDS object.}
\item{gene_sets}{Nonempty list of gene vectors. Supply unique names for stable score column names.}
\item{method}{Scoring backend: UCell, GSVA, or AUCell.}
\item{assay}{Assay to use. NULL uses the default assay; SCENIC defaults to RNA.}
\item{prefix}{Prefix for score metadata columns; NULL or an empty string omits it.}
\item{return_type}{Return the updated Seurat object, a metadata data frame, or a cells-by-signatures matrix.}
\item{layer}{Expression layer to use; data denotes normalized expression. Split Seurat v5 layers are joined on a local copy.}
}
\description{Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.}
\value{An updated Seurat object, data frame, or numeric matrix as selected by return_type. Rows in score outputs follow the original cell order.}
