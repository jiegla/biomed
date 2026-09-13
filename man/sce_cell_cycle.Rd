\name{sce_cell_cycle}
\alias{sce_cell_cycle}
\title{Assign cell-cycle scores and phases}
\usage{
sce_cell_cycle(seurat_obj, species = c("human", "mouse"), assay = NULL, ...)
}
\arguments{
\item{seurat_obj}{A Seurat object.}
\item{species}{Human or mouse gene-symbol conventions.}
\item{assay}{Assay to use. NULL uses the default assay; SCENIC defaults to RNA.}
\item{...}{Additional arguments to Seurat::CellCycleScoring, such as nbin or ctrl.}
}
\description{Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.}
\value{The Seurat object with S.Score, G2M.Score and Phase metadata. Cell identities are preserved.}
