#' Assign cell-cycle scores and phases
#'
#' Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.
#' @param seurat_obj A Seurat object.
#' @param species Human or mouse gene-symbol conventions.
#' @param assay Assay to use. NULL uses the default assay; SCENIC defaults to RNA.
#' @param ... Additional arguments to Seurat::CellCycleScoring, such as nbin or ctrl.
#' @return The Seurat object with S.Score, G2M.Score and Phase metadata. Cell identities are preserved.
#' @export
sce_cell_cycle <- function(seurat_obj, species = c("human", "mouse"), assay = NULL, ...) {
  species <- match.arg(species)
  .biomed_require("Seurat")
  if (!inherits(seurat_obj, "Seurat")) stop("`seurat_obj` must be a Seurat object.")
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(seurat_obj)
  if (!assay %in% names(seurat_obj@assays)) stop("Assay not found: ", assay)
  SeuratObject::DefaultAssay(seurat_obj) <- assay
  genes <- Seurat::cc.genes.updated.2019
  if (species == "mouse") {
    title_gene <- function(x) paste0(toupper(substr(x, 1L, 1L)), tolower(substring(x, 2L)))
    genes <- lapply(genes, title_gene)
  }
  s <- intersect(genes$s.genes, rownames(seurat_obj[[assay]]))
  g2m <- intersect(genes$g2m.genes, rownames(seurat_obj[[assay]]))
  if (!length(s) || !length(g2m)) stop("Selected assay must contain both S and G2M genes.")
  Seurat::CellCycleScoring(seurat_obj, s.features = s, g2m.features = g2m,
                          assay = assay, set.ident = FALSE, ...)
}
