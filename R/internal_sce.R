.biomed_sce_expression <- function(sce, assay = NULL, layer = "data") {
  .biomed_require("SeuratObject")
  if (!inherits(sce, "Seurat")) stop("`sce` must be a Seurat object.")
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(sce)
  if (!assay %in% names(sce@assays)) stop("Assay not found: ", assay)
  if (length(layer) != 1L || is.na(layer) || !nzchar(layer)) stop("Supply one layer name.")
  a <- sce[[assay]]
  if (inherits(a, "Assay5")) {
    available <- SeuratObject::Layers(a)
    matched <- available[available == layer | startsWith(available, paste0(layer, "."))]
    if (!length(matched)) stop("Layer not found: ", layer)
    if (length(matched) > 1L) {
      a <- SeuratObject::JoinLayers(a, layers = layer, new = layer)
      matched <- layer
    }
    x <- SeuratObject::LayerData(a, layer = matched)
  } else {
    x <- SeuratObject::GetAssayData(sce, assay = assay, slot = layer)
  }
  cells <- colnames(sce)
  if (!all(cells %in% colnames(x))) stop("Selected layer does not contain every Seurat cell.")
  x[, cells, drop = FALSE]
}
