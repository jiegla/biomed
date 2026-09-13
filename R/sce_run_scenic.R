#' Run SCENIC and add regulon activity metadata
#'
#' Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.
#' @param sce A Seurat object; for SCENIC, also accepts the path to an RDS object.
#' @param db_dir Existing directory containing species-compatible cisTarget ranking databases.
#' @param dbs Database filenames within db_dir; NULL discovers feather files. Supply only compatible ranking databases.
#' @param species Human or mouse gene-symbol conventions.
#' @param assay Assay to use. NULL uses the default assay; SCENIC defaults to RNA.
#' @param layer Expression layer to use; data denotes normalized expression. Split Seurat v5 layers are joined on a local copy.
#' @param output_dir Directory for generated tables, figures, and analysis files.
#' @param n_cores Number of SCENIC worker cores.
#' @param min_counts_per_gene Minimum summed expression for geneFiltering; NULL uses 3 times 1 percent of the cell count.
#' @param min_samples Minimum cells expressing a gene; NULL uses 1 percent of the cell count.
#' @param top_n_regulons Number of most variable regulons added as metadata.
#' @param top_n_plot Maximum genes or regulons selected for feature plots.
#' @param celltype Metadata column containing cell-type annotations.
#' @param reduction Existing dimensional reduction used for plots.
#' @param seed Random seed for reproducibility; the caller RNG state is restored.
#' @param save_pdf Save PDF figures.
#' @param save_png Save PNG figures.
#' @param overwrite Permit reuse of a nonempty output directory. Existing SCENIC output files can be replaced.
#' @return Invisibly, a list with the updated Seurat object, SCENIC options, regulon AUC, retained genes, metadata-name mapping, mean AUC and output directory.
#' @export
sce_run_scenic <- function(sce, db_dir, dbs = NULL, species = c("human", "mouse"),
                           assay = "RNA", layer = "data", output_dir = "scenic_output",
                           n_cores = 8L, min_counts_per_gene = NULL, min_samples = NULL,
                           top_n_regulons = 50L, top_n_plot = 6L, celltype = "celltype",
                           reduction = "umap", seed = 123L, save_pdf = TRUE,
                           save_png = TRUE, overwrite = FALSE) {
  species <- match.arg(species)
  if (!dir.exists(db_dir)) stop("cisTarget database directory does not exist: ", db_dir)
  db_dir <- normalizePath(db_dir, winslash = "/", mustWork = TRUE)
  if (is.null(dbs)) dbs <- list.files(db_dir, pattern = "\\.feather$")
  if (!length(dbs) || !all(file.exists(file.path(db_dir, dbs)))) {
    stop("Supply existing cisTarget .feather database filenames in `dbs`.")
  }
  for (x in list(n_cores, top_n_regulons, top_n_plot)) {
    if (length(x) != 1L || is.na(x) || x < 1 || x != as.integer(x)) stop("Counts must be positive integers.")
  }
  if (is.character(sce) && length(sce) == 1L) sce <- readRDS(sce)
  expr <- as.matrix(.biomed_sce_expression(sce, assay, layer))
  for (pkg in c("SCENIC", "AUCell", "Seurat", "pheatmap")) .biomed_require(pkg)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) && !overwrite) {
    stop("Output directory is not empty. Choose a new directory or set overwrite = TRUE.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  setwd(output_dir)
  withr::local_seed(seed)
  cell_info <- sce[[]]
  exprMat <- expr
  cellInfo <- cell_info
  save(exprMat, file = "exprMat_scenic.rda")
  save(cellInfo, file = "cellInfo_scenic.rda")
  saveRDS(expr, "exprMat_scenic.rds")
  saveRDS(cell_info, "cellInfo_scenic.rds")
  options <- SCENIC::initializeScenic(org = if (species == "human") "hgnc" else "mgi",
    dbDir = db_dir, dbs = dbs, datasetTitle = "biomed SCENIC", nCores = n_cores)
  options@settings$seed <- seed
  options@inputDatasetInfo$cellInfo <- "cellInfo_scenic.rds"
  if (is.null(min_counts_per_gene)) min_counts_per_gene <- 3 * 0.01 * ncol(expr)
  if (is.null(min_samples)) min_samples <- 0.01 * ncol(expr)
  genes <- SCENIC::geneFiltering(expr, options, minCountsPerGene = min_counts_per_gene,
                                minSamples = min_samples)
  if (!length(genes)) stop("No genes passed SCENIC filtering.")
  filtered <- expr[genes, , drop = FALSE]
  SCENIC::runCorrelation(filtered, options)
  SCENIC::runGenie3(filtered, options, resumePreviousRun = FALSE)
  options <- SCENIC::runSCENIC_1_coexNetwork2modules(options)
  options <- SCENIC::runSCENIC_2_createRegulons(options)
  options <- SCENIC::runSCENIC_3_scoreCells(options, exprMat = expr)
  saveRDS(options, "scenicOptions.rds")
  auc <- AUCell::getAUC(SCENIC::loadInt(options, "aucell_regulonAUC"))
  if (!all(colnames(sce) %in% colnames(auc))) stop("SCENIC AUC is missing cells.")
  auc <- auc[, colnames(sce), drop = FALSE]
  regulonAUC <- auc
  save(regulonAUC, file = "regulonAUC.rda")
  saveRDS(auc, "regulonAUC.rds")
  utils::write.csv(auc, "regulonAUC.csv")
  variance <- apply(auc, 1L, stats::var)
  selected <- utils::head(order(variance, decreasing = TRUE, na.last = NA), top_n_regulons)
  if (!length(selected)) stop("No regulons available for metadata export.")
  meta <- as.data.frame(t(auc[selected, , drop = FALSE]))
  original_names <- rownames(auc)[selected]
  existing <- colnames(sce[[]])
  new_names <- utils::tail(make.unique(c(existing, make.names(original_names))), ncol(meta))
  colnames(meta) <- new_names
  sce <- SeuratObject::AddMetaData(sce, meta)
  mapping <- data.frame(regulon = original_names, metadata_column = new_names)
  utils::write.csv(mapping, "regulon_metadata_mapping.csv", row.names = FALSE)
  saveRDS(sce, "sce_with_regulonAUC.rds")
  if ((save_pdf || save_png) && reduction %in% names(sce@reductions)) {
    p <- Seurat::FeaturePlot(sce, features = utils::head(new_names, top_n_plot), reduction = reduction)
    if (save_pdf) ggplot2::ggsave("regulon_FeaturePlot.pdf", p, width = 12, height = 8)
    if (save_png) ggplot2::ggsave("regulon_FeaturePlot.png", p, width = 12, height = 8, dpi = 300)
  }
  means <- NULL
  if (!is.null(celltype) && celltype %in% names(cell_info)) {
    groups <- as.character(cell_info[[celltype]])
    levels <- unique(groups[!is.na(groups)])
    if (length(levels)) {
      means <- vapply(levels, function(g) rowMeans(auc[selected, which(groups == g), drop = FALSE]),
                      numeric(length(selected)))
      rownames(means) <- original_names
      utils::write.csv(means, "regulonAUC_by_celltype.csv")
      draw <- function() pheatmap::pheatmap(means, cluster_rows = nrow(means) > 1L,
        cluster_cols = ncol(means) > 1L, main = "Mean regulon AUC by cell type")
      if (save_pdf) {
        grDevices::pdf("regulonAUC_heatmap.pdf", width = 10, height = 10)
        tryCatch(draw(), finally = grDevices::dev.off())
      }
      if (save_png) {
        grDevices::png("regulonAUC_heatmap.png", width = 3000, height = 3000, res = 300)
        tryCatch(draw(), finally = grDevices::dev.off())
      }
    }
  }
  invisible(list(sce = sce, scenic_options = options, regulon_auc = auc,
    genes_kept = genes, metadata_mapping = mapping, mean_auc_by_celltype = means,
    output_dir = output_dir))
}
