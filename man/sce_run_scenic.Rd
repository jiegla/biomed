\name{sce_run_scenic}
\alias{sce_run_scenic}
\title{Run SCENIC and add regulon activity metadata}
\usage{
sce_run_scenic(sce, db_dir, dbs = NULL, species = c("human", "mouse"),
                           assay = "RNA", layer = "data", output_dir = "scenic_output",
                           n_cores = 8L, min_counts_per_gene = NULL, min_samples = NULL,
                           top_n_regulons = 50L, top_n_plot = 6L, celltype = "celltype",
                           reduction = "umap", seed = 123L, save_pdf = TRUE,
                           save_png = TRUE, overwrite = FALSE)
}
\arguments{
\item{sce}{A Seurat object; for SCENIC, also accepts the path to an RDS object.}
\item{db_dir}{Existing directory containing species-compatible cisTarget ranking databases.}
\item{dbs}{Database filenames within db_dir; NULL discovers feather files. Supply only compatible ranking databases.}
\item{species}{Human or mouse gene-symbol conventions.}
\item{assay}{Assay to use. NULL uses the default assay; SCENIC defaults to RNA.}
\item{layer}{Expression layer to use; data denotes normalized expression. Split Seurat v5 layers are joined on a local copy.}
\item{output_dir}{Directory for generated tables, figures, and analysis files.}
\item{n_cores}{Number of SCENIC worker cores.}
\item{min_counts_per_gene}{Minimum summed expression for geneFiltering; NULL uses 3 times 1 percent of the cell count.}
\item{min_samples}{Minimum cells expressing a gene; NULL uses 1 percent of the cell count.}
\item{top_n_regulons}{Number of most variable regulons added as metadata.}
\item{top_n_plot}{Maximum genes or regulons selected for feature plots.}
\item{celltype}{Metadata column containing cell-type annotations.}
\item{reduction}{Existing dimensional reduction used for plots.}
\item{seed}{Random seed for reproducibility; the caller RNG state is restored.}
\item{save_pdf}{Save PDF figures.}
\item{save_png}{Save PNG figures.}
\item{overwrite}{Permit reuse of a nonempty output directory. Existing SCENIC output files can be replaced.}
}
\description{Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.}
\value{Invisibly, a list with the updated Seurat object, SCENIC options, regulon AUC, retained genes, metadata-name mapping, mean AUC and output directory.}
\examples{
\dontrun{
sce <- readRDS("annotated_seurat.rds")
scenic <- sce_run_scenic(
  sce,
  db_dir = "cisTarget_databases",
  species = "human",
  celltype = "celltype",
  output_dir = "scenic_output",
  n_cores = 8
)
}
}
