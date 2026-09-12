\name{sce_analyze_annotated_seurat}
\alias{sce_analyze_annotated_seurat}
\title{Analyze an annotated Seurat object}
\usage{
sce_analyze_annotated_seurat(
    seurat,
    clinical_var = NULL,
    celltype = "cell_type__custom",
    statistic_level = c("samples"),
    sample_col = "orig.ident",
    assay = NULL,
    layer = "data",
    reduction = "umap",
    output_dir = "annotated_cell_analysis",
    file_prefix = "annotated",
    export_analysis_data = NULL,
    save_pdf = TRUE,
    save_png = TRUE,
    width = NULL,
    height = NULL,
    group_order = NULL,
    run_markers = TRUE,
    marker_only_pos = TRUE,
    marker_logfc_threshold = 0.25,
    marker_min_pct = 0.10,
    marker_padj_cutoff = 0.05,
    top_n_marker = 50,
    top_n_plot = 5,
    marker_test = "wilcox",
    run_celltype_correlation = TRUE,
    run_sample_celltype_correlation = TRUE,
    correlation_method = "spearman",
    correlation_features = c("variable", "all"),
    min_cells_per_sample_celltype = 10,
    run_diversity = TRUE,
    run_composition_correlation = TRUE,
    custom_colors = biomed_colors(),
    max_panels_per_page = 12,
    max_cells_heatmap = 3000,
    seed = 1234,
    verbose = TRUE
)
}
\arguments{
\item{seurat}{An annotated Seurat object.}
\item{clinical_var}{Clinical metadata columns for stratification. Each must be constant within a sample.}
\item{celltype}{Metadata column containing cell-type annotations.}
\item{statistic_level}{Use samples as independent units by default; cells requests exploratory cell-level comparisons.}
\item{sample_col}{Sample identifier column in Seurat metadata.}
\item{assay}{Assay to use. NULL uses the default assay; SCENIC defaults to RNA.}
\item{layer}{Expression layer to use; data denotes normalized expression. Split Seurat v5 layers are joined on a local copy.}
\item{reduction}{Existing dimensional reduction used for plots.}
\item{output_dir}{Directory for generated tables, figures, and analysis files.}
\item{file_prefix}{Prefix for generated filenames.}
\item{export_analysis_data}{Export underlying analysis data. NULL uses the original sample-level default.}
\item{save_pdf}{Save PDF figures.}
\item{save_png}{Save PNG figures.}
\item{width}{Plot width in inches; NULL selects automatic dimensions where supported.}
\item{height}{Plot height in inches; NULL selects automatic dimensions where supported.}
\item{group_order}{Optional group order; a named list provides separate orders for different grouping columns. Observed unlisted levels are retained.}
\item{run_markers}{Run marker discovery and export marker tables and plots.}
\item{marker_only_pos}{Restrict marker discovery to positive markers.}
\item{marker_logfc_threshold}{Log fold-change threshold passed to FindAllMarkers.}
\item{marker_min_pct}{Minimum detection fraction passed to FindAllMarkers.}
\item{marker_padj_cutoff}{Adjusted P-value cutoff for selecting markers; fallback ranking is reported when none pass.}
\item{top_n_marker}{Maximum markers retained per cell type in summary tables.}
\item{top_n_plot}{Maximum genes or regulons selected for feature plots.}
\item{marker_test}{Test passed to Seurat FindAllMarkers.}
\item{run_celltype_correlation}{Compute correlations between average cell-type expression profiles.}
\item{run_sample_celltype_correlation}{Compute sample correlations within each cell type.}
\item{correlation_method}{Correlation method passed to stats::cor.}
\item{correlation_features}{Use variable features when available, or all informative features.}
\item{min_cells_per_sample_celltype}{Minimum cells required for a sample/cell-type expression profile.}
\item{run_diversity}{Compute sample-level Shannon and Simpson diversity and comparisons.}
\item{run_composition_correlation}{Compute sample-level cell-type fraction correlations.}
\item{custom_colors}{Discrete palette; defaults to the 14 built-in biomed colors.}
\item{max_panels_per_page}{Maximum split panels on a plot page.}
\item{max_cells_heatmap}{Maximum number of cells sampled for marker heatmaps.}
\item{seed}{Random seed for reproducibility; the caller RNG state is restored.}
\item{verbose}{Print progress messages.}
}
\description{Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.}
\value{Invisibly, a list containing composition, marker, correlation, diversity, statistical tables, parameters and output directories.}
