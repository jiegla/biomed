\name{sce_compare_gene_expression_groups}
\alias{sce_compare_gene_expression_groups}
\title{Compare gene expression across groups}
\usage{
sce_compare_gene_expression_groups(
    sce,
    genes,
    statistic_group,
    split_by_celltype = NULL,
    statistic_level = c("samples", "cells"),
    sample_col = "orig.ident",
    Drop_negative = FALSE,
    assay = NULL,
    layer = "data",
    sample_aggregation = c("mean", "median"),
    group_order = NULL,
    reduction = "umap",
    output_dir = "gene_expression_comparison",
    file_prefix = "GeneExpression",
    export_analysis_data = NULL,
    save_pdf = TRUE,
    save_png = TRUE,
    width = 10,
    height = 7
)
}
\arguments{
\item{sce}{A Seurat object; for SCENIC, also accepts the path to an RDS object.}
\item{genes}{Character vector of genes to compare; missing genes are reported.}
\item{statistic_group}{One or more metadata columns defining comparison groups.}
\item{split_by_celltype}{Optional cell-type metadata column; NULL analyzes all cells together.}
\item{statistic_level}{Use samples as independent units by default; cells requests exploratory cell-level comparisons.}
\item{sample_col}{Sample identifier column in Seurat metadata.}
\item{Drop_negative}{Exclude expression values less than or equal to zero before testing/aggregation. Detection summaries still use all cells.}
\item{assay}{Assay to use. NULL uses the default assay; SCENIC defaults to RNA.}
\item{layer}{Expression layer to use; data denotes normalized expression. Split Seurat v5 layers are joined on a local copy.}
\item{sample_aggregation}{Sample-level expression summary: mean or median.}
\item{group_order}{Optional group order; a named list provides separate orders for different grouping columns. Observed unlisted levels are retained.}
\item{reduction}{Existing dimensional reduction used for plots.}
\item{output_dir}{Directory for generated tables, figures, and analysis files.}
\item{file_prefix}{Prefix for generated filenames.}
\item{export_analysis_data}{Export underlying analysis data. NULL uses the original sample-level default.}
\item{save_pdf}{Save PDF figures.}
\item{save_png}{Save PNG figures.}
\item{width}{Plot width in inches; NULL selects automatic dimensions where supported.}
\item{height}{Plot height in inches; NULL selects automatic dimensions where supported.}
}
\description{Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.}
\value{A named list per grouping variable containing parameters, summaries, tests, detection rates, analysis data, plots and output files.}
