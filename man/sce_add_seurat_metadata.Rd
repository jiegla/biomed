\name{sce_add_seurat_metadata}
\alias{sce_add_seurat_metadata}
\title{Add sample metadata to a Seurat object}
\usage{
sce_add_seurat_metadata(
    sce,
    metadata_file,
    cols = -1,
    sample_col = "orig.ident",
    replace_sample = FALSE,
    sheet = 1,
    file_encoding = NULL,
    strict_match = FALSE,
    allow_missing = FALSE
)
}
\arguments{
\item{sce}{A Seurat object; for SCENIC, also accepts the path to an RDS object.}
\item{metadata_file}{Path to CSV, XLSX, or XLS metadata. The first column holds sample IDs.}
\item{cols}{Tidyselect expression selecting incoming metadata columns; the first ID column is excluded.}
\item{sample_col}{Sample identifier column in Seurat metadata.}
\item{replace_sample}{Replace existing metadata columns when TRUE; otherwise create .new, .new2, etc.}
\item{sheet}{Excel sheet name or index.}
\item{file_encoding}{CSV encoding, or NULL to try detected encoding, UTF-8, GB18030 and GBK.}
\item{strict_match}{Require exactly the same sample ID sets in Seurat and the metadata file.}
\item{allow_missing}{Explicitly permit unmatched Seurat samples to receive NA. Default FALSE requires all Seurat samples to match.}
}
\description{Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.}
\value{The updated Seurat object, with cells in their original order.}
