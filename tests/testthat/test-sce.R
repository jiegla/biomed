sce_fixture <- function() {
  skip_if_not_installed("Seurat")
  set.seed(71)
  genes <- unique(c("G1", "G2", Seurat::cc.genes.updated.2019$s.genes,
                    Seurat::cc.genes.updated.2019$g2m.genes, paste0("BG", 1:600)))
  counts <- Matrix::Matrix(
    matrix(rpois(length(genes) * 40, 4), length(genes), 40,
           dimnames = list(genes, paste0("c", 1:40))),
    sparse = TRUE
  )
  x <- SeuratObject::CreateSeuratObject(counts)
  x$orig.ident <- rep(c("001", "002", "003", "004"), each = 10)
  x$group <- rep(c("A", "B"), each = 20)
  x$celltype <- rep(c("T", "B"), 20)
  x <- Seurat::NormalizeData(x, verbose = FALSE)
  x[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(rnorm(80), 40, 2, dimnames = list(colnames(x), c("UMAP_1", "UMAP_2"))),
    key = "UMAP_", assay = "RNA")
  x
}

test_that("sample metadata preserves leading zeros, aligns cells and handles collisions", {
  x <- sce_fixture()
  path <- tempfile(fileext = ".csv")
  writeLines(c("sample,group,value", "004,D,4", "002,B,2", "001,A,1", "003,C,3", "005,E,5"), path)
  y <- sce_add_seurat_metadata(x, path)
  expect_identical(colnames(y), colnames(x))
  expect_identical(y$group, x$group)
  expect_equal(unname(y$value), rep(1:4, each = 10))
  expect_equal(unname(y$group.new), rep(LETTERS[1:4], each = 10))
  expect_error(sce_add_seurat_metadata(x, path, strict_match = TRUE), "exactly")
  writeLines(c("sample,value", "001,1"), path)
  expect_error(sce_add_seurat_metadata(x, path), "Metadata missing")
  y <- sce_add_seurat_metadata(x, path, allow_missing = TRUE)
  expect_equal(sum(is.na(y$value)), 30)
  writeLines(c("sample,value", "001,1", "001,2"), path)
  expect_error(sce_add_seurat_metadata(x, path))
})

test_that("UTF-8 and GB18030 metadata are decoded without losing sample IDs", {
  x <- sce_fixture()
  text <- paste(c("sample,label", paste0(sprintf("%03d", 1:4), ",\u4e34\u5e8a")), collapse = "\n")
  for (encoding in c("UTF-8", "GB18030")) {
    bytes <- iconv(text, from = "UTF-8", to = encoding)
    expect_false(is.na(bytes))
    path <- tempfile(fileext = ".csv")
    writeBin(charToRaw(bytes), path)
    y <- sce_add_seurat_metadata(x, path, file_encoding = encoding)
    expect_true(all(y$label == "\u4e34\u5e8a"))
  }
})

test_that("expression helper selects the requested assay and rejoins split layers", {
  x <- sce_fixture()
  expected <- biomed:::.biomed_sce_expression(x)
  alt <- methods::as(expected + 2, "dgCMatrix")
  x[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = alt, data = alt)
  expect_equal(as.matrix(biomed:::.biomed_sce_expression(x, "ALT")), as.matrix(expected + 2))
  x[["RNA"]] <- split(x[["RNA"]], f = x$group)
  expect_equal(as.matrix(biomed:::.biomed_sce_expression(x)), as.matrix(expected))
})

test_that("gene comparisons aggregate the selected assay at sample level", {
  x <- sce_fixture()
  mat <- biomed:::.biomed_sce_expression(x)
  alt <- methods::as(mat + 2, "dgCMatrix")
  x[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = alt, data = alt)
  result <- sce_compare_gene_expression_groups(x, "G1", "group", assay = "ALT",
    output_dir = tempfile(), save_pdf = FALSE, save_png = FALSE)
  dat <- result$group$analysis_data
  expected <- tapply(as.numeric(mat["G1", ]) + 2, x$orig.ident, mean)
  expect_equal(dat$expression, as.numeric(expected[dat$sample]))
  expect_equal(nrow(dat), 4)
  expect_true(file.exists(result$group$files$xlsx))
})

test_that("annotation workflow preserves composition totals and exports results", {
  x <- sce_fixture()
  y <- sce_analyze_annotated_seurat(x, clinical_var = "group", celltype = "celltype",
    output_dir = tempfile(), run_markers = FALSE, run_celltype_correlation = FALSE,
    run_sample_celltype_correlation = FALSE, run_composition_correlation = FALSE,
    save_pdf = FALSE, save_png = FALSE, verbose = FALSE)
  fractions <- y$sample_composition_long
  expect_equal(as.numeric(tapply(fractions$fraction, fractions$orig.ident, sum)), rep(1, 4))
  expect_equal(sum(y$overall_composition$n_cells), ncol(x))
  expect_equal(nrow(y$diversity), 4)
  expect_equal(y$parameters$custom_colors, biomed_colors())
})

test_that("all scoring backends return cell-aligned signature scores", {
  x <- sce_fixture()
  sets <- list(one = rownames(x)[1:25], two = rownames(x)[26:50])
  for (method in c("UCell", "GSVA", "AUCell")) {
    expect_true(suppressWarnings(requireNamespace(method, quietly = TRUE)))
    scores <- suppressWarnings(
      sce_run_scoring(x, sets, method = method, return_type = "matrix")
    )
    expect_equal(dim(scores), c(40L, 2L))
    expect_identical(rownames(scores), colnames(x))
    expect_true(all(is.finite(scores)))
  }
})

test_that("cell-cycle scoring preserves identities and adds metadata", {
  x <- sce_fixture()
  y <- sce_cell_cycle(x, nbin = 5, ctrl = 10)
  expect_identical(SeuratObject::Idents(x), SeuratObject::Idents(y))
  expect_true(all(c("S.Score", "G2M.Score", "Phase") %in% names(y[[]])))
  expect_true(all(is.finite(y$S.Score)))
})

test_that("SCENIC requires an explicit existing database before running", {
  expect_error(sce_run_scenic(NULL, tempfile()), "database directory")
  folder <- tempfile()
  dir.create(folder)
  expect_error(sce_run_scenic(NULL, folder), "database filenames")
})
