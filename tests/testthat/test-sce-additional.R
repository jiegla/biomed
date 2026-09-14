additional_sce_fixture <- function(ncells = 40L) {
  skip_if_not_installed("Seurat")
  counts <- matrix(rep(seq_len(30), ncells), nrow = 30,
                   dimnames = list(paste0("G", 1:30), paste0("cell", seq_len(ncells))))
  counts[1, ] <- rep(c(0, 1, 3, 7), length.out = ncells)
  x <- SeuratObject::CreateSeuratObject(counts)
  x$sampleid <- rep(paste0("s", 1:4), length.out = ncells)
  x$condition <- ifelse(x$sampleid %in% c("s1", "s2"), "A", "B")
  x$celltype <- rep(rep(c("T", "B"), each = 4), length.out = ncells)
  Seurat::NormalizeData(x, verbose = FALSE)
}

test_that("cell fractions preserve the original schema and Excel export", {
  skip_if_not_installed("openxlsx")
  x <- additional_sce_fixture()
  file <- file.path(tempdir(), "fractions", "fractions.xlsx")
  out <- sce_calculate_cell_fraction(x, file, "sampleid", "celltype")
  expect_named(out, c("orig.ident", "celltype_major", "cell_number", "fraction"))
  expect_equal(sum(out$cell_number), ncol(x))
  expect_equal(unname(tapply(out$fraction, out$orig.ident, sum)), rep(1, 4))
  expect_true(file.exists(file))
  expect_equal(nrow(openxlsx::read.xlsx(file)), nrow(out))
})

test_that("expression statistics use sample summaries and preserve reporting modes", {
  x <- additional_sce_fixture()
  path <- tempfile("expression-report-")
  args <- list(object = x, genes = "G1", sample_col = "sampleid",
               celltype_col = "celltype", group_col = "condition",
               min_cells = 2, make_plot = FALSE)
  r <- do.call(sce_gene_expression_statistic, c(args, list(output_dir = path)))
  mat <- as.matrix(biomed:::.biomed_sce_expression(x))
  expected <- aggregate(as.numeric(mat["G1", ]),
                        list(sample_id = x$sampleid, celltype = x$celltype), mean)
  joined <- merge(as.data.frame(r$sample_summary), expected, by = c("sample_id", "celltype"))
  expect_equal(joined$expression_value, joined$x)
  expect_equal(nrow(r$analysis_data), 8L)
  expect_equal(r$overall_statistics$p_adjust_global,
               p.adjust(r$overall_statistics$p_value, "BH"))
  expect_true(any(grepl("\\.rds$", list.files(path), ignore.case = TRUE)))
  expect_true(any(grepl("\\.(xlsx|csv)$", list.files(path), ignore.case = TRUE)))
  for (i in seq_len(nrow(r$pairwise_statistics))) {
    z <- r$pairwise_statistics[i, ]
    d <- r$analysis_data[r$analysis_data$celltype == z$celltype, ]
    a <- d$expression_value[d$group == z$group1]
    b <- d$expression_value[d$group == z$group2]
    expect_equal(z$p_value, suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value))
  }
  cell <- do.call(sce_gene_expression_statistic, c(args, list(statistic_by = "cell_level")))
  expect_equal(nrow(cell$analysis_data), ncol(x))
  positive <- do.call(sce_gene_expression_statistic, c(args, list(remove_zero = TRUE)))
  expect_true(all(positive$expression_long$expression > 0))
  expect_true(all(positive$sample_summary$positive_ratio == 1))
  expect_equal(ncol(x), 40L)
})

test_that("DEG runs without GMT and retains configuration and exports", {
  skip_if_not_installed("openxlsx")
  x <- additional_sce_fixture()
  path <- tempfile("deg-report-")
  set.seed(456)
  old_seed <- .Random.seed
  r <- sce_run_two_group_deg_enrichment(
    x, "condition", "A", "B", output_dir = path,
    config = list(de = list(min_cells_per_group = 5, join_layers = FALSE),
                  plot = list(run = FALSE)))
  expect_identical(.Random.seed, old_seed)
  expect_false(r$config$ora$run)
  expect_false(r$config$gsea$run)
  expect_equal(r$status$status[r$status$step == "DEG"], "success")
  expect_true(nrow(r$deg$all) > 0)
  expect_true(file.exists(file.path(path, "01_DEG", "DEG_by_celltype.xlsx")))
  expect_true(file.exists(file.path(path, "two_group_deg_enrichment_results.rda")))
  env <- new.env()
  load(file.path(path, "two_group_deg_enrichment_results.rda"), env)
  expect_equal(env$results$deg, r$deg)
})

test_that("KO QC arguments adapt without silently dropping original settings", {
  args <- list(countMatrix = matrix(1, 10, 100), qc = TRUE)
  old <- biomed:::.biomed_knk_qc_arguments(args, "qc_minLSize", 0.2, 500, 25)
  new <- biomed:::.biomed_knk_qc_arguments(args, "qc_minLibSize", 0.2, 500, 25)
  expect_equal(old$qc_minCells, 25)
  expect_equal(old$qc_mtThreshold, 0.2)
  expect_equal(new$qc_minPCT, 0.25)
  expect_equal(new$qc_maxMTratio, 0.2)
  expect_equal(new$qc_minLibSize, 500)
  expect_true(is.null(new$qc_minCells))
})

test_that("KO wrapper preserves reports with a deterministic backend fixture", {
  skip_if_not_installed("scTenifoldKnk")
  skip_if_not_installed("openxlsx")
  local_mocked_bindings(scTenifoldKnk = function(countMatrix, ...) {
    genes <- rownames(countMatrix)
    list(diffRegulation = data.frame(
      gene = genes, distance = seq_along(genes), FC = seq_along(genes),
      p.value = rep(0.001, length(genes)), p.adj = rep(0.01, length(genes))))
  }, .package = "scTenifoldKnk")
  counts <- matrix(2, 12, 100, dimnames = list(paste0("G", 1:12), paste0("c", 1:100)))
  path <- tempfile("ko-report-")
  r <- sce_run_scTenifoldKnk_KO(counts, "G1", output_dir = path,
                               save_network_plot = FALSE, top_n = 3, label_n = 0)
  expect_named(r, "G1")
  expect_equal(nrow(r$G1$diffRegulation), 12L)
  expect_equal(nrow(r$G1$top), 3L)
  expect_true(file.exists(file.path(r$G1$run_dir, "analysis_summary.txt")))
  expect_true(file.exists(file.path(r$G1$run_dir, "06_analysis_bundle.rds")))
  expect_true(file.exists(file.path(r$G1$run_dir, "02_diffRegulation_all.xlsx")))
  expect_length(list.files(file.path(r$G1$run_dir, "figures"), pattern = "\\.png$"), 3)
})

additional_spatial_fixture <- function() {
  x <- additional_sce_fixture(16L)
  xy <- expand.grid(imagecol = seq(10, 40, 10), imagerow = seq(10, 40, 10))
  coords <- data.frame(tissue = 1, row = rep(1:4, 4), col = rep(1:4, each = 4),
                       imagerow = xy$imagerow, imagecol = xy$imagecol,
                       row.names = colnames(x))
  img <- methods::new("VisiumV1", image = array(1, c(50, 50, 3)),
                      scale.factors = Seurat::scalefactors(spot = 1, fiducial = 1, hires = 1, lowres = 1),
                      coordinates = coords, spot.radius = 0.02,
                      assay = "RNA", key = "slice_")
  x[["slice"]] <- img
  x$signal <- seq_len(ncol(x))
  x
}

test_that("spatial surfaces align spot barcodes and reject invalid smoothing", {
  skip_if_not_installed("FNN")
  x <- additional_spatial_fixture()
  r <- sce_plot_spatial_continuous(x, "signal", grid_n = 12, smooth_k = 2,
                                   cap_quantile = NULL, show_he = FALSE, return_data = TRUE)
  expect_s3_class(r$plot, "ggplot")
  expect_equal(r$spot_data$value, unname(x$signal[match(r$spot_data$spot, colnames(x))]))
  expect_true(any(is.finite(r$surface_data$z)))
  expect_equal(r$parameters$spot_spacing, 10)
  expect_error(sce_plot_spatial_continuous(x, "signal", sigma_factor = 0), "positive")
  x$constant <- 1
  expect_error(sce_plot_spatial_continuous(x, "constant", transform = "zscore"), "constant")
})

test_that("pathway RDS orientation and literal case-insensitive selection are retained", {
  skip_if_not_installed("FNN")
  x <- additional_spatial_fixture()
  score_dir <- tempfile("scores-")
  dir.create(score_dir)
  scores <- matrix(seq_len(16), 1, dimnames = list("HALLMARK_TEST", rev(colnames(x))))
  saveRDS(scores, file.path(score_dir, "scores.rds"))
  out <- tempfile("pathway-report-")
  r <- sce_plot_pathway_rds_continuous(x, score_dir = score_dir,
    selected_pathways = "hallmark", pathway_match = "contains", output_dir = out,
    grid_n = 12, show_he = FALSE, save_pdf = FALSE, save_png = FALSE, keep_plots = TRUE,
    verbose = FALSE)
  expect_equal(r$summary$pathway, "HALLMARK_TEST")
  expect_equal(r$summary$valid_spots, 16L)
  expect_equal(r$log$status, "success")
  expect_length(r$plots, 1L)
  expect_true(file.exists(file.path(out, "pathway_plot_summary.csv")))
  expect_false(any(grepl("TEMP_UCELL", names(x[[]]))))
})
