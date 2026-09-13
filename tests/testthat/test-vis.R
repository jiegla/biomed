vis_fixture <- function() {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  set.seed(913)
  counts <- matrix(rpois(6 * 48, 4), 6, 48,
                   dimnames = list(c("ADRA1A", "ADRB1", "G1", "G2", "G3", "G4"), paste0("spot", 1:48)))
  x <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  x$sample_id <- rep(sprintf("%03d", 1:6), each = 8)
  x$condition <- rep(c("A", "B"), each = 24)
  x$niche <- rep(rep(c("N1", "N2"), each = 4), 6)
  x$project <- rep(c("study1", "study2"), each = 24)
  x$CT1 <- as.numeric(counts["G1", ])
  x$CT2 <- as.numeric(counts["G2", ])
  x
}

test_that("spatial fractions retain all denominator types and sample pairing", {
  skip_if_not_installed("nlme")
  x <- vis_fixture()
  out <- vis_analyze_celltype_composition(x, "niche", celltype_cols = "CT1",
    all_celltype_cols = c("CT1", "CT2"), outdir = tempfile(), save_individual = FALSE,
    save_multipage_pdf = FALSE, save_combined = FALSE, export_csv = FALSE,
    export_xlsx = FALSE, verbose = FALSE)
  d <- out$sample_group_data
  expect_equal(nrow(d), 12)
  md <- x[[]]
  expected <- aggregate(md[, c("CT1", "CT2")], md[, c("sample_id", "niche")], sum)
  idx <- match(paste(d$sample_id, d$group), paste(expected$sample_id, expected$niche))
  expect_equal(d$analysis_value, 100 * expected$CT1[idx] / (expected$CT1[idx] + expected$CT2[idx]))
  expect_equal(out$pairwise_stats$n_paired_samples, 6)
  md$CT1[1] <- -1
  expect_error(vis_analyze_celltype_composition(md, "niche", all_celltype_cols = c("CT1", "CT2"), verbose = FALSE), "nonnegative")
})

test_that("niche score files align by spot ID and pair only finite samples", {
  skip_if_not_installed("nlme")
  x <- vis_fixture()
  folder <- tempfile(); dir.create(folder)
  scores <- data.frame(P1 = seq_len(48) / 48, P2 = sin(seq_len(48)), row.names = colnames(x))
  saveRDS(scores[48:1, , drop = FALSE], file.path(folder, "scores.rds"))
  y <- vis_analyze_niche_pathway_scores(x, "niche", score_dir = folder,
    outdir = tempfile(), save_individual_png = FALSE, save_individual_pdf = FALSE,
    save_multipage_pdf = FALSE)
  expect_equal(nrow(y$sample_niche_score), 24)
  d <- y$sample_niche_score
  expect_equal(d$mean_score[d$pathway == "P1" & d$sample_id == "001" & d$niche == "N1"], mean((1:4)/48))
  expect_true(all(y$pairwise_stats$n_paired_samples == 6))
  scores$P1[1:4] <- Inf
  saveRDS(scores, file.path(folder, "scores.rds"))
  expect_warning(z <- vis_analyze_niche_pathway_scores(x, "niche", score_dir = folder,
    outdir = tempfile(), save_individual_png = FALSE, save_individual_pdf = FALSE,
    save_multipage_pdf = FALSE), "non-finite")
  expect_equal(z$pairwise_stats$n_paired_samples[z$pairwise_stats$pathway == "P1"], 5)
})

test_that("spatial gene summaries join split layers and use independent samples", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("openxlsx")
  x <- vis_fixture()
  expected <- tapply(as.numeric(SeuratObject::LayerData(x, layer = "counts")["G1", ]), x$sample_id, mean)
  alt <- SeuratObject::LayerData(x, layer = "counts") + 2
  x[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = methods::as(alt, "CsparseMatrix"))
  x[["ALT"]] <- split(x[["ALT"]], f = x$condition)
  expected <- expected + 2
  original_assay <- SeuratObject::DefaultAssay(x)
  state <- .Random.seed
  y <- vis_compare_spatial_gene_expression_groups(x, "G1", "condition", sample_col = "sample_id",
    assay = "ALT", layer = "counts", sample_aggregation = "mean", run_spatial_plots = FALSE, run_moran = FALSE,
    output_dir = tempfile(), save_pdf = FALSE, save_png = FALSE,
    cell2location_metadata_cols = "CT2", min_spots_for_correlation = 3)
  expect_identical(.Random.seed, state)
  expect_identical(SeuratObject::DefaultAssay(x), original_assay)
  expect_equal(nrow(y$sample_gene_summary), 6)
  expect_equal(y$sample_gene_summary$expression, as.numeric(expected[y$sample_gene_summary$.sample]))
  expect_equal(nrow(y$cell2location$correlations), 6)
  expect_true(file.exists(y$files$xlsx))
  x$condition[1] <- "other"
  expect_error(vis_compare_spatial_gene_expression_groups(x, "G1", "condition", sample_col = "sample_id",
    layer = "counts", output_dir = tempfile(), run_moran = FALSE), "not constant")
})

test_that("neural scoring matches gene z-scores and handles a single condition", {
  x <- vis_fixture()
  expr <- as.matrix(SeuratObject::LayerData(x, layer = "counts"))
  y <- vis_analyze_neural_signaling_spatial(x, "condition", layer = "counts",
    region_col = "niche", pathways = list(neural = c("ADRA1A", "ADRB1")))
  expected <- colMeans(t(scale(t(expr[c("ADRA1A", "ADRB1"), ]))))
  expect_equal(as.numeric(y$score_matrix["neural", ]), as.numeric(expected))
  expect_equal(nrow(y$score_sample), 6)
  expect_equal(nrow(y$score_sample_region), 12)
  expect_null(y$output_dir)
  expect_s3_class(y$plots$heatmap, "ggplot")
  z <- vis_analyze_neural_signaling_spatial(x, "condition", layer = "counts",
    project_col = "project", project = "study1", pathways = list(neural = "ADRA1A"))
  expect_equal(nrow(z$score_sample), 3)
  expect_equal(nrow(z$condition_tests), 0)
  x$condition[1] <- "other"
  expect_error(vis_analyze_neural_signaling_spatial(x, "condition", layer = "counts"), "constant within")
})

test_that("spatial images yield matched coordinates, Moran tests and figure exports", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ape")
  skip_if_not_installed("openxlsx")
  x <- vis_fixture()
  for (id in unique(x$sample_id)) {
    cells <- colnames(x)[x$sample_id == id]
    coords <- data.frame(tissue = 1L, row = rep(1:2, 4), col = rep(1:4, each = 2),
                         imagerow = rep(1:2, 4), imagecol = rep(1:4, each = 2), row.names = cells)
    image <- methods::new("VisiumV1", image = array(1, c(8, 8, 3)),
      scale.factors = Seurat::scalefactors(), coordinates = coords, spot.radius = 0.05,
      assay = "RNA", key = paste0("image", id, "_"))
    x[[paste0("image", id)]] <- image
  }
  SeuratObject::LayerData(x, layer = "custom") <- SeuratObject::LayerData(x, layer = "counts")
  y <- vis_compare_spatial_gene_expression_groups(x, "G1", sample_col = "sample_id",
    layer = "custom", run_spatial_plots = TRUE, run_moran = TRUE, moran_max_spots = 6,
    moran_large_sample = "subsample", output_dir = tempfile(), save_pdf = TRUE, save_png = FALSE)
  expect_equal(nrow(y$coordinates), 48)
  expect_equal(nrow(y$moran), 6)
  expect_true(all(y$moran$n_spots_moran == 6))
  expect_equal(nrow(y$image_sample_map), 6)
  expect_true(length(y$files$spatial_feature_plots) > 0)
  expect_true(all(file.exists(y$files$spatial_feature_plots)))
})

test_that("neural reports export workbook and reusable PDF figures", {
  skip_if_not_installed("openxlsx")
  x <- vis_fixture()
  folder <- tempfile()
  y <- vis_analyze_neural_signaling_spatial(x, "condition", layer = "counts",
    pathways = list(neural = c("ADRA1A", "ADRB1")), output_dir = folder, save_png = FALSE)
  expect_true(file.exists(file.path(folder, "neural_signaling_results.xlsx")))
  expect_length(list.files(folder, pattern = "[.]pdf$"), 4)
  expect_equal(nrow(y$condition_tests), 1)
})
