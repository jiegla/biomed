#' Analyze spatial neural-signaling gene expression and signatures
#'
#' Spatial transcriptomics workflow adapted from contributed analysis code.
#' Sample-level summaries are used for group comparisons. Spot-level correlations
#' and Moran statistics describe spatial association and do not establish causality.
#' @param vis A spatial Seurat object; niche-score analysis also accepts a spot-indexed metadata data frame.
#' @param condition_col Sample-level condition metadata column; required and constant within each sample.
#' @param sample_col Metadata column identifying independent samples. Leading-zero character IDs are preserved.
#' @param group_col Metadata grouping column. For neural analysis, NULL uses condition_col; this must be constant within each sample.
#' @param region_col Optional spot-level region column; NULL uses one All_spots region.
#' @param project_col Optional cohort metadata column used with project.
#' @param project Optional cohort values to retain before scoring and testing. NULL retains all cohorts; select an appropriate cohort explicitly when study and condition are confounded.
#' @param assay Expression assay. NULL uses the default assay.
#' @param layer Expression layer; split Seurat v5 layers are joined on a local copy and aligned by spot names.
#' @param pathways Named list of gene vectors. NULL uses the 12 contributed neural-signaling signatures. Scores average variable-gene z-scores across retained spots; they are relative signatures, not GSEA enrichment scores.
#' @param output_dir Report output directory. Neural analysis defaults to NULL (no files); gene comparison writes its Excel report by default.
#' @param save_pdf Save PDF figures.
#' @param save_png Save PNG figures.
#' @return Invisibly, a list with signature coverage, score_matrix, spot/sample/region summaries, condition and multigroup tests, four ggplot objects, parameters and output_dir.
#' @details Existing report files are overwritten in the selected output directory.
#' Cell-type fractions use the complete denominator cell-type set. Niche-score
#' global tests use sample random intercepts; pairwise tests match samples across niches.
#' Neural tests use one summary per sample; region summaries are descriptive.
#' @export
vis_analyze_neural_signaling_spatial <- function(
  vis, condition_col, sample_col = "sample_id", group_col = NULL,
  region_col = NULL, project_col = NULL, project = NULL,
  assay = NULL, layer = "data", pathways = NULL,
  output_dir = NULL, save_pdf = TRUE, save_png = TRUE
) {
  # Local bindings for columns evaluated in dplyr/ggplot data masks.
  analysis_group <- condition <- expression <- gene <- mean_expression <- mean_pct_positive <- median_expression <- median_score <- niche_annot_v2 <- pathway <- pct_positive <- sample_id <- score <- NULL


save_plot <- function(p, stem, w, h) {
  if (is.null(output_dir)) return(invisible(NULL))
  if (save_pdf) ggplot2::ggsave(paste0(stem, ".pdf"), p, width = w, height = h, limitsize = FALSE)
  if (save_png) ggplot2::ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300, limitsize = FALSE)
}
mean_na <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x[is.finite(x)])
median_na <- function(x) if (all(!is.finite(x))) NA_real_ else stats::median(x[is.finite(x)])

for (pkg in c("SeuratObject", "dplyr", "tidyr")) .biomed_require(pkg)
if (!inherits(vis, "Seurat")) stop("vis must be a Seurat object.")
meta <- vis[[]]
.biomed_required_columns(meta, c(sample_col, condition_col, group_col, region_col, project_col))
if (!is.null(project)) {
  if (is.null(project_col)) stop("project_col is required when project is supplied.")
  cells_keep <- rownames(meta)[!is.na(meta[[project_col]]) & meta[[project_col]] %in% project]
  if (!length(cells_keep)) stop("No spots match the requested project.")
  vis <- vis[, cells_keep]
  meta <- vis[[]]
}
if (is.null(pathways)) pathways <- .biomed_vis_neural_pathways()
if (!is.list(pathways) || !length(pathways) || is.null(names(pathways)) ||
    anyNA(names(pathways)) || any(!nzchar(names(pathways))) || anyDuplicated(names(pathways)) ||
    any(!vapply(pathways, function(x) is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x)), logical(1))))
  stop("pathways must be a named list of nonempty gene vectors with unique names.")
pathways <- lapply(pathways, unique)
reserved <- c("sample_id", "condition", "analysis_group", "niche_annot_v2", "gene", "pathway", "expression", "score")
if (any(c(names(pathways), unlist(pathways)) %in% reserved)) stop("Signature or gene names collide with reserved analysis columns.")
meta <- data.frame(sample_id = as.character(meta[[sample_col]]),
                   condition = as.character(meta[[condition_col]]),
                   analysis_group = as.character(meta[[if (is.null(group_col)) condition_col else group_col]]),
                   niche_annot_v2 = if (is.null(region_col)) "All_spots" else as.character(meta[[region_col]]),
                   row.names = rownames(meta))
keep <- !is.na(meta$sample_id) & nzchar(meta$sample_id)
meta <- meta[keep, , drop = FALSE]
if (!nrow(meta)) stop("No spots with valid sample IDs.")
for (column in c("condition", "analysis_group")) {
  bad <- vapply(split(meta[[column]], meta$sample_id), function(x) length(unique(x[!is.na(x)])) > 1L, logical(1))
  if (any(bad)) stop(column, " must be constant within each sample: ", paste(names(bad)[bad], collapse = ", "))
}
expr <- .biomed_sce_expression(vis, assay, layer)[, rownames(meta), drop = FALSE]
if (!any(unique(unlist(pathways)) %in% rownames(expr))) stop("None of the pathway genes were found in the selected layer.")
if (any(!is.finite(expr))) stop("Expression contains non-finite values.")
if (!is.null(output_dir)) {
  .biomed_require("openxlsx")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
out_dir <- output_dir

coverage <- dplyr::bind_rows(lapply(names(pathways), function(pw) {
  requested <- pathways[[pw]]
  found <- intersect(requested, rownames(expr))
  data.frame(pathway = pw, n_requested = length(requested), n_detected = length(found),
             detected_genes = paste(found, collapse = ";"), missing_genes = paste(setdiff(requested, found), collapse = ";"))
}))
if (!is.null(output_dir)) utils::write.csv(coverage, file.path(out_dir, "pathway_gene_coverage.csv"), row.names = FALSE)

all_genes <- unique(unlist(pathways))
genes_found <- intersect(all_genes, rownames(expr))
gene_mat <- expr[genes_found, , drop = FALSE]

# Gene-level spot summaries and sample-level medians.
gene_long <- as.data.frame(t(as.matrix(gene_mat)))
gene_long$sample_id <- meta[rownames(gene_long), "sample_id"]
gene_long$condition <- meta[rownames(gene_long), "condition"]
gene_long$analysis_group <- meta[rownames(gene_long), "analysis_group"]
gene_long$niche_annot_v2 <- meta[rownames(gene_long), "niche_annot_v2"]
gene_long <- tidyr::pivot_longer(gene_long, cols = dplyr::all_of(genes_found), names_to = "gene", values_to = "expression")

gene_sample <- gene_long |>
  dplyr::group_by(sample_id, condition, analysis_group, niche_annot_v2, gene) |>
  dplyr::summarise(n_spots = dplyr::n(), pct_positive = mean(expression > 0) * 100,
            median_expression = median_na(expression), mean_expression = mean_na(expression), .groups = "drop")

# Signature score: per-gene z-score across all spots in this cohort, averaged
# within each pathway. It is a relative pathway activity score, not GSEA NES.
score_mat <- matrix(NA_real_, nrow = length(pathways), ncol = ncol(expr),
                    dimnames = list(names(pathways), colnames(expr)))
for (pw in names(pathways)) {
  g <- intersect(pathways[[pw]], rownames(expr))
  if (length(g) == 0) next
  x <- as.matrix(expr[g, , drop = FALSE])
  mu <- rowMeans(x)
  s <- apply(x, 1, stats::sd)
  keep <- is.finite(s) & s > 0
  if (any(keep)) score_mat[pw, ] <- colMeans((x[keep, , drop = FALSE] - mu[keep]) / s[keep])
}

score_long <- as.data.frame(t(score_mat))
score_long$sample_id <- meta[rownames(score_long), "sample_id"]
score_long$condition <- meta[rownames(score_long), "condition"]
score_long$analysis_group <- meta[rownames(score_long), "analysis_group"]
score_long$niche_annot_v2 <- meta[rownames(score_long), "niche_annot_v2"]
score_long <- tidyr::pivot_longer(score_long, cols = dplyr::all_of(names(pathways)), names_to = "pathway", values_to = "score")

score_sample_region <- score_long |>
  dplyr::group_by(sample_id, condition, analysis_group, niche_annot_v2, pathway) |>
  dplyr::summarise(n_spots = dplyr::n(), median_score = median_na(score), mean_score = mean_na(score), .groups = "drop")

# Independent inferential unit: one spatial sample, never one spot or region.
score_sample <- score_long |>
  dplyr::group_by(sample_id, condition, analysis_group, pathway) |>
  dplyr::summarise(n_spots = dplyr::n(), median_score = median_na(score), mean_score = mean_na(score), .groups = "drop")

gene_sample_region <- gene_sample
gene_sample <- gene_long |>
  dplyr::group_by(sample_id, condition, analysis_group, gene) |>
  dplyr::summarise(n_spots = dplyr::n(), pct_positive = mean(expression > 0) * 100,
            median_expression = median_na(expression), mean_expression = mean_na(expression), .groups = "drop")

run_tests <- function(dat, group_col, value_col, feature_col) {
  dat <- dat[!is.na(dat[[group_col]]) & is.finite(dat[[value_col]]), , drop = FALSE]
  rows <- lapply(split(dat, dat[[feature_col]]), function(d) {
    lev <- unique(as.character(d[[group_col]]))
    if (length(lev) < 2L) return(NULL)
    z <- tryCatch(if (length(lev) == 2L) stats::wilcox.test(d[[value_col]] ~ d[[group_col]], exact = FALSE)
      else stats::kruskal.test(d[[value_col]] ~ d[[group_col]]), error = identity)
    data.frame(feature = d[[feature_col]][1], test = if (length(lev) == 2L) "Wilcoxon rank-sum" else "Kruskal-Wallis",
               n_groups = length(lev), statistic = if (inherits(z, "error")) NA_real_ else unname(z$statistic),
               p_value = if (inherits(z, "error")) NA_real_ else z$p.value,
               error = if (inherits(z, "error")) conditionMessage(z) else NA_character_)
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) out <- data.frame(feature = character(), test = character(), n_groups = integer(),
                                   statistic = numeric(), p_value = numeric(), error = character())
  out$p_adj_BH <- stats::p.adjust(out$p_value, method = "BH")
  out
}

condition_scores <- score_sample |> dplyr::filter(!is.na(condition))
condition_tests <- run_tests(condition_scores, "condition", "median_score", "pathway")
multigroup_scores <- score_sample |> dplyr::filter(!is.na(analysis_group))
multigroup_tests <- run_tests(multigroup_scores, "analysis_group", "median_score", "pathway")
condition_gene_tests <- run_tests(gene_sample |> dplyr::filter(!is.na(condition)), "condition", "median_expression", "gene")
multigroup_gene_tests <- run_tests(gene_sample |> dplyr::filter(!is.na(analysis_group)), "analysis_group", "median_expression", "gene")

theme_neural <- ggplot2::theme_classic(base_size = 10) + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
  strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70"),
  panel.spacing = grid::unit(0.35, "lines")
)

# Figure contract: hero panel = sample × pathway score heatmap; supporting panels
# = condition and clinically available multi-group distribution comparisons.
p_heat <- ggplot2::ggplot(score_sample |> dplyr::distinct(sample_id, condition, pathway, median_score),
                 ggplot2::aes(pathway, sample_id, fill = median_score)) +
  ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  ggplot2::facet_grid(condition ~ ., scales = "free_y", space = "free_y") + theme_neural +
  ggplot2::labs(title = "Neural-signaling signature scores by sample", x = NULL, y = "Sample", fill = "Median\nsignature score")
save_plot(p_heat, file.path(out_dir, "01_pathway_signature_heatmap"), 11, 8)

p_condition <- ggplot2::ggplot(condition_scores, ggplot2::aes(condition, median_score, fill = condition)) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.35, na.rm = TRUE) + ggplot2::geom_boxplot(width = .16, outlier.shape = NA) +
  ggplot2::geom_jitter(width = .12, size = 1.8, alpha = .85) + ggplot2::facet_wrap(~pathway, scales = "free_y", ncol = 4) +
  ggplot2::scale_fill_manual(values = biomed_colors(length(unique(condition_scores$condition)))) + theme_neural + ggplot2::theme(legend.position = "none") +
  ggplot2::labs(title = "Pathway signature scores by condition", x = NULL, y = "Sample median signature score")
save_plot(p_condition, file.path(out_dir, "02_condition_pathway_scores"), 12, 9)

p_multi <- ggplot2::ggplot(multigroup_scores, ggplot2::aes(analysis_group, median_score, fill = analysis_group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = .7) + ggplot2::geom_jitter(width = .13, size = 1.6, alpha = .85) +
  ggplot2::facet_wrap(~pathway, scales = "free_y", ncol = 4) + theme_neural + ggplot2::theme(legend.position = "none") +
  ggplot2::labs(title = "Pathway signature scores by analysis group", x = NULL, y = "Sample median signature score")
save_plot(p_multi, file.path(out_dir, "03_MG_status_pathway_scores"), 13, 9)

gene_group <- gene_sample |> dplyr::filter(!is.na(condition)) |>
  dplyr::group_by(condition, gene) |>
  dplyr::summarise(mean_expression = mean(median_expression, na.rm = TRUE), mean_pct_positive = mean(pct_positive, na.rm = TRUE), .groups = "drop")
p_gene <- ggplot2::ggplot(gene_group, ggplot2::aes(gene, condition, size = mean_pct_positive, fill = mean_expression)) +
  ggplot2::geom_point(shape = 21, colour = "grey20") + ggplot2::scale_size(range = c(1, 10), limits = c(0, 100)) +
  ggplot2::scale_fill_viridis_c() + theme_neural +
  ggplot2::labs(title = "Neural-signaling genes: sample-level expression", x = NULL, y = NULL, size = "% positive spots", fill = "Mean sample\nmedian expression")
save_plot(p_gene, file.path(out_dir, "04_condition_gene_dotplot"), 16, 4.5)

if (!is.null(output_dir)) {
wb <- openxlsx::createWorkbook()
add <- function(name, x) { openxlsx::addWorksheet(wb, name); openxlsx::writeData(wb, name, x); openxlsx::freezePane(wb, name, firstRow = TRUE) }
add("Pathway_gene_coverage", coverage)
add("Pathway_scores_spot", score_long)
add("Pathway_scores_sample", score_sample)
add("Pathway_scores_region", score_sample_region)
add("Pathway_test_condition", condition_tests)
add("Pathway_test_multigroup", multigroup_tests)
add("Gene_scores_sample", gene_sample)
add("Gene_scores_region", gene_sample_region)
add("Gene_test_condition", condition_gene_tests)
add("Gene_test_multigroup", multigroup_gene_tests)
add("Analysis_note", data.frame(note = c(
  paste0("Cohort filter: ", if (is.null(project)) "none" else paste(project, collapse = ", ")),
  paste0("Pathway score: mean of variable-gene z-scores in layer ", layer, "."),
  paste0("Analysis grouping column: ", if (is.null(group_col)) condition_col else group_col),
  "Sample is the inferential unit; region-level values are descriptive only and are not used as independent observations in tests."
)))
openxlsx::saveWorkbook(wb, file.path(out_dir, "neural_signaling_results.xlsx"), overwrite = TRUE)

utils::write.csv(condition_tests, file.path(out_dir, "pathway_condition_tests.csv"), row.names = FALSE)
utils::write.csv(multigroup_tests, file.path(out_dir, "pathway_multigroup_tests.csv"), row.names = FALSE)
utils::write.csv(condition_gene_tests, file.path(out_dir, "gene_condition_tests.csv"), row.names = FALSE)
message("Completed neural signaling analysis: ", normalizePath(out_dir))

}

invisible(list(coverage = coverage, score_matrix = score_mat, gene_sample = gene_sample,
  gene_sample_region = gene_sample_region, score_spot = score_long, score_sample = score_sample,
  score_sample_region = score_sample_region, condition_tests = condition_tests,
  multigroup_tests = multigroup_tests, condition_gene_tests = condition_gene_tests,
  multigroup_gene_tests = multigroup_gene_tests,
  plots = list(heatmap = p_heat, condition = p_condition, group = p_multi, genes = p_gene),
  parameters = list(sample_col = sample_col, condition_col = condition_col, group_col = group_col,
    region_col = region_col, project_col = project_col, project = project, assay = assay, layer = layer),
  output_dir = output_dir))
}
