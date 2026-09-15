# Supplementary examples for exported biomed functions.
#
# These blocks are intentionally kept separate from implementation code so that
# examples can be improved without touching analysis logic. Running
# `devtools::document()` merges them into the corresponding help topics.

#' @rdname analysis_utilities
#' @examples
#' clinical <- data.frame(
#'   group = factor(rep(c("A", "B"), each = 10)),
#'   response = factor(rep(c(0, 1), 10)),
#'   time = seq_len(20),
#'   status = rep(c(1, 0, 1, 1, 0), 4),
#'   marker1 = c(1:10, 5:14),
#'   marker2 = c(10:1, 8:17)
#' )
#'
#' # Create Low/High groups from a continuous biomarker
#' grouped <- makegroup(clinical, "marker1", method = "median")
#'
#' # Batch ANOVA and categorical association
#' anova_batch(clinical, target = "group", feature = c("marker1", "marker2"))
#' chi_square_batch_df(clinical, "group", response = "response")
#'
#' # Univariable Cox regression
#' cox_batch(c("marker1", "marker2"), clinical, time = "time", status = "status")
#'
#' \dontrun{
#' # ROC summaries and publication-style ROC plots require pROC/ggplot2
#' roc_batch(c("marker1", "marker2"), clinical, outcome = "response")
#' plot_roc_batch(
#'   clinical,
#'   response = "response",
#'   variables = c("marker1", "marker2"),
#'   fig.path = "ROC"
#' )
#' }
NULL

#' @rdname biomed_standard_api
#' @examples
#' clinical <- data.frame(
#'   arm = factor(rep(c("TKI", "ICI"), each = 10)),
#'   response = factor(rep(c("No", "Yes"), 10)),
#'   PFS = seq_len(20),
#'   PFS_status = rep(c(1, 0, 1, 1, 0), 4),
#'   PD_L1 = c(1:10, 5:14),
#'   TMB = c(10:1, 8:17)
#' )
#'
#' make_group(clinical, "PD_L1", method = "median")
#' batch_anova(clinical, variables = c("PD_L1", "TMB"), group = "arm")
#' batch_chi_square(clinical, variables = "arm", response = "response")
#' batch_cox(clinical, variables = c("PD_L1", "TMB"),
#'           time = "PFS", status = "PFS_status")
#'
#' \dontrun{
#' batch_roc(clinical, variables = c("PD_L1", "TMB"), response = "response")
#' plot_roc(clinical, variables = c("PD_L1", "TMB"), response = "response")
#' plot_donut(clinical, variables = c("arm", "response"))
#' plot_stacked_bar(clinical, variables = "response", group = "arm")
#' plot_venn(list(TKI = c("EGFR", "TP53"), ICI = c("TP53", "LRP1B")))
#' }
NULL

#' @rdname best_cutoff
#' @examples
#' \dontrun{
#' surv_data <- data.frame(
#'   time = c(6, 8, 12, 15, 18, 21, 24, 30, 36, 42),
#'   status = c(1, 1, 0, 1, 0, 1, 0, 1, 0, 0),
#'   biomarker = c(2.1, 3.0, 1.7, 5.2, 4.1, 6.3, 2.8, 7.1, 3.5, 5.9)
#' )
#' best_cutoff(surv_data, "biomarker", time = "time", status = "status")
#' }
NULL

#' @rdname find_survival_cutoff
#' @examples
#' \dontrun{
#' surv_data <- data.frame(
#'   time = c(6, 8, 12, 15, 18, 21, 24, 30, 36, 42),
#'   status = c(1, 1, 0, 1, 0, 1, 0, 1, 0, 0),
#'   biomarker = c(2.1, 3.0, 1.7, 5.2, 4.1, 6.3, 2.8, 7.1, 3.5, 5.9)
#' )
#' cutoff <- find_survival_cutoff(
#'   surv_data, "biomarker", time = "time", status = "status",
#'   fallback = "median"
#' )
#' cutoff$cutoff
#' cutoff$data
#' }
NULL

#' @rdname batch_surv
#' @examples
#' \dontrun{
#' surv_data <- data.frame(
#'   time = seq(3, 90, length.out = 40),
#'   status = rep(c(1, 0, 1, 1), 10),
#'   marker1 = seq_len(40),
#'   marker2 = rev(seq_len(40))
#' )
#' batch_surv(
#'   surv_data,
#'   variable = c("marker1", "marker2"),
#'   time = "time",
#'   status = "status",
#'   min_sample = 20,
#'   verbose = FALSE
#' )
#' }
NULL

#' @rdname batch_survival
#' @examples
#' \dontrun{
#' surv_data <- data.frame(
#'   PFS = seq(3, 90, length.out = 40),
#'   PFS_status = rep(c(1, 0, 1, 1), 10),
#'   PD_L1 = seq_len(40),
#'   TMB = rev(seq_len(40))
#' )
#' result <- batch_survival(
#'   surv_data,
#'   variables = c("PD_L1", "TMB"),
#'   time = "PFS",
#'   status = "PFS_status"
#' )
#' result$results
#' result$failed
#' }
NULL

#' @rdname pairwise_survival
#' @examples
#' \dontrun{
#' surv_data <- data.frame(
#'   arm = factor(rep(c("Chemo", "ICI", "AK112"), each = 12)),
#'   PFS = seq(2, 72, length.out = 36),
#'   PFS_status = rep(c(1, 0, 1), 12)
#' )
#' pairwise_survival(
#'   surv_data,
#'   group = "arm",
#'   time = "PFS",
#'   status = "PFS_status",
#'   reference = "Chemo"
#' )
#' }
NULL

#' @rdname plot_survival
#' @examples
#' \dontrun{
#' surv_data <- data.frame(
#'   arm = factor(rep(c("Control", "Treatment"), each = 20)),
#'   PFS = seq(2, 80, length.out = 40),
#'   PFS_status = rep(c(1, 0, 1, 1), 10)
#' )
#' p <- plot_survival(
#'   surv_data,
#'   group = "arm",
#'   time = "PFS",
#'   status = "PFS_status",
#'   max_time = 60,
#'   show_hr = TRUE,
#'   risk_table = TRUE
#' )
#' p
#' }
NULL

#' @rdname surv_fig_hr
#' @examples
#' \dontrun{
#' surv_data <- data.frame(
#'   Treatment = factor(rep(c("Control", "Treatment"), each = 20)),
#'   PFS = seq(2, 80, length.out = 40),
#'   PFS_status = rep(c(1, 0, 1, 1), 10)
#' )
#' surv_fig_hr(
#'   group_var = "Treatment",
#'   df = surv_data,
#'   time = "PFS",
#'   status = "PFS_status",
#'   max_time = 60,
#'   output_dir = NULL
#' )
#' }
NULL

#' @rdname get_cor
#' @examples
#' \dontrun{
#' biomarker <- data.frame(
#'   PD_L1 = c(1, 5, 10, 20, 30, 50, 70, 90),
#'   TMB = c(2.1, 2.4, 3.2, 4.0, 5.1, 6.2, 7.0, 8.4),
#'   arm = rep(c("A", "B"), 4)
#' )
#' get_cor(
#'   biomarker,
#'   var1 = "PD_L1",
#'   var2 = "TMB",
#'   subtype = "arm",
#'   scale = FALSE,
#'   show_plot = FALSE,
#'   save_plot = FALSE
#' )
#' }
NULL

#' @rdname plot_correlation
#' @examples
#' \dontrun{
#' biomarker <- data.frame(
#'   PD_L1 = c(1, 5, 10, 20, 30, 50, 70, 90),
#'   TMB = c(2.1, 2.4, 3.2, 4.0, 5.1, 6.2, 7.0, 8.4),
#'   arm = rep(c("A", "B"), 4)
#' )
#' plot_correlation(
#'   biomarker,
#'   x = "PD_L1",
#'   y = "TMB",
#'   group = "arm",
#'   method = "spearman",
#'   scale = FALSE
#' )
#' }
NULL

#' @rdname sce_add_seurat_metadata
#' @examples
#' \dontrun{
#' sce <- readRDS("annotated_seurat.rds")
#' sce <- sce_add_seurat_metadata(
#'   sce,
#'   metadata_file = "sample_metadata.xlsx",
#'   sample_col = "orig.ident"
#' )
#' }
NULL

#' @rdname sce_analyze_annotated_seurat
#' @examples
#' \dontrun{
#' sce <- readRDS("annotated_seurat.rds")
#' result <- sce_analyze_annotated_seurat(
#'   sce,
#'   clinical_var = c("response", "stage"),
#'   celltype = "celltype",
#'   sample_col = "orig.ident",
#'   output_dir = "annotated_cell_analysis"
#' )
#' result$composition
#' }
NULL

#' @rdname sce_calculate_cell_fraction
#' @examples
#' \dontrun{
#' sce <- readRDS("annotated_seurat.rds")
#' fraction <- sce_calculate_cell_fraction(
#'   sce,
#'   sample_col = "orig.ident",
#'   celltype_col = "celltype",
#'   output_file = "celltype_fraction_by_sample.xlsx"
#' )
#' head(fraction)
#' }
NULL

#' @rdname sce_cell_cycle
#' @examples
#' \dontrun{
#' sce <- readRDS("seurat.rds")
#' sce <- sce_cell_cycle(sce, species = "human", assay = "RNA")
#' table(sce$Phase)
#' }
NULL

#' @rdname sce_compare_gene_expression_groups
#' @examples
#' \dontrun{
#' sce <- readRDS("annotated_seurat.rds")
#' result <- sce_compare_gene_expression_groups(
#'   sce,
#'   genes = c("PDCD1", "LAG3", "CXCL13"),
#'   statistic_group = "response",
#'   split_by_celltype = "celltype",
#'   sample_col = "orig.ident",
#'   statistic_level = "samples"
#' )
#' }
NULL

#' @rdname sce_gene_expression_statistic
#' @examples
#' \dontrun{
#' sce <- readRDS("annotated_seurat.rds")
#' result <- sce_gene_expression_statistic(
#'   sce,
#'   genes = c("PDCD1", "LAG3"),
#'   sample_col = "orig.ident",
#'   celltype_col = "celltype",
#'   group_col = "response",
#'   statistic_by = "sample_level",
#'   output_dir = "gene_expression_statistics"
#' )
#' result$overall_statistics
#' }
NULL

#' @rdname sce_run_scoring
#' @examples
#' \dontrun{
#' sce <- readRDS("annotated_seurat.rds")
#' gene_sets <- list(
#'   Exhaustion = c("PDCD1", "LAG3", "TIGIT", "HAVCR2"),
#'   Cytotoxicity = c("NKG7", "GNLY", "GZMB", "PRF1")
#' )
#' sce <- sce_run_scoring(
#'   sce,
#'   gene_sets = gene_sets,
#'   method = "UCell",
#'   prefix = "Pathway"
#' )
#' }
NULL

#' @rdname sce_run_scenic
#' @examples
#' \dontrun{
#' sce <- readRDS("annotated_seurat.rds")
#' scenic <- sce_run_scenic(
#'   sce,
#'   db_dir = "cisTarget_databases",
#'   species = "human",
#'   celltype = "celltype",
#'   output_dir = "scenic_output",
#'   n_cores = 8
#' )
#' }
NULL

#' @rdname sce_run_two_group_deg_enrichment
#' @examples
#' \dontrun{
#' sce <- readRDS("annotated_seurat.rds")
#' result <- sce_run_two_group_deg_enrichment(
#'   input_file = sce,
#'   group_col = "response",
#'   ident_1 = "Responder",
#'   ident_2 = "Non-responder",
#'   celltype_col = "celltype",
#'   gmt_dir = "genesets",
#'   output_dir = "DEG_enrichment"
#' )
#' result$status
#' }
NULL

#' @rdname sce_run_scTenifoldKnk_KO
#' @examples
#' \dontrun{
#' sce <- readRDS("annotated_seurat.rds")
#' ko <- sce_run_scTenifoldKnk_KO(
#'   object = sce,
#'   gKO = c("RAMP1", "GABBR1"),
#'   subset_col = "celltype",
#'   subset_values = "Macrophage",
#'   output_dir = "scTenifoldKnk_KO",
#'   nCores = 4
#' )
#' }
NULL

#' @rdname sce_plot_spatial_continuous
#' @examples
#' \dontrun{
#' spatial <- readRDS("spatial_seurat.rds")
#' p <- sce_plot_spatial_continuous(
#'   spatial,
#'   feature = "GABBR1",
#'   image_scale = "lowres",
#'   colors = "magma",
#'   show_he = TRUE,
#'   show_spots = TRUE
#' )
#' p
#' }
NULL

#' @rdname sce_plot_pathway_rds_continuous
#' @examples
#' \dontrun{
#' spatial <- readRDS("spatial_seurat.rds")
#' result <- sce_plot_pathway_rds_continuous(
#'   spatial,
#'   score_dir = "gmt_pathway_scores",
#'   selected_pathways = c("GABAergic", "Neuropeptide_Pain"),
#'   output_dir = "Spatial_pathway_continuous"
#' )
#' result$summary
#' }
NULL
