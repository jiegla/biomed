#' Compare spatial gene expression across sample groups
#'
#' Spatial transcriptomics workflow adapted from contributed analysis code.
#' Sample-level summaries are used for group comparisons. Spot-level correlations
#' and Moran statistics describe spatial association and do not establish causality.
#' @param vis A spatial Seurat object; niche-score analysis also accepts a spot-indexed metadata data frame.
#' @param genes Character vector of gene names to compare.
#' @param statistic_group Optional sample-level clinical grouping columns, each constant within a sample.
#' @param sample_col Metadata column identifying independent samples. Leading-zero character IDs are preserved.
#' @param assay Expression assay. NULL uses the default assay.
#' @param layer Expression layer; split Seurat v5 layers are joined on a local copy and aligned by spot names.
#' @param sample_aggregation Use the median or mean expression per sample and region for comparisons.
#' @param spot_positive_cutoff Expression threshold above which a spot is positive.
#' @param Drop_negative Exclude values at or below spot_positive_cutoff from expression summaries; detection percentages still use all finite observations.
#' @param group_order Optional group order, or a named list of orders for multiple clinical columns.
#' @param split_by_region Optional metadata column defining spatial regions; NULL pools all spots within each sample.
#' @param run_spatial_plots Attempt tissue-backed spatial plots, with coordinate-plot fallback.
#' @param spatial_plot_keep_scale Scale mode for spatial expression panels: feature, all or sample.
#' @param spatial_pt_size Spatial point-size factor.
#' @param spatial_alpha Spatial expression point alpha range.
#' @param spatial_min_cutoff Lower expression cutoff passed to Seurat spatial plotting.
#' @param spatial_max_cutoff Upper expression cutoff passed to Seurat spatial plotting.
#' @param image_alpha Histology image opacity.
#' @param image_scale Image resolution passed to Seurat spatial plotting.
#' @param image_sample_map Optional data frame with image and sample columns overriding automatic image-to-sample mapping.
#' @param run_moran Compute within-image Moran spatial autocorrelation when coordinates are available.
#' @param moran_max_spots Maximum spots per image for Moran analysis, at least three.
#' @param moran_large_sample Skip or randomly subsample images above moran_max_spots.
#' @param moran_seed Random seed; the caller RNG state is restored.
#' @param cell2location_assay Optional assay with cell types as features and spots as columns.
#' @param cell2location_layer Layer in the cell2location assay.
#' @param cell2location_metadata_cols Optional metadata columns containing cell-type abundances.
#' @param cell2location_pattern Optional regular expression for abundance metadata column names.
#' @param cell2location_features Optional subset of cell types for abundance correlations.
#' @param run_cell2location_correlation Compute per-sample gene-abundance correlations when abundance data are available.
#' @param correlation_method Spearman or Pearson correlation.
#' @param min_spots_for_correlation Minimum matched spots per sample for a correlation, at least three.
#' @param output_dir Report output directory. Neural analysis defaults to NULL (no files); gene comparison writes its Excel report by default.
#' @param file_prefix Prefix for report filenames.
#' @param export_spot_data Return and export spot-level expression data where supported by Excel row limits.
#' @param save_pdf Save PDF figures.
#' @param save_png Save PNG figures.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @return Invisibly, a list with parameters, gene coverage, image mapping, coordinates, sample summaries, clinical comparisons, Moran results, cell2location correlations and report paths.
#' @details Existing report files are overwritten in the selected output directory.
#' Cell-type fractions use the complete denominator cell-type set. Niche-score
#' global tests use sample random intercepts; pairwise tests match samples across niches.
#' Neural tests use one summary per sample; region summaries are descriptive.
#' @export
vis_compare_spatial_gene_expression_groups <- function(
    vis,
    genes,
    statistic_group = NULL,
    sample_col = "orig.ident",
    assay = NULL,
    layer = "data",
    sample_aggregation = c("median", "mean"),
    spot_positive_cutoff = 0,
    Drop_negative = FALSE,
    group_order = NULL,
    split_by_region = NULL,

    # spatial plotting
    run_spatial_plots = TRUE,
    spatial_plot_keep_scale = c("feature", "all", "sample"),
    spatial_pt_size = 1.6,
    spatial_alpha = c(0.1, 1),
    spatial_min_cutoff = NA,
    spatial_max_cutoff = "q95",
    image_alpha = 1,
    image_scale = "lowres",
    image_sample_map = NULL,

    # Moran's I
    run_moran = TRUE,
    moran_max_spots = 8000,
    moran_large_sample = c("skip", "subsample"),
    moran_seed = 1234,

    # cell2location
    cell2location_assay = NULL,
    cell2location_layer = "data",
    cell2location_metadata_cols = NULL,
    cell2location_pattern = NULL,
    cell2location_features = NULL,
    run_cell2location_correlation = TRUE,
    correlation_method = c("spearman", "pearson"),
    min_spots_for_correlation = 20,

    # output
    output_dir = "spatial_gene_comparison",
    file_prefix = "SpatialGeneExpression",
    export_spot_data = FALSE,
    save_pdf = TRUE,
    save_png = TRUE,
    width = 12,
    height = 8
) {
  # Local bindings for columns evaluated in dplyr/ggplot data masks.
  .region <- gene <- NULL


  # ============================================================
  # 0. Packages and arguments
  # ============================================================
  required_pkgs <- c(
    "Seurat", "SeuratObject", "dplyr", "tidyr",
    "ggplot2", "openxlsx", "withr"
  )

  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      "Please install missing package(s): ",
      paste(missing_pkgs, collapse = ", ")
    )
  }

  .biomed_vis_positive_integer(moran_max_spots, "moran_max_spots", minimum = 3L)
  .biomed_vis_positive_integer(min_spots_for_correlation, "min_spots_for_correlation", minimum = 3L)
  withr::local_seed(moran_seed)
  if (!inherits(vis, "Seurat")) {
    stop("vis must be a Seurat object.")
  }

  if (length(genes) == 0) {
    stop("genes cannot be empty.")
  }

  genes <- unique(as.character(genes))
  .biomed_validate_columns(genes)
  statistic_group <- unique(as.character(statistic_group))
  statistic_group <- statistic_group[!is.na(statistic_group) & nzchar(statistic_group)]

  sample_aggregation <- match.arg(sample_aggregation)
  spatial_plot_keep_scale <- match.arg(spatial_plot_keep_scale)
  moran_large_sample <- match.arg(moran_large_sample)
  correlation_method <- match.arg(correlation_method)

  if (length(spot_positive_cutoff) != 1 || !is.finite(spot_positive_cutoff)) {
    stop("spot_positive_cutoff must be one finite numeric value.")
  }

  Drop_negative <- isTRUE(Drop_negative)
  run_spatial_plots <- isTRUE(run_spatial_plots)
  run_moran <- isTRUE(run_moran)
  run_cell2location_correlation <- isTRUE(run_cell2location_correlation)
  export_spot_data <- isTRUE(export_spot_data)

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(vis)
  }

  if (!assay %in% names(vis@assays)) {
    stop("Assay '", assay, "' was not found in the Seurat object.")
  }

  # Keep tissue plots on the same assay/layer as the statistical summaries.
  SeuratObject::DefaultAssay(vis) <- assay
  if (inherits(vis[[assay]], "Assay5")) {
    available_layers <- SeuratObject::Layers(vis[[assay]])
    matching_layers <- available_layers[available_layers == layer | startsWith(available_layers, paste0(layer, "."))]
    if (length(matching_layers) > 1L) {
      vis[[assay]] <- SeuratObject::JoinLayers(vis[[assay]], layers = layer, new = layer)
    }
  }

  meta <- vis[[]]
  reserved <- c(".cell", ".sample", ".region", "gene", "raw_expression", "expression_for_summary", "detected")
  if (any(genes %in% c(reserved, statistic_group)) || any(statistic_group %in% reserved)) stop("Gene or clinical column names collide with reserved analysis columns.")

  if (!sample_col %in% colnames(meta)) {
    stop("sample_col '", sample_col, "' was not found in metadata.")
  }

  if (length(statistic_group) > 0) {
    missing_groups <- setdiff(statistic_group, colnames(meta))
    if (length(missing_groups) > 0) {
      stop(
        "The following statistic_group column(s) are missing: ",
        paste(missing_groups, collapse = ", ")
      )
    }
  }

  use_region <- !is.null(split_by_region) &&
    length(split_by_region) == 1 &&
    !is.na(split_by_region) &&
    nzchar(as.character(split_by_region))

  if (!is.null(split_by_region) && length(split_by_region) != 1) {
    stop("split_by_region must be NULL or one metadata column name.")
  }

  if (use_region && !split_by_region %in% colnames(meta)) {
    stop("split_by_region '", split_by_region, "' was not found in metadata.")
  }

  if (!is.null(group_order) &&
      length(statistic_group) > 1 &&
      !is.list(group_order)) {
    stop(
      "When statistic_group contains multiple variables, group_order must ",
      "be NULL or a named list, e.g. ",
      "list(response=c('response','non_response'))."
    )
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ============================================================
  # 1. Helpers
  # ============================================================
  sanitize_filename <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x[x == ""] <- "unnamed"
    x
  }

  mean_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else mean(x)
  }

  median_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else stats::median(x)
  }

  sd_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1) NA_real_ else stats::sd(x)
  }

  sem_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }

  save_plot <- function(p, filename_no_ext, w = width, h = height) {
    if (isTRUE(save_pdf)) {
      ggplot2::ggsave(
        filename = paste0(filename_no_ext, ".pdf"),
        plot = p,
        width = w,
        height = h,
        limitsize = FALSE
      )
    }

    if (isTRUE(save_png)) {
      ggplot2::ggsave(
        filename = paste0(filename_no_ext, ".png"),
        plot = p,
        width = w,
        height = h,
        dpi = 300,
        limitsize = FALSE
      )
    }
  }

  plot_paths <- function(stem) {
    c(if (isTRUE(save_pdf)) paste0(stem, ".pdf"),
      if (isTRUE(save_png)) paste0(stem, ".png"))
  }

  default_palette <- biomed_colors()

  get_group_levels <- function(group_var, observed_values) {
    observed_values <- as.character(observed_values)
    observed_values <- observed_values[!is.na(observed_values)]
    observed_values <- unique(observed_values)

    requested_order <- NULL

    if (is.list(group_order)) {
      if (!is.null(names(group_order)) && group_var %in% names(group_order)) {
        requested_order <- as.character(group_order[[group_var]])
      }
    } else if (!is.null(group_order) && length(statistic_group) == 1) {
      requested_order <- as.character(group_order)
    }

    if (!is.null(requested_order)) {
      return(c(
        intersect(requested_order, observed_values),
        setdiff(observed_values, requested_order)
      ))
    }

    original_var <- meta[[group_var]]

    if (is.factor(original_var)) {
      return(intersect(levels(original_var), observed_values))
    }

    observed_values
  }

  # Generic layer extractor; supports Seurat v5 LayerData and older fallback.
  get_assay_matrix <- function(object, assay_name, layer_name, features = NULL, cells = NULL) {
    mat <- .biomed_sce_expression(object, assay_name, layer_name)

    if (!is.null(features)) {
      features <- intersect(features, rownames(mat))
      mat <- mat[features, , drop = FALSE]
    }

    if (!is.null(cells)) {
      cells <- intersect(cells, colnames(mat))
      mat <- mat[, cells, drop = FALSE]
    }

    mat
  }

  # Try to obtain cells belonging to one image.
  get_image_cells <- function(object, image_name) {
    out <- tryCatch(
      SeuratObject::Cells(object[[image_name]]),
      error = function(e) NULL
    )

    if (is.null(out)) {
      out <- tryCatch(
        colnames(object[[image_name]]),
        error = function(e) NULL
      )
    }

    unique(intersect(as.character(out), colnames(object)))
  }

  # Get coordinates from one image and normalize coordinate names to x/y.
  get_image_coordinates <- function(object, image_name) {
    coords <- tryCatch(
      SeuratObject::GetTissueCoordinates(object, image = image_name),
      error = function(e) NULL
    )

    if (is.null(coords)) {
      coords <- tryCatch(
        SeuratObject::GetTissueCoordinates(object[[image_name]]),
        error = function(e) NULL
      )
    }

    if (is.null(coords) || nrow(coords) == 0) {
      return(NULL)
    }

    coords <- as.data.frame(coords)

    if (is.null(rownames(coords)) || any(rownames(coords) == "")) {
      candidate_id_cols <- intersect(
        c("cell", "barcode", "key", "name"),
        colnames(coords)
      )
      if (length(candidate_id_cols) > 0) {
        rownames(coords) <- as.character(coords[[candidate_id_cols[1]]])
      }
    }

    x_candidates <- intersect(
      c("x", "imagecol", "col", "pxl_col_in_fullres"),
      colnames(coords)
    )
    y_candidates <- intersect(
      c("y", "imagerow", "row", "pxl_row_in_fullres"),
      colnames(coords)
    )

    if (length(x_candidates) > 0 && length(y_candidates) > 0) {
      xcol <- x_candidates[1]
      ycol <- y_candidates[1]
    } else {
      numeric_cols <- colnames(coords)[
        vapply(coords, is.numeric, logical(1))
      ]
      if (length(numeric_cols) < 2) {
        return(NULL)
      }
      xcol <- numeric_cols[1]
      ycol <- numeric_cols[2]
    }

    data.frame(
      .cell = rownames(coords),
      x = as.numeric(coords[[xcol]]),
      y = as.numeric(coords[[ycol]]),
      stringsAsFactors = FALSE
    )
  }

  # Infer image -> sample mapping by cell overlap with sample_col.
  infer_image_sample_map <- function(object, sample_col_name) {
    image_names <- names(object@images)

    if (length(image_names) == 0) {
      return(data.frame(
        image = character(),
        sample = character(),
        n_cells = integer(),
        status = character(),
        stringsAsFactors = FALSE
      ))
    }

    md <- object[[]]
    rows <- vector("list", length(image_names))

    for (i in seq_along(image_names)) {
      img <- image_names[i]
      cells_img <- get_image_cells(object, img)

      if (length(cells_img) == 0) {
        rows[[i]] <- data.frame(
          image = img,
          sample = NA_character_,
          n_cells = 0L,
          status = "no_cells",
          stringsAsFactors = FALSE
        )
        next
      }

      vals <- unique(as.character(md[cells_img, sample_col_name]))
      vals <- vals[!is.na(vals)]

      if (length(vals) == 1) {
        rows[[i]] <- data.frame(
          image = img,
          sample = vals,
          n_cells = length(cells_img),
          status = "auto_unique",
          stringsAsFactors = FALSE
        )
      } else {
        tb <- sort(table(as.character(md[cells_img, sample_col_name])), decreasing = TRUE)
        sample_top <- if (length(tb) > 0) names(tb)[1] else NA_character_

        rows[[i]] <- data.frame(
          image = img,
          sample = sample_top,
          n_cells = length(cells_img),
          status = if (length(vals) > 1) "multiple_samples_top_used" else "no_sample",
          stringsAsFactors = FALSE
        )
      }
    }

    dplyr::bind_rows(rows)
  }

  make_custom_spatial_plot <- function(plot_dat, genes_now, sample_now, keep_scale_mode) {
    dd <- plot_dat |>
      dplyr::filter(.data$.sample == sample_now) |>
      dplyr::select(".cell", "x", "y", dplyr::all_of(genes_now))

    if (nrow(dd) == 0) {
      return(NULL)
    }

    long <- tidyr::pivot_longer(
      dd,
      cols = dplyr::all_of(genes_now),
      names_to = "gene",
      values_to = "expression"
    )

    p <- ggplot2::ggplot(
      long,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$expression)
    ) +
      ggplot2::geom_point(size = max(0.1, spatial_pt_size * 0.7)) +
      ggplot2::coord_fixed() +
      ggplot2::scale_y_reverse() +
      ggplot2::facet_wrap(ggplot2::vars(gene), scales = "fixed") +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::labs(
        color = "Expression",
        title = paste0(sample_now, " | ", assay, ":", layer)
      )

    if (identical(keep_scale_mode, "all")) {
      rng <- range(long$expression, finite = TRUE)
      if (all(is.finite(rng)) && diff(rng) > 0) {
        p <- p + ggplot2::scale_color_viridis_c(limits = rng)
      } else {
        p <- p + ggplot2::scale_color_viridis_c()
      }
    } else {
      p <- p + ggplot2::scale_color_viridis_c()
    }

    p
  }

  # ============================================================
  # 2. Validate genes and extract expression
  # ============================================================
  assay_features <- rownames(vis[[assay]])
  genes_found <- intersect(genes, assay_features)
  genes_missing <- setdiff(genes, genes_found)

  if (length(genes_missing) > 0) {
    warning(
      "Gene(s) not found in assay '", assay, "' and will be skipped: ",
      paste(genes_missing, collapse = ", ")
    )
  }

  if (length(genes_found) == 0) {
    stop("None of the requested genes were found in assay '", assay, "'.")
  }

  message(
    "Extracting ", length(genes_found), " gene(s) from assay='",
    assay, "', layer='", layer, "' ..."
  )

  expr_mat <- get_assay_matrix(
    object = vis,
    assay_name = assay,
    layer_name = layer,
    features = genes_found
  )

  if (is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    stop(
      "The requested assay/layer contains no expression values: assay='",
      assay, "', layer='", layer,
      "'. Choose a populated layer (for example 'counts' or a normalized layer), ",
      "or normalize the object before calling this function."
    )
  }

  common_cells <- intersect(colnames(expr_mat), rownames(meta))
  expr_mat <- expr_mat[, common_cells, drop = FALSE]
  meta_use <- meta[common_cells, , drop = FALSE]

  expr_df <- as.data.frame(t(as.matrix(expr_mat)))
  expr_df$.cell <- rownames(expr_df)

  common_meta <- data.frame(
    .cell = rownames(meta_use),
    .sample = as.character(meta_use[[sample_col]]),
    stringsAsFactors = FALSE
  )

  if (use_region) {
    common_meta$.region <- as.character(meta_use[[split_by_region]])
  } else {
    common_meta$.region <- "All_spots"
  }

  for (g in statistic_group) {
    common_meta[[g]] <- as.character(meta_use[[g]])
  }

  spot_data <- dplyr::left_join(
    expr_df,
    common_meta,
    by = ".cell"
  )

  spot_data <- spot_data |>
    dplyr::filter(!is.na(.data$.sample), !is.na(.data$.region))

  if (nrow(spot_data) == 0) {
    stop("No usable spots remained after matching expression and metadata.")
  }

  # ============================================================
  # 3. Check sample-level clinical metadata consistency
  # ============================================================
  if (length(statistic_group) > 0) {
    for (group_var in statistic_group) {
      ck <- spot_data |>
        dplyr::distinct(.data$.sample, .data[[group_var]]) |>
        dplyr::filter(!is.na(.data[[group_var]])) |>
        dplyr::group_by(.data$.sample) |>
        dplyr::summarise(
          n_group = dplyr::n_distinct(.data[[group_var]]),
          values = paste(unique(.data[[group_var]]), collapse = ";"),
          .groups = "drop"
        )

      bad <- ck |>
        dplyr::filter(.data$n_group > 1)

      if (nrow(bad) > 0) {
        stop(
          "Clinical variable '", group_var,
          "' is not constant within sample_col='", sample_col, "'.\n",
          "Problematic sample(s): ",
          paste(bad$.sample, collapse = ", ")
        )
      }
    }
  }

  # ============================================================
  # 4. Spot-level long data and per-sample pseudobulk summaries
  # ============================================================
  id_cols <- c(".cell", ".sample", ".region", statistic_group)

  long_spot <- tidyr::pivot_longer(
    spot_data,
    cols = dplyr::all_of(genes_found),
    names_to = "gene",
    values_to = "raw_expression"
  ) |>
    dplyr::filter(is.finite(.data$raw_expression)) |>
    dplyr::mutate(
      detected = .data$raw_expression > spot_positive_cutoff,
      expression_for_summary = dplyr::if_else(
        Drop_negative & .data$raw_expression <= spot_positive_cutoff,
        NA_real_,
        .data$raw_expression
      )
    )

  group_cols_sample <- c(".sample", ".region", statistic_group, "gene")

  sample_gene_summary <- long_spot |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols_sample))) |>
    dplyr::summarise(
      n_spots = dplyr::n(),
      n_positive_spots = sum(.data$detected, na.rm = TRUE),
      pct_positive_spots = mean(.data$detected, na.rm = TRUE) * 100,
      mean_expression = mean_or_na(.data$expression_for_summary),
      median_expression = median_or_na(.data$expression_for_summary),
      sd_spot_expression = sd_or_na(.data$expression_for_summary),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      expression = if (identical(sample_aggregation, "mean")) {
        .data$mean_expression
      } else {
        .data$median_expression
      }
    )

  sample_gene_summary$sample <- sample_gene_summary$.sample
  sample_gene_summary$region <- sample_gene_summary$.region

  # ============================================================
  # 5. All-sample descriptive statistics
  # ============================================================
  all_sample_summary <- sample_gene_summary |>
    dplyr::group_by(.data$.region, .data$gene) |>
    dplyr::summarise(
      n_samples = dplyr::n_distinct(.data$.sample),
      total_spots = sum(.data$n_spots, na.rm = TRUE),
      mean_sample_expression = mean_or_na(.data$expression),
      median_sample_expression = median_or_na(.data$expression),
      sd_sample_expression = sd_or_na(.data$expression),
      sem_sample_expression = sem_or_na(.data$expression),
      mean_pct_positive_spots = mean_or_na(.data$pct_positive_spots),
      median_pct_positive_spots = median_or_na(.data$pct_positive_spots),
      .groups = "drop"
    )

  # Sample x gene heatmap-like tile plot.
  heatmap_dir <- file.path(output_dir, "00_all_samples", "SampleHeatmap")
  dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

  for (region_now in unique(as.character(sample_gene_summary$.region))) {
    dd <- sample_gene_summary |>
      dplyr::filter(.data$.region == region_now)

    p_heat <- ggplot2::ggplot(
      dd,
      ggplot2::aes(
        x = .data$gene,
        y = .data$.sample,
        fill = .data$expression
      )
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        x = NULL,
        y = sample_col,
        fill = paste0(sample_aggregation, "\nexpression"),
        title = paste0("All samples | ", region_now)
      )

    save_plot(
      p_heat,
      file.path(
        heatmap_dir,
        paste0(
          sanitize_filename(file_prefix), "_",
          sanitize_filename(region_now),
          "_SampleHeatmap"
        )
      ),
      w = max(width, length(genes_found) * 0.6 + 5),
      h = max(height, length(unique(dd$.sample)) * 0.25 + 3)
    )
  }

  # ============================================================
  # 6. Clinical-group comparisons
  # ============================================================
  group_results <- list()

  if (length(statistic_group) > 0) {

    for (group_var in statistic_group) {

      message("")
      message("============================================================")
      message("Clinical grouping: ", group_var)
      message("============================================================")

      group_dir <- file.path(
        output_dir,
        paste0("clinical_", sanitize_filename(group_var))
      )
      dot_dir <- file.path(group_dir, "DotPlot")
      box_dir <- file.path(group_dir, "SampleBoxplot")

      dir.create(dot_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(box_dir, recursive = TRUE, showWarnings = FALSE)

      dd <- sample_gene_summary |>
        dplyr::filter(!is.na(.data[[group_var]]))

      group_levels <- get_group_levels(group_var, dd[[group_var]])
      dd$.group <- factor(as.character(dd[[group_var]]), levels = group_levels)

      group_summary <- dd |>
        dplyr::group_by(.data$.region, .data$gene, .data$.group) |>
        dplyr::summarise(
          n_samples = dplyr::n_distinct(.data$.sample),
          total_spots = sum(.data$n_spots, na.rm = TRUE),
          mean_expression = mean_or_na(.data$expression),
          median_expression = median_or_na(.data$expression),
          sd_expression = sd_or_na(.data$expression),
          sem_expression = sem_or_na(.data$expression),
          mean_pct_positive_spots = mean_or_na(.data$pct_positive_spots),
          median_pct_positive_spots = median_or_na(.data$pct_positive_spots),
          .groups = "drop"
        )

      # --------------------------------------------------------
      # 6.1 Tests
      # --------------------------------------------------------
      overall_rows <- list()
      pair_rows <- list()
      oi <- 1L
      pi <- 1L

      strata <- unique(dd[, c(".region", "gene"), drop = FALSE])

      if (nrow(strata) > 0) {
        for (ii in seq_len(nrow(strata))) {
          region_now <- strata$.region[ii]
          gene_now <- strata$gene[ii]

          d0 <- dd |>
            dplyr::filter(
              .data$.region == region_now,
              .data$gene == gene_now,
              is.finite(.data$expression)
            )

          present_groups <- group_levels[
            group_levels %in% as.character(unique(d0$.group))
          ]

          if (length(present_groups) == 2) {
            g1 <- present_groups[1]
            g2 <- present_groups[2]

            x <- d0$expression[as.character(d0$.group) == g1]
            y <- d0$expression[as.character(d0$.group) == g2]

            tst <- tryCatch(
              suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)),
              error = function(e) NULL
            )

            overall_rows[[oi]] <- data.frame(
              region = region_now,
              gene = gene_now,
              n_groups = 2L,
              test = "Wilcoxon rank-sum",
              statistic = if (is.null(tst)) NA_real_ else unname(tst$statistic),
              p_value = if (is.null(tst)) NA_real_ else tst$p.value,
              stringsAsFactors = FALSE
            )
            oi <- oi + 1L

          } else if (length(present_groups) >= 3) {
            dtest <- d0 |>
              dplyr::filter(as.character(.data$.group) %in% present_groups)
            dtest$.group <- droplevels(dtest$.group)

            tst <- tryCatch(
              stats::kruskal.test(expression ~ .group, data = dtest),
              error = function(e) NULL
            )

            overall_rows[[oi]] <- data.frame(
              region = region_now,
              gene = gene_now,
              n_groups = length(present_groups),
              test = "Kruskal-Wallis",
              statistic = if (is.null(tst)) NA_real_ else unname(tst$statistic),
              p_value = if (is.null(tst)) NA_real_ else tst$p.value,
              stringsAsFactors = FALSE
            )
            oi <- oi + 1L
          }

          if (length(present_groups) >= 2) {
            pairs <- utils::combn(present_groups, 2, simplify = FALSE)

            for (pair in pairs) {
              g1 <- pair[1]
              g2 <- pair[2]

              x <- d0$expression[as.character(d0$.group) == g1]
              y <- d0$expression[as.character(d0$.group) == g2]

              tst <- tryCatch(
                suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)),
                error = function(e) NULL
              )

              pair_rows[[pi]] <- data.frame(
                region = region_now,
                gene = gene_now,
                group1 = g1,
                group2 = g2,
                contrast = paste0(g2, " vs ", g1),
                n_group1 = length(x),
                n_group2 = length(y),
                mean_group1 = mean_or_na(x),
                mean_group2 = mean_or_na(y),
                median_group1 = median_or_na(x),
                median_group2 = median_or_na(y),
                delta_mean = mean_or_na(y) - mean_or_na(x),
                delta_median = median_or_na(y) - median_or_na(x),
                statistic = if (is.null(tst)) NA_real_ else unname(tst$statistic),
                p_value = if (is.null(tst)) NA_real_ else tst$p.value,
                stringsAsFactors = FALSE
              )
              pi <- pi + 1L
            }
          }
        }
      }

      if (length(overall_rows) > 0) {
        overall_test <- dplyr::bind_rows(overall_rows)
        overall_test$p_adj_BH_global <- stats::p.adjust(
          overall_test$p_value,
          method = "BH"
        )
      } else {
        overall_test <- data.frame()
      }

      if (length(pair_rows) > 0) {
        pairwise_test <- dplyr::bind_rows(pair_rows) |>
          dplyr::group_by(.data$region, .data$gene) |>
          dplyr::mutate(
            p_adj_BH_within_gene_region = stats::p.adjust(
              .data$p_value,
              method = "BH"
            )
          ) |>
          dplyr::ungroup()

        pairwise_test$p_adj_BH_global <- stats::p.adjust(
          pairwise_test$p_value,
          method = "BH"
        )
      } else {
        pairwise_test <- data.frame()
      }

      # --------------------------------------------------------
      # 6.2 DotPlot
      # size = positive-spot percentage
      # fill = sample-level expression
      # --------------------------------------------------------
      p_dot <- ggplot2::ggplot(
        group_summary,
        ggplot2::aes(
          x = .data$gene,
          y = .data$.group,
          size = .data$mean_pct_positive_spots,
          fill = .data$mean_expression
        )
      ) +
        ggplot2::geom_point(shape = 21, na.rm = TRUE) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(.region),
          scales = "free_y",
          space = "free_y"
        ) +
        ggplot2::scale_size(range = c(1, 9), limits = c(0, 100)) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          strip.text.y = ggplot2::element_text(angle = 0)
        ) +
        ggplot2::labs(
          x = NULL,
          y = group_var,
          size = "% positive\nspots",
          fill = paste0("Mean sample\n", sample_aggregation),
          title = paste0(group_var, " | sample-level spatial expression")
        )

      save_plot(
        p_dot,
        file.path(
          dot_dir,
          paste0(
            sanitize_filename(file_prefix), "_",
            sanitize_filename(group_var),
            "_DotPlot"
          )
        ),
        w = max(width, length(genes_found) * 0.65 + 5),
        h = max(height, length(unique(group_summary$.region)) * 2.0 + 3)
      )

      # --------------------------------------------------------
      # 6.3 Sample-level box/violin plots
      # --------------------------------------------------------
      region_plots <- list()

      for (region_now in unique(as.character(dd$.region))) {
        dplot <- dd |>
          dplyr::filter(
            .data$.region == region_now,
            is.finite(.data$expression)
          )

        if (nrow(dplot) == 0) next

        dplot$.group <- factor(dplot$.group, levels = group_levels)

        p_box <- ggplot2::ggplot(
          dplot,
          ggplot2::aes(
            x = .data$.group,
            y = .data$expression,
            fill = .data$.group
          )
        ) +
          ggplot2::geom_violin(
            trim = FALSE,
            scale = "width",
            alpha = 0.35,
            na.rm = TRUE
          ) +
          ggplot2::geom_boxplot(
            width = 0.15,
            outlier.shape = NA,
            alpha = 0.65,
            na.rm = TRUE
          ) +
          ggplot2::geom_jitter(
            width = 0.12,
            height = 0,
            size = 1.6,
            alpha = 0.8
          ) +
          ggplot2::facet_wrap(
            ggplot2::vars(gene),
            scales = "free_y"
          ) +
          ggplot2::theme_bw(base_size = 12) +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "none"
          ) +
          ggplot2::labs(
            x = group_var,
            y = paste0("Sample ", sample_aggregation, " expression"),
            title = paste0(region_now, " | ", group_var)
          )

        if (length(group_levels) <= length(default_palette)) {
          p_box <- p_box +
            ggplot2::scale_fill_manual(
              values = default_palette[seq_along(group_levels)],
              drop = FALSE
            )
        }

        region_plots[[region_now]] <- p_box

        save_plot(
          p_box,
          file.path(
            box_dir,
            paste0(
              sanitize_filename(file_prefix), "_",
              sanitize_filename(group_var), "_",
              sanitize_filename(region_now),
              "_SampleBoxplot"
            )
          ),
          w = width,
          h = height
        )
      }

      group_results[[group_var]] <- list(
        group_levels = group_levels,
        summary = group_summary,
        overall_test = overall_test,
        pairwise_test = pairwise_test,
        plots = list(
          dotplot = p_dot,
          sample_boxplots = region_plots
        )
      )
    }
  }

  # ============================================================
  # 7. Image <-> sample map and coordinates
  # ============================================================
  image_map <- infer_image_sample_map(vis, sample_col)

  if (!is.null(image_sample_map)) {
    # Accept named character vector: c(image1="sample1", image2="sample2")
    if (is.vector(image_sample_map) && !is.null(names(image_sample_map))) {
      manual_map <- data.frame(
        image = names(image_sample_map),
        sample = as.character(image_sample_map),
        stringsAsFactors = FALSE
      )
    } else if (
      is.data.frame(image_sample_map) &&
      all(c("image", "sample") %in% colnames(image_sample_map))
    ) {
      manual_map <- image_sample_map[, c("image", "sample"), drop = FALSE]
    } else {
      stop(
        "image_sample_map must be a named vector c(image='sample') ",
        "or a data.frame with columns image and sample."
      )
    }

    if (nrow(image_map) == 0) {
      image_map <- manual_map |>
        dplyr::mutate(
          n_cells = NA_integer_,
          status = "manual"
        )
    } else {
      manual_map2 <- manual_map |>
        dplyr::rename(sample_manual = "sample")

      image_map <- image_map |>
        dplyr::rename(sample_auto = "sample") |>
        dplyr::left_join(manual_map2, by = "image") |>
        dplyr::mutate(
          sample = dplyr::coalesce(
            as.character(.data$sample_manual),
            as.character(.data$sample_auto)
          ),
          status = dplyr::if_else(
            !is.na(.data$sample_manual),
            "manual",
            .data$status
          )
        ) |>
        dplyr::select(
          .data$image, .data$sample, .data$n_cells, .data$status
        )
    }
  }

  coord_rows <- list()
  ci <- 1L

  if (nrow(image_map) > 0) {
    for (ii in seq_len(nrow(image_map))) {
      img <- image_map$image[ii]
      sample_now <- image_map$sample[ii]

      if (is.na(sample_now) || !nzchar(sample_now)) next

      coords <- get_image_coordinates(vis, img)
      if (is.null(coords) || nrow(coords) == 0) next

      coords$.image <- img
      coords$.sample <- sample_now
      coord_rows[[ci]] <- coords
      ci <- ci + 1L
    }
  }

  coordinates <- if (length(coord_rows) > 0) {
    dplyr::bind_rows(coord_rows) |>
      dplyr::distinct(.data$.cell, .keep_all = TRUE)
  } else {
    data.frame()
  }

  # ============================================================
  # 8. Spatial expression plots per sample
  # ============================================================
  spatial_plot_files <- character()
  spatial_plot_dir <- file.path(output_dir, "SpatialFeaturePlot")
  dir.create(spatial_plot_dir, recursive = TRUE, showWarnings = FALSE)

  if (run_spatial_plots) {

    if (nrow(image_map) == 0) {
      warning(
        "No spatial images were found. Spatial plots were skipped, ",
        "but statistical summaries are still available."
      )
    } else {

      standard_slots <- c("counts", "data", "scale.data")
      can_use_seurat_spatialplot <- layer %in% standard_slots

      for (ii in seq_len(nrow(image_map))) {
        img <- image_map$image[ii]
        sample_now <- image_map$sample[ii]

        if (is.na(sample_now) || !nzchar(sample_now)) next

        cells_sample <- rownames(meta)[
          as.character(meta[[sample_col]]) == sample_now
        ]
        cells_img <- get_image_cells(vis, img)
        cells_use <- intersect(cells_sample, cells_img)

        if (length(cells_use) == 0) next

        p_sp <- NULL

        if (can_use_seurat_spatialplot) {
          p_sp <- tryCatch(
            Seurat::SpatialFeaturePlot(
              object = vis,
              features = genes_found,
              images = img,
              slot = layer,
              keep.scale = if (identical(spatial_plot_keep_scale, "sample")) {
                NULL
              } else {
                spatial_plot_keep_scale
              },
              min.cutoff = spatial_min_cutoff,
              max.cutoff = spatial_max_cutoff,
              ncol = min(4, length(genes_found)),
              pt.size.factor = spatial_pt_size,
              alpha = spatial_alpha,
              image.alpha = image_alpha,
              image.scale = image_scale,
              combine = TRUE
            ),
            error = function(e) {
              warning(
                "Seurat::SpatialFeaturePlot failed for image '", img,
                "'. Falling back to coordinate plot. Error: ",
                conditionMessage(e)
              )
              NULL
            }
          )
        }

        # Fallback for custom layers (e.g. scvi_normalized) or plot failure.
        if (is.null(p_sp) && nrow(coordinates) > 0) {
          coords_now <- coordinates |>
            dplyr::filter(.data$.image == img)

          expr_now <- expr_df |>
            dplyr::filter(.data$.cell %in% coords_now$.cell)

          custom_dat <- dplyr::left_join(
            coords_now,
            expr_now,
            by = ".cell"
          )

          p_sp <- make_custom_spatial_plot(
            plot_dat = custom_dat,
            genes_now = genes_found,
            sample_now = sample_now,
            keep_scale_mode = spatial_plot_keep_scale
          )
        }

        if (!is.null(p_sp)) {
          basename <- file.path(
            spatial_plot_dir,
            paste0(
              sanitize_filename(file_prefix), "_",
              sanitize_filename(sample_now), "_",
              sanitize_filename(img),
              "_SpatialFeaturePlot"
            )
          )

          save_plot(
            p_sp,
            basename,
            w = max(width, min(4, length(genes_found)) * 3.2),
            h = max(height, ceiling(length(genes_found) / 4) * 3.2)
          )

          spatial_plot_files <- c(spatial_plot_files, plot_paths(basename))
        }
      }
    }
  }

  # Optional region/domain spatial plots.
  region_plot_files <- character()
  region_plot_dir <- file.path(output_dir, "SpatialRegionPlot")
  dir.create(region_plot_dir, recursive = TRUE, showWarnings = FALSE)

  if (run_spatial_plots && use_region && nrow(image_map) > 0) {
    for (ii in seq_len(nrow(image_map))) {
      img <- image_map$image[ii]
      sample_now <- image_map$sample[ii]

      if (is.na(sample_now) || !nzchar(sample_now)) next

      p_region <- tryCatch(
        Seurat::SpatialDimPlot(
          object = vis,
          group.by = split_by_region,
          images = img,
          cols = if (length(unique(stats::na.omit(meta[[split_by_region]]))) <=
                     length(default_palette)) default_palette else NULL,
          pt.size.factor = spatial_pt_size,
          image.alpha = image_alpha,
          image.scale = image_scale,
          combine = TRUE
        ),
        error = function(e) NULL
      )

      if (!is.null(p_region)) {
        basename <- file.path(
          region_plot_dir,
          paste0(
            sanitize_filename(file_prefix), "_",
            sanitize_filename(sample_now), "_",
            sanitize_filename(img),
            "_", sanitize_filename(split_by_region),
            "_SpatialRegion"
          )
        )

        save_plot(p_region, basename, w = width, h = height)
        region_plot_files <- c(region_plot_files, plot_paths(basename))
      }
    }
  }

  # ============================================================
  # 9. Moran's I per sample x gene
  # ============================================================
  moran_results <- data.frame()

  if (run_moran) {
    moran_dir <- file.path(output_dir, "MoranI")
    dir.create(moran_dir, recursive = TRUE, showWarnings = FALSE)

    if (nrow(coordinates) == 0) {
      warning("Spatial coordinates were unavailable; Moran's I was skipped.")
    } else {

      moran_rows <- list()
      mi <- 1L

      sample_image_pairs <- image_map |>
        dplyr::filter(!is.na(.data$sample))

      for (ii in seq_len(nrow(sample_image_pairs))) {
        sample_now <- sample_image_pairs$sample[ii]
        img <- sample_image_pairs$image[ii]

        coords_now <- coordinates |>
          dplyr::filter(.data$.image == img)

        cells_now <- intersect(coords_now$.cell, colnames(expr_mat))

        if (length(cells_now) < 3) next

        if (length(cells_now) > moran_max_spots) {
          if (identical(moran_large_sample, "skip")) {
            warning(
              "Skipping Moran's I for sample '", sample_now,
              "' because n_spots=", length(cells_now),
              " > moran_max_spots=", moran_max_spots, "."
            )
            next
          } else {
            cells_now <- sample(
              cells_now,
              size = moran_max_spots,
              replace = FALSE
            )
          }
        }

        coords_now <- coords_now[
          match(cells_now, coords_now$.cell),
          ,
          drop = FALSE
        ]

        pos <- as.matrix(coords_now[, c("x", "y"), drop = FALSE])
        rownames(pos) <- coords_now$.cell

        dat <- as.matrix(expr_mat[genes_found, cells_now, drop = FALSE])

        # Remove genes with zero variance in this sample.
        keep_gene <- apply(dat, 1, function(x) {
          x <- x[is.finite(x)]
          length(x) >= 3 && stats::sd(x) > 0
        })

        dat <- dat[keep_gene, , drop = FALSE]
        if (nrow(dat) == 0) next

        mr <- tryCatch(
          Seurat::RunMoransI(
            data = dat,
            pos = pos,
            verbose = FALSE
          ),
          error = function(e) {
            warning(
              "RunMoransI failed for sample '", sample_now,
              "', image '", img, "': ",
              conditionMessage(e)
            )
            NULL
          }
        )

        if (is.null(mr)) next

        mr <- as.data.frame(mr)

        # Seurat returns rownames as features; preserve generically.
        if (is.null(rownames(mr))) {
          mr$gene <- rownames(dat)[seq_len(nrow(mr))]
        } else {
          mr$gene <- rownames(mr)
        }

        mr$sample <- sample_now
        mr$image <- img
        mr$n_spots_moran <- length(cells_now)

        moran_rows[[mi]] <- mr
        mi <- mi + 1L
      }

      if (length(moran_rows) > 0) {
        moran_results <- dplyr::bind_rows(moran_rows)

        pcols <- grep("p.*value|p_value|pval", colnames(moran_results),
                      ignore.case = TRUE, value = TRUE)

        if (length(pcols) > 0) {
          pcol <- pcols[1]
          moran_results$p_adj_BH_global <- stats::p.adjust(
            moran_results[[pcol]],
            method = "BH"
          )
        }
      }
    }
  }

  # ============================================================
  # 10. cell2location abundance extraction
  # ============================================================
  abundance_data <- NULL
  abundance_source <- NULL

  # Option A: metadata columns.
  if (!is.null(cell2location_metadata_cols) || !is.null(cell2location_pattern)) {

    abundance_cols <- character()

    if (!is.null(cell2location_metadata_cols)) {
      abundance_cols <- unique(as.character(cell2location_metadata_cols))
    }

    if (!is.null(cell2location_pattern)) {
      abundance_cols <- unique(c(
        abundance_cols,
        grep(
          cell2location_pattern,
          colnames(meta_use),
          value = TRUE
        )
      ))
    }

    missing_abundance_cols <- setdiff(abundance_cols, colnames(meta_use))
    if (length(missing_abundance_cols) > 0) {
      warning(
        "cell2location metadata column(s) not found and skipped: ",
        paste(missing_abundance_cols, collapse = ", ")
      )
    }

    abundance_cols <- intersect(abundance_cols, colnames(meta_use))

    if (!is.null(cell2location_features)) {
      abundance_cols <- intersect(
        abundance_cols,
        as.character(cell2location_features)
      )
    }

    if (length(abundance_cols) > 0) {
      abundance_data <- meta_use[, abundance_cols, drop = FALSE]
      abundance_data$.cell <- rownames(abundance_data)
      abundance_source <- "metadata"
    }
  }

  # Option B: assay where features are cell types and columns are spots.
  if (is.null(abundance_data) && !is.null(cell2location_assay)) {
    if (!cell2location_assay %in% names(vis@assays)) {
      warning(
        "cell2location_assay '", cell2location_assay,
        "' was not found; cell2location correlation will be skipped."
      )
    } else {
      abundance_mat <- get_assay_matrix(
        object = vis,
        assay_name = cell2location_assay,
        layer_name = cell2location_layer,
        cells = common_cells
      )

      abundance_features_found <- rownames(abundance_mat)

      if (!is.null(cell2location_features)) {
        abundance_features_found <- intersect(
          as.character(cell2location_features),
          abundance_features_found
        )
        abundance_mat <- abundance_mat[
          abundance_features_found,
          ,
          drop = FALSE
        ]
      }

      if (nrow(abundance_mat) > 0) {
        abundance_data <- as.data.frame(t(as.matrix(abundance_mat)))
        abundance_data$.cell <- rownames(abundance_data)
        abundance_source <- paste0(
          "assay:", cell2location_assay, ":", cell2location_layer
        )
      }
    }
  }

  # ============================================================
  # 11. Gene x cell2location correlations within each sample
  # ============================================================
  cell2location_correlations <- data.frame()
  cell2location_group_tests <- list()
  c2l_plot_files <- character()

  if (run_cell2location_correlation && !is.null(abundance_data)) {

    message("")
    message("Running gene x cell2location abundance correlations ...")

    abundance_features <- setdiff(colnames(abundance_data), ".cell")

    cdat <- expr_df |>
      dplyr::left_join(abundance_data, by = ".cell") |>
      dplyr::left_join(
        common_meta[, c(".cell", ".sample", statistic_group), drop = FALSE],
        by = ".cell"
      )

    cor_rows <- list()
    cri <- 1L

    samples_present <- unique(as.character(cdat$.sample))
    samples_present <- samples_present[!is.na(samples_present)]

    for (sample_now in samples_present) {
      ds <- cdat |>
        dplyr::filter(.data$.sample == sample_now)

      for (gene_now in genes_found) {
        gx <- ds[[gene_now]]

        for (celltype_now in abundance_features) {
          ay <- suppressWarnings(as.numeric(as.character(ds[[celltype_now]])))

          keep <- is.finite(gx) & is.finite(ay)
          n_used <- sum(keep)

          if (n_used < min_spots_for_correlation) next
          if (stats::sd(gx[keep]) == 0 || stats::sd(ay[keep]) == 0) next

          ct <- tryCatch(
            suppressWarnings(
              stats::cor.test(
                gx[keep],
                ay[keep],
                method = correlation_method,
                exact = FALSE
              )
            ),
            error = function(e) NULL
          )

          if (is.null(ct)) next

          row_now <- data.frame(
            sample = sample_now,
            gene = gene_now,
            celltype = celltype_now,
            n_spots = n_used,
            method = correlation_method,
            correlation = unname(ct$estimate),
            p_value = ct$p.value,
            stringsAsFactors = FALSE
          )

          for (group_var in statistic_group) {
            vals <- unique(as.character(ds[[group_var]]))
            vals <- vals[!is.na(vals)]
            row_now[[group_var]] <- if (length(vals) == 1) vals else NA_character_
          }

          cor_rows[[cri]] <- row_now
          cri <- cri + 1L
        }
      }
    }

    if (length(cor_rows) > 0) {
      cell2location_correlations <- dplyr::bind_rows(cor_rows)
      cell2location_correlations$p_adj_BH_global <- stats::p.adjust(
        cell2location_correlations$p_value,
        method = "BH"
      )

      # --------------------------------------------------------
      # 11.1 Correlation heatmaps: sample x gene/celltype
      # --------------------------------------------------------
      c2l_dir <- file.path(output_dir, "Cell2locationCorrelation")
      dir.create(c2l_dir, recursive = TRUE, showWarnings = FALSE)

      for (gene_now in unique(cell2location_correlations$gene)) {
        dh <- cell2location_correlations |>
          dplyr::filter(.data$gene == gene_now)

        p_cor <- ggplot2::ggplot(
          dh,
          ggplot2::aes(
            x = .data$celltype,
            y = .data$sample,
            fill = .data$correlation
          )
        ) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient2(
            midpoint = 0,
            limits = c(-1, 1),
            na.value = "grey90"
          ) +
          ggplot2::theme_bw(base_size = 10) +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
          ) +
          ggplot2::labs(
            x = "cell2location cell type",
            y = sample_col,
            fill = paste0(correlation_method, "\nr"),
            title = paste0(
              gene_now,
              " vs cell2location abundance"
            )
          )

        basename <- file.path(
          c2l_dir,
          paste0(
            sanitize_filename(file_prefix), "_",
            sanitize_filename(gene_now),
            "_Cell2locationCorrelation"
          )
        )

        save_plot(
          p_cor,
          basename,
          w = max(width, length(unique(dh$celltype)) * 0.35 + 5),
          h = max(height, length(unique(dh$sample)) * 0.25 + 3)
        )

        c2l_plot_files <- c(c2l_plot_files, plot_paths(basename))
      }

      # --------------------------------------------------------
      # 11.2 Compare per-sample correlations across clinical groups
      # --------------------------------------------------------
      if (length(statistic_group) > 0) {

        for (group_var in statistic_group) {

          dg <- cell2location_correlations |>
            dplyr::filter(!is.na(.data[[group_var]]))

          test_rows <- list()
          ti <- 1L

          strata <- unique(dg[, c("gene", "celltype"), drop = FALSE])

          for (ii in seq_len(nrow(strata))) {
            gene_now <- strata$gene[ii]
            celltype_now <- strata$celltype[ii]

            d0 <- dg |>
              dplyr::filter(
                .data$gene == gene_now,
                .data$celltype == celltype_now,
                is.finite(.data$correlation)
              )

            lev <- get_group_levels(group_var, d0[[group_var]])
            lev <- lev[lev %in% unique(as.character(d0[[group_var]]))]

            if (length(lev) == 2) {
              x <- d0$correlation[as.character(d0[[group_var]]) == lev[1]]
              y <- d0$correlation[as.character(d0[[group_var]]) == lev[2]]

              tst <- tryCatch(
                suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)),
                error = function(e) NULL
              )

              test_rows[[ti]] <- data.frame(
                gene = gene_now,
                celltype = celltype_now,
                group1 = lev[1],
                group2 = lev[2],
                n_group1 = length(x),
                n_group2 = length(y),
                median_r_group1 = median_or_na(x),
                median_r_group2 = median_or_na(y),
                delta_median_r = median_or_na(y) - median_or_na(x),
                test = "Wilcoxon rank-sum",
                p_value = if (is.null(tst)) NA_real_ else tst$p.value,
                stringsAsFactors = FALSE
              )
              ti <- ti + 1L

            } else if (length(lev) >= 3) {
              dtest <- d0
              dtest$.tmp_group <- factor(
                as.character(dtest[[group_var]]),
                levels = lev
              )

              tst <- tryCatch(
                stats::kruskal.test(correlation ~ .tmp_group, data = dtest),
                error = function(e) NULL
              )

              test_rows[[ti]] <- data.frame(
                gene = gene_now,
                celltype = celltype_now,
                group1 = NA_character_,
                group2 = NA_character_,
                n_group1 = NA_integer_,
                n_group2 = NA_integer_,
                median_r_group1 = NA_real_,
                median_r_group2 = NA_real_,
                delta_median_r = NA_real_,
                test = "Kruskal-Wallis",
                p_value = if (is.null(tst)) NA_real_ else tst$p.value,
                stringsAsFactors = FALSE
              )
              ti <- ti + 1L
            }
          }

          if (length(test_rows) > 0) {
            tt <- dplyr::bind_rows(test_rows)
            tt$p_adj_BH_global <- stats::p.adjust(tt$p_value, method = "BH")
            cell2location_group_tests[[group_var]] <- tt
          }
        }
      }
    }
  }

  # ============================================================
  # 12. Excel workbook
  # ============================================================
  xlsx_file <- file.path(
    output_dir,
    paste0(sanitize_filename(file_prefix), "_results.xlsx")
  )

  parameter_table <- data.frame(
    parameter = c(
      "genes_requested",
      "genes_found",
      "genes_missing",
      "statistic_group",
      "sample_col",
      "assay",
      "layer",
      "sample_aggregation",
      "spot_positive_cutoff",
      "Drop_negative",
      "split_by_region",
      "run_spatial_plots",
      "run_moran",
      "moran_max_spots",
      "cell2location_source",
      "correlation_method",
      "min_spots_for_correlation"
    ),
    value = c(
      paste(genes, collapse = ";"),
      paste(genes_found, collapse = ";"),
      paste(genes_missing, collapse = ";"),
      if (length(statistic_group) == 0) "NULL" else paste(statistic_group, collapse = ";"),
      sample_col,
      assay,
      layer,
      sample_aggregation,
      as.character(spot_positive_cutoff),
      as.character(Drop_negative),
      if (use_region) split_by_region else "NULL",
      as.character(run_spatial_plots),
      as.character(run_moran),
      as.character(moran_max_spots),
      if (is.null(abundance_source)) "NULL" else abundance_source,
      correlation_method,
      as.character(min_spots_for_correlation)
    ),
    stringsAsFactors = FALSE
  )

  wb <- openxlsx::createWorkbook()

  add_sheet <- function(wb, name, dat) {
    name <- substr(name, 1, 31)
    existing <- openxlsx::sheets(wb)

    if (name %in% existing) {
      base <- substr(name, 1, 27)
      k <- 2L
      repeat {
        candidate <- paste0(base, "_", k)
        if (!candidate %in% existing) {
          name <- candidate
          break
        }
        k <- k + 1L
      }
    }

    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, dat)
    openxlsx::freezePane(wb, sheet = name, firstRow = TRUE)
    name
  }

  add_sheet(wb, "Parameters", parameter_table)
  add_sheet(wb, "Image_Sample_Map", image_map)
  add_sheet(wb, "AllSample_Summary", all_sample_summary)
  add_sheet(wb, "Sample_Gene_Summary", sample_gene_summary)

  if (nrow(moran_results) > 0) {
    add_sheet(wb, "MoranI", moran_results)
  }

  if (nrow(cell2location_correlations) > 0) {
    add_sheet(wb, "C2L_Correlations", cell2location_correlations)
  }

  if (length(group_results) > 0) {
    for (group_var in names(group_results)) {
      gr <- group_results[[group_var]]

      add_sheet(
        wb,
        paste0("GrpSum_", sanitize_filename(group_var)),
        gr$summary
      )

      if (nrow(gr$overall_test) > 0) {
        add_sheet(
          wb,
          paste0("Overall_", sanitize_filename(group_var)),
          gr$overall_test
        )
      }

      if (nrow(gr$pairwise_test) > 0) {
        add_sheet(
          wb,
          paste0("Pairwise_", sanitize_filename(group_var)),
          gr$pairwise_test
        )
      }
    }
  }

  if (length(cell2location_group_tests) > 0) {
    for (group_var in names(cell2location_group_tests)) {
      add_sheet(
        wb,
        paste0("C2Ltest_", sanitize_filename(group_var)),
        cell2location_group_tests[[group_var]]
      )
    }
  }

  if (export_spot_data) {
    if (nrow(long_spot) <= 1048575) {
      add_sheet(wb, "Spot_Data", long_spot)
    } else {
      warning(
        "Spot_Data has more rows than Excel supports and was not exported."
      )
    }
  }

  openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)

  # ============================================================
  # 13. Return
  # ============================================================
  message("")
  message("============================================================")
  message("Spatial gene comparison completed.")
  message("Samples: ", dplyr::n_distinct(sample_gene_summary$.sample))
  message("Genes found: ", length(genes_found))
  message("Output: ", normalizePath(output_dir, mustWork = FALSE))
  message("Excel: ", normalizePath(xlsx_file, mustWork = FALSE))
  message("============================================================")

  invisible(list(
    parameters = parameter_table,
    genes_found = genes_found,
    genes_missing = genes_missing,
    image_sample_map = image_map,
    coordinates = coordinates,
    spot_data = if (export_spot_data) long_spot else NULL,
    sample_gene_summary = sample_gene_summary,
    all_sample_summary = all_sample_summary,
    clinical_group_results = group_results,
    moran = moran_results,
    cell2location = list(
      source = abundance_source,
      correlations = cell2location_correlations,
      clinical_group_tests = cell2location_group_tests
    ),
    files = list(
      xlsx = xlsx_file,
      output_dir = output_dir,
      spatial_feature_plots = spatial_plot_files,
      spatial_region_plots = region_plot_files,
      cell2location_correlation_plots = c2l_plot_files
    )
  ))
}
