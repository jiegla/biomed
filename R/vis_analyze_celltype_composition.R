#' Analyze spatial cell-type composition
#'
#' Spatial transcriptomics workflow adapted from contributed analysis code.
#' Sample-level summaries are used for group comparisons. Spot-level correlations
#' and Moran statistics describe spatial association and do not establish causality.
#' @param object A metadata data frame or Seurat object containing spot-level inferred cell-type abundances.
#' @param group_cols One or more spot metadata columns combined to define analysis groups.
#' @param sample_col Metadata column identifying independent samples. Leading-zero character IDs are preserved.
#' @param celltype_cols Abundance columns to analyze. NULL selects all denominator cell types.
#' @param all_celltype_cols Full abundance column set used for fraction denominators, even when analyzing a subset.
#' @param celltype_range Inclusive start/end metadata column names used to detect denominator cell types when all_celltype_cols is NULL.
#' @param metric Fraction of total inferred abundance or absolute abundance.
#' @param abundance_stat For abundance, aggregate spots using the mean, median or total.
#' @param fraction_scale Positive scale for fractions; 100 reports percentages.
#' @param group_levels Optional ordering of groups; observed groups omitted from this order are retained.
#' @param group_sep Separator used when combining multiple grouping columns.
#' @param group_colors Optional named or unnamed group colors; defaults to biomed_colors().
#' @param paired Use sample-matched paired Wilcoxon comparisons; FALSE uses rank-sum comparisons.
#' @param min_paired_samples Minimum matched samples required for a paired comparison.
#' @param min_samples_global Minimum distinct samples required for the global mixed-effects test.
#' @param p_adjust_method Multiple-testing adjustment method accepted by stats::p.adjust.
#' @param alpha Adjusted P-value threshold used to flag significant results.
#' @param na_abundance_as_zero Replace missing abundance values by zero before aggregation.
#' @param facet_ncol Maximum columns in the combined faceted plot.
#' @param individual_width Individual figure width in inches.
#' @param individual_height Individual figure height in inches.
#' @param combined_width Combined figure width in inches.
#' @param combined_panel_height Height per combined facet row in inches.
#' @param dpi Raster output resolution in dots per inch.
#' @param outdir Output directory. NULL constructs a directory from the grouping column(s); these report workflows export files by default.
#' @param save_individual Save each cell type as PDF and PNG.
#' @param save_multipage_pdf Save all panels in a multipage PDF.
#' @param save_combined Save the combined faceted figure.
#' @param export_csv Export CSV result tables.
#' @param export_xlsx Export an Excel workbook when openxlsx is available.
#' @param verbose Print progress messages.
#' @return A list with sample_group_data, group_summary, global and pairwise statistics, significant-result tables, spot counts, a combined plot and output paths.
#' @details Existing report files are overwritten in the selected output directory.
#' Cell-type fractions use the complete denominator cell-type set. Niche-score
#' global tests use sample random intercepts; pairwise tests match samples across niches.
#' Neural tests use one summary per sample; region summaries are descriptive.
#' @export
vis_analyze_celltype_composition <- function(
    object,
    group_cols,
    sample_col = "sample_id",

    # Which cell types to ANALYZE/PLOT.
    # NULL = analyze all detected denominator cell types.
    celltype_cols = NULL,

    # Columns used as the denominator for fraction.
    # NULL = automatically detect using celltype_range.
    all_celltype_cols = NULL,
    celltype_range = c("ASC", "pDC"),

    # Analysis type
    metric = c("fraction", "abundance"),

    # Only used when metric = "abundance"
    abundance_stat = c("mean", "median", "total"),

    # Fraction is returned as percentage by default
    fraction_scale = 100,

    # Optional final order of analysis groups.
    # For multiple group_cols these should be the COMBINED labels.
    group_levels = NULL,
    group_sep = " | ",

    # Optional named or unnamed color vector
    group_colors = NULL,

    # Statistics
    paired = TRUE,
    min_paired_samples = 3,
    min_samples_global = 3,
    p_adjust_method = "BH",
    alpha = 0.05,

    # Missing abundance values
    na_abundance_as_zero = TRUE,

    # Plot settings
    facet_ncol = 5,
    individual_width = 6,
    individual_height = 5,
    combined_width = 18,
    combined_panel_height = 2.5,
    dpi = 300,

    # Output
    outdir = NULL,
    save_individual = TRUE,
    save_multipage_pdf = TRUE,
    save_combined = TRUE,
    export_csv = TRUE,
    export_xlsx = TRUE,

    verbose = TRUE
) {
  # Local bindings for columns evaluated in dplyr/ggplot data masks.
  abundance <- analysis_group <- analysis_value <- celltype <- FDR_within_celltype <- global_FDR <- global_p <- group <- sample_id <- sample_id_internal <- total_abundance <- total_all_celltype_abundance <- NULL


  # ==========================================================
  # 0. Packages
  # ==========================================================

  required_pkgs <- c(
    "dplyr",
    "tidyr",
    "tibble",
    "ggplot2",
    "nlme"
  )

  missing_pkgs <- required_pkgs[
    !vapply(
      required_pkgs,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      "Please install required package(s): ",
      paste(missing_pkgs, collapse = ", ")
    )
  }

  .biomed_vis_positive_integer(min_paired_samples, "min_paired_samples")
  .biomed_vis_positive_integer(min_samples_global, "min_samples_global")
  .biomed_vis_positive_integer(facet_ncol, "facet_ncol")
  if (length(fraction_scale) != 1L || !is.finite(fraction_scale) || fraction_scale <= 0) stop("fraction_scale must be positive.")
  .biomed_probability(alpha, "alpha")
  p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)
  metric <- match.arg(metric)
  abundance_stat <- match.arg(abundance_stat)

  # ==========================================================
  # 1. Extract metadata
  # ==========================================================

  if (is.data.frame(object)) {

    md <- object

  } else {

    md <- tryCatch(
      object@meta.data,
      error = function(e) NULL
    )

    if (is.null(md)) {
      stop(
        "`object` must be a data.frame or an object containing `@meta.data`."
      )
    }
  }

  if (is.null(rownames(md))) {
    rownames(md) <- paste0("row_", seq_len(nrow(md)))
  }

  # ==========================================================
  # 2. Basic checks
  # ==========================================================

  if (missing(group_cols) || length(group_cols) < 1) {
    stop("Please provide at least one `group_cols` column.")
  }

  group_cols <- unique(as.character(group_cols))

  missing_group_cols <- setdiff(group_cols, colnames(md))

  if (length(missing_group_cols) > 0) {
    stop(
      "Cannot find group column(s): ",
      paste(missing_group_cols, collapse = ", ")
    )
  }

  if (!sample_col %in% colnames(md)) {
    stop("Cannot find sample column: ", sample_col)
  }

  # ==========================================================
  # 3. Detect all cell-type abundance columns
  #
  # all_celltype_cols:
  #   denominator for fraction
  #   and default analysis cell types
  # ==========================================================

  if (is.null(all_celltype_cols)) {

    if (length(celltype_range) != 2) {
      stop(
        "`celltype_range` must contain exactly 2 column names, ",
        "for example c('ASC', 'pDC')."
      )
    }

    start_idx <- match(celltype_range[1], colnames(md))
    end_idx   <- match(celltype_range[2], colnames(md))

    if (is.na(start_idx) || is.na(end_idx)) {

      if (metric == "abundance" && !is.null(celltype_cols)) {

        # For abundance no denominator is needed.
        all_celltype_cols <- celltype_cols

        if (verbose) {
          message(
            "Could not detect celltype_range; for abundance analysis ",
            "using `celltype_cols` directly."
          )
        }

      } else {

        stop(
          "Cannot automatically detect all cell-type columns from `celltype_range` = ",
          paste(celltype_range, collapse = " -> "),
          ".\n",
          "Please provide `all_celltype_cols` explicitly."
        )
      }

    } else {

      if (end_idx < start_idx) {
        stop(
          "`celltype_range[2]` appears before `celltype_range[1]` in metadata."
        )
      }

      all_celltype_cols <- colnames(md)[start_idx:end_idx]
    }
  }

  all_celltype_cols <- unique(as.character(all_celltype_cols))
  .biomed_validate_columns(all_celltype_cols)

  missing_all_celltypes <- setdiff(
    all_celltype_cols,
    colnames(md)
  )

  if (length(missing_all_celltypes) > 0) {
    stop(
      "Cannot find all_celltype_cols: ",
      paste(missing_all_celltypes, collapse = ", ")
    )
  }

  if (is.null(celltype_cols)) {
    celltype_cols <- all_celltype_cols
  }

  celltype_cols <- unique(as.character(celltype_cols))
  .biomed_validate_columns(celltype_cols)

  missing_celltypes <- setdiff(
    celltype_cols,
    colnames(md)
  )

  if (length(missing_celltypes) > 0) {
    stop(
      "Cannot find requested celltype_cols: ",
      paste(missing_celltypes, collapse = ", ")
    )
  }

  # For fraction, analyzed cell types should normally be in denominator.
  not_in_denominator <- setdiff(
    celltype_cols,
    all_celltype_cols
  )

  if (
    metric == "fraction" &&
    length(not_in_denominator) > 0
  ) {
    stop(
      "For fraction analysis, the following requested cell type(s) are ",
      "not included in `all_celltype_cols`: ",
      paste(not_in_denominator, collapse = ", ")
    )
  }

  if (verbose) {
    message("========================================")
    message("vis_analyze_celltype_composition")
    message("metric              = ", metric)
    message("sample_col          = ", sample_col)
    message("group_cols          = ", paste(group_cols, collapse = ", "))
    message("cell types analyzed = ", length(celltype_cols))

    if (metric == "fraction") {
      message(
        "fraction denominator cell types = ",
        length(all_celltype_cols)
      )
    } else {
      message("abundance_stat      = ", abundance_stat)
    }

    message("========================================")
  }

  # ==========================================================
  # 4. Prepare metadata
  # ==========================================================

  convert_numeric <- function(x) {

    if (is.factor(x)) {
      x <- as.character(x)
    }

    out <- .biomed_as_numeric(x, "Cell-type abundance")
    if (any(out < 0, na.rm = TRUE)) stop("Cell-type abundance must be nonnegative.")
    out
  }

  md2 <- md |>
    tibble::rownames_to_column("spot") |>
    dplyr::select(
      dplyr::all_of(
        unique(
          c(
            "spot",
            sample_col,
            group_cols,
            all_celltype_cols
          )
        )
      )
    )

  md2[[sample_col]] <- as.character(
    md2[[sample_col]]
  )

  for (g in group_cols) {
    md2[[g]] <- as.character(md2[[g]])
  }

  for (ct in all_celltype_cols) {
    md2[[ct]] <- convert_numeric(md2[[ct]])
  }

  # Remove rows with missing sample/group information
  keep <- !is.na(md2[[sample_col]]) &
    nzchar(md2[[sample_col]])

  for (g in group_cols) {
    keep <- keep &
      !is.na(md2[[g]]) &
      nzchar(md2[[g]])
  }

  md2 <- md2[keep, , drop = FALSE]

  if (nrow(md2) == 0) {
    stop("No rows remain after removing missing sample/group values.")
  }

  # Standard internal sample name
  md2$sample_id_internal <- md2[[sample_col]]

  # Create one analysis-group label
  if (length(group_cols) == 1) {

    md2$analysis_group <- md2[[group_cols]]

  } else {

    md2$analysis_group <- do.call(
      paste,
      c(
        md2[group_cols],
        sep = group_sep
      )
    )
  }

  # ==========================================================
  # 5. Group order
  # ==========================================================

  observed_groups <- unique(
    as.character(md2$analysis_group)
  )

  if (is.null(group_levels)) {

    # Preserve factor order for a single grouping column if possible
    if (
      length(group_cols) == 1 &&
      is.factor(md[[group_cols]])
    ) {

      group_levels <- levels(md[[group_cols]])

      group_levels <- group_levels[
        group_levels %in% observed_groups
      ]

    } else {

      group_levels <- observed_groups
    }

  } else {

    group_levels <- as.character(group_levels)

    absent_requested_groups <- setdiff(
      group_levels,
      observed_groups
    )

    if (
      length(absent_requested_groups) > 0 &&
      verbose
    ) {
      message(
        "Ignoring group_levels not present in data: ",
        paste(absent_requested_groups, collapse = ", ")
      )
    }

    group_levels <- group_levels[
      group_levels %in% observed_groups
    ]

    # Add any observed groups omitted from user-supplied group_levels
    omitted_observed_groups <- setdiff(
      observed_groups,
      group_levels
    )

    group_levels <- c(
      group_levels,
      omitted_observed_groups
    )
  }

  md2$analysis_group <- factor(
    md2$analysis_group,
    levels = group_levels
  )

  # ==========================================================
  # 6. Output directory
  # ==========================================================

  safe_text <- function(x) {
    gsub(
      "[^A-Za-z0-9_.-]",
      "_",
      x
    )
  }

  group_tag <- paste(
    safe_text(group_cols),
    collapse = "_"
  )

  metric_tag <- if (metric == "fraction") {
    "fraction"
  } else {
    paste0("abundance_", abundance_stat)
  }

  if (is.null(outdir)) {
    outdir <- paste0(
      "Celltype_",
      metric_tag,
      "_by_",
      group_tag
    )
  }

  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (save_individual) {
    dir.create(
      file.path(outdir, "individual_plots"),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  # ==========================================================
  # 7. Spot number
  # ==========================================================

  spot_number <- md2 |>
    dplyr::count(
      sample_id_internal,
      analysis_group,
      name = "n_spots"
    ) |>
    dplyr::rename(
      sample_id = sample_id_internal,
      group = analysis_group
    )

  # ==========================================================
  # 8. Long format
  # ==========================================================

  abundance_long <- md2 |>
    dplyr::select(
      sample_id_internal,
      analysis_group,
      dplyr::all_of(all_celltype_cols)
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(all_celltype_cols),
      names_to = "celltype",
      values_to = "abundance"
    )

  if (na_abundance_as_zero) {
    abundance_long$abundance[
      is.na(abundance_long$abundance)
    ] <- 0
  }

  # ==========================================================
  # 9. Aggregate each sample × group × cell type
  # ==========================================================

  sample_group_all_celltypes <- abundance_long |>
    dplyr::group_by(
      sample_id_internal,
      analysis_group,
      celltype
    ) |>
    dplyr::summarise(
      n_spots_nonmissing = sum(!is.na(abundance)),
      total_abundance = sum(
        abundance,
        na.rm = TRUE
      ),
      mean_abundance = if (
        all(is.na(abundance))
      ) {
        NA_real_
      } else {
        mean(
          abundance,
          na.rm = TRUE
        )
      },
      median_abundance = if (
        all(is.na(abundance))
      ) {
        NA_real_
      } else {
        stats::median(
          abundance,
          na.rm = TRUE
        )
      },
      sd_abundance = if (
        sum(!is.na(abundance)) < 2
      ) {
        NA_real_
      } else {
        stats::sd(
          abundance,
          na.rm = TRUE
        )
      },
      .groups = "drop"
    )

  # ==========================================================
  # 10. Fraction denominator
  # ==========================================================

  denominator_table <- sample_group_all_celltypes |>
    dplyr::group_by(
      sample_id_internal,
      analysis_group
    ) |>
    dplyr::summarise(
      total_all_celltype_abundance = sum(
        total_abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  # ==========================================================
  # 11. Build final sample-level analysis table
  # ==========================================================

  analysis_data <- sample_group_all_celltypes |>
    dplyr::filter(
      celltype %in% celltype_cols
    ) |>
    dplyr::left_join(
      denominator_table,
      by = c(
        "sample_id_internal",
        "analysis_group"
      )
    ) |>
    dplyr::left_join(
      spot_number |>
        dplyr::rename(
          sample_id_internal = sample_id,
          analysis_group = group
        ),
      by = c(
        "sample_id_internal",
        "analysis_group"
      )
    ) |>
    dplyr::mutate(
      fraction = ifelse(
        total_all_celltype_abundance > 0,
        fraction_scale *
          total_abundance /
          total_all_celltype_abundance,
        NA_real_
      )
    )

  if (metric == "fraction") {

    analysis_data$analysis_value <- analysis_data$fraction

    y_label <- if (fraction_scale == 100) {
      "Cell-type fraction (%)"
    } else {
      "Cell-type fraction"
    }

  } else {

    analysis_data$analysis_value <- switch(
      abundance_stat,
      mean = analysis_data$mean_abundance,
      median = analysis_data$median_abundance,
      total = analysis_data$total_abundance
    )

    y_label <- switch(
      abundance_stat,
      mean = "Mean inferred abundance per spot",
      median = "Median inferred abundance per spot",
      total = "Total inferred abundance"
    )
  }

  analysis_data <- analysis_data |>
    dplyr::filter(
      is.finite(analysis_value)
    ) |>
    dplyr::mutate(
      analysis_group = factor(
        analysis_group,
        levels = group_levels
      )
    ) |>
    dplyr::rename(
      sample_id = sample_id_internal,
      group = analysis_group
    ) |>
    dplyr::arrange(
      celltype,
      group,
      sample_id
    )

  # ==========================================================
  # 12. Descriptive group summary
  # ==========================================================

  group_summary <- analysis_data |>
    dplyr::group_by(
      celltype,
      group
    ) |>
    dplyr::summarise(
      n_samples = dplyr::n_distinct(sample_id),
      mean = mean(
        analysis_value,
        na.rm = TRUE
      ),
      sd = stats::sd(
        analysis_value,
        na.rm = TRUE
      ),
      median = stats::median(
        analysis_value,
        na.rm = TRUE
      ),
      Q1 = stats::quantile(
        analysis_value,
        0.25,
        na.rm = TRUE,
        names = FALSE
      ),
      Q3 = stats::quantile(
        analysis_value,
        0.75,
        na.rm = TRUE,
        names = FALSE
      ),
      min = min(
        analysis_value,
        na.rm = TRUE
      ),
      max = max(
        analysis_value,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  # ==========================================================
  # 13. Global mixed-effects model
  #
  # fraction:
  #   asin(sqrt(fraction)) ~ group + (1 | sample)
  #
  # abundance:
  #   log1p(abundance) ~ group + (1 | sample)
  # ==========================================================

  global_stat_list <- lapply(
    celltype_cols,
    function(ct) {

      d <- analysis_data |>
        dplyr::filter(
          celltype == ct,
          is.finite(analysis_value)
        )

      n_groups <- dplyr::n_distinct(d$group)
      n_samples <- dplyr::n_distinct(d$sample_id)

      if (
        n_groups < 2 ||
        n_samples < min_samples_global
      ) {
        return(
          data.frame(
            celltype = ct,
            n_samples = n_samples,
            n_groups = n_groups,
            global_p = NA_real_
          )
        )
      }

      if (metric == "fraction") {

        fraction01 <- d$analysis_value /
          fraction_scale

        fraction01 <- pmin(
          pmax(
            fraction01,
            0
          ),
          1
        )

        d$y <- asin(
          sqrt(fraction01)
        )

      } else {

        # log1p requires non-negative abundance
        if (any(d$analysis_value < 0, na.rm = TRUE)) {
          return(
            data.frame(
              celltype = ct,
              n_samples = n_samples,
              n_groups = n_groups,
              global_p = NA_real_
            )
          )
        }

        d$y <- log1p(
          d$analysis_value
        )
      }

      fit <- try(
        nlme::lme(
          fixed = y ~ group,
          random = ~1 | sample_id,
          data = d,
          method = "ML",
          na.action = stats::na.omit,
          control = nlme::lmeControl(
            returnObject = FALSE
          )
        ),
        silent = TRUE
      )

      if (inherits(fit, "try-error")) {
        return(
          data.frame(
            celltype = ct,
            n_samples = n_samples,
            n_groups = n_groups,
            global_p = NA_real_,
            model_message = as.character(fit)
          )
        )
      }

      a <- try(
        stats::anova(fit),
        silent = TRUE
      )

      if (inherits(a, "try-error")) {

        p_value <- NA_real_

      } else {

        p_value <- tryCatch(
          a["group", "p-value"],
          error = function(e) NA_real_
        )
      }

      data.frame(
        celltype = ct,
        n_samples = n_samples,
        n_groups = n_groups,
        global_p = as.numeric(p_value),
        model_message = if (inherits(a, "try-error")) as.character(a) else NA_character_
      )
    }
  )

  global_stats <- dplyr::bind_rows(
    global_stat_list
  )

  global_stats$global_FDR <- stats::p.adjust(
    global_stats$global_p,
    method = p_adjust_method
  )

  global_stats <- global_stats |>
    dplyr::arrange(
      global_FDR,
      global_p
    )

  # ==========================================================
  # 14. Pairwise tests
  # ==========================================================

  pairwise_template <- data.frame(
    celltype = character(0),
    group1 = character(0),
    group2 = character(0),
    paired = logical(0),
    n_group1 = integer(0),
    n_group2 = integer(0),
    n_paired_samples = integer(0),
    median_group1 = numeric(0),
    median_group2 = numeric(0),
    median_difference_group2_minus_group1 = numeric(0),
    p_value = numeric(0),
    stringsAsFactors = FALSE
  )

  pairwise_list <- list()
  counter <- 1

  for (ct in celltype_cols) {

    d <- analysis_data |>
      dplyr::filter(
        celltype == ct
      )

    groups_present <- group_levels[
      group_levels %in%
        as.character(unique(d$group))
    ]

    if (length(groups_present) < 2) {
      next
    }

    pairs <- utils::combn(
      groups_present,
      2,
      simplify = FALSE
    )

    for (pair in pairs) {

      g1 <- pair[1]
      g2 <- pair[2]

      d1 <- d |>
        dplyr::filter(
          group == g1
        ) |>
        dplyr::select(
          sample_id,
          analysis_value
        ) |>
        dplyr::rename(
          value_g1 = analysis_value
        )

      d2 <- d |>
        dplyr::filter(
          group == g2
        ) |>
        dplyr::select(
          sample_id,
          analysis_value
        ) |>
        dplyr::rename(
          value_g2 = analysis_value
        )

      n_group1 <- nrow(d1)
      n_group2 <- nrow(d2)

      if (paired) {

        paired_data <- dplyr::inner_join(
          d1,
          d2,
          by = "sample_id"
        )

        n_pair <- nrow(
          paired_data
        )

        if (n_pair >= min_paired_samples) {

          test <- try(
            stats::wilcox.test(
              paired_data$value_g1,
              paired_data$value_g2,
              paired = TRUE,
              exact = FALSE
            ),
            silent = TRUE
          )

          p_value <- if (
            inherits(test, "try-error")
          ) {
            NA_real_
          } else {
            test$p.value
          }

        } else {

          p_value <- NA_real_
        }

        median_g1 <- if (n_pair > 0) {
          stats::median(
            paired_data$value_g1,
            na.rm = TRUE
          )
        } else {
          NA_real_
        }

        median_g2 <- if (n_pair > 0) {
          stats::median(
            paired_data$value_g2,
            na.rm = TRUE
          )
        } else {
          NA_real_
        }

        median_diff <- if (n_pair > 0) {
          stats::median(
            paired_data$value_g2 -
              paired_data$value_g1,
            na.rm = TRUE
          )
        } else {
          NA_real_
        }

      } else {

        n_pair <- NA_integer_

        x <- d1$value_g1[
          is.finite(d1$value_g1)
        ]

        y <- d2$value_g2[
          is.finite(d2$value_g2)
        ]

        if (
          length(x) >= min_paired_samples &&
          length(y) >= min_paired_samples
        ) {

          test <- try(
            stats::wilcox.test(
              x,
              y,
              paired = FALSE,
              exact = FALSE
            ),
            silent = TRUE
          )

          p_value <- if (
            inherits(test, "try-error")
          ) {
            NA_real_
          } else {
            test$p.value
          }

        } else {

          p_value <- NA_real_
        }

        median_g1 <- if (length(x) > 0) {
          stats::median(
            x,
            na.rm = TRUE
          )
        } else {
          NA_real_
        }

        median_g2 <- if (length(y) > 0) {
          stats::median(
            y,
            na.rm = TRUE
          )
        } else {
          NA_real_
        }

        median_diff <- median_g2 - median_g1
      }

      pairwise_list[[counter]] <- data.frame(
        celltype = ct,
        group1 = g1,
        group2 = g2,
        paired = paired,
        n_group1 = n_group1,
        n_group2 = n_group2,
        n_paired_samples = n_pair,
        median_group1 = median_g1,
        median_group2 = median_g2,
        median_difference_group2_minus_group1 =
          median_diff,
        p_value = as.numeric(p_value),
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }
  }

  if (length(pairwise_list) == 0) {

    pairwise_stats <- pairwise_template

  } else {

    pairwise_stats <- dplyr::bind_rows(
      pairwise_list
    )

    pairwise_stats <- pairwise_stats |>
      dplyr::group_by(
        celltype
      ) |>
      dplyr::mutate(
        FDR_within_celltype = stats::p.adjust(
          p_value,
          method = p_adjust_method
        )
      ) |>
      dplyr::ungroup()

    pairwise_stats$FDR_all_tests <- stats::p.adjust(
      pairwise_stats$p_value,
      method = p_adjust_method
    )

    pairwise_stats <- pairwise_stats |>
      dplyr::arrange(
        FDR_within_celltype,
        p_value
      )
  }

  # ==========================================================
  # 15. Significant-result tables
  # ==========================================================

  global_significant <- global_stats |>
    dplyr::filter(
      !is.na(global_FDR),
      global_FDR < alpha
    )

  if (nrow(pairwise_stats) > 0) {

    pairwise_significant <- pairwise_stats |>
      dplyr::filter(
        !is.na(FDR_within_celltype),
        FDR_within_celltype < alpha
      )

  } else {

    pairwise_significant <- pairwise_stats
  }

  # ==========================================================
  # 16. Colors
  # ==========================================================

  default_palette <- biomed_colors(length(group_levels))

  if (is.null(group_colors)) {

    group_colors <- stats::setNames(
      rep(
        default_palette,
        length.out = length(group_levels)
      ),
      group_levels
    )

  } else {

    if (is.null(names(group_colors))) {

      if (length(group_colors) < length(group_levels)) {
        stop(
          "Unnamed `group_colors` has fewer colors than group levels."
        )
      }

      group_colors <- stats::setNames(
        group_colors[
          seq_along(group_levels)
        ],
        group_levels
      )

    } else {

      missing_color_groups <- setdiff(
        group_levels,
        names(group_colors)
      )

      if (length(missing_color_groups) > 0) {

        extra_colors <- rep(
          default_palette,
          length.out = length(missing_color_groups)
        )

        names(extra_colors) <- missing_color_groups

        group_colors <- c(
          group_colors,
          extra_colors
        )
      }

      group_colors <- group_colors[
        group_levels
      ]
    }
  }

  # ==========================================================
  # 17. Individual plot function
  # ==========================================================

  make_plot <- function(ct) {

    d <- analysis_data |>
      dplyr::filter(
        celltype == ct
      )

    fdr <- global_stats$global_FDR[
      match(
        ct,
        global_stats$celltype
      )
    ]

    subtitle_text <- if (
      length(fdr) == 1 &&
      !is.na(fdr)
    ) {

      paste0(
        "Mixed-effects global FDR = ",
        format(
          fdr,
          digits = 3,
          scientific = TRUE
        )
      )

    } else {

      "Mixed-effects global FDR = NA"
    }

    ggplot2::ggplot(
      d,
      ggplot2::aes(
        x = group,
        y = analysis_value,
        fill = group
      )
    ) +
      ggplot2::geom_boxplot(
        width = 0.62,
        outlier.shape = NA,
        alpha = 0.65,
        linewidth = 0.55
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(
          color = group
        ),
        width = 0.12,
        height = 0,
        size = 2.2,
        alpha = 0.85
      ) +
      ggplot2::scale_fill_manual(
        values = group_colors,
        drop = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = group_colors,
        drop = FALSE
      ) +
      ggplot2::scale_x_discrete(
        drop = FALSE
      ) +
      ggplot2::labs(
        title = ct,
        subtitle = subtitle_text,
        x = NULL,
        y = y_label
      ) +
      ggplot2::theme_classic(
        base_size = 12
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 14,
          hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          size = 9,
          hjust = 0.5
        ),
        axis.text.x = ggplot2::element_text(
          angle = 35,
          hjust = 1,
          vjust = 1
        ),
        legend.position = "none"
      )
  }

  # ==========================================================
  # 18. Save multipage PDF
  # ==========================================================

  if (save_multipage_pdf) {

    grDevices::pdf(
      file.path(
        outdir,
        paste0(
          metric_tag,
          "_all_celltypes_multipage.pdf"
        )
      ),
      width = individual_width,
      height = individual_height
    )

    opened_device <- grDevices::dev.cur()
    on.exit(if (opened_device %in% grDevices::dev.list()) grDevices::dev.off(opened_device), add = TRUE)
    for (ct in celltype_cols) {
      print(
        make_plot(ct)
      )
    }

    grDevices::dev.off()
  }

  # ==========================================================
  # 19. Save individual plots
  # ==========================================================

  if (save_individual) {

    for (ct in celltype_cols) {

      p <- make_plot(ct)

      safe_ct <- safe_text(ct)

      ggplot2::ggsave(
        filename = file.path(
          outdir,
          "individual_plots",
          paste0(
            safe_ct,
            "_",
            metric_tag,
            ".png"
          )
        ),
        plot = p,
        width = individual_width,
        height = individual_height,
        dpi = dpi
      )

      ggplot2::ggsave(
        filename = file.path(
          outdir,
          "individual_plots",
          paste0(
            safe_ct,
            "_",
            metric_tag,
            ".pdf"
          )
        ),
        plot = p,
        width = individual_width,
        height = individual_height
      )
    }
  }

  # ==========================================================
  # 20. Combined facet plot
  # ==========================================================

  n_celltypes <- length(
    celltype_cols
  )

  actual_facet_ncol <- min(
    facet_ncol,
    max(1, n_celltypes)
  )

  facet_nrow <- ceiling(
    n_celltypes /
      actual_facet_ncol
  )

  p_all <- ggplot2::ggplot(
    analysis_data,
    ggplot2::aes(
      x = group,
      y = analysis_value,
      fill = group
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.62,
      outlier.shape = NA,
      alpha = 0.65,
      linewidth = 0.35
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(
        color = group
      ),
      width = 0.11,
      height = 0,
      size = 0.9,
      alpha = 0.70
    ) +
    ggplot2::facet_wrap(
      ~celltype,
      scales = "free_y",
      ncol = actual_facet_ncol
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors,
      drop = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = group_colors,
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = y_label
    ) +
    ggplot2::theme_classic(
      base_size = 9
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 8
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 6
      ),
      legend.position = "none"
    )

  combined_height <- max(
    6,
    facet_nrow *
      combined_panel_height
  )

  if (save_combined) {

    ggplot2::ggsave(
      filename = file.path(
        outdir,
        paste0(
          metric_tag,
          "_all_celltypes_combined.png"
        )
      ),
      plot = p_all,
      width = combined_width,
      height = combined_height,
      dpi = dpi,
      limitsize = FALSE
    )

    ggplot2::ggsave(
      filename = file.path(
        outdir,
        paste0(
          metric_tag,
          "_all_celltypes_combined.pdf"
        )
      ),
      plot = p_all,
      width = combined_width,
      height = combined_height,
      limitsize = FALSE
    )
  }

  # ==========================================================
  # 21. Parameter table
  # ==========================================================

  parameter_table <- data.frame(
    parameter = c(
      "metric",
      "abundance_stat",
      "sample_col",
      "group_cols",
      "group_levels",
      "celltype_cols_analyzed",
      "all_celltype_cols",
      "fraction_scale",
      "paired",
      "min_paired_samples",
      "min_samples_global",
      "p_adjust_method",
      "alpha",
      "na_abundance_as_zero"
    ),
    value = c(
      metric,
      abundance_stat,
      sample_col,
      paste(group_cols, collapse = ";"),
      paste(group_levels, collapse = ";"),
      paste(celltype_cols, collapse = ";"),
      paste(all_celltype_cols, collapse = ";"),
      as.character(fraction_scale),
      as.character(paired),
      as.character(min_paired_samples),
      as.character(min_samples_global),
      p_adjust_method,
      as.character(alpha),
      as.character(na_abundance_as_zero)
    ),
    stringsAsFactors = FALSE
  )

  # ==========================================================
  # 22. CSV output
  # ==========================================================

  if (export_csv) {

    utils::write.csv(
      analysis_data,
      file.path(
        outdir,
        paste0(
          "sample_group_celltype_",
          metric_tag,
          ".csv"
        )
      ),
      row.names = FALSE
    )

    utils::write.csv(
      group_summary,
      file.path(
        outdir,
        "group_summary.csv"
      ),
      row.names = FALSE
    )

    utils::write.csv(
      global_stats,
      file.path(
        outdir,
        "global_statistics.csv"
      ),
      row.names = FALSE
    )

    utils::write.csv(
      pairwise_stats,
      file.path(
        outdir,
        "pairwise_statistics.csv"
      ),
      row.names = FALSE
    )

    utils::write.csv(
      spot_number,
      file.path(
        outdir,
        "spot_number.csv"
      ),
      row.names = FALSE
    )
  }

  # ==========================================================
  # 23. Excel output
  # ==========================================================

  xlsx_file <- NULL

  if (export_xlsx) {

    if (!requireNamespace("openxlsx", quietly = TRUE)) {

      warning(
        "Package `openxlsx` is not installed; XLSX was not exported."
      )

    } else {

      xlsx_file <- file.path(
        outdir,
        paste0(
          "celltype_",
          metric_tag,
          "_analysis.xlsx"
        )
      )

      wb <- openxlsx::createWorkbook()

      sheets <- list(
        Parameters = parameter_table,
        Sample_group_data = analysis_data,
        Group_summary = group_summary,
        Global_statistics = global_stats,
        Global_significant = global_significant,
        Pairwise_statistics = pairwise_stats,
        Pairwise_significant = pairwise_significant,
        Spot_number = spot_number
      )

      for (nm in names(sheets)) {

        dat <- sheets[[nm]]

        openxlsx::addWorksheet(
          wb,
          nm
        )

        openxlsx::writeData(
          wb,
          nm,
          dat,
          withFilter = ncol(dat) > 0
        )

        openxlsx::freezePane(
          wb,
          nm,
          firstRow = TRUE
        )

        if (ncol(dat) > 0) {
          openxlsx::setColWidths(
            wb,
            nm,
            cols = seq_len(ncol(dat)),
            widths = "auto"
          )
        }
      }

      openxlsx::saveWorkbook(
        wb,
        xlsx_file,
        overwrite = TRUE
      )
    }
  }

  # ==========================================================
  # 24. Finish
  # ==========================================================

  if (verbose) {

    message("----------------------------------------")
    message("Analysis finished.")
    message("Output directory: ", normalizePath(
      outdir,
      winslash = "/",
      mustWork = FALSE
    ))
    message("Analyzed cell types: ", length(celltype_cols))
    message(
      "Global significant (FDR < ",
      alpha,
      "): ",
      nrow(global_significant)
    )
    message(
      "Pairwise significant (within-celltype FDR < ",
      alpha,
      "): ",
      nrow(pairwise_significant)
    )
    message("----------------------------------------")
  }

  return(
    list(
      metric = metric,
      abundance_stat = abundance_stat,
      celltypes = celltype_cols,
      all_celltype_cols = all_celltype_cols,
      group_levels = group_levels,
      sample_group_data = analysis_data,
      group_summary = group_summary,
      global_stats = global_stats,
      global_significant = global_significant,
      pairwise_stats = pairwise_stats,
      pairwise_significant = pairwise_significant,
      spot_number = spot_number,
      combined_plot = p_all,
      xlsx_file = xlsx_file,
      outdir = outdir
    )
  )
}
