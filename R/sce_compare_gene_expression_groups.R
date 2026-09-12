#' Compare gene expression across groups
#'
#' Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.
#' @param sce A Seurat object; for SCENIC, also accepts the path to an RDS object.
#' @param genes Character vector of genes to compare; missing genes are reported.
#' @param statistic_group One or more metadata columns defining comparison groups.
#' @param split_by_celltype Optional cell-type metadata column; NULL analyzes all cells together.
#' @param statistic_level Use samples as independent units by default; cells requests exploratory cell-level comparisons.
#' @param sample_col Sample identifier column in Seurat metadata.
#' @param Drop_negative Exclude expression values less than or equal to zero before testing/aggregation. Detection summaries still use all cells.
#' @param assay Assay to use. NULL uses the default assay; SCENIC defaults to RNA.
#' @param layer Expression layer to use; data denotes normalized expression. Split Seurat v5 layers are joined on a local copy.
#' @param sample_aggregation Sample-level expression summary: mean or median.
#' @param group_order Optional group order; a named list provides separate orders for different grouping columns. Observed unlisted levels are retained.
#' @param reduction Existing dimensional reduction used for plots.
#' @param output_dir Directory for generated tables, figures, and analysis files.
#' @param file_prefix Prefix for generated filenames.
#' @param export_analysis_data Export underlying analysis data. NULL uses the original sample-level default.
#' @param save_pdf Save PDF figures.
#' @param save_png Save PNG figures.
#' @param width Plot width in inches; NULL selects automatic dimensions where supported.
#' @param height Plot height in inches; NULL selects automatic dimensions where supported.
#' @return A named list per grouping variable containing parameters, summaries, tests, detection rates, analysis data, plots and output files.
#' @export
sce_compare_gene_expression_groups <- function(
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
) {

  required_pkgs <- c("Seurat", "dplyr", "tidyr", "ggplot2", "openxlsx")
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0) {
    stop(
      "Please install missing package(s): ",
      paste(missing_pkgs, collapse = ", ")
    )
  }

  if (!inherits(sce, "Seurat")) {
    stop("sce must be a Seurat object.")
  }

  statistic_level <- match.arg(statistic_level)
  sample_aggregation <- match.arg(sample_aggregation)

  if (length(genes) == 0) {
    stop("genes cannot be empty.")
  }
  genes <- unique(as.character(genes))

  if (length(statistic_group) == 0) {
    stop("statistic_group cannot be empty.")
  }
  statistic_group <- unique(as.character(statistic_group))

  if (length(Drop_negative) != 1 || is.na(Drop_negative)) {
    stop("Drop_negative must be a single TRUE or FALSE.")
  }
  Drop_negative <- isTRUE(Drop_negative)

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(sce)
  }
  if (!assay %in% names(sce@assays)) {
    stop("Assay '", assay, "' was not found in the Seurat object.")
  }

  if (is.null(export_analysis_data)) {
    export_analysis_data <- identical(statistic_level, "samples")
  }

  meta <- sce[[]]

  missing_group_vars <- setdiff(statistic_group, colnames(meta))
  if (length(missing_group_vars) > 0) {
    stop(
      "The following statistic_group variable(s) were not found in metadata: ",
      paste(missing_group_vars, collapse = ", ")
    )
  }

  use_celltype <- !is.null(split_by_celltype) &&
    length(split_by_celltype) == 1 &&
    !is.na(split_by_celltype) &&
    nzchar(as.character(split_by_celltype))

  if (!is.null(split_by_celltype) && length(split_by_celltype) != 1) {
    stop("split_by_celltype must be NULL or a single metadata column name.")
  }

  if (use_celltype && !split_by_celltype %in% colnames(meta)) {
    stop(
      "split_by_celltype '", split_by_celltype,
      "' was not found in Seurat metadata."
    )
  }

  if (identical(statistic_level, "samples")) {
    if (length(sample_col) != 1 || !sample_col %in% colnames(meta)) {
      stop(
        "For statistic_level = 'samples', sample_col must be one valid metadata column. ",
        "Current sample_col: ", paste(sample_col, collapse = ", ")
      )
    }
  }

  if (!is.null(group_order) &&
      length(statistic_group) > 1 &&
      !is.list(group_order)) {
    stop(
      "When statistic_group contains multiple metadata variables, ",
      "group_order must be NULL or a named list, e.g. ",
      "list(response = c('MPR','non-MPR'), therapy = c('TKI','TKI+chemo'))."
    )
  }

  assay_features <- rownames(sce[[assay]])
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

  reduction_available <- reduction %in% names(sce@reductions)
  if (!reduction_available) {
    warning(
      "Reduction '", reduction,
      "' was not found. Statistical analysis will run, but FeaturePlot will be skipped."
    )
  }

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

  message(
    "Extracting expression once from assay='", assay,
    "', layer='", layer, "' ..."
  )

  expr <- as.data.frame(t(as.matrix(
    .biomed_sce_expression(sce, assay, layer)[genes_found, , drop = FALSE]
  )), check.names = FALSE)
  if (any(genes_found %in% c(".cell", ".celltype", ".sample", ".group"))) {
    stop("Gene names conflict with reserved internal columns.")
  }
  expr$.cell <- rownames(expr)

  common_meta <- data.frame(
    .cell = rownames(meta),
    stringsAsFactors = FALSE
  )

  if (use_celltype) {
    common_meta$.celltype <- as.character(meta[[split_by_celltype]])
  } else {
    common_meta$.celltype <- "All_cells"
  }

  if (identical(statistic_level, "samples")) {
    common_meta$.sample <- as.character(meta[[sample_col]])
  }



  base_dat <- dplyr::left_join(expr, common_meta, by = ".cell")

  all_results <- vector("list", length(statistic_group))
  names(all_results) <- statistic_group

  for (group_var in statistic_group) {

    message("")
    message("============================================================")
    message("Analyzing statistic_group: ", group_var)
    message("============================================================")

    group_dir <- file.path(output_dir, sanitize_filename(group_var))
    dot_dir <- file.path(group_dir, "DotPlot")
    violin_dir <- file.path(group_dir, "ViolinPlot")
    feature_dir <- file.path(group_dir, "FeaturePlot")
    dir.create(dot_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(violin_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)

    current_dat <- base_dat
    current_dat$.group <- as.character(meta[[group_var]][match(current_dat$.cell, rownames(meta))])

    keep <- !is.na(current_dat$.group) & !is.na(current_dat$.celltype)
    if (identical(statistic_level, "samples")) {
      keep <- keep & !is.na(current_dat$.sample)
    }
    current_dat <- current_dat[keep, , drop = FALSE]

    if (nrow(current_dat) == 0) {
      warning("No usable cells for statistic_group='", group_var, "'. Skipping.")
      next
    }

    group_levels <- get_group_levels(group_var, current_dat$.group)
    current_dat$.group <- factor(current_dat$.group, levels = group_levels)

    if (length(group_levels) < 2) {
      warning(
        "statistic_group='", group_var,
        "' has fewer than 2 non-missing groups. Tests will be empty, but summaries/plots will still be produced."
      )
    }

    if (identical(statistic_level, "samples")) {
      sample_group_check <- current_dat |>
        dplyr::distinct(.data$.sample, .data$.group) |>
        dplyr::group_by(.data$.sample) |>
        dplyr::summarise(
          n_group = dplyr::n_distinct(.data$.group),
          groups = paste(unique(as.character(.data$.group)), collapse = ";"),
          .groups = "drop"
        )

      bad_samples <- sample_group_check |>
        dplyr::filter(.data$n_group > 1)

      if (nrow(bad_samples) > 0) {
        stop(
          "For statistic_group='", group_var,
          "', some sample IDs map to more than one group. ",
          "sample_col='", sample_col, "' must uniquely identify a biological sample within this comparison.\n",
          "Problematic samples: ",
          paste(bad_samples$.sample, collapse = ", ")
        )
      }
    }

    long_raw <- tidyr::pivot_longer(
      current_dat,
      cols = dplyr::all_of(genes_found),
      names_to = "gene",
      values_to = "raw_expression"
    ) |>
      dplyr::filter(is.finite(.data$raw_expression))

    cell_detection <- long_raw |>
      dplyr::group_by(.data$.celltype, .data$gene, .data$.group) |>
      dplyr::summarise(
        n_cells = dplyr::n(),
        n_positive_cells = sum(.data$raw_expression > 0),
        pct_positive_cells = mean(.data$raw_expression > 0) * 100,
        mean_all_cells = mean(.data$raw_expression),
        median_all_cells = stats::median(.data$raw_expression),
        mean_positive_cells = {
          x <- .data$raw_expression[.data$raw_expression > 0]
          if (length(x) == 0) NA_real_ else mean(x)
        },
        .groups = "drop"
      )

    if (identical(statistic_level, "cells")) {

      stat_data_all <- long_raw |>
        dplyr::transmute(
          statistical_unit = .data$.cell,
          cell = .data$.cell,
          celltype = .data$.celltype,
          group = .data$.group,
          gene = .data$gene,
          raw_expression = .data$raw_expression,
          detected = .data$raw_expression > 0,
          expression = dplyr::if_else(
            Drop_negative & .data$raw_expression <= 0,
            NA_real_,
            .data$raw_expression
          )
        )

    } else {

      stat_data_all <- long_raw |>
        dplyr::group_by(.data$.sample, .data$.celltype, .data$.group, .data$gene) |>
        dplyr::summarise(
          n_cells = dplyr::n(),
          n_positive_cells = sum(.data$raw_expression > 0),
          pct_positive_cells = mean(.data$raw_expression > 0) * 100,
          detected = any(.data$raw_expression > 0),
          expression = {
            x <- .data$raw_expression
            if (Drop_negative) {
              x <- x[x > 0]
            }
            if (length(x) == 0) {
              NA_real_
            } else if (identical(sample_aggregation, "mean")) {
              mean(x)
            } else {
              stats::median(x)
            }
          },
          .groups = "drop"
        ) |>
        dplyr::transmute(
          statistical_unit = .data$.sample,
          sample = .data$.sample,
          celltype = .data$.celltype,
          group = .data$.group,
          gene = .data$gene,
          n_cells = .data$n_cells,
          n_positive_cells = .data$n_positive_cells,
          pct_positive_cells = .data$pct_positive_cells,
          detected = .data$detected,
          expression = .data$expression
        )
    }

    stat_data_all$group <- factor(stat_data_all$group, levels = group_levels)
    stat_data <- stat_data_all |>
      dplyr::filter(is.finite(.data$expression))

    summary_table <- stat_data_all |>
      dplyr::group_by(.data$celltype, .data$gene, .data$group) |>
      dplyr::summarise(
        n_units_total = dplyr::n(),
        n_units_used = sum(is.finite(.data$expression)),
        n_positive_units = sum(.data$detected, na.rm = TRUE),
        pct_positive_units = mean(.data$detected, na.rm = TRUE) * 100,
        mean_expression = mean_or_na(.data$expression),
        median_expression = median_or_na(.data$expression),
        sd_expression = sd_or_na(.data$expression),
        sem_expression = {
          n_used <- sum(is.finite(.data$expression))
          sd_val <- sd_or_na(.data$expression)
          if (n_used <= 1 || is.na(sd_val)) NA_real_ else sd_val / sqrt(n_used)
        },
        min_expression = {
          x <- .data$expression[is.finite(.data$expression)]
          if (length(x) == 0) NA_real_ else min(x)
        },
        max_expression = {
          x <- .data$expression[is.finite(.data$expression)]
          if (length(x) == 0) NA_real_ else max(x)
        },
        .groups = "drop"
      )

    overall_rows <- list()
    pairwise_rows <- list()
    overall_i <- 1L
    pair_i <- 1L

    strata <- unique(stat_data[, c("celltype", "gene"), drop = FALSE])

    if (nrow(strata) > 0) {
      for (ii in seq_len(nrow(strata))) {
        ct_now <- strata$celltype[ii]
        gene_now <- strata$gene[ii]

        dd <- stat_data |>
          dplyr::filter(.data$celltype == ct_now, .data$gene == gene_now)

        present_groups <- group_levels[
          group_levels %in% as.character(unique(dd$group))
        ]
        present_groups <- present_groups[
          vapply(
            present_groups,
            function(g) sum(as.character(dd$group) == g) > 0,
            logical(1)
          )
        ]

        if (length(present_groups) == 2) {
          g1 <- present_groups[1]
          g2 <- present_groups[2]
          x <- dd$expression[as.character(dd$group) == g1]
          y <- dd$expression[as.character(dd$group) == g2]

          test_obj <- tryCatch(
            suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)),
            error = function(e) NULL
          )

          overall_rows[[overall_i]] <- data.frame(
            celltype = ct_now,
            gene = gene_now,
            n_groups = 2L,
            test = "Wilcoxon rank-sum",
            statistic = if (is.null(test_obj)) NA_real_ else unname(test_obj$statistic),
            p_value = if (is.null(test_obj)) NA_real_ else test_obj$p.value,
            stringsAsFactors = FALSE
          )
          overall_i <- overall_i + 1L

        } else if (length(present_groups) >= 3) {
          dd_test <- dd |>
            dplyr::filter(as.character(.data$group) %in% present_groups)
          dd_test$group <- droplevels(dd_test$group)

          test_obj <- tryCatch(
            stats::kruskal.test(expression ~ group, data = dd_test),
            error = function(e) NULL
          )

          overall_rows[[overall_i]] <- data.frame(
            celltype = ct_now,
            gene = gene_now,
            n_groups = length(present_groups),
            test = "Kruskal-Wallis",
            statistic = if (is.null(test_obj)) NA_real_ else unname(test_obj$statistic),
            p_value = if (is.null(test_obj)) NA_real_ else test_obj$p.value,
            stringsAsFactors = FALSE
          )
          overall_i <- overall_i + 1L
        }

        if (length(present_groups) >= 2) {
          group_pairs <- utils::combn(present_groups, 2, simplify = FALSE)

          for (pair in group_pairs) {
            g1 <- pair[1]
            g2 <- pair[2]
            x <- dd$expression[as.character(dd$group) == g1]
            y <- dd$expression[as.character(dd$group) == g2]

            test_obj <- tryCatch(
              suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)),
              error = function(e) NULL
            )

            pairwise_rows[[pair_i]] <- data.frame(
              celltype = ct_now,
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
              statistic = if (is.null(test_obj)) NA_real_ else unname(test_obj$statistic),
              p_value = if (is.null(test_obj)) NA_real_ else test_obj$p.value,
              stringsAsFactors = FALSE
            )
            pair_i <- pair_i + 1L
          }
        }
      }
    }

    if (length(overall_rows) > 0) {
      overall_test <- dplyr::bind_rows(overall_rows)
      overall_test$p_adj_BH_global <- stats::p.adjust(overall_test$p_value, method = "BH")
    } else {
      overall_test <- data.frame(
        celltype = character(),
        gene = character(),
        n_groups = integer(),
        test = character(),
        statistic = numeric(),
        p_value = numeric(),
        p_adj_BH_global = numeric(),
        stringsAsFactors = FALSE
      )
    }

    if (length(pairwise_rows) > 0) {
      pairwise_test <- dplyr::bind_rows(pairwise_rows) |>
        dplyr::group_by(.data$celltype, .data$gene) |>
        dplyr::mutate(
          p_adj_BH_within_gene_celltype = stats::p.adjust(.data$p_value, method = "BH")
        ) |>
        dplyr::ungroup()
      pairwise_test$p_adj_BH_global <- stats::p.adjust(pairwise_test$p_value, method = "BH")
    } else {
      pairwise_test <- data.frame(
        celltype = character(),
        gene = character(),
        group1 = character(),
        group2 = character(),
        contrast = character(),
        n_group1 = integer(),
        n_group2 = integer(),
        mean_group1 = numeric(),
        mean_group2 = numeric(),
        median_group1 = numeric(),
        median_group2 = numeric(),
        delta_mean = numeric(),
        delta_median = numeric(),
        statistic = numeric(),
        p_value = numeric(),
        p_adj_BH_within_gene_celltype = numeric(),
        p_adj_BH_global = numeric(),
        stringsAsFactors = FALSE
      )
    }

    dotplot_data <- stat_data_all |>
      dplyr::group_by(.data$celltype, .data$gene, .data$group) |>
      dplyr::summarise(
        n_units = dplyr::n(),
        n_units_used = sum(is.finite(.data$expression)),
        avg_expression = mean_or_na(.data$expression),
        pct_positive = mean(.data$detected, na.rm = TRUE) * 100,
        .groups = "drop"
      )

    dotplot_data$group <- factor(dotplot_data$group, levels = group_levels)

    p_dot <- ggplot2::ggplot(
      dotplot_data,
      ggplot2::aes(
        x = .data$gene,
        y = .data$group,
        size = .data$pct_positive,
        fill = .data$avg_expression
      )
    ) +
      ggplot2::geom_point(shape = 21, na.rm = TRUE) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$celltype),
        scales = "free_y",
        space = "free_y"
      ) +
      ggplot2::scale_size(range = c(1, 8), limits = c(0, 100)) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        strip.text.y = ggplot2::element_text(angle = 0)
      ) +
      ggplot2::labs(
        x = NULL,
        y = group_var,
        size = if (identical(statistic_level, "samples")) {
          "% positive samples"
        } else {
          "% positive cells"
        },
        fill = "Average\nexpression",
        title = paste0(group_var, " | ", statistic_level, "-level DotPlot")
      )

    dot_height <- max(height, 2.5 + length(unique(dotplot_data$celltype)) * 1.3)
    save_plot(
      p_dot,
      file.path(
        dot_dir,
        paste0(sanitize_filename(file_prefix), "_", sanitize_filename(group_var), "_DotPlot")
      ),
      w = width,
      h = dot_height
    )

    violin_plots <- list()
    celltypes_present <- unique(as.character(stat_data_all$celltype))

    for (ct in celltypes_present) {
      dd_plot <- stat_data_all |>
        dplyr::filter(.data$celltype == ct, is.finite(.data$expression))

      if (nrow(dd_plot) == 0) {
        next
      }

      dd_plot$group <- factor(dd_plot$group, levels = group_levels)

      p_vln <- ggplot2::ggplot(
        dd_plot,
        ggplot2::aes(x = .data$group, y = .data$expression)
      ) +
        ggplot2::geom_violin(trim = FALSE, scale = "width", na.rm = TRUE) +
        ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, na.rm = TRUE)

      if (identical(statistic_level, "samples")) {
        p_vln <- p_vln +
          ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.5, alpha = 0.75)
      } else if (nrow(dd_plot) <= 5000) {
        p_vln <- p_vln +
          ggplot2::geom_jitter(width = 0.12, height = 0, size = 0.35, alpha = 0.25)
      }

      p_vln <- p_vln +
        ggplot2::facet_wrap(ggplot2::vars(.data$gene), scales = "free_y") +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        ) +
        ggplot2::labs(
          x = group_var,
          y = if (identical(statistic_level, "samples")) {
            paste0("Sample ", sample_aggregation, " expression")
          } else {
            "Cell expression"
          },
          title = paste0(ct, " | ", group_var)
        )

      violin_plots[[ct]] <- p_vln

      save_plot(
        p_vln,
        file.path(
          violin_dir,
          paste0(
            sanitize_filename(file_prefix), "_",
            sanitize_filename(group_var), "_",
            sanitize_filename(ct), "_Violin"
          )
        ),
        w = width,
        h = height
      )
    }

    featureplot_files <- character()

    if (reduction_available) {
      for (ct in unique(as.character(long_raw$.celltype))) {
        cells_ct <- unique(long_raw$.cell[long_raw$.celltype == ct])
        cells_ct <- intersect(cells_ct, colnames(sce))

        if (length(cells_ct) == 0) {
          next
        }

        obj_ct <- subset(sce, cells = cells_ct)
        Seurat::DefaultAssay(obj_ct) <- assay

        group_values_ct <- obj_ct[[]][[group_var]]
        keep_cells <- colnames(obj_ct)[!is.na(group_values_ct)]

        if (length(keep_cells) == 0) {
          next
        }

        obj_ct <- subset(obj_ct, cells = keep_cells)

        groups_ct <- unique(as.character(obj_ct[[]][[group_var]]))
        groups_ct <- group_levels[group_levels %in% groups_ct]

        if (length(groups_ct) == 0) {
          next
        }

        obj_ct[[group_var]] <- factor(
          obj_ct[[]][[group_var]],
          levels = group_levels
        )

        p_feature <- tryCatch(
          Seurat::FeaturePlot(
            object = obj_ct,
            features = genes_found,
            reduction = reduction,
            slot = layer,
            split.by = group_var,
            order = TRUE,
            ncol = max(1, length(groups_ct)),
            combine = TRUE
          ),
          error = function(e) {
            warning(
              "FeaturePlot failed for statistic_group='", group_var,
              "', celltype='", ct, "': ", conditionMessage(e)
            )
            NULL
          }
        )

        if (!is.null(p_feature)) {
          feature_basename <- file.path(
            feature_dir,
            paste0(
              sanitize_filename(file_prefix), "_",
              sanitize_filename(group_var), "_",
              sanitize_filename(ct), "_FeaturePlot"
            )
          )

          feature_w <- max(width, length(groups_ct) * 3.5)
          feature_h <- max(height, length(genes_found) * 3.0)

          save_plot(
            p_feature,
            feature_basename,
            w = feature_w,
            h = feature_h
          )

          featureplot_files <- c(featureplot_files, feature_basename)
        }
      }
    }

    parameter_table <- data.frame(
      parameter = c(
        "genes_requested",
        "genes_found",
        "genes_missing",
        "statistic_group",
        "group_levels",
        "split_by_celltype",
        "statistic_level",
        "sample_col",
        "Drop_negative",
        "assay",
        "layer",
        "sample_aggregation",
        "reduction"
      ),
      value = c(
        paste(genes, collapse = ";"),
        paste(genes_found, collapse = ";"),
        paste(genes_missing, collapse = ";"),
        group_var,
        paste(group_levels, collapse = ";"),
        if (use_celltype) split_by_celltype else "All_cells",
        statistic_level,
        if (identical(statistic_level, "samples")) sample_col else "NA",
        as.character(Drop_negative),
        assay,
        layer,
        if (identical(statistic_level, "samples")) sample_aggregation else "NA",
        reduction
      ),
      stringsAsFactors = FALSE
    )

    xlsx_file <- file.path(
      group_dir,
      paste0(
        sanitize_filename(file_prefix), "_",
        sanitize_filename(group_var),
        "_statistics.xlsx"
      )
    )

    wb <- openxlsx::createWorkbook()

    openxlsx::addWorksheet(wb, "Parameters")
    openxlsx::writeData(wb, "Parameters", parameter_table)

    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", summary_table)

    openxlsx::addWorksheet(wb, "Cell_Detection")
    openxlsx::writeData(wb, "Cell_Detection", cell_detection)

    openxlsx::addWorksheet(wb, "Overall_Test")
    openxlsx::writeData(wb, "Overall_Test", overall_test)

    openxlsx::addWorksheet(wb, "Pairwise_Test")
    openxlsx::writeData(wb, "Pairwise_Test", pairwise_test)

    openxlsx::addWorksheet(wb, "DotPlot_Data")
    openxlsx::writeData(wb, "DotPlot_Data", dotplot_data)

    if (isTRUE(export_analysis_data)) {
      if (nrow(stat_data_all) <= 1048575) {
        openxlsx::addWorksheet(wb, "Analysis_Data")
        openxlsx::writeData(wb, "Analysis_Data", stat_data_all)
      } else {
        warning(
          "Analysis_Data for statistic_group='", group_var,
          "' has more rows than Excel can store and was not written to the workbook."
        )
      }
    }

    for (sheet_name in openxlsx::sheets(wb)) {
      openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
    }

    openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)

    all_results[[group_var]] <- list(
      parameters = parameter_table,
      summary = summary_table,
      cell_detection = cell_detection,
      overall_test = overall_test,
      pairwise_test = pairwise_test,
      dotplot_data = dotplot_data,
      analysis_data = stat_data_all,
      plots = list(
        dotplot = p_dot,
        violin = violin_plots
      ),
      files = list(
        xlsx = xlsx_file,
        output_dir = group_dir,
        dotplot_dir = dot_dir,
        violin_dir = violin_dir,
        featureplot_dir = feature_dir,
        featureplot_files = featureplot_files
      )
    )

    message("Finished: ", group_var)
    message("Excel: ", xlsx_file)
  }

  message("")
  message("============================================================")
  message("All requested statistic_group analyses completed.")
  message("Output root: ", normalizePath(output_dir, mustWork = FALSE))
  message("============================================================")

  invisible(all_results)
}
