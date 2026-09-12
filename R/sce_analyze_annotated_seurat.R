#' Analyze an annotated Seurat object
#'
#' Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.
#' @param seurat An annotated Seurat object.
#' @param clinical_var Clinical metadata columns for stratification. Each must be constant within a sample.
#' @param celltype Metadata column containing cell-type annotations.
#' @param statistic_level Use samples as independent units by default; cells requests exploratory cell-level comparisons.
#' @param sample_col Sample identifier column in Seurat metadata.
#' @param assay Assay to use. NULL uses the default assay; SCENIC defaults to RNA.
#' @param layer Expression layer to use; data denotes normalized expression. Split Seurat v5 layers are joined on a local copy.
#' @param reduction Existing dimensional reduction used for plots.
#' @param output_dir Directory for generated tables, figures, and analysis files.
#' @param file_prefix Prefix for generated filenames.
#' @param export_analysis_data Export underlying analysis data. NULL uses the original sample-level default.
#' @param save_pdf Save PDF figures.
#' @param save_png Save PNG figures.
#' @param width Plot width in inches; NULL selects automatic dimensions where supported.
#' @param height Plot height in inches; NULL selects automatic dimensions where supported.
#' @param group_order Optional group order; a named list provides separate orders for different grouping columns. Observed unlisted levels are retained.
#' @param run_markers Run marker discovery and export marker tables and plots.
#' @param marker_only_pos Restrict marker discovery to positive markers.
#' @param marker_logfc_threshold Log fold-change threshold passed to FindAllMarkers.
#' @param marker_min_pct Minimum detection fraction passed to FindAllMarkers.
#' @param marker_padj_cutoff Adjusted P-value cutoff for selecting markers; fallback ranking is reported when none pass.
#' @param top_n_marker Maximum markers retained per cell type in summary tables.
#' @param top_n_plot Maximum genes or regulons selected for feature plots.
#' @param marker_test Test passed to Seurat FindAllMarkers.
#' @param run_celltype_correlation Compute correlations between average cell-type expression profiles.
#' @param run_sample_celltype_correlation Compute sample correlations within each cell type.
#' @param correlation_method Correlation method passed to stats::cor.
#' @param correlation_features Use variable features when available, or all informative features.
#' @param min_cells_per_sample_celltype Minimum cells required for a sample/cell-type expression profile.
#' @param run_diversity Compute sample-level Shannon and Simpson diversity and comparisons.
#' @param run_composition_correlation Compute sample-level cell-type fraction correlations.
#' @param custom_colors Discrete palette; defaults to the 14 built-in biomed colors.
#' @param max_panels_per_page Maximum split panels on a plot page.
#' @param max_cells_heatmap Maximum number of cells sampled for marker heatmaps.
#' @param seed Random seed for reproducibility; the caller RNG state is restored.
#' @param verbose Print progress messages.
#' @return Invisibly, a list containing composition, marker, correlation, diversity, statistical tables, parameters and output directories.
#' @export
sce_analyze_annotated_seurat <- function(
    seurat,
    clinical_var = NULL,
    celltype = "cell_type__custom",
    statistic_level = c("samples"),
    sample_col = "orig.ident",
    assay = NULL,
    layer = "data",
    reduction = "umap",
    output_dir = "annotated_cell_analysis",
    file_prefix = "annotated",
    export_analysis_data = NULL,
    save_pdf = TRUE,
    save_png = TRUE,
    width = NULL,
    height = NULL,
    group_order = NULL,
    run_markers = TRUE,
    marker_only_pos = TRUE,
    marker_logfc_threshold = 0.25,
    marker_min_pct = 0.10,
    marker_padj_cutoff = 0.05,
    top_n_marker = 50,
    top_n_plot = 5,
    marker_test = "wilcox",
    run_celltype_correlation = TRUE,
    run_sample_celltype_correlation = TRUE,
    correlation_method = "spearman",
    correlation_features = c("variable", "all"),
    min_cells_per_sample_celltype = 10,
    run_diversity = TRUE,
    run_composition_correlation = TRUE,
    custom_colors = biomed_colors(),
    max_panels_per_page = 12,
    max_cells_heatmap = 3000,
    seed = 1234,
    verbose = TRUE
) {

  pkgs <- c(
    "Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "tidyr",
    "tibble", "patchwork", "openxlsx", "pheatmap", "gplots"
  )
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      "\nInstall them before running this function."
    )
  }

  if (!inherits(seurat, "Seurat")) stop("`seurat` must be a Seurat object.")

  statistic_level[statistic_level == "cells"] <- "cell"
  statistic_level <- unique(match.arg(
    statistic_level,
    choices = c("samples", "cell"),
    several.ok = TRUE
  ))
  correlation_features <- match.arg(correlation_features)
  if (!correlation_method %in% c("pearson", "spearman", "kendall")) {
    stop("`correlation_method` must be pearson, spearman, or kendall.")
  }
  if (top_n_marker < 1 || top_n_plot < 1) stop("top_n_marker/top_n_plot must be >= 1.")
  if (top_n_plot > top_n_marker) top_n_plot <- top_n_marker
  if (max_panels_per_page < 1) max_panels_per_page <- 12
  if (max_cells_heatmap < 100) max_cells_heatmap <- 100

  if (is.null(export_analysis_data)) export_analysis_data <- TRUE
  export_analysis_data <- isTRUE(export_analysis_data)

  msg <- function(...) {
    if (isTRUE(verbose)) message(sprintf(...))
  }

  clean_filename <- function(x) {
    x <- gsub("[\\\\/:*?\"<>|]", "_", as.character(x))
    x <- gsub("\\s+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    ifelse(nchar(x) == 0, "unnamed", x)
  }

  clean_sheet <- function(x, used = character()) {
    x <- gsub("[\\[\\]:*?/\\\\]", "_", as.character(x))
    x <- substr(x, 1, 31)
    if (!x %in% used) return(x)
    base <- substr(x, 1, 27)
    i <- 1
    candidate <- paste0(base, "_", i)
    while (candidate %in% used) {
      i <- i + 1
      candidate <- paste0(base, "_", i)
    }
    substr(candidate, 1, 31)
  }

  ensure_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
  }

  clamp <- function(x, lo, hi) max(lo, min(hi, x))

  auto_dim <- function(type = c("generic", "dim_split", "sample_bar", "dot", "heatmap"),
                       n1 = 1, n2 = 1, user_w = width, user_h = height) {
    type <- match.arg(type)
    if (!is.null(user_w) && !is.null(user_h)) return(c(user_w, user_h))

    if (type == "dim_split") {
      w <- clamp(4.2 * min(n1, 4), 8, 22)
      h <- clamp(4.0 * ceiling(n1 / 4), 5, 24)
    } else if (type == "sample_bar") {
      w <- clamp(5 + 0.28 * n1, 10, 30)
      h <- clamp(5 + 0.18 * n2, 6, 14)
    } else if (type == "dot") {
      w <- clamp(6 + 0.18 * n1, 10, 30)
      h <- clamp(4 + 0.35 * n2, 6, 20)
    } else if (type == "heatmap") {
      w <- clamp(5 + 0.30 * n1, 7, 28)
      h <- clamp(5 + 0.18 * n2, 7, 28)
    } else {
      w <- 10
      h <- 7
    }
    if (!is.null(user_w)) w <- user_w
    if (!is.null(user_h)) h <- user_h
    c(w, h)
  }

  save_plot <- function(plot, filename_no_ext, w = 10, h = 7, dpi = 220) {
    if (isTRUE(save_pdf)) {
      ggplot2::ggsave(
        filename = paste0(filename_no_ext, ".pdf"), plot = plot,
        width = w, height = h, units = "in", limitsize = FALSE
      )
    }
    if (isTRUE(save_png)) {
      ggplot2::ggsave(
        filename = paste0(filename_no_ext, ".png"), plot = plot,
        width = w, height = h, units = "in", dpi = dpi, limitsize = FALSE
      )
    }
  }

  save_pheatmap <- function(mat, filename_no_ext, w = 9, h = 8,
                            annotation_row = NULL, annotation_col = NULL,
                            main = NULL, cluster_rows = TRUE, cluster_cols = TRUE,
                            show_rownames = TRUE, show_colnames = TRUE) {
    if (nrow(mat) < 2 || ncol(mat) < 2) return(invisible(NULL))
    draw_one <- function() {
      pheatmap::pheatmap(
        mat,
        annotation_row = annotation_row,
        annotation_col = annotation_col,
        main = main,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_rownames = show_rownames,
        show_colnames = show_colnames,
        border_color = NA,
        silent = FALSE
      )
    }
    if (isTRUE(save_pdf)) {
      grDevices::pdf(paste0(filename_no_ext, ".pdf"), width = w, height = h)
      tryCatch(draw_one(), finally = grDevices::dev.off())
    }
    if (isTRUE(save_png)) {
      grDevices::png(
        paste0(filename_no_ext, ".png"),
        width = w, height = h, units = "in", res = 220
      )
      tryCatch(draw_one(), finally = grDevices::dev.off())
    }
  }

  write_xlsx_safe <- function(x, file, freeze_first_row = TRUE) {
    if (!is.list(x) || is.data.frame(x) || is.matrix(x)) x <- list(Sheet1 = x)
    if (is.null(names(x)) || any(!nzchar(names(x)))) {
      names(x) <- paste0("Sheet", seq_along(x))
    }
    wb <- openxlsx::createWorkbook()
    used <- character()
    for (nm in names(x)) {
      sh <- clean_sheet(nm, used)
      used <- c(used, sh)
      openxlsx::addWorksheet(wb, sh)
      openxlsx::writeData(wb, sh, x[[nm]], rowNames = FALSE)
      if (isTRUE(freeze_first_row)) {
        openxlsx::freezePane(wb, sh, firstRow = TRUE)
      }
    }
    openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
  }

  write_csv_safe <- function(x, file) {
    utils::write.csv(x, file = file, row.names = FALSE)
  }

  write_rds_csv <- function(x, rds_file, csv_file, rowname_col = NULL) {
    saveRDS(x, rds_file)
    if (is.matrix(x)) {
      out <- as.data.frame(x, check.names = FALSE)
      if (!is.null(rowname_col)) {
        out <- cbind(stats::setNames(data.frame(rownames(x), stringsAsFactors = FALSE), rowname_col), out)
      }
    } else {
      out <- x
    }
    write_csv_safe(out, csv_file)
  }

  discrete_palette <- function(values, palette = custom_colors) {
    vals <- unique(as.character(values))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) == 0 || is.null(palette) || length(palette) < length(vals)) return(NULL)
    stats::setNames(palette[seq_along(vals)], vals)
  }

  add_manual_fill <- function(p, pal) {
    if (is.null(pal)) return(p)
    p + ggplot2::scale_fill_manual(values = pal, drop = FALSE)
  }

  balloon_dim <- function(n_rows, n_cols) {
    w <- clamp(5 + 0.45 * n_cols, 7, 24)
    h <- clamp(4 + 0.32 * n_rows, 6, 30)
    if (!is.null(width)) w <- width
    if (!is.null(height)) h <- height
    c(w, h)
  }

  save_balloonplot <- function(tab, filename_no_ext, title, xlab, ylab) {
    if (length(tab) == 0 || nrow(tab) == 0 || ncol(tab) == 0) return(invisible(NULL))
    d <- balloon_dim(nrow(tab), ncol(tab))
    draw_one <- function() {
      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(oldpar), add = TRUE)
      graphics::par(mar = c(8, 10, 4, 2) + 0.1)
      gplots::balloonplot(
        tab,
        main = title,
        xlab = xlab,
        ylab = ylab,
        dotcolor = if (!is.null(custom_colors) && length(custom_colors) >= 1) custom_colors[1] else "skyblue",
        label = FALSE,
        show.margins = TRUE
      )
    }
    if (isTRUE(save_pdf)) {
      grDevices::pdf(paste0(filename_no_ext, ".pdf"), width = d[1], height = d[2])
      tryCatch(draw_one(), finally = grDevices::dev.off())
    }
    if (isTRUE(save_png)) {
      grDevices::png(paste0(filename_no_ext, ".png"), width = d[1], height = d[2], units = "in", res = 220)
      tryCatch(draw_one(), finally = grDevices::dev.off())
    }
    invisible(NULL)
  }

  get_expr <- function(obj, assay_use, layer_use) {
    .biomed_sce_expression(obj, assay_use, layer_use)
  }

  avg_by_group <- function(expr, groups, min_cells = 1) {
    groups <- as.character(groups)
    ok <- !is.na(groups) & nzchar(groups)
    expr <- expr[, ok, drop = FALSE]
    groups <- groups[ok]
    lev <- unique(groups)
    counts <- table(groups)
    lev <- lev[counts[lev] >= min_cells]
    if (length(lev) == 0) return(NULL)
    out <- lapply(lev, function(g) {
      idx <- which(groups == g)
      Matrix::rowMeans(expr[, idx, drop = FALSE])
    })
    out <- do.call(cbind, out)
    if (is.null(dim(out))) out <- matrix(out, ncol = 1)
    rownames(out) <- rownames(expr)
    colnames(out) <- lev
    out
  }

  choose_cor_features <- function(mat, obj, assay_use, mode = "variable") {
    feats <- rownames(mat)
    if (mode == "variable") {
      vf <- tryCatch(
        SeuratObject::VariableFeatures(obj[[assay_use]]),
        error = function(e) character()
      )
      vf <- intersect(vf, feats)
      if (length(vf) >= 50) feats <- vf
    }
    if (ncol(mat) >= 2) {
      vars <- apply(mat[feats, , drop = FALSE], 1, stats::var, na.rm = TRUE)
      feats <- feats[is.finite(vars) & vars > 0]
    }
    feats
  }

  safe_cor <- function(mat, method = "spearman") {
    if (is.null(mat) || ncol(mat) < 2 || nrow(mat) < 2) return(NULL)
    stats::cor(mat, method = method, use = "pairwise.complete.obs")
  }

  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(seurat)
  if (!assay %in% names(seurat@assays)) stop("Assay not found: ", assay)
  if (!reduction %in% names(seurat@reductions)) {
    stop("Reduction `", reduction, "` not found. Available: ",
         paste(names(seurat@reductions), collapse = ", "))
  }

  md <- seurat[[]]
  required_meta <- unique(c(celltype, sample_col, clinical_var))
  missing_meta <- setdiff(required_meta, colnames(md))
  if (length(missing_meta) > 0) {
    stop("Metadata columns not found: ", paste(missing_meta, collapse = ", "))
  }

  split_vars <- unique(c(sample_col, clinical_var))
  split_vars <- split_vars[!is.na(split_vars) & nzchar(split_vars)]

  if (any(is.na(md[[celltype]]) | !nzchar(as.character(md[[celltype]])))) {
    n_bad <- sum(is.na(md[[celltype]]) | !nzchar(as.character(md[[celltype]])))
    stop("`", celltype, "` contains ", n_bad,
         " NA/empty cells. Please finish annotation or remove these cells first.")
  }
  if (any(is.na(md[[sample_col]]) | !nzchar(as.character(md[[sample_col]])))) {
    stop("`", sample_col, "` contains NA/empty values.")
  }

  if (is.factor(md[[celltype]])) {
    celltype_levels <- levels(droplevels(md[[celltype]]))
  } else {
    celltype_levels <- sort(unique(as.character(md[[celltype]])))
  }

  if (!is.null(group_order)) {
    if (!is.list(group_order) || is.null(names(group_order))) {
      stop("`group_order` must be a named list, e.g. list(group=c('A','B')).")
    }
    for (v in intersect(names(group_order), colnames(md))) {
      md[[v]] <- factor(md[[v]], levels = unique(c(group_order[[v]], as.character(md[[v]]))))
      group_order[[v]] <- levels(md[[v]])
      seurat[[v]] <- md[[v]]
    }
  }

  samples <- unique(as.character(md[[sample_col]]))
  sample_meta <- data.frame(sample_id = samples, stringsAsFactors = FALSE)
  names(sample_meta)[1] <- sample_col

  for (v in setdiff(split_vars, sample_col)) {
    vals_out <- vector("character", length(samples))
    for (i in seq_along(samples)) {
      z <- md[[v]][as.character(md[[sample_col]]) == samples[i]]
      z <- unique(as.character(z[!is.na(z)]))
      if (length(z) > 1) {
        stop(
          "Clinical variable `", v, "` is not constant within sample `",
          samples[i], "`: ", paste(z, collapse = ", "),
          ". Sample-level statistics would be invalid."
        )
      }
      vals_out[i] <- if (length(z) == 0) NA_character_ else z[1]
    }
    sample_meta[[v]] <- vals_out
    if (!is.null(group_order) && v %in% names(group_order)) {
      sample_meta[[v]] <- factor(sample_meta[[v]], levels = levels(md[[v]]))
    }
  }

  obj <- seurat
  if (inherits(obj[[assay]], "Assay5")) {
    assay_layers <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character())
    matched_layers <- assay_layers[grepl(paste0("^", layer, "($|\\.)"), assay_layers)]
    if (length(matched_layers) > 1) {
      msg("Joining %d Seurat v5 `%s` layers before downstream analysis.",
          length(matched_layers), layer)
      obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]])
    }
  }

  expr <- get_expr(obj, assay, layer)
  if (ncol(expr) != ncol(obj)) stop("Expression matrix and Seurat cells are not aligned.")

  SeuratObject::Idents(obj) <- factor(as.character(obj[[celltype]][, 1]), levels = celltype_levels)

  output_dir <- ensure_dir(output_dir)
  dirs <- list(
    composition = ensure_dir(file.path(output_dir, "01_composition")),
    dimplot = ensure_dir(file.path(output_dir, "02_dimplot")),
    markers = ensure_dir(file.path(output_dir, "03_markers")),
    marker_plots = ensure_dir(file.path(output_dir, "04_marker_plots")),
    correlation = ensure_dir(file.path(output_dir, "05_expression_correlation")),
    sample_cor = ensure_dir(file.path(output_dir, "06_sample_celltype_correlation")),
    qc = ensure_dir(file.path(output_dir, "07_diversity_QC")),
    runinfo = ensure_dir(file.path(output_dir, "00_run_info"))
  )

  prefix <- clean_filename(file_prefix)
  withr::local_seed(seed)

  celltype_palette <- discrete_palette(celltype_levels)
  split_palettes <- lapply(split_vars, function(v) discrete_palette(md[[v]]))
  names(split_palettes) <- split_vars

  utils::capture.output(utils::sessionInfo(), file = file.path(dirs$runinfo, paste0(prefix, "_sessionInfo.txt")))
  writeLines(
    c(
      paste0("Date: ", Sys.time()),
      paste0("Cells: ", ncol(obj)),
      paste0("Features: ", nrow(obj)),
      paste0("Assay: ", assay),
      paste0("Layer: ", layer),
      paste0("Reduction: ", reduction),
      paste0("Celltype column: ", celltype),
      paste0("Sample column: ", sample_col),
      paste0("Clinical variables: ", paste(split_vars, collapse = ", ")),
      paste0("Statistic level(s): ", paste(statistic_level, collapse = ", ")),
      paste0("Custom colors available: ", if (is.null(custom_colors)) 0 else length(custom_colors)),
      paste0("Correlation: ", correlation_method, " / features=", correlation_features)
    ),
    con = file.path(dirs$runinfo, paste0(prefix, "_run_summary.txt"))
  )

  msg("Starting analysis: %d cells, %d cell types, %d samples.",
      ncol(obj), length(celltype_levels), length(samples))

  comp_cell_overall <- data.frame(
    celltype = factor(as.character(md[[celltype]]), levels = celltype_levels),
    stringsAsFactors = FALSE
  ) |>
    dplyr::count(.data$celltype, name = "n_cells", .drop = FALSE) |>
    dplyr::mutate(
      fraction = .data$n_cells / sum(.data$n_cells),
      percent = 100 * .data$fraction
    )
  names(comp_cell_overall)[1] <- celltype

  sample_vec <- factor(as.character(md[[sample_col]]), levels = samples)
  ct_vec <- factor(as.character(md[[celltype]]), levels = celltype_levels)
  sample_comp <- as.data.frame(table(sample_vec, ct_vec), stringsAsFactors = FALSE)
  names(sample_comp) <- c(sample_col, celltype, "n_cells")
  sample_comp[[sample_col]] <- as.character(sample_comp[[sample_col]])
  sample_comp[[celltype]] <- as.character(sample_comp[[celltype]])
  sample_totals <- stats::aggregate(sample_comp$n_cells,
                             by = list(sample_comp[[sample_col]]), sum)
  names(sample_totals) <- c(sample_col, "sample_total_cells")
  sample_comp <- merge(sample_comp, sample_totals, by = sample_col, all.x = TRUE)
  sample_comp$fraction <- sample_comp$n_cells / sample_comp$sample_total_cells
  sample_comp$percent <- 100 * sample_comp$fraction
  sample_comp <- merge(sample_comp, sample_meta, by = sample_col, all.x = TRUE, sort = FALSE)

  sample_count_wide <- sample_comp |>
    dplyr::select(dplyr::all_of(c(sample_col, celltype, "n_cells"))) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(celltype), values_from = "n_cells", values_fill = 0)
  sample_fraction_wide <- sample_comp |>
    dplyr::select(dplyr::all_of(c(sample_col, celltype, "fraction"))) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(celltype), values_from = "fraction", values_fill = 0)

  write_xlsx_safe(
    list(
      overall_celltype = comp_cell_overall,
      sample_celltype_long = sample_comp,
      sample_counts_wide = sample_count_wide,
      sample_fraction_wide = sample_fraction_wide,
      sample_metadata = sample_meta
    ),
    file.path(dirs$composition, paste0(prefix, "_celltype_composition.xlsx"))
  )

  low_count <- sample_comp[sample_comp$n_cells < min_cells_per_sample_celltype, , drop = FALSE]
  write_xlsx_safe(
    list(
      low_count_sample_celltype = low_count,
      sample_total_cells = sample_totals
    ),
    file.path(dirs$qc, paste0(prefix, "_low_cell_count_QC.xlsx"))
  )

  p_overall <- ggplot2::ggplot(
    comp_cell_overall,
    ggplot2::aes(x = .data[[celltype]], y = .data$percent, fill = .data[[celltype]])
  ) +
    ggplot2::geom_col(width = 0.8, show.legend = FALSE) +
    ggplot2::labs(x = celltype, y = "Cells (%)", title = "Overall cell-type composition") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  p_overall <- add_manual_fill(p_overall, celltype_palette)
  d <- auto_dim("sample_bar", n1 = length(celltype_levels), n2 = length(celltype_levels))
  save_plot(p_overall, file.path(dirs$composition, paste0(prefix, "_composition_overall")), d[1], d[2])

  p_sample <- ggplot2::ggplot(
    sample_comp,
    ggplot2::aes(x = .data[[sample_col]], y = .data$fraction, fill = .data[[celltype]])
  ) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(x = sample_col, y = "Cell fraction", fill = celltype,
                  title = "Cell-type composition by sample") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  p_sample <- add_manual_fill(p_sample, celltype_palette)
  d <- auto_dim("sample_bar", n1 = length(samples), n2 = length(celltype_levels))
  save_plot(p_sample, file.path(dirs$composition, paste0(prefix, "_composition_by_sample")), d[1], d[2])

  clinical_composition_tables <- list()
  for (v in split_vars) {
    v_clean <- clean_filename(v)
    dfv <- data.frame(
      group = md[[v]],
      celltype_tmp = factor(as.character(md[[celltype]]), levels = celltype_levels),
      stringsAsFactors = FALSE
    )
    dfv <- dfv[!is.na(dfv$group), , drop = FALSE]
    pooled <- dfv |>
      dplyr::count(.data$group, .data$celltype_tmp, name = "n_cells", .drop = FALSE) |>
      dplyr::group_by(.data$group) |>
      dplyr::mutate(
        fraction = .data$n_cells / sum(.data$n_cells),
        percent = 100 * .data$fraction
      ) |>
      dplyr::ungroup()
    names(pooled)[2] <- celltype

    clinical_composition_tables[[paste0(v, "_pooled")]] <- pooled

    p_pool <- ggplot2::ggplot(
      pooled,
      ggplot2::aes(x = .data$group, y = .data$fraction, fill = .data[[celltype]])
    ) +
      ggplot2::geom_col(width = 0.85) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
      ggplot2::labs(x = v, y = "Cell fraction", fill = celltype,
                    title = paste0("Pooled-cell composition by ", v)) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    p_pool <- add_manual_fill(p_pool, celltype_palette)
    ng <- length(unique(pooled$group))
    d <- auto_dim("sample_bar", n1 = ng, n2 = length(celltype_levels))
    save_plot(p_pool, file.path(dirs$composition, paste0(prefix, "_composition_pooled_by_", v_clean)), d[1], d[2])

    if (v != sample_col && v %in% colnames(sample_comp)) {
      tmp <- sample_comp[!is.na(sample_comp[[v]]), , drop = FALSE]
      mean_comp <- tmp |>
        dplyr::group_by(.data[[v]], .data[[celltype]]) |>
        dplyr::summarise(
          n_samples = dplyr::n(),
          mean_fraction = mean(.data$fraction, na.rm = TRUE),
          median_fraction = stats::median(.data$fraction, na.rm = TRUE),
          sd_fraction = stats::sd(.data$fraction, na.rm = TRUE),
          sem_fraction = .data$sd_fraction / sqrt(.data$n_samples),
          .groups = "drop"
        )
      clinical_composition_tables[[paste0(v, "_sample_summary")]] <- mean_comp

      p_mean <- ggplot2::ggplot(
        mean_comp,
        ggplot2::aes(x = .data[[v]], y = .data$mean_fraction, fill = .data[[celltype]])
      ) +
        ggplot2::geom_col(width = 0.85) +
        ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
        ggplot2::labs(x = v, y = "Mean sample-level fraction", fill = celltype,
                      title = paste0("Mean sample composition by ", v)) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      p_mean <- add_manual_fill(p_mean, celltype_palette)
      ng <- length(unique(tmp[[v]]))
      d <- auto_dim("sample_bar", n1 = ng, n2 = length(celltype_levels))
      save_plot(p_mean, file.path(dirs$composition, paste0(prefix, "_composition_sample_mean_by_", v_clean)), d[1], d[2])

      p_box <- ggplot2::ggplot(
        tmp,
        ggplot2::aes(x = .data[[v]], y = .data$fraction, fill = .data[[v]])
      ) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75) +
        ggplot2::geom_jitter(width = 0.15, size = 1.4, alpha = 0.75) +
        ggplot2::facet_wrap(stats::as.formula(paste("~", celltype)), scales = "free_y") +
        ggplot2::labs(x = v, y = "Sample-level cell fraction",
                      title = paste0("Sample-level cell fractions by ", v)) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
      p_box <- add_manual_fill(p_box, split_palettes[[v]])
      nct <- length(celltype_levels)
      d <- c(clamp(4 * min(4, ceiling(sqrt(nct))), 12, 24),
             clamp(3.2 * ceiling(nct / min(4, ceiling(sqrt(nct)))), 7, 26))
      if (!is.null(width)) d[1] <- width
      if (!is.null(height)) d[2] <- height
      save_plot(p_box, file.path(dirs$composition, paste0(prefix, "_composition_sample_boxplot_by_", v_clean)), d[1], d[2])
    }
  }

  if (length(clinical_composition_tables) > 0) {
    wb_list <- list()
    used <- character()
    for (nm in names(clinical_composition_tables)) {
      sh <- clean_sheet(nm, used)
      used <- c(used, sh)
      wb_list[[sh]] <- clinical_composition_tables[[nm]]
    }
    write_xlsx_safe(
      wb_list,
      file.path(dirs$composition, paste0(prefix, "_composition_by_clinical_variables.xlsx"))
    )
  }

  msg("Generating balloon plots for cell counts...")
  for (v in split_vars) {
    v_clean <- clean_filename(v)
    vv <- as.character(md[[v]])
    cc <- factor(as.character(md[[celltype]]), levels = celltype_levels)
    keep <- !is.na(vv) & nzchar(vv) & !is.na(cc)
    if (!any(keep)) next
    tb <- table(vv[keep], cc[keep])
    if (nrow(tb) == 0 || ncol(tb) == 0) next
    save_balloonplot(
      tb,
      file.path(dirs$composition, paste0(prefix, "_balloonplot_", v_clean, "_by_", clean_filename(celltype))),
      title = paste0("Cell counts: ", v, " \u00d7 ", celltype),
      xlab = celltype,
      ylab = v
    )
  }

  sample_test_results <- list()
  pairwise_test_results <- list()
  cell_test_results <- list()

  if ("samples" %in% statistic_level) {
    for (v in setdiff(split_vars, sample_col)) {
      tmp <- sample_comp[!is.na(sample_comp[[v]]), , drop = FALSE]
      groups <- unique(as.character(tmp[[v]]))
      groups <- groups[!is.na(groups)]
      if (length(groups) < 2) next

      res_v <- list()
      pair_v <- list()
      for (ct in celltype_levels) {
        dct <- tmp[tmp[[celltype]] == ct, , drop = FALSE]
        dct <- dct[!is.na(dct[[v]]), , drop = FALSE]
        g <- as.character(dct[[v]])
        lev <- if (!is.null(group_order) && v %in% names(group_order)) {
          intersect(group_order[[v]], unique(g))
        } else {
          unique(g)
        }
        if (length(lev) < 2) next

        medians <- tapply(dct$fraction, g, stats::median, na.rm = TRUE)
        means <- tapply(dct$fraction, g, mean, na.rm = TRUE)
        ns <- table(g)

        if (length(lev) == 2) {
          x <- dct$fraction[g == lev[1]]
          y <- dct$fraction[g == lev[2]]
          p <- if (length(x) >= 2 && length(y) >= 2) {
            tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value,
                     error = function(e) NA_real_)
          } else NA_real_
          res_v[[ct]] <- data.frame(
            clinical_var = v,
            celltype = ct,
            test = "Wilcoxon rank-sum",
            group1 = lev[1], group2 = lev[2],
            n_group1 = unname(ns[lev[1]]), n_group2 = unname(ns[lev[2]]),
            median_group1 = unname(medians[lev[1]]),
            median_group2 = unname(medians[lev[2]]),
            mean_group1 = unname(means[lev[1]]),
            mean_group2 = unname(means[lev[2]]),
            delta_median_group2_minus_group1 = unname(medians[lev[2]] - medians[lev[1]]),
            p_value = p,
            stringsAsFactors = FALSE
          )
        } else {
          p <- tryCatch(
            stats::kruskal.test(dct$fraction, as.factor(g))$p.value,
            error = function(e) NA_real_
          )
          res_v[[ct]] <- data.frame(
            clinical_var = v,
            celltype = ct,
            test = "Kruskal-Wallis",
            group1 = NA_character_, group2 = NA_character_,
            n_group1 = NA_integer_, n_group2 = NA_integer_,
            median_group1 = NA_real_, median_group2 = NA_real_,
            mean_group1 = NA_real_, mean_group2 = NA_real_,
            delta_median_group2_minus_group1 = NA_real_,
            p_value = p,
            stringsAsFactors = FALSE
          )

          cmb <- utils::combn(lev, 2, simplify = FALSE)
          pair_rows <- lapply(cmb, function(z) {
            x <- dct$fraction[g == z[1]]
            y <- dct$fraction[g == z[2]]
            pv <- if (length(x) >= 2 && length(y) >= 2) {
              tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value,
                       error = function(e) NA_real_)
            } else NA_real_
            data.frame(
              clinical_var = v, celltype = ct,
              group1 = z[1], group2 = z[2],
              n_group1 = length(x), n_group2 = length(y),
              median_group1 = stats::median(x, na.rm = TRUE),
              median_group2 = stats::median(y, na.rm = TRUE),
              delta_median_group2_minus_group1 = stats::median(y, na.rm = TRUE) - stats::median(x, na.rm = TRUE),
              p_value = pv,
              stringsAsFactors = FALSE
            )
          })
          pair_v[[ct]] <- dplyr::bind_rows(pair_rows)
        }
      }
      res_v <- dplyr::bind_rows(res_v)
      if (nrow(res_v) > 0) {
        res_v$p_adj_BH <- stats::p.adjust(res_v$p_value, method = "BH")
        sample_test_results[[v]] <- res_v
      }
      pair_v <- dplyr::bind_rows(pair_v)
      if (nrow(pair_v) > 0) {
        pair_v$p_adj_BH <- stats::p.adjust(pair_v$p_value, method = "BH")
        pairwise_test_results[[v]] <- pair_v
      }
    }
  }

  if ("cell" %in% statistic_level) {
    for (v in setdiff(split_vars, sample_col)) {
      z <- data.frame(group = md[[v]], ct = as.character(md[[celltype]]), stringsAsFactors = FALSE)
      z <- z[!is.na(z$group), , drop = FALSE]
      if (length(unique(z$group)) < 2) next
      rows <- list()
      for (ct in celltype_levels) {
        tab <- table(z$group, z$ct == ct)
        if (nrow(tab) < 2 || ncol(tab) < 2) next
        chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
        use_fisher <- any(chi$expected < 5) && prod(dim(tab)) <= 20
        if (use_fisher) {
          tst <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
          p <- if (is.null(tst)) chi$p.value else tst$p.value
          test_name <- if (is.null(tst)) "Chi-square" else "Fisher exact"
        } else {
          p <- chi$p.value
          test_name <- "Chi-square"
        }
        rows[[ct]] <- data.frame(
          clinical_var = v,
          celltype = ct,
          test = test_name,
          p_value = p,
          note = "Cell-level exploratory test; cells are not independent biological replicates.",
          stringsAsFactors = FALSE
        )
      }
      rows <- dplyr::bind_rows(rows)
      if (nrow(rows) > 0) {
        rows$p_adj_BH <- stats::p.adjust(rows$p_value, method = "BH")
        cell_test_results[[v]] <- rows
      }
    }
  }

  stat_wb <- list()
  if (length(sample_test_results)) {
    stat_wb[["sample_level_global"]] <- dplyr::bind_rows(sample_test_results)
  }
  if (length(pairwise_test_results)) {
    stat_wb[["sample_level_pairwise"]] <- dplyr::bind_rows(pairwise_test_results)
  }
  if (length(cell_test_results)) {
    stat_wb[["cell_level_exploratory"]] <- dplyr::bind_rows(cell_test_results)
  }
  if (length(stat_wb) > 0) {
    write_xlsx_safe(
      stat_wb,
      file.path(dirs$composition, paste0(prefix, "_celltype_composition_statistics.xlsx"))
    )
  }

  msg("Generating %s plots...", reduction)

  p_dim_ct <- Seurat::DimPlot(
    obj, reduction = reduction, group.by = celltype,
    cols = celltype_palette,
    label = length(celltype_levels) <= 30, repel = TRUE
  ) + ggplot2::ggtitle(paste0(toupper(reduction), " by ", celltype))
  save_plot(p_dim_ct, file.path(dirs$dimplot, paste0(prefix, "_", reduction, "_celltype")),
            ifelse(is.null(width), 10, width), ifelse(is.null(height), 8, height))

  for (v in split_vars) {
    v_clean <- clean_filename(v)
    vals <- unique(as.character(md[[v]]))
    vals <- vals[!is.na(vals)]

    p_group <- Seurat::DimPlot(
      obj, reduction = reduction, group.by = v,
      cols = split_palettes[[v]],
      label = length(vals) <= 20, repel = TRUE
    ) + ggplot2::ggtitle(paste0(toupper(reduction), " grouped by ", v))
    d <- c(clamp(8 + 0.15 * length(vals), 9, 16), 8)
    if (!is.null(width)) d[1] <- width
    if (!is.null(height)) d[2] <- height
    save_plot(p_group, file.path(dirs$dimplot, paste0(prefix, "_", reduction, "_groupby_", v_clean)), d[1], d[2])

    pages <- split(vals, ceiling(seq_along(vals) / max_panels_per_page))
    for (pi in seq_along(pages)) {
      keep_vals <- pages[[pi]]
      keep_cells <- rownames(md)[as.character(md[[v]]) %in% keep_vals]
      n_pan <- length(keep_vals)
      p_split <- Seurat::DimPlot(
        obj,
        cells = keep_cells,
        reduction = reduction,
        group.by = celltype,
        split.by = v,
        cols = celltype_palette,
        ncol = min(4, n_pan),
        label = FALSE,
        raster = NULL
      )
      d <- auto_dim("dim_split", n1 = n_pan)
      suffix <- if (length(pages) > 1) paste0("_page", pi) else ""
      save_plot(
        p_split,
        file.path(dirs$dimplot, paste0(prefix, "_", reduction, "_celltype_splitby_", v_clean, suffix)),
        d[1], d[2]
      )
    }
  }

  markers <- NULL
  top50 <- NULL
  top_plot_list <- NULL

  if (isTRUE(run_markers)) {
    msg("Running FindAllMarkers on `%s`...", celltype)
    markers <- Seurat::FindAllMarkers(
      object = obj,
      assay = assay,
      group.by = celltype,
      only.pos = marker_only_pos,
      logfc.threshold = marker_logfc_threshold,
      min.pct = marker_min_pct,
      test.use = marker_test,
      slot = layer,
      return.thresh = 1,
      verbose = verbose,
      random.seed = seed
    )

    if (nrow(markers) > 0) {
      write_rds_csv(
        markers,
        file.path(dirs$markers, paste0(prefix, "_FindAllMarkers_all.rds")),
        file.path(dirs$markers, paste0(prefix, "_FindAllMarkers_all.csv"))
      )

      fc_candidates <- c("avg_log2FC", "avg_logFC", "avg_diff")
      fc_col <- fc_candidates[fc_candidates %in% colnames(markers)][1]
      if (is.na(fc_col) || length(fc_col) == 0) {
        stop("Cannot identify fold-change column in FindAllMarkers output.")
      }

      m2 <- markers
      if ("p_val_adj" %in% colnames(m2)) {
        m_sig <- m2[!is.na(m2$p_val_adj) & m2$p_val_adj <= marker_padj_cutoff, , drop = FALSE]
        if (nrow(m_sig) == 0) {
          warning("No markers pass p_val_adj cutoff; ranking all returned markers instead.")
          m_sig <- m2
        }
      } else {
        m_sig <- m2
      }

      if (any(m_sig[[fc_col]] > 0, na.rm = TRUE)) {
        m_pos <- m_sig[m_sig[[fc_col]] > 0 & !is.na(m_sig[[fc_col]]), , drop = FALSE]
      } else {
        m_pos <- m_sig
      }

      top50 <- m_pos |>
        dplyr::group_by(.data$cluster) |>
        dplyr::slice_max(order_by = .data[[fc_col]], n = top_n_marker, with_ties = FALSE) |>
        dplyr::ungroup()

      marker_wb <- list(Combined = top50)
      used <- "Combined"
      for (ct in unique(as.character(top50$cluster))) {
        sh <- clean_sheet(ct, used)
        used <- c(used, sh)
        marker_wb[[sh]] <- top50[as.character(top50$cluster) == ct, , drop = FALSE]
      }
      write_xlsx_safe(
        marker_wb,
        file.path(dirs$markers, paste0(prefix, "_top", top_n_marker, "_markers_by_celltype.xlsx"))
      )

      top_plot <- top50 |>
        dplyr::group_by(.data$cluster) |>
        dplyr::slice_max(order_by = .data[[fc_col]], n = top_n_plot, with_ties = FALSE) |>
        dplyr::ungroup()

      top_plot_list <- split(top_plot$gene, as.character(top_plot$cluster))
      top_plot_list <- lapply(top_plot_list, unique)
      top_genes <- unique(unlist(top_plot_list, use.names = FALSE))
      top_genes <- intersect(top_genes, rownames(obj[[assay]]))

      if (length(top_genes) > 0) {
        dp <- Seurat::DotPlot(
          obj,
          features = top_genes,
          assay = assay,
          group.by = celltype
        ) +
          ggplot2::coord_flip() +
          ggplot2::labs(title = paste0("Top ", top_n_plot, " markers per cell type")) +
          ggplot2::theme_bw(base_size = 10)
        d <- auto_dim("dot", n1 = length(top_genes), n2 = length(celltype_levels))
        save_plot(dp, file.path(dirs$marker_plots, paste0(prefix, "_top", top_n_plot, "_DotPlot")), d[1], d[2])

        obj_heat <- tryCatch(
          Seurat::ScaleData(obj, assay = assay, features = top_genes, verbose = FALSE),
          error = function(e) {
            warning("ScaleData failed for marker heatmap: ", conditionMessage(e))
            NULL
          }
        )
        if (!is.null(obj_heat)) {
          cells_by_ct <- split(colnames(obj_heat), as.character(obj_heat[[celltype]][, 1]))
          per_ct <- max(10, floor(max_cells_heatmap / max(1, length(cells_by_ct))))
          heat_cells <- unlist(lapply(cells_by_ct, function(z) {
            if (length(z) <= per_ct) z else sample(z, per_ct)
          }), use.names = FALSE)
          hp <- Seurat::DoHeatmap(
            obj_heat,
            features = top_genes,
            cells = heat_cells,
            group.by = celltype,
            group.colors = celltype_palette,
            assay = assay,
            slot = "scale.data",
            raster = TRUE
          ) + ggplot2::ggtitle(paste0("Top ", top_n_plot, " markers per cell type"))
          d <- auto_dim("heatmap", n1 = length(celltype_levels), n2 = length(top_genes))
          save_plot(hp, file.path(dirs$marker_plots, paste0(prefix, "_top", top_n_plot, "_DoHeatmap")), d[1], d[2])
        }

        fp_dir <- ensure_dir(file.path(dirs$marker_plots, "FeaturePlot_top_markers"))
        for (ct in names(top_plot_list)) {
          gs <- intersect(top_plot_list[[ct]], rownames(obj[[assay]]))
          if (length(gs) == 0) next
          fp <- Seurat::FeaturePlot(
            obj,
            features = gs,
            reduction = reduction,
            ncol = min(3, length(gs)),
            order = TRUE,
            raster = NULL
          ) + patchwork::plot_annotation(title = paste0(ct, ": top marker genes"))
          nr <- ceiling(length(gs) / min(3, length(gs)))
          fw <- ifelse(is.null(width), 4.2 * min(3, length(gs)), width)
          fh <- ifelse(is.null(height), 4.0 * nr, height)
          save_plot(fp, file.path(fp_dir, paste0(prefix, "_", clean_filename(ct), "_FeaturePlot")), fw, fh)
        }
      }
    } else {
      warning("FindAllMarkers returned zero rows.")
    }
  }

  avg_celltype <- NULL
  celltype_cor <- NULL

  if (isTRUE(run_celltype_correlation)) {
    msg("Calculating average expression and cell-type correlation...")
    avg_celltype <- avg_by_group(expr, as.character(md[[celltype]]), min_cells = 1)
    if (!is.null(avg_celltype)) {
      write_rds_csv(
        avg_celltype,
        file.path(dirs$correlation, paste0(prefix, "_average_expression_by_celltype.rds")),
        file.path(dirs$correlation, paste0(prefix, "_average_expression_by_celltype.csv")),
        rowname_col = "gene"
      )

      cor_feats <- choose_cor_features(avg_celltype, obj, assay, correlation_features)
      if (length(cor_feats) >= 2 && ncol(avg_celltype) >= 2) {
        celltype_cor <- safe_cor(avg_celltype[cor_feats, , drop = FALSE], correlation_method)
        if (!is.null(celltype_cor)) {
          write_xlsx_safe(
            list(correlation = data.frame(celltype = rownames(celltype_cor), celltype_cor, check.names = FALSE)),
            file.path(dirs$correlation, paste0(prefix, "_celltype_expression_correlation.xlsx"))
          )
          d <- auto_dim("heatmap", n1 = ncol(celltype_cor), n2 = nrow(celltype_cor))
          save_pheatmap(
            celltype_cor,
            file.path(dirs$correlation, paste0(prefix, "_celltype_expression_correlation_heatmap")),
            w = d[1], h = d[2],
            main = paste0("Cell-type expression correlation (", correlation_method, ")")
          )
        }
      }
    }
  }

  composition_cor <- NULL
  if (isTRUE(run_composition_correlation) && nrow(sample_fraction_wide) >= 3) {
    mat_comp <- as.matrix(sample_fraction_wide[, setdiff(colnames(sample_fraction_wide), sample_col), drop = FALSE])
    rownames(mat_comp) <- sample_fraction_wide[[sample_col]]
    keep <- apply(mat_comp, 2, stats::sd, na.rm = TRUE) > 0
    mat_comp <- mat_comp[, keep, drop = FALSE]
    if (ncol(mat_comp) >= 2) {
      composition_cor <- stats::cor(mat_comp, method = correlation_method, use = "pairwise.complete.obs")
      write_xlsx_safe(
        list(correlation = data.frame(celltype = rownames(composition_cor), composition_cor, check.names = FALSE)),
        file.path(dirs$correlation, paste0(prefix, "_sample_composition_celltype_correlation.xlsx"))
      )
      d <- auto_dim("heatmap", n1 = ncol(composition_cor), n2 = nrow(composition_cor))
      save_pheatmap(
        composition_cor,
        file.path(dirs$correlation, paste0(prefix, "_sample_composition_celltype_correlation_heatmap")),
        w = d[1], h = d[2],
        main = paste0("Cell-type co-variation across samples (", correlation_method, ")")
      )
    }
  }

  sample_celltype_cor_list <- list()
  sample_celltype_avg_paths <- list()

  if (isTRUE(run_sample_celltype_correlation)) {
    msg("Calculating sample-by-sample expression correlation within each cell type...")
    avg_dir <- ensure_dir(file.path(dirs$sample_cor, "average_expression_matrices"))
    heat_dir <- ensure_dir(file.path(dirs$sample_cor, "correlation_heatmaps"))
    xlsx_dir <- ensure_dir(file.path(dirs$sample_cor, "correlation_tables"))

    for (ct in celltype_levels) {
      ct_cells_idx <- which(as.character(md[[celltype]]) == ct)
      if (length(ct_cells_idx) == 0) next
      g <- as.character(md[[sample_col]][ct_cells_idx])
      ct_expr <- expr[, ct_cells_idx, drop = FALSE]
      avg_s <- avg_by_group(ct_expr, g, min_cells = min_cells_per_sample_celltype)
      if (is.null(avg_s) || ncol(avg_s) < 2) next

      ct_clean <- clean_filename(ct)
      rds_path <- file.path(avg_dir, paste0(prefix, "_", ct_clean, "_average_expression_by_sample.rds"))
      csv_path <- file.path(avg_dir, paste0(prefix, "_", ct_clean, "_average_expression_by_sample.csv"))
      write_rds_csv(
        avg_s,
        rds_path,
        csv_path,
        rowname_col = "gene"
      )
      sample_celltype_avg_paths[[ct]] <- c(rds = rds_path, csv = csv_path)

      cor_feats <- choose_cor_features(avg_s, obj, assay, correlation_features)
      if (length(cor_feats) < 2) next
      cs <- safe_cor(avg_s[cor_feats, , drop = FALSE], correlation_method)
      if (is.null(cs)) next
      sample_celltype_cor_list[[ct]] <- cs

      write_xlsx_safe(
        list(correlation = data.frame(sample = rownames(cs), cs, check.names = FALSE)),
        file.path(xlsx_dir, paste0(prefix, "_", ct_clean, "_sample_expression_correlation.xlsx"))
      )

      ann_vars <- setdiff(colnames(sample_meta), sample_col)
      ann <- NULL
      if (length(ann_vars) > 0) {
        ann <- sample_meta[match(colnames(cs), sample_meta[[sample_col]]), ann_vars, drop = FALSE]
        rownames(ann) <- colnames(cs)
        ann[] <- lapply(ann, function(z) factor(as.character(z)))
      }
      d <- auto_dim("heatmap", n1 = ncol(cs), n2 = nrow(cs))
      save_pheatmap(
        cs,
        file.path(heat_dir, paste0(prefix, "_", ct_clean, "_sample_expression_correlation_heatmap")),
        w = d[1], h = d[2],
        annotation_row = ann,
        annotation_col = ann,
        main = paste0(ct, ": sample expression correlation")
      )
    }
  }

  diversity <- NULL
  diversity_tests <- list()
  if (isTRUE(run_diversity)) {
    diversity <- sample_comp |>
      dplyr::group_by(.data[[sample_col]]) |>
      dplyr::summarise(
        total_cells = max(.data$sample_total_cells),
        detected_celltypes = sum(.data$n_cells > 0),
        shannon = -sum(ifelse(.data$fraction > 0, .data$fraction * log(.data$fraction), 0)),
        simpson = 1 - sum(.data$fraction^2),
        .groups = "drop"
      ) |>
      dplyr::left_join(sample_meta, by = sample_col)

    for (v in setdiff(split_vars, sample_col)) {
      dv <- diversity[!is.na(diversity[[v]]), , drop = FALSE]
      ng <- length(unique(dv[[v]]))
      if (ng < 2) next
      for (metric in c("detected_celltypes", "shannon", "simpson")) {
        if (ng == 2) {
          lev <- unique(as.character(dv[[v]]))
          x <- dv[[metric]][as.character(dv[[v]]) == lev[1]]
          y <- dv[[metric]][as.character(dv[[v]]) == lev[2]]
          pv <- if (length(x) >= 2 && length(y) >= 2) {
            tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value,
                     error = function(e) NA_real_)
          } else NA_real_
          rr <- data.frame(clinical_var = v, metric = metric,
                           test = "Wilcoxon rank-sum", p_value = pv)
        } else {
          pv <- tryCatch(stats::kruskal.test(dv[[metric]], as.factor(dv[[v]]))$p.value,
                         error = function(e) NA_real_)
          rr <- data.frame(clinical_var = v, metric = metric,
                           test = "Kruskal-Wallis", p_value = pv)
        }
        diversity_tests[[paste(v, metric, sep = "__")]] <- rr
      }
    }
    diversity_test_df <- dplyr::bind_rows(diversity_tests)
    if (nrow(diversity_test_df) > 0) {
      diversity_test_df$p_adj_BH <- stats::p.adjust(diversity_test_df$p_value, method = "BH")
    }
    write_xlsx_safe(
      list(
        diversity_by_sample = diversity,
        diversity_tests = diversity_test_df
      ),
      file.path(dirs$qc, paste0(prefix, "_sample_celltype_diversity.xlsx"))
    )
  }

  result <- list(
    parameters = list(
      celltype = celltype,
      sample_col = sample_col,
      clinical_var = clinical_var,
      split_vars = split_vars,
      statistic_level = statistic_level,
      assay = assay,
      layer = layer,
      reduction = reduction,
      correlation_method = correlation_method,
      correlation_features = correlation_features,
      min_cells_per_sample_celltype = min_cells_per_sample_celltype,
      custom_colors = custom_colors
    ),
    sample_metadata = sample_meta,
    overall_composition = comp_cell_overall,
    sample_composition_long = sample_comp,
    sample_count_matrix = sample_count_wide,
    sample_fraction_matrix = sample_fraction_wide,
    sample_level_statistics = sample_test_results,
    sample_level_pairwise_statistics = pairwise_test_results,
    cell_level_statistics = cell_test_results,
    markers = markers,
    top_markers = top50,
    average_expression_by_celltype = avg_celltype,
    celltype_expression_correlation = celltype_cor,
    composition_correlation = composition_cor,
    sample_within_celltype_correlations = sample_celltype_cor_list,
    diversity = diversity,
    low_count_QC = low_count,
    output_dirs = dirs
  )

  if (export_analysis_data) {
    saveRDS(
      result,
      file.path(output_dir, paste0(prefix, "_analysis_results.rds"))
    )
  }

  msg("Analysis finished. Output: %s", normalizePath(output_dir, mustWork = FALSE))
  invisible(result)
}

