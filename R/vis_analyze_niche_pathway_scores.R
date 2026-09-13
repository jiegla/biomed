#' Compare spatial niche pathway scores across samples
#'
#' Spatial transcriptomics workflow adapted from contributed analysis code.
#' Sample-level summaries are used for group comparisons. Spot-level correlations
#' and Moran statistics describe spatial association and do not establish causality.
#' @param vis A spatial Seurat object; niche-score analysis also accepts a spot-indexed metadata data frame.
#' @param group_col Metadata grouping column. For neural analysis, NULL uses condition_col; this must be constant within each sample.
#' @param sample_col Metadata column identifying independent samples. Leading-zero character IDs are preserved.
#' @param score_dir Directory containing spot-by-pathway score RDS files.
#' @param selected_rds Score RDS filenames relative to score_dir. NULL discovers all files with an RDS extension. Numeric columns must contain pathway scores only.
#' @param outdir Output directory. NULL constructs a directory from the grouping column(s); these report workflows export files by default.
#' @param min_spots_per_sample_niche Minimum spots retained per sample and niche.
#' @param save_individual_png Save individual pathway PNG figures.
#' @param save_individual_pdf Save individual pathway PDF figures.
#' @param save_multipage_pdf Save all panels in a multipage PDF.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi Raster output resolution in dots per inch.
#' @return Invisibly, a list with sample_niche_score, group_summary, global_stats, pairwise_stats, pathway_index and pathway_manifest.
#' @details Existing report files are overwritten in the selected output directory.
#' Cell-type fractions use the complete denominator cell-type set. Niche-score
#' global tests use sample random intercepts; pairwise tests match samples across niches.
#' Neural tests use one summary per sample; region summaries are descriptive.
#' @export
vis_analyze_niche_pathway_scores <- function(
  vis,
  group_col,
  sample_col = "sample_id",
  score_dir = "gmt_pathway_scores",
  selected_rds = NULL,
  outdir = NULL,
  min_spots_per_sample_niche = 1,
  save_individual_png = TRUE,
  save_individual_pdf = TRUE,
  save_multipage_pdf = TRUE,
  width = 6,
  height = 5,
  dpi = 300
) {
  # Local bindings for columns evaluated in dplyr/ggplot data masks.
  FDR_all_tests <- global_FDR_all <- mean_score <- median_score <- n_spots <- niche <- pathway <- sample_id <- sd_score <- source_file <- NULL


  for (pkg in c("dplyr", "tidyr", "nlme")) .biomed_require(pkg)
  .biomed_vis_positive_integer(min_spots_per_sample_niche, "min_spots_per_sample_niche")
  md <- .biomed_vis_score_dataframe(vis)
  if (is.null(selected_rds)) selected_rds <- list.files(score_dir, pattern = "[.]rds$", ignore.case = TRUE)
  if (!length(selected_rds) || anyNA(selected_rds) || anyDuplicated(selected_rds)) stop("Supply distinct score RDS filenames in selected_rds.")
  vis_spots <- rownames(md)

  if (is.null(vis_spots)) stop("vis@meta.data has no rownames.")
  if (!group_col %in% colnames(md)) stop("Cannot find group column: ", group_col)
  if (!sample_col %in% colnames(md)) stop("Cannot find sample column: ", sample_col)

  rds_paths <- file.path(score_dir, selected_rds)
  missing_files <- rds_paths[!file.exists(rds_paths)]
  if (length(missing_files) > 0) {
    stop("Cannot find these RDS files:\n", paste(missing_files, collapse = "\n"))
  }

  if (is.null(outdir)) outdir <- paste0("Niche_pathway_score_", group_col)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "individual_plots"), recursive = TRUE, showWarnings = FALSE)

  group_levels <- if (is.factor(md[[group_col]])) levels(droplevels(md[[group_col]])) else unique(as.character(md[[group_col]]))
  group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]
  niche_colors <- stats::setNames(biomed_colors(length(group_levels)), group_levels)

  present_groups <- unique(as.character(md[[group_col]]))
  group_levels <- group_levels[group_levels %in% present_groups]
  if (length(group_levels) < 2) {
    stop("Fewer than two valid groups were found in ", group_col, ".")
  }
  niche_colors <- niche_colors[group_levels]

  niche_info <- data.frame(
    spot = vis_spots,
    sample_id = as.character(md[[sample_col]]),
    niche = as.character(md[[group_col]]),
    stringsAsFactors = FALSE,
    row.names = vis_spots
  )

  sample_score_list <- vector("list", length(selected_rds))
  pathway_manifest_list <- vector("list", length(selected_rds))

  for (i in seq_along(selected_rds)) {
    rds_file <- selected_rds[i]
    rds_path <- file.path(score_dir, rds_file)

    message("")
    message("========================================")
    message("Reading: ", rds_file)
    message("========================================")

    score_obj <- readRDS(rds_path)
    score_df <- .biomed_vis_score_dataframe(score_obj)

    if (is.null(rownames(score_df))) {
      stop("RDS has no rownames/spot names: ", rds_file)
    }

    common_spots <- intersect(vis_spots, rownames(score_df))
    if (length(common_spots) == 0) {
      stop("No spot names in ", rds_file, " match rownames(vis@meta.data).")
    }

    numeric_cols <- colnames(score_df)[vapply(score_df, is.numeric, logical(1))]
    if (length(numeric_cols) == 0) {
      stop("No numeric pathway score columns found in: ", rds_file)
    }

    score_df <- score_df[common_spots, numeric_cols, drop = FALSE]

    # UCell RDS created by return_type = "meta" should contain pathway scores only.
    # Warn if any numeric column contains non-finite values.
    nonfinite_n <- sum(!is.finite(as.matrix(score_df)))
    if (nonfinite_n > 0) {
      warning(rds_file, " contains ", nonfinite_n, " non-finite score values; treated as missing.")
      score_df[] <- lapply(score_df, function(x) { x[!is.finite(x)] <- NA_real_; x })
    }

    message("Matched spots: ", length(common_spots), " / ", length(vis_spots))
    message("Pathways: ", ncol(score_df))

    if (any(numeric_cols %in% c("sample_id", "niche", "n_spots"))) stop("Score columns collide with reserved metadata names.")
    meta_part <- niche_info[common_spots, c("sample_id", "niche"), drop = FALSE]
    tmp <- cbind(meta_part, score_df)
    tmp <- tmp[!is.na(tmp$sample_id) & !is.na(tmp$niche) & tmp$niche %in% group_levels, , drop = FALSE]

    if (nrow(tmp) == 0) {
      stop("No valid sample/niche spots remain for: ", rds_file)
    }

    # IMPORTANT:
    # Aggregate while data are still WIDE.
    # This avoids creating ~15 million spot \u00d7 pathway long-format rows.
    mean_wide <- tmp |>
      dplyr::group_by(sample_id, niche) |>
      dplyr::summarise(
        n_spots = dplyr::n(),
        dplyr::across(dplyr::all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )

    median_wide <- tmp |>
      dplyr::group_by(sample_id, niche) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(numeric_cols), ~ stats::median(.x, na.rm = TRUE)),
        .groups = "drop"
      )

    sd_wide <- tmp |>
      dplyr::group_by(sample_id, niche) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(numeric_cols), ~ stats::sd(.x, na.rm = TRUE)),
        .groups = "drop"
      )

    mean_long <- mean_wide |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(numeric_cols),
        names_to = "pathway",
        values_to = "mean_score"
      )

    median_long <- median_wide |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(numeric_cols),
        names_to = "pathway",
        values_to = "median_score"
      )

    sd_long <- sd_wide |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(numeric_cols),
        names_to = "pathway",
        values_to = "sd_score"
      )

    sample_score_list[[i]] <- mean_long |>
      dplyr::left_join(median_long, by = c("sample_id", "niche", "pathway")) |>
      dplyr::left_join(sd_long, by = c("sample_id", "niche", "pathway")) |>
      dplyr::mutate(source_file = rds_file) |>
      dplyr::select(source_file, pathway, sample_id, niche, n_spots, mean_score, median_score, sd_score) |>
      dplyr::filter(
        n_spots >= min_spots_per_sample_niche,
        is.finite(mean_score)
      )

    pathway_manifest_list[[i]] <- data.frame(
      source_file = rds_file,
      pathway = numeric_cols,
      stringsAsFactors = FALSE
    )
  }

  sample_niche_score <- dplyr::bind_rows(sample_score_list)
  pathway_manifest <- dplyr::bind_rows(pathway_manifest_list) |> dplyr::distinct()

  sample_niche_score$niche <- factor(sample_niche_score$niche, levels = group_levels)

  if (nrow(sample_niche_score) == 0) {
    stop("No sample \u00d7 niche pathway scores remain after filtering.")
  }

  # Group-level descriptive summary across samples.
  # Different names are intentionally used here so dplyr does not overwrite mean_score.
  group_summary <- sample_niche_score |>
    dplyr::group_by(source_file, pathway, niche) |>
    dplyr::summarise(
      n_samples = dplyr::n_distinct(sample_id),
      group_mean_score = mean(mean_score, na.rm = TRUE),
      group_sd_score = stats::sd(mean_score, na.rm = TRUE),
      group_median_score = stats::median(mean_score, na.rm = TRUE),
      Q1 = stats::quantile(mean_score, 0.25, na.rm = TRUE),
      Q3 = stats::quantile(mean_score, 0.75, na.rm = TRUE),
      min_score = min(mean_score, na.rm = TRUE),
      max_score = max(mean_score, na.rm = TRUE),
      .groups = "drop"
    )

  pathway_index <- sample_niche_score |>
    dplyr::distinct(source_file, pathway)

  # Global comparison:
  # mean pathway score ~ niche + random intercept for sample_id
  global_list <- vector("list", nrow(pathway_index))

  for (i in seq_len(nrow(pathway_index))) {
    src <- pathway_index$source_file[i]
    pw <- pathway_index$pathway[i]

    d <- sample_niche_score |>
      dplyr::filter(source_file == src, pathway == pw, is.finite(mean_score)) |>
      droplevels()

    if (dplyr::n_distinct(d$niche) < 2 || dplyr::n_distinct(d$sample_id) < 3) {
      global_list[[i]] <- data.frame(
        source_file = src,
        pathway = pw,
        n_samples = dplyr::n_distinct(d$sample_id),
        n_niches = dplyr::n_distinct(d$niche),
        global_p = NA_real_,
        model_message = "Too few samples or groups"
      )
      next
    }

    fit <- try(
      nlme::lme(
        fixed = mean_score ~ niche,
        random = ~1 | sample_id,
        data = d,
        method = "ML",
        na.action = stats::na.omit,
        control = nlme::lmeControl(returnObject = FALSE)
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      global_p <- NA_real_
    } else {
      a <- stats::anova(fit)
      global_p <- tryCatch(
        as.numeric(a["niche", "p-value"]),
        error = function(e) NA_real_
      )
    }

    global_list[[i]] <- data.frame(
      source_file = src,
      pathway = pw,
      n_samples = dplyr::n_distinct(d$sample_id),
      n_niches = dplyr::n_distinct(d$niche),
      global_p = global_p,
      model_message = if (inherits(fit, "try-error")) as.character(fit) else NA_character_
    )
  }

  global_stats <- dplyr::bind_rows(global_list) |>
    dplyr::mutate(
      global_FDR_all = stats::p.adjust(global_p, method = "BH")
    ) |>
    dplyr::group_by(source_file) |>
    dplyr::mutate(
      global_FDR_within_source = stats::p.adjust(global_p, method = "BH")
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      significance = dplyr::case_when(
        global_FDR_all < 0.001 ~ "***",
        global_FDR_all < 0.01 ~ "**",
        global_FDR_all < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) |>
    dplyr::arrange(global_FDR_all, global_p)

  # Pairwise paired Wilcoxon:
  # only samples represented in both niches contribute to each comparison.
  pairwise_list <- list()
  k <- 1L

  for (i in seq_len(nrow(pathway_index))) {
    src <- pathway_index$source_file[i]
    pw <- pathway_index$pathway[i]

    d <- sample_niche_score |>
      dplyr::filter(source_file == src, pathway == pw)

    groups_present <- group_levels[group_levels %in% as.character(unique(d$niche))]
    if (length(groups_present) < 2) next

    comparisons <- utils::combn(groups_present, 2, simplify = FALSE)

    for (pair in comparisons) {
      g1 <- pair[1]
      g2 <- pair[2]

      d1 <- d |>
        dplyr::filter(niche == g1) |>
        dplyr::select(sample_id, mean_score) |>
        dplyr::rename(score_g1 = mean_score)

      d2 <- d |>
        dplyr::filter(niche == g2) |>
        dplyr::select(sample_id, mean_score) |>
        dplyr::rename(score_g2 = mean_score)

      paired_data <- dplyr::inner_join(d1, d2, by = "sample_id")
      paired_data <- paired_data[is.finite(paired_data$score_g1) & is.finite(paired_data$score_g2), , drop = FALSE]
      n_pair <- nrow(paired_data)

      if (n_pair >= 3) {
        wt <- try(
          stats::wilcox.test(
            paired_data$score_g1,
            paired_data$score_g2,
            paired = TRUE,
            exact = FALSE
          ),
          silent = TRUE
        )
        p_value <- if (inherits(wt, "try-error")) NA_real_ else wt$p.value
      } else {
        p_value <- NA_real_
      }

      pairwise_list[[k]] <- data.frame(
        source_file = src,
        pathway = pw,
        group1 = g1,
        group2 = g2,
        n_paired_samples = n_pair,
        median_group1 = if (n_pair > 0) stats::median(paired_data$score_g1, na.rm = TRUE) else NA_real_,
        median_group2 = if (n_pair > 0) stats::median(paired_data$score_g2, na.rm = TRUE) else NA_real_,
        median_difference_group2_minus_group1 =
          if (n_pair > 0) stats::median(paired_data$score_g2 - paired_data$score_g1, na.rm = TRUE) else NA_real_,
        p_value = p_value,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  if (length(pairwise_list) > 0) {
    pairwise_stats <- dplyr::bind_rows(pairwise_list) |>
      dplyr::group_by(source_file, pathway) |>
      dplyr::mutate(
        FDR_within_pathway = stats::p.adjust(p_value, method = "BH")
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        FDR_all_tests = stats::p.adjust(p_value, method = "BH"),
        significance = dplyr::case_when(
          FDR_within_pathway < 0.001 ~ "***",
          FDR_within_pathway < 0.01 ~ "**",
          FDR_within_pathway < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      ) |>
      dplyr::arrange(FDR_all_tests, p_value)
  } else {
    pairwise_stats <- data.frame(
      source_file = character(),
      pathway = character(),
      group1 = character(),
      group2 = character(),
      n_paired_samples = integer(),
      median_group1 = numeric(),
      median_group2 = numeric(),
      median_difference_group2_minus_group1 = numeric(),
      p_value = numeric(),
      FDR_within_pathway = numeric(),
      FDR_all_tests = numeric(),
      significance = character(),
      stringsAsFactors = FALSE
    )
  }

  make_plot <- function(src, pw) {
    d <- sample_niche_score |>
      dplyr::filter(source_file == src, pathway == pw)

    stat_row <- global_stats |>
      dplyr::filter(source_file == src, pathway == pw)

    subtitle_text <- NULL
    if (nrow(stat_row) == 1 && !is.na(stat_row$global_FDR_all)) {
      subtitle_text <- paste0(
        "Mixed-effects global FDR = ",
        format(stat_row$global_FDR_all, digits = 3, scientific = TRUE),
        "  ", stat_row$significance
      )
    }

    ggplot2::ggplot(d, ggplot2::aes(x = niche, y = mean_score, fill = niche)) +
      ggplot2::geom_boxplot(
        width = 0.62,
        outlier.shape = NA,
        alpha = 0.65,
        linewidth = 0.55
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(color = niche),
        width = 0.12,
        height = 0,
        size = 2.2,
        alpha = 0.85
      ) +
      ggplot2::scale_fill_manual(values = niche_colors, drop = FALSE) +
      ggplot2::scale_color_manual(values = niche_colors, drop = FALSE) +
      ggplot2::labs(
        title = pw,
        subtitle = subtitle_text,
        x = NULL,
        y = "Mean UCell score"
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5),
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1),
        legend.position = "none"
      )
  }

  opened_device <- NULL
  on.exit(if (!is.null(opened_device) && opened_device %in% grDevices::dev.list()) grDevices::dev.off(opened_device), add = TRUE)
  source_files <- unique(pathway_index$source_file)

  for (src in source_files) {
    source_short <- tools::file_path_sans_ext(basename(src))
    source_dir <- file.path(outdir, "individual_plots", sanitize_filename(source_short))
    dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

    current_pathways <- pathway_index |>
      dplyr::filter(source_file == src) |>
      dplyr::pull(pathway)

    if (save_multipage_pdf) {
      grDevices::pdf(
        file.path(
          outdir,
          paste0(sanitize_filename(source_short), "_", group_col, "_all_pathways.pdf")
        ),
        width = width,
        height = height
      )
    }

    if (save_multipage_pdf) opened_device <- grDevices::dev.cur()
    for (j in seq_along(current_pathways)) {
      pw <- current_pathways[j]
      p <- make_plot(src, pw)

      if (save_multipage_pdf) print(p)

      safe_pw <- substr(sanitize_filename(pw), 1, 140)
      file_prefix <- paste0(sprintf("%03d", j), "_", safe_pw)

      if (save_individual_png) {
        ggplot2::ggsave(
          filename = file.path(source_dir, paste0(file_prefix, ".png")),
          plot = p,
          width = width,
          height = height,
          dpi = dpi,
          bg = "white"
        )
      }

      if (save_individual_pdf) {
        ggplot2::ggsave(
          filename = file.path(source_dir, paste0(file_prefix, ".pdf")),
          plot = p,
          width = width,
          height = height,
          bg = "white"
        )
      }
    }

    if (save_multipage_pdf) { grDevices::dev.off(opened_device); opened_device <- NULL }
  }

  utils::write.csv(
    sample_niche_score,
    file.path(outdir, "sample_niche_pathway_mean_score.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    group_summary,
    file.path(outdir, "pathway_group_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    global_stats,
    file.path(outdir, "pathway_global_statistics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    pairwise_stats,
    file.path(outdir, "pathway_pairwise_statistics.csv"),
    row.names = FALSE
  )

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()

    sheets <- list(
      Sample_niche_score = sample_niche_score,
      Group_summary = group_summary,
      Global_statistics = global_stats,
      Pairwise_statistics = pairwise_stats,
      Pathway_sources = pathway_manifest
    )

    for (nm in names(sheets)) {
      openxlsx::addWorksheet(wb, nm)
      openxlsx::writeData(wb, nm, sheets[[nm]], withFilter = TRUE)
      openxlsx::freezePane(wb, nm, firstRow = TRUE)

      if (ncol(sheets[[nm]]) > 0) {
        openxlsx::setColWidths(
          wb,
          nm,
          cols = seq_len(ncol(sheets[[nm]])),
          widths = "auto"
        )
      }
    }

    openxlsx::saveWorkbook(
      wb,
      file.path(outdir, paste0(group_col, "_pathway_score_analysis.xlsx")),
      overwrite = TRUE
    )
  }

  message("")
  message("========================================")
  message("Finished: ", group_col)
  message("RDS files: ", length(unique(sample_niche_score$source_file)))
  message("Pathways: ", nrow(pathway_index))
  message("Samples: ", dplyr::n_distinct(sample_niche_score$sample_id))
  message("Output: ", normalizePath(outdir, mustWork = FALSE))
  message("========================================")

  invisible(list(
    sample_niche_score = sample_niche_score,
    group_summary = group_summary,
    global_stats = global_stats,
    pairwise_stats = pairwise_stats,
    pathway_index = pathway_index,
    pathway_manifest = pathway_manifest
  ))
}
