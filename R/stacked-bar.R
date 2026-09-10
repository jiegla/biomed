#' Draw proportional stacked bar charts
#'
#' Creates one proportional stacked bar chart and one association test per
#' outcome variable. File output is opt-in through `save`.
#'
#' @param dataframe A data frame.
#' @param x_var Grouping column shown on the x-axis.
#' @param y_vars Categorical outcome columns.
#' @param out_dir Output directory.
#' @param color Optional colour vector, preferably named by outcome level.
#' @param test_method One of `"auto"`, `"chi_square"`, or `"fisher"`.
#' @param show_label,show_n,show_p Display slice labels, sample counts, and test
#'   results.
#' @param percent_accuracy Accuracy passed to [scales::percent()].
#' @param palette_name Retained for backward compatibility; the package palette
#'   is used when `color` is `NULL`.
#' @param width,height Figure dimensions in inches.
#' @param dpi PNG resolution.
#' @param horizontal Flip the coordinates.
#' @param remove_na Remove incomplete observations.
#' @param return_plot Return plots together with statistics.
#' @param theme_size Base theme font size.
#' @param save Save PDF and PNG plots plus an XLSX or CSV statistics file.
#' @param statistics_format Statistics file type when `save = TRUE`.
#'
#' @return A statistics data frame, or a list with plots and statistics when
#'   `return_plot = TRUE`.
#' @export
draw_stack_barplot <- function(
    dataframe,
    x_var,
    y_vars,
    out_dir = "StackBar",
    color = NULL,
    test_method = c("auto", "chi_square", "fisher"),
    show_label = TRUE,
    show_n = TRUE,
    show_p = TRUE,
    percent_accuracy = 0.1,
    palette_name = NULL,
    width = 4.5,
    height = 5.5,
    dpi = 600,
    horizontal = FALSE,
    remove_na = TRUE,
    return_plot = FALSE,
    theme_size = 16,
    save = FALSE,
    statistics_format = c("xlsx", "csv")) {
  test_method <- match.arg(test_method)
  statistics_format <- match.arg(statistics_format)
  dataframe <- as.data.frame(dataframe)
  if (length(x_var) != 1L || !is.character(x_var) ||
      !is.character(y_vars) || !length(y_vars)) {
    cli::cli_abort("{.arg x_var} must name one column and {.arg y_vars} must name one or more columns.")
  }
  .biomed_required_columns(dataframe, c(x_var, y_vars))
  if (!is.null(palette_name)) {
    cli::cli_warn("{.arg palette_name} is deprecated; use {.arg color} instead.")
  }
  if (isTRUE(save)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (statistics_format == "xlsx") {
      .biomed_require("openxlsx", "It is used to save the statistics workbook.")
    }
  }

  plots <- list()
  stats_rows <- list()
  for (y in y_vars) {
    dat <- data.frame(.x = dataframe[[x_var]], .y = dataframe[[y]])
    if (remove_na) {
      dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    } else {
      dat$.x <- as.character(dat$.x)
      dat$.y <- as.character(dat$.y)
      dat$.x[is.na(dat$.x)] <- "NA"
      dat$.y[is.na(dat$.y)] <- "NA"
    }
    dat$.x <- factor(dat$.x)
    dat$.y <- factor(dat$.y)
    dat$.x <- droplevels(dat$.x)
    dat$.y <- droplevels(dat$.y)
    tab <- table(dat$.x, dat$.y)
    if (any(dim(tab) < 2L)) {
      cli::cli_warn("{.val {y}} has fewer than two observed levels and was skipped.")
      next
    }

    chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
    use_fisher <- test_method == "fisher" ||
      (test_method == "auto" && all(dim(tab) == c(2L, 2L)) && any(chi$expected < 5))
    test <- if (use_fisher) stats::fisher.test(tab) else chi
    method <- if (use_fisher) "Fisher exact test" else "Chi-square test"
    p_value <- unname(test$p.value)
    stats_rows[[y]] <- data.frame(
      variable = y, method = method, p_value = p_value,
      significance = .biomed_significance(p_value), N = sum(tab),
      stringsAsFactors = FALSE
    )

    plot_data <- as.data.frame(tab, stringsAsFactors = FALSE)
    names(plot_data) <- c(".x", ".y", "n")
    totals <- stats::ave(plot_data$n, plot_data$.x, FUN = sum)
    plot_data$prop <- plot_data$n / totals
    plot_data$label <- paste0(
      plot_data$n, " (", scales::percent(plot_data$prop, accuracy = percent_accuracy), ")"
    )
    n_data <- unique(data.frame(.x = plot_data$.x, n = totals))
    n_data$label <- paste0("n=", n_data$n)
    fill_values <- .biomed_colour_values(levels(dat$.y), color)

    g <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data[[".x"]], y = .data[["prop"]], fill = .data[[".y"]])
    ) +
      ggplot2::geom_col(width = 0.65, colour = "black", linewidth = 0.6) +
      ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
      ggplot2::scale_y_continuous(
        labels = scales::label_percent(), breaks = seq(0, 1, 0.2),
        limits = c(0, 1.18), expand = c(0, 0)
      ) +
      ggplot2::labs(x = x_var, y = paste0(y, " proportion"), fill = y) +
      ggplot2::theme_classic(base_size = theme_size) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(colour = "black"),
        legend.position = "top"
      )
    if (show_label) {
      g <- g + ggplot2::geom_text(
        ggplot2::aes(label = .data[["label"]]),
        position = ggplot2::position_stack(vjust = 0.5), size = 4
      )
    }
    if (show_n) {
      g <- g + ggplot2::geom_text(
        data = n_data,
        ggplot2::aes(x = .data[[".x"]], y = 1.03, label = .data[["label"]]),
        inherit.aes = FALSE, size = 4
      )
    }
    if (show_p) {
      g <- g + ggplot2::annotate(
        "text", x = (nlevels(dat$.x) + 1) / 2, y = 1.12,
        label = paste0(method, "\nP=", signif(p_value, 3), " ",
                       .biomed_significance(p_value)), size = 4
      )
    }
    if (horizontal) g <- g + ggplot2::coord_flip()
    plots[[y]] <- g

    if (isTRUE(save)) {
      stem <- file.path(
        out_dir,
        paste0(.biomed_safe_name(x_var), "__", .biomed_safe_name(y), "_stack_bar")
      )
      ggplot2::ggsave(paste0(stem, ".pdf"), g, width = width, height = height, bg = "white")
      ggplot2::ggsave(
        paste0(stem, ".png"), g, width = width, height = height, dpi = dpi, bg = "white"
      )
    }
  }

  statistics <- if (length(stats_rows)) {
    do.call(rbind, unname(stats_rows))
  } else {
    data.frame()
  }
  rownames(statistics) <- NULL
  if (isTRUE(save) && nrow(statistics)) {
    path <- file.path(out_dir, paste0(.biomed_safe_name(x_var), "_statistics.", statistics_format))
    if (statistics_format == "xlsx") {
      openxlsx::write.xlsx(statistics, path, overwrite = TRUE)
    } else {
      utils::write.csv(statistics, path, row.names = FALSE)
    }
  }
  if (return_plot) list(plot = plots, statistics = statistics) else statistics
}
