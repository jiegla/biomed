#' Draw donut charts for categorical variables
#'
#' Creates one donut chart per categorical variable. Missing values can either
#' be removed or displayed as a separate category.
#'
#' @param df A data frame.
#' @param vars Categorical columns to plot.
#' @param blank Inner radius control between `0` and `2.5`.
#' @param linewidth Width of the white borders between slices.
#' @param showtext Show labels inside slices.
#' @param text_type Label content: count and percent, percent, or count.
#' @param show_center_n Show the number of analysed observations in the centre.
#' @param removeNA Remove missing values when `TRUE`.
#' @param na_label Label used for missing values when `removeNA = FALSE`.
#' @param color_list A colour vector shared by all variables, or a named list
#'   containing one vector per variable. Named vectors map colours to levels.
#' @param legend Display a legend.
#'
#' @return A named list of ggplot objects.
#' @export
donut_pie <- function(
    df,
    vars,
    blank = 0.8,
    linewidth = 2.5,
    showtext = TRUE,
    text_type = c("both", "percent", "count"),
    show_center_n = TRUE,
    removeNA = TRUE,
    na_label = "NA",
    color_list = NULL,
    legend = TRUE) {
  text_type <- match.arg(text_type)
  if (!is.data.frame(df)) cli::cli_abort("{.arg df} must be a data frame.")
  if (!is.character(vars) || !length(vars)) {
    cli::cli_abort("{.arg vars} must be a non-empty character vector.")
  }
  .biomed_required_columns(df, vars)
  if (!is.numeric(blank) || length(blank) != 1L || is.na(blank) ||
      blank < 0 || blank >= 2.5) {
    cli::cli_abort("{.arg blank} must be one number in `[0, 2.5)`. ")
  }
  if (!is.numeric(linewidth) || length(linewidth) != 1L || is.na(linewidth) ||
      linewidth < 0) {
    cli::cli_abort("{.arg linewidth} must be one non-negative number.")
  }
  logical_args <- base::list(showtext, show_center_n, removeNA, legend)
  valid_logical <- vapply(
    logical_args,
    function(x) is.logical(x) && length(x) == 1L && !is.na(x),
    logical(1)
  )
  if (!all(valid_logical)) {
    cli::cli_abort("Display and missing-value options must be single logical values.")
  }

  make_plot <- function(varname) {
    source <- df[[varname]]
    values <- as.character(source)
    if (removeNA) {
      values <- values[!is.na(source)]
    } else {
      values[is.na(source)] <- na_label
    }
    if (!length(values)) {
      return(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0, y = 0, label = "No valid data") +
          ggplot2::labs(title = varname) +
          ggplot2::theme_void()
      )
    }

    groups <- if (is.factor(source)) {
      factor_levels <- levels(source)
      if (!removeNA && any(is.na(source))) factor_levels <- c(factor_levels, na_label)
      factor_levels[factor_levels %in% values]
    } else {
      sort(unique(values))
    }
    counts <- table(factor(values, levels = groups))
    plot_data <- data.frame(
      Group = factor(names(counts), levels = groups),
      n = as.integer(counts), stringsAsFactors = FALSE
    )
    plot_data$fraction <- plot_data$n / sum(plot_data$n)
    plot_data$percent <- scales::percent(plot_data$fraction, accuracy = 0.1)
    plot_data$label <- switch(
      text_type,
      both = paste0(plot_data$n, "\n", plot_data$percent),
      percent = plot_data$percent,
      count = as.character(plot_data$n)
    )

    colours <- color_list
    if (is.list(color_list)) {
      colours <- color_list[[varname]]
      if (is.null(colours)) {
        cli::cli_abort("{.arg color_list} has no entry for {.val {varname}}.")
      }
    }
    fill_values <- .biomed_colour_values(as.character(plot_data$Group), colours)
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = 2, y = .data[["fraction"]], fill = .data[["Group"]])
    ) +
      ggplot2::geom_col(width = 1, colour = "white", linewidth = linewidth) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::xlim(blank, 2.5) +
      ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
      ggplot2::labs(title = varname, fill = varname) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = if (legend) "right" else "none"
      )
    if (showtext) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data[["label"]]),
        position = ggplot2::position_stack(vjust = 0.5), size = 4
      )
    }
    if (show_center_n) {
      p <- p + ggplot2::annotate(
        "text", x = blank + (2 - blank) / 2, y = 0,
        label = paste0("N=", sum(plot_data$n)), size = 5, fontface = "bold"
      )
    }
    p
  }

  plots <- lapply(vars, make_plot)
  names(plots) <- vars
  plots
}
