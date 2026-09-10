#' Legacy correlation plot
#'
#' Compute correlation and plot the same filtered observations. Original zero values and missing groups are filtered before scaling and testing. Prefer plot_correlation for standardized names.
#' @param df Data frame or matrix for correlation; data frame for survival.
#' @param var1 X-axis column name or numeric position.
#' @param var2 Y-axis column name or numeric position.
#' @param is.matrix Compatibility argument; matrix input is converted to a data frame.
#' @param scale Standardize x and y after filtering.
#' @param subtype Optional grouping column for colors and group-specific regression lines.
#' @param na.subtype.rm Drop missing group values before correlation calculation and plotting.
#' @param color_subtype Optional group colors, optionally named.
#' @param palette Built-in biomed (default), or a ggsci palette: jama, npg, aaas, nejm, lancet, jco or d3; used when group colors are not supplied.
#' @param index Optional output filename prefix; defaults to 1.
#' @param method Correlation method: spearman, pearson or kendall.
#' @param show_cor_result Print the cor.test result.
#' @param col_line Regression-line color; default depends on overall correlation sign.
#' @param id Optional text-label column.
#' @param show_label Show text labels using id.
#' @param point_size Point size.
#' @param title Plot title; defaults to the compared column names.
#' @param alpha Point opacity.
#' @param title_size Title size relative to text_size.
#' @param text_size Base text size.
#' @param axis_angle X-axis text angle.
#' @param hjust Horizontal justification for title and x-axis text.
#' @param show_plot Print the plot.
#' @param save_plot Save the plot when path is supplied.
#' @param path Optional plot output directory.
#' @param fig.format ggsave output format, e.g. png or pdf.
#' @param fig.width Figure width in inches.
#' @param fig.height Figure height in inches.
#' @param add.hdr.line Add highest-density-region contours using optional ggdensity.
#' @param remove.x.zero Remove original x values equal to zero before scaling.
#' @param remove.y.zero Remove original y values equal to zero before scaling.
#' @param point.fill Ungrouped point fill.
#' @param point.color Ungrouped point outline color.
#' @param point.shape Ungrouped point shape.
#' @param stroke Ungrouped point outline width.
#' @param add.regress Add linear regression line(s) with confidence bands, regardless of correlation method; grouped plots fit a line within each group.
#' @param save_data Also export analyzed plot data as RData. Legacy get_cor defaults to TRUE when saving a plot; plot_correlation defaults to FALSE.
#' @return A ggplot invisibly, with a correlation attribute containing x, y, method, n, estimate and p_value. Insufficient or constant data warns and returns NULL. NULL df returns NULL.
#' @export
get_cor <- function(
    df = NULL,
    var1,
    var2,
    is.matrix = FALSE,

    scale = TRUE,
    subtype = NULL,
    na.subtype.rm = FALSE,
    color_subtype = NULL,
    palette = "biomed",
    index = NULL,
    method = c("spearman", "pearson", "kendall"),
    show_cor_result = FALSE,
    col_line = NULL,
    id = NULL,
    show_label = FALSE,
    point_size = 4,
    title = NULL,
    alpha = 0.5,
    title_size = 1.5,
    text_size = 10,
    axis_angle = 0,
    hjust = 0,
    show_plot = TRUE,
    save_plot = FALSE,
    path = NULL,
    fig.format = "png",
    fig.width = 7,
    fig.height = 7.3,
    add.hdr.line = FALSE,

    remove.x.zero = FALSE,
    remove.y.zero = FALSE,

    point.fill = "lightblue",
    point.color = "black",
    point.shape = 21,
    stroke = 0.5,

    add.regress = TRUE,
    save_data = TRUE
) {

  if (is.null(df)) return(NULL)

  method <- rlang::arg_match(method)
  if (is.null(index)) index <- 1

  if (missing(var1) || missing(var2)) {
    cli::cli_abort("Both {.arg var1} and {.arg var2} must be specified.")
  }

  df <- as.data.frame(df)

  if (is.numeric(var1)) {
    var1 <- colnames(df)[var1]
  }

  if (is.numeric(var2)) {
    var2 <- colnames(df)[var2]
  }

  if (length(var1) != 1L || is.na(var1) || !var1 %in% colnames(df)) {
    cli::cli_abort("Variable {.val {var1}} not found in df.")
  }

  if (length(var2) != 1L || is.na(var2) || !var2 %in% colnames(df)) {
    cli::cli_abort("Variable {.val {var2}} not found in df.")
  }

  if (isTRUE(is.matrix)) {
    df <- as.data.frame(df)
  }

  data <- df

  data[[var1]] <- .biomed_as_numeric(data[[var1]], var1)
  data[[var2]] <- .biomed_as_numeric(data[[var2]], var2)
  if (!is.null(subtype)) {
    .biomed_required_columns(data, subtype)
    if (na.subtype.rm) data <- data[!is.na(data[[subtype]]), , drop = FALSE]
  }
  data <- data[
    !is.na(data[[var1]]) & !is.na(data[[var2]]),
    ,
    drop = FALSE
  ]

  if (isTRUE(remove.x.zero)) {
    data <- data[data[[var1]] != 0, , drop = FALSE]
  }

  if (isTRUE(remove.y.zero)) {
    data <- data[data[[var2]] != 0, , drop = FALSE]
  }

  if (nrow(data) < 3) {
    cli::cli_warn("Insufficient data after filtering (need >= 3): {.val {var1}} vs {.val {var2}}.")
    return(NULL)
  }

  if (length(unique(data[[var1]])) < 2 || length(unique(data[[var2]])) < 2) {
    cli::cli_warn("Variable {.val {var1}} or {.val {var2}} has no variation after filtering.")
    return(NULL)
  }

  if (isTRUE(scale)) {
    data[[var1]] <- as.numeric(base::scale(data[[var1]]))
    data[[var2]] <- as.numeric(base::scale(data[[var2]]))
  }

  cli::cli_alert_info("Calculating {method} correlation (n = {nrow(data)})")

  cor_result <- stats::cor.test(
    data[[var1]],
    data[[var2]],
    method = method,
    exact = FALSE
  )

  if (show_cor_result) print(cor_result)

  pvalue <- cor_result$p.value

  cli::cli_alert_info(
    "P-value: {format(pvalue, digits = 2, scientific = TRUE)}"
  )

  if (is.null(col_line)) {
    col_line <- if (as.numeric(cor_result$estimate) > 0) {
      "darkred"
    } else {
      "steelblue"
    }
  }

  if (!is.null(subtype)) {
    if (!subtype %in% colnames(data)) {
      cli::cli_warn("Subtype {.val {subtype}} not found. Ignoring.")
      subtype <- NULL
    }
  }

  if (is.null(color_subtype) && !is.null(subtype)) {
    groups <- unique(as.character(data[[subtype]]))
    groups[is.na(groups)] <- "Not_available"
    n_colors <- length(groups)
    if (identical(palette, "biomed")) {
      color_subtype <- biomed_colors(n_colors)
    } else if (palette %in% c("jama", "npg", "aaas", "nejm", "lancet", "jco", "d3")) {
      .biomed_require("ggsci")
      palette_fun <- getExportedValue("ggsci", paste0("pal_", palette))
      n_base <- switch(palette, jama = 7L, nejm = 8L, lancet = 9L, 10L)
      base_colors <- palette_fun()(n_base)
      base_colors <- base_colors[!is.na(base_colors)]
      color_subtype <- grDevices::colorRampPalette(base_colors)(n_colors)
    } else {
      stop("Unknown palette; provide color_subtype explicitly or select a supported ggsci palette.")
    }
  }

  cor_label <- switch(
    method,
    spearman = "Spearman rho",
    pearson  = "Pearson R",
    kendall  = "Kendall tau"
  )

  label_text <- paste0(
    cor_label, " = ", round(as.numeric(cor_result$estimate), 3),
    "\nP ", ifelse(
      pvalue < 0.001,
      "< 0.001",
      paste0("= ", format(pvalue, digits = 2, scientific = TRUE))
    )
  )

  if (is.null(title)) {
    title <- paste0(var1, " vs ", var2)
  }

  if (!is.null(subtype)) {

    group_col <- .biomed_unique_name(".biomed_group", names(data))
    data[[group_col]] <- data[[subtype]]

    if (na.subtype.rm) {
      data <- data[!is.na(data[[group_col]]), , drop = FALSE]
    }

    data[[group_col]] <- as.character(data[[group_col]])
    data[[group_col]][is.na(data[[group_col]])] <- "Not_available"
    data[[group_col]] <- factor(data[[group_col]])

    color_subtype <- .biomed_colour_values(levels(data[[group_col]]), color_subtype)
    cli::cli_alert_info("Groups: {.val {levels(data[[group_col]])}}")

    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data[[var1]],
        y = .data[[var2]],
        colour = .data[[group_col]]
      )
    ) +
      ggplot2::geom_point(
        size = point_size,
        alpha = alpha
      ) +
      ggplot2::scale_color_manual(values = color_subtype)

  } else {

    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data[[var1]],
        y = .data[[var2]]
      )
    ) +
      ggplot2::geom_point(
        shape = point.shape,
        fill = point.fill,
        colour = point.color,
        size = point_size,
        alpha = alpha,
        stroke = stroke
      )
  }

  if (isTRUE(add.regress)) {
    p <- p +
      ggplot2::geom_smooth(
        method = "lm",
        se = TRUE,
        color = col_line,
        linewidth = 0.8
      )
  }

  p <- p +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = label_text,
      hjust = -0.05,
      vjust = 1.2,
      size = 5
    ) +
    ggplot2::labs(
      x = var1,
      y = var2,
      title = title
    ) +
    ggplot2::theme_classic(base_size = text_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = hjust,
        face = "bold",
        size = text_size * title_size
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = axis_angle,
        hjust = hjust
      )
    )

  if (show_label) {
    if (is.null(id) || !id %in% colnames(data)) {
      cli::cli_warn("Label column not found. Set {.arg id} to show labels.")
    } else {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = .data[[id]]),
          nudge_x = 0.25,
          nudge_y = 0.25,
          check_overlap = TRUE,
          size = 3
        )
    }
  }

  if (add.hdr.line) {
    rlang::check_installed("ggdensity")
    p <- p + ggdensity::geom_hdr_lines()
  }

  if (show_plot) print(p)

  if (save_plot) {

    if (is.null(path)) {
      cli::cli_warn("{.arg path} is NULL; plot will not be saved.")
    } else {

      if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
      }

      filename <- paste0(sanitize_filename(paste(index, var2, var1, "correlation", sep = "-")), ".", fig.format)

      ggplot2::ggsave(
        filename = file.path(path, filename),
        plot = p,
        width = fig.width,
        height = fig.height
      )

      if (isTRUE(save_data)) save(
        data,
        file = file.path(path, paste0("0-input-data-", sanitize_filename(paste(var1, var2, sep = "-")), ".RData"))
      )

      cli::cli_alert_success("Plot saved to {.path {path}}")
    }
  }

  attr(p, "correlation") <- data.frame(
    x = var1, y = var2, method = method, n = nrow(data),
    estimate = unname(cor_result$estimate), p_value = cor_result$p.value)
  invisible(p)
}
