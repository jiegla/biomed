#' Draw a Venn diagram and extract its partitions
#'
#' Draws a Venn diagram for two to six sets and optionally saves publication
#' figures and a workbook of intersection members. Two to five sets use Venn
#' geometry; six sets use an exact intersection membership matrix.
#'
#' @param list A named list of vectors representing sets.
#' @param venn_name Output filename stem.
#' @param col,fill_col,cat.col Outline, fill, and category-label colours.
#' @param palette Colour palette retained for compatibility.
#' @param showFigure Draw on the current graphics device.
#' @param file_type Any of `"tif"`, `"png"`, and `"pdf"`.
#' @param width,height Device dimensions.
#' @param units Raster device units.
#' @param res Raster resolution.
#' @param alpha,cex,cat.cex,lwd,lty,main,main.cex,sub,sub.cex,rotation.degree,reverse
#'   Graphical arguments passed to [VennDiagram::venn.diagram()].
#' @param out_dir Output directory.
#' @param save Save requested figures.
#' @param write_xlsx Save the partition table as XLSX.
#'
#' @return Invisibly, a list containing the grob, partition table, and output
#'   paths.
#' @export
vennjgl <- function(
    list,
    venn_name = "test",
    col = NULL,
    palette = c("npg", "aaas", "nejm", "lancet", "jama"),
    fill_col = NULL,
    cat.col = NULL,
    showFigure = TRUE,
    file_type = c("tif", "png", "pdf"),
    width = 6,
    height = 5,
    units = "in",
    res = 300,
    alpha = 0.3,
    cex = 1.5,
    cat.cex = 2,
    lwd = 2,
    lty = 1,
    main = NULL,
    main.cex = 2,
    sub = NULL,
    sub.cex = 1,
    rotation.degree = 0,
    reverse = TRUE,
    out_dir = ".",
    save = TRUE,
    write_xlsx = TRUE) {
  .biomed_require("VennDiagram")
  if (!is.list(list)) cli::cli_abort("{.arg list} must be a list of vectors.")
  n_set <- length(list)
  if (n_set < 2L || n_set > 6L) {
    cli::cli_abort("{.arg list} must contain between two and six sets.")
  }
  nm <- names(list)
  if (is.null(nm)) nm <- rep("", n_set)
  missing_name <- is.na(nm) | !nzchar(nm)
  nm[missing_name] <- paste0("vector", which(missing_name))
  names(list) <- make.unique(nm)
  list <- lapply(list, function(x) {
    x <- unique(as.character(x))
    x[!is.na(x) & nzchar(x)]
  })
  palette <- match.arg(palette)
  file_type <- unique(tolower(file_type))
  if (!length(file_type) || !all(file_type %in% c("tif", "png", "pdf"))) {
    cli::cli_abort("{.arg file_type} supports only {.val tif}, {.val png}, and {.val pdf}.")
  }

  auto_colours <- if (requireNamespace("ggsci", quietly = TRUE)) {
    fn <- switch(
      palette,
      npg = ggsci::pal_npg("nrc"), aaas = ggsci::pal_aaas(),
      nejm = ggsci::pal_nejm(), lancet = ggsci::pal_lancet(),
      jama = ggsci::pal_jama()
    )
    fn(n_set)
  } else {
    .biomed_palette(n_set)
  }
  col <- if (is.null(col)) auto_colours else col
  fill_col <- if (is.null(fill_col)) auto_colours else fill_col
  cat.col <- if (is.null(cat.col)) auto_colours else cat.col
  if (any(lengths(base::list(col, fill_col, cat.col)) < n_set)) {
    cli::cli_abort("Each colour vector must contain at least {n_set} colours.")
  }
  if (is.null(main)) main <- venn_name

  venn_plot <- if (n_set == 6L) {
    .biomed_intersection_matrix(list, fill_col, main, cex)
  } else VennDiagram::venn.diagram(
    x = list, filename = NULL, lwd = lwd, lty = lty, col = col,
    fill = fill_col, cat.col = cat.col, cat.cex = cat.cex,
    rotation.degree = rotation.degree, main = main, main.cex = main.cex,
    sub = sub, sub.cex = sub.cex, cex = cex, alpha = alpha,
    reverse = reverse, disable.logging = TRUE
  )
  if (isTRUE(showFigure)) {
    grid::grid.newpage()
    grid::grid.draw(venn_plot)
  }

  image_files <- character()
  if (isTRUE(save)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stem <- file.path(out_dir, .biomed_safe_name(venn_name))
    draw_to_device <- function(type) {
      path <- paste0(stem, ".", type)
      if (type == "tif") {
        grDevices::tiff(path, width = width, height = height, units = units,
                        res = res, compression = "lzw")
      } else if (type == "png") {
        grDevices::png(path, width = width, height = height, units = units, res = res)
      } else {
        divisor <- switch(units, "in" = 1, "cm" = 2.54, "mm" = 25.4, "px" = res,
                          cli::cli_abort("Unsupported graphics units."))
        grDevices::pdf(path, width = width / divisor, height = height / divisor)
      }
      on.exit(grDevices::dev.off())
      grid::grid.newpage()
      grid::grid.draw(venn_plot)
      path
    }
    image_files <- vapply(file_type, draw_to_device, character(1))
  }

  partition <- VennDiagram::get.venn.partitions(list)
  partition$values <- vapply(
    partition[["..values.."]], paste, collapse = ", ", FUN.VALUE = character(1)
  )
  partition[["..values.."]] <- NULL
  names(partition)[names(partition) == "..set.."] <- "set"
  names(partition)[names(partition) == "..count.."] <- "count"
  excel_file <- NULL
  if (isTRUE(write_xlsx)) {
    .biomed_require("openxlsx", "It is used to save Venn partitions.")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    excel_file <- file.path(out_dir, paste0(.biomed_safe_name(venn_name), ".xlsx"))
    openxlsx::write.xlsx(partition, file = excel_file, overwrite = TRUE)
  }
  invisible(base::list(
    venn_plot = venn_plot, partition = partition,
    image_files = image_files, excel_file = excel_file
  ))
}

# VennDiagram supports at most five sets. Six sets use an exact membership
# matrix, avoiding a misleading arrangement of overlapping circles.
.biomed_intersection_matrix <- function(sets, colors, title, cex) {
  universe <- unique(unlist(sets, use.names = FALSE))
  if (!length(universe)) cli::cli_abort("The sets contain no members.")
  membership <- vapply(sets, function(x) universe %in% x, logical(length(universe)))
  if (is.null(dim(membership))) membership <- matrix(membership, nrow = length(universe))
  keys <- apply(membership, 1L, paste0, collapse = "")
  counts <- sort(table(keys), decreasing = TRUE)
  rows <- match(names(counts), keys)
  m <- membership[rows, , drop = FALSE]
  d <- expand.grid(intersection = seq_len(nrow(m)), set = seq_len(ncol(m)))
  d$present <- as.vector(m)
  d$set_name <- factor(names(sets)[d$set], levels = rev(names(sets)))
  d$fill <- ifelse(d$present, colors[d$set], "grey92")
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[["intersection"]], y = .data[["set_name"]])) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data[["fill"]]), width = 0.8, height = 0.8) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(breaks = seq_along(counts), labels = as.integer(counts)) +
    ggplot2::labs(title = title, subtitle = "Exact intersections (six sets)",
                  x = "Number of members in each intersection", y = NULL) +
    ggplot2::theme_minimal(base_size = 10 * cex) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  ggplot2::ggplotGrob(p)
}
