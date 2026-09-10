.biomed_require <- function(package, reason = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    detail <- if (is.null(reason)) "" else paste0(" ", reason)
    cli::cli_abort(
      "Package {.pkg {package}} is required.{detail} Install it with {.code install.packages(\"{package}\")} ."
    )
  }
  invisible(TRUE)
}

.biomed_colour_values <- function(groups, colour = NULL) {
  if (is.null(colour)) {
    return(stats::setNames(.biomed_palette(length(groups)), groups))
  }
  if (!is.atomic(colour)) {
    cli::cli_abort("{.arg colour} must be a colour vector or a named list.")
  }
  if (is.null(names(colour))) {
    if (length(colour) < length(groups)) {
      cli::cli_abort("At least {length(groups)} colours are required.")
    }
    return(stats::setNames(colour[seq_along(groups)], groups))
  }
  missing <- setdiff(groups, names(colour))
  if (length(missing)) {
    cli::cli_abort("Colours are missing for: {missing}.")
  }
  colour[groups]
}

.biomed_required_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    cli::cli_abort(c(
      "Required columns are missing.",
      "x" = "Missing: {missing}."
    ))
  }
  invisible(TRUE)
}

.biomed_as_numeric <- function(x, name) {
  if (is.numeric(x)) {
    if (any(is.infinite(x))) cli::cli_abort("{.val {name}} contains infinite values.")
    return(as.numeric(x))
  }

  raw <- as.character(x)
  out <- suppressWarnings(as.numeric(raw))
  if (any(is.infinite(out))) cli::cli_abort("{.val {name}} contains infinite values.")
  bad <- !is.na(raw) & nzchar(raw) & is.na(out)

  if (any(bad)) {
    examples <- unique(raw[bad])
    examples <- utils::head(examples, 3L)
    cli::cli_abort(
      "{.val {name}} must be numeric; values such as {examples} cannot be converted safely."
    )
  }

  out
}

.biomed_validate_columns <- function(columns) {
  if (!is.character(columns) || !length(columns) || anyNA(columns) ||
      any(!nzchar(columns)) || anyDuplicated(columns)) {
    cli::cli_abort("Columns must be a non-empty character vector of unique, non-missing names.")
  }
  invisible(columns)
}

.biomed_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    cli::cli_abort("{.val {name}} must be a single number strictly between 0 and 1.")
  }
}

.biomed_standard_table <- function(x, mapping) {
  for (old in names(mapping)) {
    names(x)[names(x) == old] <- unname(mapping[[old]])
  }
  rownames(x) <- NULL
  x
}

.biomed_save_table <- function(x, output_dir, filename) {
  if (!is.null(output_dir)) {
    .biomed_require("openxlsx")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    openxlsx::write.xlsx(x, file.path(output_dir, filename), overwrite = TRUE)
  }
  x
}

.biomed_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "output")
}

.biomed_palette <- function(n) {
  base <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
  )
  if (n <= length(base)) {
    return(base[seq_len(n)])
  }
  grDevices::colorRampPalette(base)(n)
}

.biomed_significance <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
  )
}
