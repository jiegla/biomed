.biomed_survival_data <- function(data, group, time, status, status_encoding, reference = NULL) {
  .biomed_required_columns(data, c(group, time, status))
  d <- data.frame(.time = .biomed_as_numeric(data[[time]], time),
                  .status = .biomed_status01(data[[status]], status_encoding),
                  .group = data[[group]])
  if (any(d$.time < 0, na.rm = TRUE)) stop("Time must be nonnegative.")
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (!nrow(d)) stop("No complete survival observations.")
  d$.group <- droplevels(factor(d$.group))
  if (!is.null(reference)) {
    if (length(reference) != 1L || !reference %in% levels(d$.group)) stop("reference must name an observed group.")
    d$.group <- stats::relevel(d$.group, reference)
  }
  d
}

.biomed_pairwise_survival <- function(d, p_adjust = "BH") {
  lv <- levels(d$.group)
  empty <- data.frame(group = character(), reference = character(), n = integer(),
                      events = integer(), hr = double(), conf_low = double(),
                      conf_high = double(), p_value = double(), logrank_p = double(),
                      warning = character(), p_adjusted = double())
  if (length(lv) < 2L) return(empty)
  pairs <- utils::combn(lv, 2L, simplify = FALSE)
  result <- lapply(pairs, function(pair) {
    sub <- d[d$.group %in% pair, , drop = FALSE]
    sub$.group <- factor(sub$.group, levels = pair)
    warnings <- character()
    capture <- function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
    fit <- withCallingHandlers(tryCatch(
      survival::coxph(survival::Surv(.time, .status) ~ .group, data = sub),
      error = function(e) { warnings <<- c(warnings, conditionMessage(e)); NULL }),
      warning = capture)
    vals <- rep(NA_real_, 4L)
    if (!is.null(fit)) {
      s <- summary(fit)
      vals <- c(s$conf.int[1L, "exp(coef)"], s$conf.int[1L, "lower .95"],
                s$conf.int[1L, "upper .95"], s$coefficients[1L, "Pr(>|z|)"])
    }
    lp <- withCallingHandlers(tryCatch({
      lr <- survival::survdiff(survival::Surv(.time, .status) ~ .group, data = sub)
      stats::pchisq(lr$chisq, 1L, lower.tail = FALSE)
    }, error = function(e) NA_real_), warning = capture)
    data.frame(group = pair[2L], reference = pair[1L], n = nrow(sub),
               events = sum(sub$.status), hr = vals[1L], conf_low = vals[2L],
               conf_high = vals[3L], p_value = vals[4L], logrank_p = lp,
               warning = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_)
  })
  result <- do.call(rbind, result)
  result$p_adjusted <- stats::p.adjust(result$logrank_p, method = p_adjust)
  rownames(result) <- NULL
  result
}

.biomed_save_survival <- function(plot, output_dir, filename, formats, width, height) {
  if (is.null(output_dir)) return(invisible(NULL))
  formats <- match.arg(formats, c("pdf", "png"), several.ok = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (ext in formats) {
    file <- file.path(output_dir, paste0(sanitize_filename(filename), ".", ext))
    if (ext == "pdf") {
      grDevices::pdf(file, width = width, height = height)
    } else {
      grDevices::png(file, width = width, height = height, units = "in", res = 150)
    }
    tryCatch(print(plot), finally = grDevices::dev.off())
  }
  stats <- attr(plot, "statistics")
  writeLines(utils::capture.output(print(stats)),
             file.path(output_dir, paste0(sanitize_filename(filename), "_statistics.txt")))
  invisible(NULL)
}

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
  if (!is.character(columns) || anyNA(columns)) {
    cli::cli_abort("Column names must be a character vector without missing values.")
  }
  empty <- !nzchar(trimws(columns))
  if (any(empty)) {
    cli::cli_abort("Column names must not be empty or contain only whitespace.")
  }
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    cli::cli_abort(c(
      "Required columns are missing.",
      "x" = "Missing: {paste(missing, collapse = ', ')}.",
      "i" = "Available columns: {paste(names(data), collapse = ', ')}."
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
  biomed_colors(n)
}

.biomed_significance <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
  )
}

.biomed_numeric_values <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

.biomed_status01 <- function(x, encoding = "auto") {
  encoding <- match.arg(encoding, c("auto", "01", "12", "labels"))
  if (is.logical(x) && encoding %in% c("01", "auto", "labels")) return(as.numeric(x))
  raw <- trimws(tolower(as.character(x)))
  numeric_x <- suppressWarnings(as.numeric(raw))
  present <- !is.na(x)
  if (encoding %in% c("01", "12")) {
    allowed <- if (encoding == "01") c(0, 1) else c(1, 2)
    if (any(present & (is.na(numeric_x) | !numeric_x %in% allowed))) {
      stop("Status must use the selected ", encoding, " encoding.")
    }
    return(if (encoding == "12") numeric_x - 1 else numeric_x)
  }
  if (encoding == "auto" && all(!present | !is.na(numeric_x))) {
    observed <- unique(numeric_x[present])
    if (all(observed %in% c(0, 1))) return(numeric_x)
    if (all(observed %in% c(1, 2))) return(numeric_x - 1)
    stop("Unsupported numeric status; use 0/1 or 1/2.")
  }
  events <- c("1", "true", "t", "yes", "y", "event", "dead", "death",
              "deceased", "progression", "progressed", "pd", "relapse",
              "recurred", "recurrence")
  censored <- c("0", "false", "f", "no", "n", "censored", "alive", "living",
                "no event", "noevent", "non-event", "non event", "no progression")
  out <- rep(NA_real_, length(x))
  out[raw %in% events] <- 1
  out[raw %in% censored] <- 0
  if (any(present & is.na(out))) stop("Unrecognized status labels.")
  out
}

.biomed_survival_options <- function(minprop, counts = numeric()) {
  if (length(minprop) != 1L || !is.numeric(minprop) || !is.finite(minprop) ||
      minprop <= 0 || minprop > 0.5) stop("minprop must be in (0, 0.5].")
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts))) stop("Sample/event thresholds must be nonnegative integers.")
}

.biomed_valid_cutoffs <- function(x, minprop = 0.1) {
  x <- x[!is.na(x)]
  n <- length(x)

  if (n == 0) return(numeric(0))

  ux <- sort(unique(x))

  if (length(ux) < 2) return(numeric(0))

  cuts <- ux[-length(ux)]
  min_n <- ceiling(n * minprop)

  cuts[vapply(
    cuts,
    function(cut) {
      n_low <- sum(x <= cut, na.rm = TRUE)
      n_high <- sum(x > cut, na.rm = TRUE)
      n_low >= min_n && n_high >= min_n
    },
    logical(1)
  )]
}

.biomed_extract_cox <- function(fit) {
  s <- summary(fit)

  co <- as.data.frame(s$coefficients, check.names = FALSE)
  ci <- as.data.frame(s$conf.int, check.names = FALSE)

  low_col <- grep("lower", colnames(ci), value = TRUE)[1]
  high_col <- grep("upper", colnames(ci), value = TRUE)[1]
  p_col <- grep("^Pr", colnames(co), value = TRUE)[1]

  data.frame(
    term = rownames(co),
    HR = ci[["exp(coef)"]],
    lower95 = ci[[low_col]],
    upper95 = ci[[high_col]],
    P = co[[p_col]],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.biomed_add_failure <- function(fail_df, ID, stage, reason) {
  rbind(
    fail_df,
    data.frame(
      ID = as.character(ID),
      stage = as.character(stage),
      reason = as.character(reason),
      stringsAsFactors = FALSE
    )
  )
}

.biomed_unique_name <- function(base, nms) {
  nm <- base
  while (nm %in% nms) {
    nm <- paste0(".", nm)
  }
  nm
}

.biomed_roc_data <- function(data, response, predictor, positive_class) {
  .biomed_required_columns(data, c(response, predictor))
  dat <- data.frame(
    .response = data[[response]],
    .predictor = .biomed_as_numeric(data[[predictor]], predictor)
  )
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  y <- as.character(dat$.response)
  observed <- if (is.factor(dat$.response)) levels(droplevels(dat$.response)) else sort(unique(y))
  if (length(observed) != 2L) {
    cli::cli_abort("{.arg response} must contain exactly two observed classes.")
  }
  positive <- if (is.null(positive_class)) observed[2L] else as.character(positive_class)
  if (length(positive) != 1L || is.na(positive)) {
    cli::cli_abort("{.arg positive_class} must identify exactly one observed class.")
  }
  if (!positive %in% observed) {
    cli::cli_abort("Positive class {.val {positive}} was not found in {.arg response}.")
  }
  negative <- setdiff(observed, positive)
  dat$.response <- factor(y, levels = c(negative, positive))
  list(data = dat, negative = negative, positive = positive)
}

.biomed_threshold_metrics <- function(roc_obj, dat, threshold) {
  positive <- levels(dat$.response)[2L]
  predicted <- if (roc_obj$direction == "<") {
    dat$.predictor >= threshold
  } else {
    dat$.predictor <= threshold
  }
  actual <- dat$.response == positive
  tp <- sum(predicted & actual)
  tn <- sum(!predicted & !actual)
  fp <- sum(predicted & !actual)
  fn <- sum(!predicted & actual)
  safe_div <- function(x, y) if (y == 0) NA_real_ else x / y
  c(
    Specificity = safe_div(tn, tn + fp),
    Sensitivity = safe_div(tp, tp + fn),
    Accuracy = safe_div(tp + tn, length(actual)),
    Precision = safe_div(tp, tp + fp),
    Recall = safe_div(tp, tp + fn)
  )
}

.biomed_partitions <- function(sets) {
  universe <- unique(unlist(sets, use.names = FALSE))
  masks <- as.matrix(expand.grid(rep(base::list(c(FALSE, TRUE)), length(sets))))
  masks <- masks[rowSums(masks) > 0L, , drop = FALSE]
  result <- as.data.frame(masks)
  names(result) <- names(sets)
  members <- lapply(seq_len(nrow(masks)), function(i) {
    keep <- rep(TRUE, length(universe))
    for (j in seq_along(sets)) {
      keep <- keep & ((universe %in% sets[[j]]) == masks[i, j])
    }
    universe[keep]
  })
  result$set <- apply(masks, 1L, function(x) paste(names(sets)[x], collapse = " & "))
  result$count <- lengths(members)
  result$values <- vapply(members, paste, collapse = ", ", FUN.VALUE = character(1))
  result
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
