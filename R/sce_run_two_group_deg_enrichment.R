# Two-group differential expression + ORA + GSEA for Seurat objects
# RStudio-friendly version with per-celltype RDA + Excel enrichment exports
# RStudio-friendly: accepts an in-memory Seurat object OR a Seurat file path.
# Also compatible with pipeline resource/functions/ usage.
# ORA/GSEA objects and ggplot objects are saved as .rda files grouped by cell type.

`%biomed_deg_or%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

.biomed_deg_merge_config <- function(default, user) {
  if (is.null(user)) return(default)
  if (!is.list(user)) stop("config must be a list, a YAML file path, or NULL.")
  utils::modifyList(default, user, keep.null = TRUE)
}

.biomed_deg_default_deg_enrichment_config <- function() {
  list(
    output_dir = "results/two_group_deg_enrichment",
    project_name = "group1_vs_group2",
    seed = 1234,
    de = list(
      assay = "RNA",
      slot = "data",
      test_use = "wilcox",
      # Keep this low/zero so GSEA receives a broad ranked gene list.
      test_logfc_threshold = 0,
      min_pct = 0.01,
      min_diff_pct = -Inf,
      min_cells_per_group = 20,
      min_cells_feature = 3,
      max_cells_per_ident = Inf,
      latent_vars = NULL,
      join_layers = TRUE,
      prep_sct = FALSE,
      verbose = FALSE
    ),
    thresholds = list(
      # Used to label DEG as Up/Down; not used to pre-filter FindMarkers.
      fc = 1,
      pvalue = 0.05,
      p_column = "p_val_adj"
    ),
    gene_id = list(
      de_id = "SYMBOL",
      gmt_id = "SYMBOL",
      orgdb = "org.Hs.eg.db",
      strip_ensembl_version = TRUE
    ),
    gmt = list(
      dir = NULL,
      files = NULL,
      pattern = "\\.gmt$",
      recursive = FALSE
    ),
    ora = list(
      run = TRUE,
      directions = c("Up", "Down", "All"),
      pvalue_cutoff = 1,
      p_adjust_method = "BH",
      qvalue_cutoff = 1,
      min_gs_size = 10,
      max_gs_size = 500,
      use_tested_genes_as_universe = TRUE
    ),
    gsea = list(
      run = TRUE,
      rank_column = "auto",
      exponent = 1,
      pvalue_cutoff = 1,
      p_adjust_method = "BH",
      min_gs_size = 10,
      max_gs_size = 500,
      eps = 1e-10,
      by = "fgsea",
      seed = TRUE
    ),
    plot = list(
      run = TRUE,
      show_n = 20,
      only_significant = FALSE,
      pvalue_cutoff = 0.05,
      width = 8,
      height = 6,
      dpi = 300,
      formats = c("png", "pdf")
    ),
    export = list(
      write_excel = TRUE,
      write_rda = TRUE,
      row_names = FALSE
    )
  )
}

.biomed_deg_read_analysis_config <- function(config = NULL) {
  default <- .biomed_deg_default_deg_enrichment_config()
  if (is.null(config)) return(default)

  if (is.character(config) && length(config) == 1L) {
    if (!file.exists(config)) stop("Config file does not exist: ", config)
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' is required when config is a YAML file.")
    }
    user <- yaml::read_yaml(config)
    return(.biomed_deg_merge_config(default, user))
  }

  .biomed_deg_merge_config(default, config)
}

.biomed_deg_read_seurat_object <- function(input) {
  # RStudio-friendly mode: accept an in-memory Seurat object directly.
  if (inherits(input, "Seurat")) {
    return(input)
  }

  # Pipeline/file mode: accept a single Seurat file path.
  if (!is.character(input) || length(input) != 1L || is.na(input) || !nzchar(input)) {
    stop(
      "The first argument must be either a Seurat object or one file path ",
      "(.qs, .rds, .rda, or .RData)."
    )
  }

  input <- path.expand(input)
  if (!file.exists(input)) stop("Seurat input file does not exist: ", input)
  ext <- tolower(tools::file_ext(input))

  if (identical(ext, "qs")) {
    if (!requireNamespace("qs", quietly = TRUE)) {
      stop("Package 'qs' is required for .qs files.")
    }
    obj <- qs::qread(input)
  } else if (identical(ext, "rds")) {
    obj <- readRDS(input)
  } else if (ext %in% c("rda", "rdata")) {
    env <- new.env(parent = emptyenv())
    loaded <- load(input, envir = env)
    candidates <- mget(loaded, envir = env, inherits = FALSE)
    is_seurat <- vapply(candidates, inherits, logical(1), what = "Seurat")
    candidates <- candidates[is_seurat]
    if (length(candidates) != 1L) {
      stop("An .rda/.RData input must contain exactly one Seurat object.")
    }
    obj <- candidates[[1]]
  } else {
    stop("Unsupported input format: .", ext, ". Use .qs, .rds, .rda, or .RData.")
  }

  if (!inherits(obj, "Seurat")) stop("The loaded object is not a Seurat object.")
  obj
}

.biomed_deg_input_source_label <- function(input) {
  if (inherits(input, "Seurat")) {
    return("in-memory Seurat object")
  }
  normalizePath(path.expand(input), winslash = "/", mustWork = FALSE)
}

.biomed_deg_has_gmt_config <- function(gmt_cfg) {
  has_dir <- !is.null(gmt_cfg$dir) &&
    length(gmt_cfg$dir) > 0L &&
    any(!is.na(gmt_cfg$dir) & nzchar(as.character(gmt_cfg$dir)))

  has_files <- !is.null(gmt_cfg$files) &&
    length(gmt_cfg$files) > 0L &&
    any(!is.na(gmt_cfg$files) & nzchar(as.character(gmt_cfg$files)))

  isTRUE(has_dir || has_files)
}

.biomed_deg_safe_filename <- function(x, max_len = 100L) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "NA"
  x <- gsub("[\\\\/:*?\"<>|]+", "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)
  substr(x, 1L, max_len)
}

.biomed_deg_safe_sheet_names <- function(x) {
  x <- as.character(x)
  for (bad in c("[", "]", ":", "*", "?", "/", "\\")) {
    x <- gsub(bad, "_", x, fixed = TRUE)
  }
  x[!nzchar(x)] <- "Sheet"
  x <- substr(x, 1L, 31L)

  out <- character(length(x))
  used <- character(0)
  for (i in seq_along(x)) {
    base <- x[[i]]
    candidate <- base
    k <- 1L
    while (candidate %in% used) {
      suffix <- paste0("_", k)
      candidate <- paste0(substr(base, 1L, 31L - nchar(suffix)), suffix)
      k <- k + 1L
    }
    out[[i]] <- candidate
    used <- c(used, candidate)
  }
  out
}

# Companion writer for the user's xlsx2dflist() reader.
# All worksheets freeze the first row for easier inspection in Excel.
.biomed_deg_dflist2xlsx <- function(df_list, filepath, row_names = FALSE, overwrite = TRUE) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for Excel export.")
  }
  if (!is.list(df_list)) stop("df_list must be a named list of data.frames.")
  if (length(df_list) == 0L) return(invisible(NULL))
  if (is.null(names(df_list))) names(df_list) <- paste0("Sheet", seq_along(df_list))

  sheet_names <- .biomed_deg_safe_sheet_names(names(df_list))
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(filepath) && !isTRUE(overwrite)) {
    stop("Excel file already exists and overwrite = FALSE: ", filepath)
  }

  wb <- openxlsx::createWorkbook()

  for (i in seq_along(df_list)) {
    dat <- df_list[[i]]
    if (is.null(dat)) dat <- data.frame()
    if (!is.data.frame(dat)) {
      dat <- tryCatch(
        as.data.frame(dat),
        error = function(e) data.frame(Note = paste0("Could not convert object to data.frame: ", conditionMessage(e)))
      )
    }

    openxlsx::addWorksheet(wb, sheetName = sheet_names[[i]])
    openxlsx::writeData(
      wb,
      sheet = sheet_names[[i]],
      x = dat,
      rowNames = row_names
    )
    openxlsx::freezePane(
      wb,
      sheet = sheet_names[[i]],
      firstRow = TRUE
    )
  }

  openxlsx::saveWorkbook(wb, filepath, overwrite = overwrite)
  invisible(filepath)
}

# Convert per-celltype ORA/GSEA object lists into Excel-ready data.frame lists.
# ORA: one worksheet per GMT x direction (Up/Down/All).
# GSEA: one worksheet per GMT collection.
.biomed_deg_enrichment_objects_to_excel_sheets <- function(object_list, analysis = c("ORA", "GSEA")) {
  analysis <- match.arg(analysis)
  if (!is.list(object_list) || length(object_list) == 0L) return(list())

  out <- list()
  gmt_names <- names(object_list)
  if (is.null(gmt_names)) gmt_names <- paste0("GMT", seq_along(object_list))

  if (identical(analysis, "ORA")) {
    for (i in seq_along(object_list)) {
      gmt_name <- gmt_names[[i]]
      gmt_obj <- object_list[[i]]

      # Current ORA structure is: GMT -> direction -> enrichResult.
      if (is.list(gmt_obj) && !inherits(gmt_obj, "enrichResult")) {
        if (length(gmt_obj) == 0L) {
          out[[gmt_name]] <- data.frame(
            Note = "No ORA result for this GMT collection.",
            stringsAsFactors = FALSE
          )
          next
        }

        direction_names <- names(gmt_obj)
        if (is.null(direction_names)) direction_names <- paste0("Result", seq_along(gmt_obj))

        for (j in seq_along(gmt_obj)) {
          direction <- direction_names[[j]]
          obj <- gmt_obj[[j]]
          df <- tryCatch(as.data.frame(obj), error = function(e) data.frame())
          if (nrow(df) == 0L) {
            df <- data.frame(
              Note = paste0("No enriched terms for direction: ", direction),
              stringsAsFactors = FALSE
            )
          }
          out[[paste0(gmt_name, "__", direction)]] <- df
        }
      } else {
        df <- tryCatch(as.data.frame(gmt_obj), error = function(e) data.frame())
        if (nrow(df) == 0L) {
          df <- data.frame(Note = "No ORA result for this GMT collection.", stringsAsFactors = FALSE)
        }
        out[[gmt_name]] <- df
      }
    }
  } else {
    for (i in seq_along(object_list)) {
      gmt_name <- gmt_names[[i]]
      obj <- object_list[[i]]
      df <- tryCatch(as.data.frame(obj), error = function(e) data.frame())
      if (nrow(df) == 0L) {
        df <- data.frame(
          Note = "No GSEA result for this GMT collection.",
          stringsAsFactors = FALSE
        )
      }
      out[[gmt_name]] <- df
    }
  }

  out
}

.biomed_deg_save_named_rda <- function(object, object_name, file, compress = TRUE) {
  if (!is.character(object_name) || length(object_name) != 1L || !nzchar(object_name)) {
    stop("object_name must be one non-empty character string.")
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  env <- new.env(parent = emptyenv())
  assign(object_name, object, envir = env)
  save(list = object_name, file = file, envir = env, compress = compress)
  invisible(file)
}

.biomed_deg_read_gmt_file <- function(file) {
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) {
    return(list(
      term2gene = data.frame(term = character(), gene = character()),
      term2name = data.frame(term = character(), name = character())
    ))
  }

  parsed <- lapply(lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]])
  parsed <- parsed[vapply(parsed, length, integer(1)) >= 3L]
  if (length(parsed) == 0L) stop("Invalid GMT file (no valid rows): ", file)

  term2name <- do.call(rbind, lapply(parsed, function(z) {
    data.frame(term = z[[1]], name = z[[2]], stringsAsFactors = FALSE)
  }))
  term2name <- unique(term2name)

  term2gene <- do.call(rbind, lapply(parsed, function(z) {
    genes <- unique(z[-c(1, 2)])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    data.frame(term = z[[1]], gene = genes, stringsAsFactors = FALSE)
  }))
  term2gene <- unique(term2gene)

  list(term2gene = term2gene, term2name = term2name)
}

.biomed_deg_collect_gmt_files <- function(gmt_cfg) {
  files <- gmt_cfg$files
  if (!is.null(files)) files <- path.expand(as.character(files))

  if (!is.null(gmt_cfg$dir)) {
    gmt_dir <- path.expand(gmt_cfg$dir)
    if (!dir.exists(gmt_dir)) stop("GMT directory does not exist: ", gmt_dir)
    from_dir <- list.files(
      path = gmt_dir,
      pattern = gmt_cfg$pattern %biomed_deg_or% "\\.gmt$",
      full.names = TRUE,
      recursive = isTRUE(gmt_cfg$recursive),
      ignore.case = TRUE
    )
    files <- c(files, from_dir)
  }

  files <- unique(files[file.exists(files)])
  if (length(files) == 0L) stop("No GMT files were found. Set config$gmt$dir or config$gmt$files.")
  files
}

.biomed_deg_get_orgdb <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("OrgDb package is required for gene ID conversion: ", pkg)
  }
  get(pkg, envir = asNamespace(pkg), inherits = FALSE)
}

.biomed_deg_clean_gene_ids <- function(x, id_type, strip_ensembl_version = TRUE) {
  x <- as.character(x)
  if (isTRUE(strip_ensembl_version) && identical(toupper(id_type), "ENSEMBL")) {
    x <- sub("\\..*$", "", x)
  }
  x
}

.biomed_deg_convert_gene_ids <- function(genes, from_type, to_type, orgdb_pkg,
                              strip_ensembl_version = TRUE) {
  from_type <- toupper(from_type)
  to_type <- toupper(to_type)
  genes <- .biomed_deg_clean_gene_ids(genes, from_type, strip_ensembl_version)
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])

  if (identical(from_type, to_type)) {
    return(data.frame(source = genes, target = genes, stringsAsFactors = FALSE))
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required for gene ID conversion.")
  }

  orgdb <- .biomed_deg_get_orgdb(orgdb_pkg)
  map <- suppressMessages(AnnotationDbi::select(
    x = orgdb,
    keys = genes,
    keytype = from_type,
    columns = to_type
  ))
  map <- map[, c(from_type, to_type), drop = FALSE]
  colnames(map) <- c("source", "target")
  map <- map[!is.na(map$target) & nzchar(map$target), , drop = FALSE]
  unique(map)
}

.biomed_deg_convert_gene_vector <- function(genes, gene_cfg) {
  map <- .biomed_deg_convert_gene_ids(
    genes = genes,
    from_type = gene_cfg$de_id,
    to_type = gene_cfg$gmt_id,
    orgdb_pkg = gene_cfg$orgdb,
    strip_ensembl_version = gene_cfg$strip_ensembl_version
  )
  unique(map$target)
}

.biomed_deg_make_ranked_gene_list <- function(deg, fc_col, gene_cfg) {
  score <- suppressWarnings(as.numeric(deg[[fc_col]]))
  genes <- .biomed_deg_clean_gene_ids(deg$gene, gene_cfg$de_id, gene_cfg$strip_ensembl_version)
  keep <- is.finite(score) & !is.na(genes) & nzchar(genes)
  score <- score[keep]
  genes <- genes[keep]

  map <- .biomed_deg_convert_gene_ids(
    genes = genes,
    from_type = gene_cfg$de_id,
    to_type = gene_cfg$gmt_id,
    orgdb_pkg = gene_cfg$orgdb,
    strip_ensembl_version = gene_cfg$strip_ensembl_version
  )
  score_df <- data.frame(source = genes, score = score, stringsAsFactors = FALSE)
  score_df <- merge(score_df, map, by = "source", all = FALSE)
  score_df <- score_df[is.finite(score_df$score) & !is.na(score_df$target), , drop = FALSE]
  if (nrow(score_df) == 0L) return(stats::setNames(numeric(), character()))

  # One target ID may map from several source IDs; keep the strongest absolute effect.
  score_df <- score_df[order(score_df$target, -abs(score_df$score)), , drop = FALSE]
  score_df <- score_df[!duplicated(score_df$target), , drop = FALSE]
  ranked <- stats::setNames(score_df$score, score_df$target)
  sort(ranked, decreasing = TRUE)
}

.biomed_deg_detect_fc_column <- function(deg, preferred = "auto") {
  if (!identical(preferred, "auto") && preferred %in% colnames(deg)) return(preferred)
  candidates <- c("avg_log2FC", "avg_logFC", "avg_log2_FC", "avg_diff")
  hit <- candidates[candidates %in% colnames(deg)]
  if (length(hit) == 0L) {
    hit <- grep("avg.*(log.*fc|diff)", colnames(deg), value = TRUE, ignore.case = TRUE)
  }
  if (length(hit) == 0L) stop("Could not detect a fold-change column in FindMarkers output.")
  hit[[1]]
}

.biomed_deg_bind_rows_safe <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, x)
  if (length(x) == 0L) return(data.frame())
  if (requireNamespace("dplyr", quietly = TRUE)) return(dplyr::bind_rows(x))

  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(all_names, names(z))
    for (nm in missing) z[[nm]] <- NA
    z[, all_names, drop = FALSE]
  })
  do.call(rbind, x)
}

.biomed_deg_add_enrichment_metadata <- function(df, celltype, gmt_name, direction = NA_character_) {
  if (nrow(df) == 0L) return(df)
  df$celltype <- celltype
  df$gmt <- gmt_name
  if (!is.na(direction)) df$direction <- direction
  df
}

.biomed_deg_enrich_result_df <- function(x) {
  if (is.null(x)) return(data.frame())
  out <- tryCatch(as.data.frame(x), error = function(e) data.frame())
  if (is.null(out) || nrow(out) == 0L) data.frame() else out
}

.biomed_deg_ratio_to_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2L) return(NA_real_)
    as.numeric(z[[1]]) / as.numeric(z[[2]])
  }, numeric(1))
}

.biomed_deg_save_plot_formats <- function(plot, stem, plot_cfg) {
  formats <- unique(tolower(plot_cfg$formats %biomed_deg_or% c("png", "pdf")))
  for (fmt in formats) {
    file <- paste0(stem, ".", fmt)
    ggplot2::ggsave(
      filename = file,
      plot = plot,
      width = plot_cfg$width,
      height = plot_cfg$height,
      dpi = plot_cfg$dpi
    )
  }
  invisible(stem)
}

.biomed_deg_plot_volcano <- function(deg, fc_col, p_col, title, show_n = 20L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for plots.")
  dat <- deg
  dat$plot_fc <- suppressWarnings(as.numeric(dat[[fc_col]]))
  dat$plot_p <- suppressWarnings(as.numeric(dat[[p_col]]))
  dat$plot_p[!is.finite(dat$plot_p) | is.na(dat$plot_p)] <- 1
  min_positive <- suppressWarnings(min(dat$plot_p[dat$plot_p > 0], na.rm = TRUE))
  if (!is.finite(min_positive)) min_positive <- 1e-300
  dat$plot_p[dat$plot_p <= 0] <- min_positive / 10
  dat$minus_log10_p <- -log10(dat$plot_p)

  sig <- dat[dat$significance != "NS", , drop = FALSE]
  if (nrow(sig) > 0L) {
    sig <- sig[order(sig$plot_p, -abs(sig$plot_fc)), , drop = FALSE]
    labels <- utils::head(sig, show_n)
  } else {
    labels <- dat[order(dat$plot_p, -abs(dat$plot_fc)), , drop = FALSE]
    labels <- utils::head(labels, show_n)
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$plot_fc, y = .data$minus_log10_p, colour = .data$significance)) +
    ggplot2::geom_point(alpha = 0.65, size = 1) +
    ggplot2::geom_vline(xintercept = c(-unique(dat$fc_threshold)[1], unique(dat$fc_threshold)[1]), linetype = 2) +
    ggplot2::geom_hline(yintercept = -log10(unique(dat$p_threshold)[1]), linetype = 2) +
    ggplot2::labs(title = title, x = fc_col, y = paste0("-log10(", p_col, ")"), colour = NULL) +
    ggplot2::theme_classic()

  if (nrow(labels) > 0L) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = labels,
        ggplot2::aes(label = .data$gene),
        size = 3,
        max.overlaps = Inf,
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = labels,
        ggplot2::aes(label = .data$gene),
        size = 2.5,
        check_overlap = TRUE,
        vjust = -0.4,
        show.legend = FALSE
      )
    }
  }
  p
}

.biomed_deg_plot_ora <- function(df, title, plot_cfg) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(df) == 0L) return(NULL)
  dat <- df
  if (isTRUE(plot_cfg$only_significant) && "p.adjust" %in% names(dat)) {
    dat <- dat[!is.na(dat$p.adjust) & dat$p.adjust <= plot_cfg$pvalue_cutoff, , drop = FALSE]
  }
  if (nrow(dat) == 0L) return(NULL)
  if (!"ID" %in% names(dat)) stop("ORA result does not contain an ID column for plotting.")

  dat$gene_ratio_numeric <- if ("GeneRatio" %in% names(dat)) .biomed_deg_ratio_to_numeric(dat$GeneRatio) else NA_real_
  dat$neg_log10_padj <- -log10(pmax(as.numeric(dat$p.adjust), 1e-300))

  dat <- dat[order(dat$p.adjust, -dat$gene_ratio_numeric), , drop = FALSE]
  dat <- utils::head(dat, plot_cfg$show_n)
  if (nrow(dat) == 0L) return(NULL)

  # Use enrichment ID, not Description, on the Y axis.
  # Keep the direction suffix because the same ID can occur in Up/Down/All ORA.
  dat$ID_label <- paste0(as.character(dat$ID), " [", dat$direction, "]")
  dat$ID_label <- stats::reorder(dat$ID_label, dat$gene_ratio_numeric)

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$gene_ratio_numeric, y = .data$ID_label)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$Count, colour = .data$neg_log10_padj)) +
    ggplot2::labs(title = title, x = "Gene ratio", y = "ID", size = "Count", colour = "-log10(adj. P)") +
    ggplot2::theme_classic()
}

.biomed_deg_plot_gsea <- function(df, title, plot_cfg) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(df) == 0L) return(NULL)
  dat <- df
  if (isTRUE(plot_cfg$only_significant) && "p.adjust" %in% names(dat)) {
    dat <- dat[!is.na(dat$p.adjust) & dat$p.adjust <= plot_cfg$pvalue_cutoff, , drop = FALSE]
  }
  if (nrow(dat) == 0L) return(NULL)
  if (!"ID" %in% names(dat)) stop("GSEA result does not contain an ID column for plotting.")

  dat$direction <- ifelse(dat$NES >= 0, "Positive", "Negative")
  dat <- dat[order(dat$p.adjust, -abs(dat$NES)), , drop = FALSE]
  dat <- utils::head(dat, plot_cfg$show_n)
  dat$ID_label <- stats::reorder(as.character(dat$ID), dat$NES)
  dat$neg_log10_padj <- -log10(pmax(as.numeric(dat$p.adjust), 1e-300))

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$NES, y = .data$ID_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(size = .data$setSize, colour = .data$neg_log10_padj)) +
    ggplot2::labs(title = title, x = "Normalized enrichment score (NES)", y = "ID",
                  size = "Set size", colour = "-log10(adj. P)") +
    ggplot2::theme_classic()
}

.biomed_deg_prepare_celltype_sets <- function(obj, celltype_col) {
  if (is.null(celltype_col) || identical(tolower(celltype_col), "all")) {
    return(list(all = colnames(obj)))
  }
  if (!celltype_col %in% colnames(obj[[]])) {
    stop("celltype_col not found in Seurat metadata: ", celltype_col)
  }
  values <- as.character(obj[[]][[celltype_col]])
  keep <- !is.na(values) & nzchar(values)
  split(colnames(obj)[keep], values[keep])
}

.biomed_deg_maybe_join_layers <- function(obj, assay, do_join = TRUE) {
  if (!isTRUE(do_join)) return(obj)
  join_fun <- NULL
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    join_fun <- tryCatch(getExportedValue("SeuratObject", "JoinLayers"), error = function(e) NULL)
  }
  if (is.null(join_fun) && requireNamespace("Seurat", quietly = TRUE)) {
    join_fun <- tryCatch(getExportedValue("Seurat", "JoinLayers"), error = function(e) NULL)
  }
  if (is.null(join_fun)) return(obj)

  tryCatch(
    join_fun(obj, assay = assay),
    error = function(e) {
      message("JoinLayers skipped: ", conditionMessage(e))
      obj
    }
  )
}

.biomed_deg_run_ora_one <- function(genes, universe, gmt, cfg) {
  if (length(genes) == 0L) return(NULL)
  clusterProfiler::enricher(
    gene = unique(genes),
    universe = if (isTRUE(cfg$use_tested_genes_as_universe)) unique(universe) else NULL,
    TERM2GENE = gmt$term2gene,
    TERM2NAME = gmt$term2name,
    pvalueCutoff = cfg$pvalue_cutoff,
    pAdjustMethod = cfg$p_adjust_method,
    qvalueCutoff = cfg$qvalue_cutoff,
    minGSSize = cfg$min_gs_size,
    maxGSSize = cfg$max_gs_size
  )
}

.biomed_deg_run_gsea_one <- function(ranked, gmt, cfg) {
  if (length(ranked) < cfg$min_gs_size) return(NULL)
  clusterProfiler::GSEA(
    geneList = ranked,
    exponent = cfg$exponent,
    minGSSize = cfg$min_gs_size,
    maxGSSize = cfg$max_gs_size,
    eps = cfg$eps,
    pvalueCutoff = cfg$pvalue_cutoff,
    pAdjustMethod = cfg$p_adjust_method,
    TERM2GENE = gmt$term2gene,
    TERM2NAME = gmt$term2name,
    verbose = FALSE,
    seed = cfg$seed,
    by = cfg$by
  )
}

#' Run two-group differential expression and GMT enrichment
#'
#' Run Seurat differential expression, optional ORA and GSEA, and export tables, plots and enrichment objects by cell type.
#'
#' @param input_file An in-memory Seurat object or a path to a QS, RDS, RDA or RData file. RDA/RData must contain exactly one Seurat object.
#' @param group_col Metadata column defining comparison groups.
#' @param ident_1 Comparison group used as numerator for differential expression.
#' @param ident_2 Reference group for differential expression.
#' @param celltype_col Metadata column identifying cell types. For DEG, NULL or 'all' combines all cells.
#' @param config Nested configuration list, YAML file path, or NULL for defaults.
#' @param output_dir Output directory; defaults and NULL behavior are shown in Usage. NULL disables exports for gene-expression statistics and selects the configuration directory for DEG.
#' @param gmt_dir Optional directory containing GMT gene sets; overrides config.
#' @param gmt_files Optional vector of GMT file paths; overrides config.
#' @param run_ora Optional logical override enabling or disabling ORA.
#' @param run_gsea Optional logical override enabling or disabling GSEA.
#' @details Positive fold changes refer to ident_1 relative to ident_2. Missing GMT configuration skips enrichment while retaining DEG analysis. Cells are the observations in FindMarkers; this workflow is not a replicate-level pseudobulk model. A status table records failed and skipped steps, so inspect it even when the function returns. Config sections are de, thresholds, gene_id, gmt, ora, gsea, plot and export; see the packaged example YAML under system.file('examples', 'sce_deg_config.yml', package = 'biomed'). Config accepts a nested list or YAML file. The default gene identifier is SYMBOL; differing identifier types require AnnotationDbi and the configured organism annotation package. Default exports include Excel tables, per-cell-type RDA enrichment objects and plots, PNG/PDF figures, the complete results RDA and the resolved configuration.
#' @return Invisibly, a list containing deg, ora, gsea, ora_objects, gsea_objects, ora_plots, gsea_plots, status, config and input/output metadata.
#' @export
sce_run_two_group_deg_enrichment <- function(
    input_file,
    group_col,
    ident_1,
    ident_2,
    celltype_col = "all",
    config = NULL,
    output_dir = NULL,
    gmt_dir = NULL,
    gmt_files = NULL,
    run_ora = NULL,
    run_gsea = NULL) {

  # Seurat is always required.
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Missing required package: Seurat")
  }

  cfg <- .biomed_deg_read_analysis_config(config)

  # Convenient RStudio overrides. These do not break the original config/YAML mode.
  if (!is.null(output_dir)) cfg$output_dir <- output_dir
  if (!is.null(gmt_dir)) cfg$gmt$dir <- gmt_dir
  if (!is.null(gmt_files)) cfg$gmt$files <- gmt_files
  if (!is.null(run_ora)) cfg$ora$run <- isTRUE(run_ora)
  if (!is.null(run_gsea)) cfg$gsea$run <- isTRUE(run_gsea)

  .biomed_require("withr")
  withr::local_seed(cfg$seed)

  obj <- .biomed_deg_read_seurat_object(input_file)
  if (inherits(input_file, "Seurat")) {
    message("[INPUT] Using in-memory Seurat object directly.")
  } else {
    message("[INPUT] Loaded Seurat object from: ", .biomed_deg_input_source_label(input_file))
  }
  meta <- obj[[]]

  if (!is.character(group_col) || length(group_col) != 1L || is.na(group_col) || !nzchar(group_col)) {
    stop("group_col must be one metadata column name.")
  }
  if (!group_col %in% colnames(meta)) stop("group_col not found in metadata: ", group_col)

  if (length(ident_1) != 1L || length(ident_2) != 1L) {
    stop("ident_1 and ident_2 must each be a single group value.")
  }
  if (identical(as.character(ident_1), as.character(ident_2))) {
    stop("ident_1 and ident_2 must be different groups.")
  }

  available_groups <- unique(as.character(meta[[group_col]]))
  available_groups <- available_groups[!is.na(available_groups)]
  missing_groups <- setdiff(c(as.character(ident_1), as.character(ident_2)), available_groups)
  if (length(missing_groups) > 0L) {
    stop("Requested group(s) not present in '", group_col, "': ", paste(missing_groups, collapse = ", "))
  }

  assay <- cfg$de$assay
  if (!assay %in% names(obj@assays)) stop("Assay not found in Seurat object: ", assay)
  Seurat::DefaultAssay(obj) <- assay
  obj <- .biomed_deg_maybe_join_layers(obj, assay, cfg$de$join_layers)

  if (isTRUE(cfg$de$prep_sct) && identical(assay, "SCT")) {
    obj <- Seurat::PrepSCTFindMarkers(obj, assay = assay, verbose = cfg$de$verbose)
  }

  enrichment_requested <- isTRUE(cfg$ora$run) || isTRUE(cfg$gsea$run)
  gmt_configured <- .biomed_deg_has_gmt_config(cfg$gmt)
  enrichment_enabled <- enrichment_requested && gmt_configured

  if (enrichment_requested && !gmt_configured) {
    message(
      "[INFO] No GMT directory/files configured. ",
      "Running DEG only; ORA and GSEA will be skipped. ",
      "To enable enrichment, set gmt_dir=... / gmt_files=... ",
      "or configure config$gmt."
    )
    cfg$ora$run <- FALSE
    cfg$gsea$run <- FALSE
  }

  if (enrichment_enabled) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      stop("Package 'clusterProfiler' is required when ORA or GSEA is enabled.")
    }
    gmt_files_resolved <- .biomed_deg_collect_gmt_files(cfg$gmt)
    gmt_list <- lapply(gmt_files_resolved, .biomed_deg_read_gmt_file)
    gmt_names <- tools::file_path_sans_ext(basename(gmt_files_resolved))
    names(gmt_list) <- make.unique(gmt_names)
  } else {
    gmt_files_resolved <- character()
    gmt_list <- list()
  }

  out_dir <- path.expand(cfg$output_dir)
  dirs <- c(
    out_dir,
    file.path(out_dir, "01_DEG"),
    file.path(out_dir, "01_DEG", "plots"),
    file.path(out_dir, "02_ORA"),
    file.path(out_dir, "02_ORA", "plots"),
    file.path(out_dir, "02_ORA", "objects"),
    file.path(out_dir, "03_GSEA"),
    file.path(out_dir, "03_GSEA", "plots"),
    file.path(out_dir, "03_GSEA", "objects")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  cell_sets <- .biomed_deg_prepare_celltype_sets(obj, celltype_col)
  deg_list <- list()
  ora_df_list <- list()
  gsea_df_list <- list()
  ora_objects <- list()
  gsea_objects <- list()
  ora_plots <- list()
  gsea_plots <- list()
  status <- list()

  comparison <- paste0(ident_1, "_vs_", ident_2)

  for (ct in names(cell_sets)) {
    message("[DEG] ", ct)

    ct_safe <- .biomed_deg_safe_filename(ct)
    ora_ct_objects <- list()
    gsea_ct_objects <- list()
    ora_ct_plots <- list()
    gsea_ct_plots <- list()

    if (isTRUE(cfg$ora$run)) {
      dir.create(file.path(out_dir, "02_ORA", "objects", ct_safe), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(out_dir, "02_ORA", "plots", ct_safe), recursive = TRUE, showWarnings = FALSE)
    }
    if (isTRUE(cfg$gsea$run)) {
      dir.create(file.path(out_dir, "03_GSEA", "objects", ct_safe), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(out_dir, "03_GSEA", "plots", ct_safe), recursive = TRUE, showWarnings = FALSE)
    }

    cells <- cell_sets[[ct]]
    group_values <- as.character(obj[[]][cells, group_col, drop = TRUE])
    keep <- !is.na(group_values) & group_values %in% c(as.character(ident_1), as.character(ident_2))
    cells <- cells[keep]
    group_values <- group_values[keep]
    n1 <- sum(group_values == as.character(ident_1))
    n2 <- sum(group_values == as.character(ident_2))

    if (n1 < cfg$de$min_cells_per_group || n2 < cfg$de$min_cells_per_group) {
      msg <- paste0("Skipped: insufficient cells (", ident_1, "=", n1, ", ", ident_2, "=", n2, ")")
      message("  ", msg)
      status[[length(status) + 1L]] <- data.frame(
        celltype = ct, step = "DEG", gmt = NA_character_, status = "skipped",
        message = msg, n_ident_1 = n1, n_ident_2 = n2, stringsAsFactors = FALSE
      )
      next
    }

    sub_obj <- subset(obj, cells = cells)
    find_args <- list(
      object = sub_obj,
      ident.1 = as.character(ident_1),
      ident.2 = as.character(ident_2),
      group.by = group_col,
      assay = assay,
      slot = cfg$de$slot,
      test.use = cfg$de$test_use,
      logfc.threshold = cfg$de$test_logfc_threshold,
      min.pct = cfg$de$min_pct,
      min.diff.pct = cfg$de$min_diff_pct,
      min.cells.feature = cfg$de$min_cells_feature,
      max.cells.per.ident = cfg$de$max_cells_per_ident,
      random.seed = cfg$seed,
      verbose = cfg$de$verbose
    )
    if (!is.null(cfg$de$latent_vars)) find_args$latent.vars <- cfg$de$latent_vars

    deg <- tryCatch(
      do.call(Seurat::FindMarkers, find_args),
      error = function(e) e
    )
    if (inherits(deg, "error")) {
      msg <- conditionMessage(deg)
      message("  Failed: ", msg)
      status[[length(status) + 1L]] <- data.frame(
        celltype = ct, step = "DEG", gmt = NA_character_, status = "failed",
        message = msg, n_ident_1 = n1, n_ident_2 = n2, stringsAsFactors = FALSE
      )
      next
    }

    deg$gene <- rownames(deg)
    rownames(deg) <- NULL
    fc_col <- .biomed_deg_detect_fc_column(deg, cfg$gsea$rank_column)
    p_col <- cfg$thresholds$p_column
    if (!p_col %in% colnames(deg)) stop("Configured p-value column not found: ", p_col)

    deg$celltype <- ct
    deg$comparison <- comparison
    deg$ident_1 <- as.character(ident_1)
    deg$ident_2 <- as.character(ident_2)
    deg$n_ident_1 <- n1
    deg$n_ident_2 <- n2
    deg$fc_threshold <- cfg$thresholds$fc
    deg$p_threshold <- cfg$thresholds$pvalue
    deg$significance <- "NS"
    deg$significance[!is.na(deg[[p_col]]) & deg[[p_col]] < cfg$thresholds$pvalue &
                       deg[[fc_col]] >= cfg$thresholds$fc] <- "Up"
    deg$significance[!is.na(deg[[p_col]]) & deg[[p_col]] < cfg$thresholds$pvalue &
                       deg[[fc_col]] <= -cfg$thresholds$fc] <- "Down"
    deg <- deg[order(deg[[p_col]], -abs(deg[[fc_col]])), , drop = FALSE]
    deg_list[[ct]] <- deg

    status[[length(status) + 1L]] <- data.frame(
      celltype = ct, step = "DEG", gmt = NA_character_, status = "success",
      message = paste0("tested_genes=", nrow(deg), "; Up=", sum(deg$significance == "Up"),
                       "; Down=", sum(deg$significance == "Down")),
      n_ident_1 = n1, n_ident_2 = n2, stringsAsFactors = FALSE
    )

    if (isTRUE(cfg$plot$run) && nrow(deg) > 0L) {
      p <- .biomed_deg_plot_volcano(
        deg = deg,
        fc_col = fc_col,
        p_col = p_col,
        title = paste0(ct, ": ", ident_1, " vs ", ident_2),
        show_n = cfg$plot$show_n
      )
      .biomed_deg_save_plot_formats(
        p,
        file.path(out_dir, "01_DEG", "plots", paste0(.biomed_deg_safe_filename(ct), "_volcano")),
        cfg$plot
      )
    }

    up_source <- deg$gene[deg$significance == "Up"]
    down_source <- deg$gene[deg$significance == "Down"]
    all_source <- unique(c(up_source, down_source))
    up_genes <- .biomed_deg_convert_gene_vector(up_source, cfg$gene_id)
    down_genes <- .biomed_deg_convert_gene_vector(down_source, cfg$gene_id)
    all_genes <- .biomed_deg_convert_gene_vector(all_source, cfg$gene_id)
    universe <- .biomed_deg_convert_gene_vector(deg$gene, cfg$gene_id)
    ranked <- .biomed_deg_make_ranked_gene_list(deg, fc_col, cfg$gene_id)

    for (gmt_name in names(gmt_list)) {
      gmt <- gmt_list[[gmt_name]]
      key <- paste(ct, gmt_name, sep = "::")

      if (isTRUE(cfg$ora$run)) {
        ora_direction_genes <- list(Up = up_genes, Down = down_genes, All = all_genes)
        ora_direction_genes <- ora_direction_genes[intersect(cfg$ora$directions, names(ora_direction_genes))]
        ora_per_gmt <- list()
        ora_obj_per_gmt <- list()

        for (direction in names(ora_direction_genes)) {
          genes <- ora_direction_genes[[direction]]
          ora_obj <- tryCatch(
            .biomed_deg_run_ora_one(genes, universe, gmt, cfg$ora),
            error = function(e) e
          )
          if (inherits(ora_obj, "error")) {
            status[[length(status) + 1L]] <- data.frame(
              celltype = ct, step = paste0("ORA_", direction), gmt = gmt_name,
              status = "failed", message = conditionMessage(ora_obj),
              n_ident_1 = n1, n_ident_2 = n2, stringsAsFactors = FALSE
            )
            next
          }
          ora_df <- .biomed_deg_enrich_result_df(ora_obj)
          ora_df <- .biomed_deg_add_enrichment_metadata(ora_df, ct, gmt_name, direction)
          ora_per_gmt[[direction]] <- ora_df
          ora_obj_per_gmt[[direction]] <- ora_obj
          status[[length(status) + 1L]] <- data.frame(
            celltype = ct, step = paste0("ORA_", direction), gmt = gmt_name,
            status = "success", message = paste0("terms=", nrow(ora_df), "; genes=", length(genes)),
            n_ident_1 = n1, n_ident_2 = n2, stringsAsFactors = FALSE
          )
        }

        ora_combined <- .biomed_deg_bind_rows_safe(ora_per_gmt)
        ora_df_list[[key]] <- ora_combined
        ora_objects[[key]] <- ora_obj_per_gmt
        ora_ct_objects[[gmt_name]] <- ora_obj_per_gmt

        if (isTRUE(cfg$plot$run) && nrow(ora_combined) > 0L) {
          p <- .biomed_deg_plot_ora(ora_combined, paste0(ct, " | ", gmt_name, " | ORA"), cfg$plot)
          if (!is.null(p)) {
            ora_ct_plots[[gmt_name]] <- p
            ora_plots[[key]] <- p
            .biomed_deg_save_plot_formats(
              p,
              file.path(
                out_dir, "02_ORA", "plots", ct_safe,
                paste0(.biomed_deg_safe_filename(gmt_name), "_ORA")
              ),
              cfg$plot
            )
          }
        }
      }

      if (isTRUE(cfg$gsea$run)) {
        gsea_obj <- tryCatch(
          .biomed_deg_run_gsea_one(ranked, gmt, cfg$gsea),
          error = function(e) e
        )
        if (inherits(gsea_obj, "error")) {
          status[[length(status) + 1L]] <- data.frame(
            celltype = ct, step = "GSEA", gmt = gmt_name, status = "failed",
            message = conditionMessage(gsea_obj), n_ident_1 = n1, n_ident_2 = n2,
            stringsAsFactors = FALSE
          )
        } else {
          gsea_df <- .biomed_deg_enrich_result_df(gsea_obj)
          gsea_df <- .biomed_deg_add_enrichment_metadata(gsea_df, ct, gmt_name)
          gsea_df_list[[key]] <- gsea_df
          gsea_objects[[key]] <- gsea_obj
          status[[length(status) + 1L]] <- data.frame(
            celltype = ct, step = "GSEA", gmt = gmt_name, status = "success",
            message = paste0("terms=", nrow(gsea_df), "; ranked_genes=", length(ranked)),
            n_ident_1 = n1, n_ident_2 = n2, stringsAsFactors = FALSE
          )
          gsea_ct_objects[[gmt_name]] <- gsea_obj

          if (isTRUE(cfg$plot$run) && nrow(gsea_df) > 0L) {
            p <- .biomed_deg_plot_gsea(gsea_df, paste0(ct, " | ", gmt_name, " | GSEA"), cfg$plot)
            if (!is.null(p)) {
              gsea_ct_plots[[gmt_name]] <- p
              gsea_plots[[key]] <- p
              .biomed_deg_save_plot_formats(
                p,
                file.path(
                  out_dir, "03_GSEA", "plots", ct_safe,
                  paste0(.biomed_deg_safe_filename(gmt_name), "_GSEA")
                ),
                cfg$plot
              )
            }
          }
        }
      }
    }

    # Save compact R workspace files grouped by cell type.
    # Each .rda contains all GMT objects/plots for the current cell type.
    if (isTRUE(cfg$export$write_rda)) {
      if (length(ora_ct_objects) > 0L) {
        .biomed_deg_save_named_rda(
          object = ora_ct_objects,
          object_name = paste0("ORA_by_", ct_safe),
          file = file.path(
            out_dir, "02_ORA", "objects", ct_safe,
            paste0("ORA_by_", ct_safe, ".rda")
          )
        )
      }

      if (length(gsea_ct_objects) > 0L) {
        .biomed_deg_save_named_rda(
          object = gsea_ct_objects,
          object_name = paste0("GSEA_by_", ct_safe),
          file = file.path(
            out_dir, "03_GSEA", "objects", ct_safe,
            paste0("GSEA_by_", ct_safe, ".rda")
          )
        )
      }

      if (length(ora_ct_plots) > 0L) {
        .biomed_deg_save_named_rda(
          object = ora_ct_plots,
          object_name = paste0("ORA_plots_by_", ct_safe),
          file = file.path(
            out_dir, "02_ORA", "plots", ct_safe,
            paste0("ORA_plots_by_", ct_safe, ".rda")
          )
        )
      }

      if (length(gsea_ct_plots) > 0L) {
        .biomed_deg_save_named_rda(
          object = gsea_ct_plots,
          object_name = paste0("GSEA_plots_by_", ct_safe),
          file = file.path(
            out_dir, "03_GSEA", "plots", ct_safe,
            paste0("GSEA_plots_by_", ct_safe, ".rda")
          )
        )
      }
    }

    # Save an Excel-readable companion next to each per-celltype enrichment .rda.
    # Every worksheet freezes the first row.
    if (isTRUE(cfg$export$write_excel)) {
      if (length(ora_ct_objects) > 0L) {
        ora_excel_sheets <- .biomed_deg_enrichment_objects_to_excel_sheets(
          ora_ct_objects,
          analysis = "ORA"
        )
        if (length(ora_excel_sheets) > 0L) {
          .biomed_deg_dflist2xlsx(
            ora_excel_sheets,
            file.path(
              out_dir, "02_ORA", "objects", ct_safe,
              paste0("ORA_by_", ct_safe, ".xlsx")
            ),
            row_names = cfg$export$row_names
          )
        }
      }

      if (length(gsea_ct_objects) > 0L) {
        gsea_excel_sheets <- .biomed_deg_enrichment_objects_to_excel_sheets(
          gsea_ct_objects,
          analysis = "GSEA"
        )
        if (length(gsea_excel_sheets) > 0L) {
          .biomed_deg_dflist2xlsx(
            gsea_excel_sheets,
            file.path(
              out_dir, "03_GSEA", "objects", ct_safe,
              paste0("GSEA_by_", ct_safe, ".xlsx")
            ),
            row_names = cfg$export$row_names
          )
        }
      }
    }
  }

  ora_all <- .biomed_deg_bind_rows_safe(ora_df_list)
  gsea_all <- .biomed_deg_bind_rows_safe(gsea_df_list)
  status_df <- .biomed_deg_bind_rows_safe(status)

  if (isTRUE(cfg$export$write_excel)) {
    if (length(deg_list) > 0L) {
      .biomed_deg_dflist2xlsx(
        deg_list,
        file.path(out_dir, "01_DEG", "DEG_by_celltype.xlsx"),
        row_names = cfg$export$row_names
      )
    }

    ora_by_celltype <- if (nrow(ora_all) > 0L) split(ora_all, ora_all$celltype) else list()
    if (length(ora_by_celltype) > 0L) {
      .biomed_deg_dflist2xlsx(
        ora_by_celltype,
        file.path(out_dir, "02_ORA", "ORA_by_celltype.xlsx"),
        row_names = cfg$export$row_names
      )
    }

    gsea_by_celltype <- if (nrow(gsea_all) > 0L) split(gsea_all, gsea_all$celltype) else list()
    if (length(gsea_by_celltype) > 0L) {
      .biomed_deg_dflist2xlsx(
        gsea_by_celltype,
        file.path(out_dir, "03_GSEA", "GSEA_by_celltype.xlsx"),
        row_names = cfg$export$row_names
      )
    }

    if (nrow(status_df) > 0L) {
      .biomed_deg_dflist2xlsx(
        list(run_status = status_df),
        file.path(out_dir, "run_status.xlsx"),
        row_names = FALSE
      )
    }
  }

  results <- list(
    deg = deg_list,
    ora = ora_all,
    gsea = gsea_all,
    ora_objects = ora_objects,
    gsea_objects = gsea_objects,
    ora_plots = ora_plots,
    gsea_plots = gsea_plots,
    status = status_df,
    config = cfg,
    input_source = .biomed_deg_input_source_label(input_file),
    input_type = if (inherits(input_file, "Seurat")) "Seurat_object" else "file",
    group_col = group_col,
    ident_1 = ident_1,
    ident_2 = ident_2,
    celltype_col = celltype_col,
    output_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  )

  if (isTRUE(cfg$export$write_rda)) {
    save(
      results,
      file = file.path(out_dir, "two_group_deg_enrichment_results.rda"),
      compress = TRUE
    )
  }

  if (requireNamespace("yaml", quietly = TRUE)) {
    yaml::write_yaml(cfg, file.path(out_dir, "config_used.yml"))
  } else {
    dput(cfg, file = file.path(out_dir, "config_used.R"))
  }

  message("Completed. Results written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  invisible(results)
}

