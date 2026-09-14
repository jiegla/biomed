# scTenifoldKnk virtual knockout wrapper for Seurat / count matrices
# Default nCores = 2
#
# Core outputs per KO gene:
#   00_run_parameters.txt
#   01_input_summary.xlsx
#   02_diffRegulation_all.xlsx
#   03_diffRegulation_significant.xlsx
#   04_top_perturbed_genes.xlsx
#   05_scTenifoldKnk_result.rds
#   06_analysis_bundle.rds
#   figures/*.pdf + *.png
#   enrichment/*.xlsx (optional)
#   sessionInfo.txt
#
# IMPORTANT:
# - scTenifoldKnk expects RAW COUNTS (genes x cells), not Seurat normalized "data".
# - scTenifoldKnk "FC" is a perturbation/statistical score derived from manifold
#   distance; it is NOT a conventional signed log2 fold-change.
# - Therefore this wrapper does not label genes as "up/down" based on FC.

#' Run virtual gene knockouts with scTenifoldKnk
#'
#' Select cells and raw-count genes, run one or more virtual knockouts, and export perturbation tables, figures, parameters and analysis bundles.
#'
#' @param object A Seurat object. For virtual knockout analysis, also accepts a genes-by-cells count matrix or an RDS/QS file containing one of these objects.
#' @param gKO One or more knockout genes; each is analyzed separately.
#' @param output_dir Output directory; defaults and NULL behavior are shown in Usage. NULL disables exports for gene-expression statistics and selects the configuration directory for DEG.
#' @param assay Assay name. For spatial plotting, NULL uses the default assay.
#' @param layer Expression layer or slot. Virtual knockout requires raw counts. Spatial plotting can choose an available layer when NULL.
#' @param subset_col Optional metadata column used to select a cell population.
#' @param subset_values Values to retain in subset_col.
#' @param cells Optional cell names to retain.
#' @param features Optional gene names to retain, in addition to KO genes.
#' @param prefilter_genes Filter genes by minimum number of expressing cells, retaining KO genes.
#' @param min_cells_gene Minimum expressing cells for the manual gene filter.
#' @param max_genes Optional cap ranked by detection frequency; KO genes are forced to remain. NULL disables the cap.
#' @param max_cells Optional maximum cell count; larger inputs are randomly downsampled.
#' @param min_ko_positive_cells Warn when a KO gene is detected in fewer cells.
#' @param ignore_gene_case Allow an unambiguous case-insensitive match of KO names.
#' @param qc Apply the backend quality-control filters.
#' @param qc_mtThreshold Maximum mitochondrial count fraction for backend QC.
#' @param qc_minLSize Minimum cell library size for backend QC.
#' @param qc_minCells Minimum expressing-cell count parameter for backend QC; see version adaptation in Details.
#' @param nc_lambda Network regularization parameter passed to scTenifoldKnk.
#' @param nc_nNet Number of subsampled networks.
#' @param nc_nCells Cells per network, capped at the number available.
#' @param nc_nComp Principal components for network construction; must exceed two.
#' @param nc_scaleScores Scale principal-component regression scores.
#' @param nc_symmetric Construct symmetric networks.
#' @param nc_q Network weight quantile threshold.
#' @param td_K Rank of the tensor decomposition.
#' @param td_maxIter Maximum tensor-decomposition iterations.
#' @param td_maxError Tensor-decomposition convergence tolerance.
#' @param td_nDecimal Decimal precision for tensor decomposition.
#' @param ma_nDim Number of manifold-alignment dimensions.
#' @param fdr_cutoff Adjusted P-value threshold for perturbed genes.
#' @param top_n Number of top genes reported and plotted.
#' @param label_n Number of genes labeled in the perturbation plot.
#' @param save_network_plot Save the backend KO-centered network plot.
#' @param network_q Edge-weight quantile for the KO network plot.
#' @param network_annotate Query enrichment services to annotate network nodes.
#' @param network_nCategories Maximum annotation categories in the network plot.
#' @param run_enrichment Run optional online Enrichr enrichment of significant perturbed genes.
#' @param enrichr_databases Enrichr database names; NULL selects available GO BP, KEGG and Reactome libraries.
#' @param seed Random seed, with the caller's random-number state restored afterwards.
#' @param nCores Positive integer number of backend worker cores.
#' @details Uses raw nonnegative counts, not normalized expression. Each KO is run separately. The FC statistic measures perturbation and is not a signed gene-expression log fold change. Gene caps always retain requested KO genes and may therefore be exceeded by those forced genes. The original QC argument names remain supported: scTenifoldKnk 1.1 and later receive qc_maxMTratio, qc_minLibSize and qc_minPCT; the latter is qc_minCells divided by the number of selected input cells, capped at one. Filtering details remain those of the installed backend version. Random-number state is restored on return. Network annotation and optional enrichR analysis access the Enrichr service and submit gene identifiers when enabled. Default outputs include XLSX tables, PNG/PDF figures, raw results and analysis-bundle RDS files, parameters and summary TXT files, and session information.
#' @return Invisibly, a named list of per-gene bundles containing result, diffRegulation, significant, top, plots, manifold_data, enrichment, parameters, input_summary and run_dir.
#' @export
sce_run_scTenifoldKnk_KO <- function(
    object,
    gKO,
    output_dir = "scTenifoldKnk_KO",
    assay = "RNA",
    layer = "counts",

    # Optional cell subset -- especially useful for merged single-cell objects
    subset_col = NULL,
    subset_values = NULL,
    cells = NULL,

    # Gene/cell preprocessing
    features = NULL,
    prefilter_genes = TRUE,
    min_cells_gene = 25,
    max_genes = 5000,
    max_cells = NULL,
    min_ko_positive_cells = 10,
    ignore_gene_case = TRUE,

    # scTenifoldKnk QC
    # For an already-QC'ed Seurat object, qc = FALSE is usually appropriate.
    # Still extract RAW counts from layer="counts".
    qc = FALSE,
    qc_mtThreshold = 0.10,
    qc_minLSize = 1000,
    qc_minCells = 25,

    # Network construction
    nc_lambda = 0,
    nc_nNet = 10,
    nc_nCells = 500,
    nc_nComp = 3,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = 0.90,

    # Tensor decomposition
    td_K = 3,
    td_maxIter = 1000,
    td_maxError = 1e-5,
    td_nDecimal = 3,

    # Manifold alignment
    ma_nDim = 2,

    # Statistics / plots
    fdr_cutoff = 0.05,
    top_n = 20,
    label_n = 10,

    # Official KO-centered network plot
    save_network_plot = TRUE,
    network_q = 0.99,
    network_annotate = FALSE,
    network_nCategories = 20,

    # Optional Enrichr enrichment of significantly perturbed genes
    run_enrichment = FALSE,
    enrichr_databases = NULL,

    # Reproducibility / compute
    seed = 1234,
    nCores = 2
) {

  # ---------------------------
  # 0. Helpers
  # ---------------------------
  `%||%` <- function(x, y) if (is.null(x)) y else x

  require_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", pkg, "' is required but not installed.\n",
        "Install it first, e.g. install.packages('", pkg, "').",
        call. = FALSE
      )
    }
  }

  write_xlsx_safe <- function(x, file, sheet_name = "Sheet1") {
    require_pkg("openxlsx")
    sheet_name <- substr(sheet_name, 1, 31)

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = sheet_name)

    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = x,
      startRow = 1,
      startCol = 1,
      rowNames = FALSE,
      colNames = TRUE,
      withFilter = TRUE
    )

    openxlsx::freezePane(
      wb,
      sheet = sheet_name,
      firstActiveRow = 2
    )

    if (ncol(x) > 0) {
      openxlsx::setColWidths(
        wb,
        sheet = sheet_name,
        cols = seq_len(ncol(x)),
        widths = "auto"
      )
    }

    openxlsx::saveWorkbook(
      wb,
      file = file,
      overwrite = TRUE
    )
  }

  sanitize <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    gsub("^_|_$", "", x)
  }

  save_gg <- function(p, stem, width = 8, height = 6) {
    require_pkg("ggplot2")
    ggplot2::ggsave(
      filename = paste0(stem, ".pdf"),
      plot = p, width = width, height = height, units = "in",
      limitsize = FALSE
    )
    ggplot2::ggsave(
      filename = paste0(stem, ".png"),
      plot = p, width = width, height = height, units = "in",
      dpi = 300, limitsize = FALSE
    )
  }

  add_labels <- function(p, dat, xvar, yvar, label_var = "gene") {
    if (nrow(dat) == 0) return(p)

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p + ggrepel::geom_text_repel(
        data = dat,
        ggplot2::aes(
          x = .data[[xvar]],
          y = .data[[yvar]],
          label = .data[[label_var]]
        ),
        size = 3.5,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.20,
        seed = seed
      )
    } else {
      p + ggplot2::geom_text(
        data = dat,
        ggplot2::aes(
          x = .data[[xvar]],
          y = .data[[yvar]],
          label = .data[[label_var]]
        ),
        size = 3.2,
        vjust = -0.7,
        check_overlap = TRUE
      )
    }
  }

  read_object_auto <- function(x) {
    if (!is.character(x) || length(x) != 1 || !file.exists(x)) {
      return(x)
    }

    ext <- tolower(tools::file_ext(x))
    message("[Input] Reading: ", x)

    if (ext == "rds") {
      return(readRDS(x))
    }
    if (ext == "qs") {
      if (!requireNamespace("qs", quietly = TRUE)) {
        stop("Legacy .qs input requires 'qs' on a compatible older R version. Read it there with qs::qread() and saveRDS() for current R; qs2 cannot read this format.")
      }
      return(qs::qread(x))
    }

    stop("Unsupported object file extension: .", ext,
         ". Please provide a Seurat object, count matrix, .rds, or .qs.")
  }

  resolve_gene_names <- function(requested, available, ignore_case = TRUE) {
    out <- character(length(requested))

    for (i in seq_along(requested)) {
      g <- requested[[i]]

      if (g %in% available) {
        out[[i]] <- g
        next
      }

      if (ignore_case) {
        idx <- which(toupper(available) == toupper(g))
        if (length(idx) == 1) {
          message("[Gene] '", g, "' matched to '", available[idx], "' by case-insensitive matching.")
          out[[i]] <- available[idx]
          next
        }
        if (length(idx) > 1) {
          stop("Ambiguous case-insensitive match for KO gene: ", g)
        }
      }

      stop("KO gene '", g, "' is not present in the count matrix.")
    }

    unique(out)
  }

  extract_counts_from_seurat <- function(obj, assay, layer) {
    require_pkg("SeuratObject")

    if (!assay %in% names(obj@assays)) {
      stop("Assay '", assay, "' not found. Available assays: ",
           paste(names(obj@assays), collapse = ", "))
    }

    lyr <- tryCatch(
      SeuratObject::Layers(obj[[assay]]),
      error = function(e) character(0)
    )

    # Seurat v5: exact requested layer exists
    if (layer %in% lyr) {
      return(SeuratObject::LayerData(obj, assay = assay, layer = layer))
    }

    # Merged Seurat v5 object may contain counts.sampleA, counts.sampleB, ...
    candidate_layers <- grep(
      paste0("^", gsub("\\.", "\\\\.", layer), "(\\.|$)"),
      lyr,
      value = TRUE
    )

    if (length(candidate_layers) > 1) {
      message(
        "[Seurat v5] Multiple '", layer,
        "' layers detected; joining them before extraction: ",
        paste(candidate_layers, collapse = ", ")
      )

      obj2 <- SeuratObject::JoinLayers(
        object = obj,
        assay = assay,
        layers = candidate_layers,
        new = layer
      )
      return(SeuratObject::LayerData(obj2, assay = assay, layer = layer))
    }

    # Backward-compatible fallback for Seurat v3/v4
    if (requireNamespace("Seurat", quietly = TRUE)) {
      out <- tryCatch(
        Seurat::GetAssayData(obj, assay = assay, slot = layer),
        error = function(e) NULL
      )
      if (!is.null(out)) return(out)
    }

    stop(
      "Could not extract layer/slot '", layer, "' from assay '", assay, "'.\n",
      "Available layers: ", paste(lyr, collapse = ", ")
    )
  }

  make_manifold_plot <- function(result, dr, fig_dir, gene, top_n) {
    require_pkg("ggplot2")

    ma <- result$manifoldAlignment
    if (is.null(ma) || ncol(ma) < 2 || is.null(rownames(ma))) return(NULL)

    rn <- rownames(ma)
    x_rows <- grep("^X_", rn)
    y_rows <- grep("^Y_", rn)

    if (length(x_rows) == 0 || length(y_rows) == 0) return(NULL)

    x_genes <- sub("^X_", "", rn[x_rows])
    y_genes <- sub("^Y_", "", rn[y_rows])
    common <- intersect(x_genes, y_genes)
    if (length(common) == 0) return(NULL)

    top_genes <- utils::head(dr$gene[order(dr$p.adj, -dr$FC, na.last = NA)], top_n)
    top_genes <- intersect(top_genes, common)
    if (length(top_genes) == 0) return(NULL)

    x_idx <- match(paste0("X_", top_genes), rn)
    y_idx <- match(paste0("Y_", top_genes), rn)

    disp <- data.frame(
      gene = top_genes,
      WT_1 = as.numeric(ma[x_idx, 1]),
      WT_2 = as.numeric(ma[x_idx, 2]),
      KO_1 = as.numeric(ma[y_idx, 1]),
      KO_2 = as.numeric(ma[y_idx, 2]),
      stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(disp) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = .data$WT_1, y = .data$WT_2,
          xend = .data$KO_1, yend = .data$KO_2
        ),
        arrow = grid::arrow(length = grid::unit(0.10, "inches")),
        alpha = 0.65
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$WT_1, y = .data$WT_2),
        shape = 1, size = 2.5
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$KO_1, y = .data$KO_2),
        shape = 16, size = 2.5
      ) +
      ggplot2::labs(
        x = "Manifold dimension 1",
        y = "Manifold dimension 2",
        title = paste0(gene, " virtual KO: manifold displacement"),
        subtitle = "Open circle = WT; filled circle = virtual KO"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5)
      )

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(x = .data$KO_1, y = .data$KO_2, label = .data$gene),
        size = 3.2, max.overlaps = Inf, seed = seed
      )
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(x = .data$KO_1, y = .data$KO_2, label = .data$gene),
        size = 3, check_overlap = TRUE, vjust = -0.6
      )
    }

    stem <- file.path(fig_dir, "04_manifold_displacement_top_genes")
    save_gg(p, stem, width = 8, height = 7)

    list(plot = p, data = disp)
  }

  save_official_network_plot <- function(result, gene, fig_dir) {
    if (!save_network_plot) return(invisible(NULL))

    draw_fun <- function() {
      scTenifoldKnk::plotKO(
        X = result,
        gKO = gene,
        q = network_q,
        annotate = network_annotate,
        nCategories = network_nCategories,
        fdrThreshold = fdr_cutoff
      )
    }

    # PDF
    tryCatch({
      grDevices::pdf(
        file.path(fig_dir, "05_KO_centered_network.pdf"),
        width = 10, height = 8, onefile = TRUE
      )
      draw_fun()
      grDevices::dev.off()
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      warning("Could not save KO network PDF: ", conditionMessage(e))
    })

    # PNG
    tryCatch({
      grDevices::png(
        file.path(fig_dir, "05_KO_centered_network.png"),
        width = 3000, height = 2400, res = 300
      )
      draw_fun()
      grDevices::dev.off()
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      warning("Could not save KO network PNG: ", conditionMessage(e))
    })

    invisible(NULL)
  }

  run_enrichr <- function(genes, outdir) {
    if (!run_enrichment) return(NULL)

    if (length(genes) < 3) {
      warning("Fewer than 3 significant perturbed genes; enrichment skipped.")
      return(NULL)
    }

    if (!requireNamespace("enrichR", quietly = TRUE)) {
      warning("Package 'enrichR' is not installed; enrichment skipped.")
      return(NULL)
    }

    dbs <- enrichr_databases

    # If databases are not specified, choose current GO BP / KEGG / Reactome
    # libraries from Enrichr when available.
    if (is.null(dbs)) {
      db_info <- tryCatch(enrichR::listEnrichrDbs(), error = function(e) NULL)
      if (is.null(db_info) || !"libraryName" %in% colnames(db_info)) {
        warning("Could not query Enrichr database list; enrichment skipped.")
        return(NULL)
      }

      libs <- db_info$libraryName

      pick_last <- function(pattern) {
        z <- grep(pattern, libs, value = TRUE)
        if (length(z) == 0) return(character(0))
        utils::tail(sort(z), 1)
      }

      dbs <- unique(c(
        pick_last("^GO_Biological_Process_"),
        pick_last("^KEGG_"),
        pick_last("^Reactome_")
      ))
    }

    if (length(dbs) == 0) {
      warning("No Enrichr databases selected; enrichment skipped.")
      return(NULL)
    }

    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    enr <- tryCatch(
      enrichR::enrichr(genes, dbs),
      error = function(e) {
        warning("Enrichr failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(enr)) return(NULL)

    for (nm in names(enr)) {
      write_xlsx_safe(
        enr[[nm]],
        file = file.path(outdir, paste0(sanitize(nm), ".xlsx")),
        sheet_name = substr(sanitize(nm), 1, 31)
      )
    }

    enr
  }

  # ---------------------------
  # 1. Package checks
  # ---------------------------
  require_pkg("scTenifoldKnk")
  require_pkg("ggplot2")
  require_pkg("Matrix")
  require_pkg("openxlsx")

  if (!is.numeric(nCores) || length(nCores) != 1 || !is.finite(nCores) || nCores < 1 || nCores != floor(nCores)) {
    stop("nCores must be a positive integer.")
  }
  nCores <- as.integer(nCores)

  .biomed_require("withr")
  withr::local_seed(seed)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------
  # 2. Load object and select cells
  # ---------------------------
  object <- read_object_auto(object)

  is_seurat <- inherits(object, "Seurat")
  is_matrix_like <- is.matrix(object) ||
    inherits(object, "Matrix") ||
    inherits(object, "sparseMatrix")

  if (!is_seurat && !is_matrix_like) {
    stop("`object` must be a Seurat object, matrix/sparse Matrix, .rds, or .qs.")
  }

  subset_tag <- "all_cells"

  if (is_seurat) {
    meta <- object@meta.data
    selected_cells <- colnames(object)

    if (!is.null(subset_col)) {
      if (!subset_col %in% colnames(meta)) {
        stop(
          "subset_col '", subset_col, "' not found in Seurat metadata.\n",
          "Available columns include: ",
          paste(utils::head(colnames(meta), 30), collapse = ", ")
        )
      }
      if (is.null(subset_values) || length(subset_values) == 0) {
        stop("When subset_col is provided, subset_values must also be provided.")
      }

      keep <- rownames(meta)[!is.na(meta[[subset_col]]) &
                              meta[[subset_col]] %in% subset_values]
      selected_cells <- intersect(selected_cells, keep)

      subset_tag <- paste0(
        sanitize(subset_col), "_",
        sanitize(paste(subset_values, collapse = "-"))
      )
    }

    if (!is.null(cells)) {
      missing_cells <- setdiff(cells, colnames(object))
      if (length(missing_cells) > 0) {
        warning(length(missing_cells), " requested cells were not found and were ignored.")
      }
      selected_cells <- intersect(selected_cells, cells)
    }

    if (length(selected_cells) == 0) {
      stop("No cells remain after subsetting.")
    }

    message(
      "[Subset] ", subset_tag,
      " | selected cells = ", length(selected_cells)
    )

    if (length(selected_cells) < 20) {
      stop(
        "Only ", length(selected_cells),
        " cells were selected for ", subset_tag,
        ". scTenifoldKnk requires a substantially larger cell population."
      )
    }

    if (!is.null(max_cells) && length(selected_cells) > max_cells) {
      selected_cells <- sample(selected_cells, max_cells)
      message("[Cells] Randomly downsampled to ", length(selected_cells), " cells.")
    }

    # Subset before joining layers to reduce memory use
    object_sub <- object[, selected_cells]

    countMatrix <- extract_counts_from_seurat(
      obj = object_sub,
      assay = assay,
      layer = layer
    )

    # Capture Seurat variable features for optional manual use/debugging
    seurat_variable_features <- tryCatch(
      SeuratObject::VariableFeatures(object_sub[[assay]]),
      error = function(e) character(0)
    )

  } else {
    countMatrix <- object
    seurat_variable_features <- character(0)

    if (!is.null(subset_col) || !is.null(subset_values)) {
      stop("subset_col/subset_values require a Seurat object with metadata.")
    }

    if (!is.null(cells)) {
      keep <- intersect(cells, colnames(countMatrix))
      if (length(keep) == 0) stop("None of the requested cells are present in the matrix.")
      countMatrix <- countMatrix[, keep, drop = FALSE]
    }

    if (!is.null(max_cells) && ncol(countMatrix) > max_cells) {
      keep <- sample(colnames(countMatrix), max_cells)
      countMatrix <- countMatrix[, keep, drop = FALSE]
      message("[Cells] Randomly downsampled to ", ncol(countMatrix), " cells.")
    }
  }

  if (is.null(rownames(countMatrix)) || is.null(colnames(countMatrix))) {
    stop("Count matrix must have gene rownames and cell colnames.")
  }

  if (anyDuplicated(rownames(countMatrix))) {
    stop("Duplicated gene names detected in count matrix. Make rownames unique first.")
  }

  # ------------------------------------------------------------------
  # Standardize extracted counts to a true 2D matrix-like object.
  # This protects against Seurat v5 / alternative layer backends that
  # can otherwise cause base::rowSums() to see a one-dimensional object.
  # ------------------------------------------------------------------
  message("[Matrix] extracted class: ", paste(class(countMatrix), collapse = ", "))
  message("[Matrix] extracted dim: ", paste(dim(countMatrix), collapse = " x "))

  if (is.null(dim(countMatrix)) || length(dim(countMatrix)) != 2L) {
    stop(
      "The extracted countMatrix is not two-dimensional. ",
      "Please check the selected subset and RNA counts layer."
    )
  }

  if (nrow(countMatrix) < 2L || ncol(countMatrix) < 2L) {
    stop(
      "The selected subset contains too few genes/cells for scTenifoldKnk: ",
      nrow(countMatrix), " genes x ", ncol(countMatrix), " cells."
    )
  }

  # Prefer a standard sparse dgCMatrix for memory efficiency.
  if (!inherits(countMatrix, "dgCMatrix")) {
    countMatrix <- tryCatch(
      methods::as(countMatrix, "dgCMatrix"),
      error = function(e1) {
        tryCatch(
          Matrix::Matrix(as.matrix(countMatrix), sparse = TRUE),
          error = function(e2) {
            stop(
              "Failed to convert countMatrix to a standard sparse matrix.\n",
              "Original class: ", paste(class(countMatrix), collapse = ", "), "\n",
              "Conversion error: ", conditionMessage(e2)
            )
          }
        )
      }
    )
  }

  message("[Matrix] standardized class: ", paste(class(countMatrix), collapse = ", "))
  message("[Matrix] standardized dim: ", nrow(countMatrix), " genes x ", ncol(countMatrix), " cells")

  # ---------------------------
  # 3. Validate raw counts
  # ---------------------------
  vals <- if (inherits(countMatrix, "sparseMatrix")) {
    countMatrix@x
  } else {
    as.numeric(countMatrix)
  }

  if (length(vals) > 100000) {
    vals <- sample(vals, 100000)
  }
  vals <- vals[is.finite(vals)]

  if (any(vals < 0)) {
    stop("Negative values detected. scTenifoldKnk should be run on raw non-negative counts.")
  }

  if (length(vals) > 0 && any(abs(vals - round(vals)) > 1e-8)) {
    warning(
      "Non-integer values detected in the selected layer. ",
      "scTenifoldKnk is designed for raw counts. ",
      "Use layer='counts', not normalized 'data' or 'scale.data'."
    )
  }

  # ---------------------------
  # 4. Resolve KO genes
  # ---------------------------
  gKO <- as.character(gKO)
  gKO <- gKO[nzchar(gKO)]

  if (length(gKO) == 0) stop("Please provide at least one KO gene.")

  gKO_resolved <- resolve_gene_names(
    requested = gKO,
    available = rownames(countMatrix),
    ignore_case = ignore_gene_case
  )

  # ---------------------------
  # 5. Manual gene filtering
  # ---------------------------
  detected_cells <- Matrix::rowSums(countMatrix != 0)

  if (prefilter_genes) {
    keep <- detected_cells >= min_cells_gene |
      rownames(countMatrix) %in% gKO_resolved

    message(
      "[Genes] min_cells_gene filtering: ",
      sum(keep), " / ", nrow(countMatrix), " genes retained."
    )

    countMatrix <- countMatrix[keep, , drop = FALSE]
    detected_cells <- detected_cells[keep]
  }

  if (!is.null(features)) {
    features <- unique(as.character(features))
    keep_features <- intersect(features, rownames(countMatrix))
    keep_features <- unique(c(keep_features, gKO_resolved))
    countMatrix <- countMatrix[
      intersect(keep_features, rownames(countMatrix)),
      ,
      drop = FALSE
    ]
    detected_cells <- detected_cells[rownames(countMatrix)]

    message("[Genes] User-specified feature set: ", nrow(countMatrix), " genes retained.")
  }

  # Practical cap for large merged datasets.
  # Ranking is by number of cells expressing each gene, while always retaining KO genes.
  if (!is.null(max_genes) && nrow(countMatrix) > max_genes) {
    if (max_genes < max(length(gKO_resolved) + 10, 100)) {
      stop("max_genes is too small for stable network analysis.")
    }

    ord <- order(detected_cells, decreasing = TRUE)
    selected <- rownames(countMatrix)[ord[seq_len(min(max_genes, length(ord)))]]
    selected <- unique(c(gKO_resolved, selected))
    selected <- intersect(selected, rownames(countMatrix))

    countMatrix <- countMatrix[selected, , drop = FALSE]
    detected_cells <- detected_cells[selected]

    message(
      "[Genes] max_genes cap applied: ",
      nrow(countMatrix), " genes retained (KO gene(s) forced to remain)."
    )
  }

  if (ncol(countMatrix) < 20) {
    stop("Fewer than 20 cells remain. This is too small for this wrapper to proceed.")
  }
  if (ncol(countMatrix) < 100) {
    warning("Fewer than 100 cells are being used; network inference may be unstable.")
  }

  # KO expression sanity check
  ko_positive <- stats::setNames(
    vapply(
      gKO_resolved,
      function(g) sum(countMatrix[g, ] > 0),
      numeric(1)
    ),
    gKO_resolved
  )

  for (g in gKO_resolved) {
    if (ko_positive[[g]] == 0) {
      stop("KO gene '", g, "' has zero detected cells in the selected population.")
    }
    if (ko_positive[[g]] < min_ko_positive_cells) {
      warning(
        "KO gene '", g, "' is detected in only ", ko_positive[[g]],
        " cells (< min_ko_positive_cells=", min_ko_positive_cells,
        "). Interpret results cautiously."
      )
    }
  }

  # scTenifoldKnk cannot sample more cells per network than available.
  nc_nCells_effective <- min(as.integer(nc_nCells), ncol(countMatrix))
  if (nc_nCells_effective != nc_nCells) {
    message(
      "[Network] nc_nCells reduced from ", nc_nCells,
      " to ", nc_nCells_effective,
      " because fewer cells are available."
    )
  }

  if (nc_nComp <= 2 || nc_nComp >= nrow(countMatrix)) {
    stop("nc_nComp must be > 2 and < number of retained genes.")
  }

  # Shared input summary
  input_summary <- data.frame(
    item = c(
      "subset", "assay", "layer", "n_cells", "n_genes",
      "min_cells_gene", "max_genes", "nCores", "seed",
      paste0("KO_positive_cells:", gKO_resolved)
    ),
    value = c(
      subset_tag, assay, layer, ncol(countMatrix), nrow(countMatrix),
      min_cells_gene,
      ifelse(is.null(max_genes), "NULL", as.character(max_genes)),
      nCores, seed,
      as.character(ko_positive[gKO_resolved])
    ),
    stringsAsFactors = FALSE
  )

  write_xlsx_safe(
    input_summary,
    file.path(output_dir, "input_summary_shared.xlsx"),
    sheet_name = "input_summary"
  )

  # ---------------------------
  # 6. Run one or more KO genes
  # ---------------------------
  all_runs <- vector("list", length(gKO_resolved))
  names(all_runs) <- gKO_resolved

  for (gene in gKO_resolved) {

    message("\n==========================================")
    message("Running virtual knockout: ", gene)
    message("Cells: ", ncol(countMatrix), " | Genes: ", nrow(countMatrix),
            " | nCores: ", nCores)
    message("==========================================")

    run_dir <- file.path(
      output_dir,
      paste0("KO_", sanitize(gene), "__", subset_tag)
    )
    fig_dir <- file.path(run_dir, "figures")
    enr_dir <- file.path(run_dir, "enrichment")

    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    set.seed(seed)

    t0 <- Sys.time()

    knk_args <- list(
      countMatrix = countMatrix,
      qc = qc,
      gKO = gene,
      nc_lambda = nc_lambda,
      nc_nNet = nc_nNet,
      nc_nCells = nc_nCells_effective,
      nc_nComp = nc_nComp,
      nc_scaleScores = nc_scaleScores,
      nc_symmetric = nc_symmetric,
      nc_q = nc_q,
      td_K = td_K,
      td_maxIter = td_maxIter,
      td_maxError = td_maxError,
      td_nDecimal = td_nDecimal,
      ma_nDim = ma_nDim,
      nCores = nCores
    )
    knk_args <- .biomed_knk_qc_arguments(
      knk_args, names(formals(scTenifoldKnk::scTenifoldKnk)),
      qc_mtThreshold, qc_minLSize, qc_minCells
    )
    result <- do.call(scTenifoldKnk::scTenifoldKnk, knk_args)

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

    if (is.null(result$diffRegulation)) {
      stop("scTenifoldKnk returned no $diffRegulation table for ", gene, ".")
    }

    dr <- as.data.frame(result$diffRegulation)

    required_cols <- c("gene", "distance", "FC", "p.value", "p.adj")
    miss <- setdiff(required_cols, colnames(dr))
    if (length(miss) > 0) {
      stop(
        "Unexpected scTenifoldKnk output. Missing columns: ",
        paste(miss, collapse = ", ")
      )
    }

    # Robust plotting/statistics columns
    pos_padj <- dr$p.adj[is.finite(dr$p.adj) & dr$p.adj > 0]
    padj_floor <- if (length(pos_padj) > 0) min(pos_padj) / 10 else 1e-300

    dr$padj_plot <- dr$p.adj
    dr$padj_plot[!is.finite(dr$padj_plot) | dr$padj_plot <= 0] <- padj_floor
    dr$neg_log10_padj <- -log10(dr$padj_plot)

    # Preserve the original FC statistic, but make a finite plotting copy.
    finite_fc <- dr$FC[is.finite(dr$FC)]
    fc_cap <- if (length(finite_fc) > 0) max(finite_fc, na.rm = TRUE) else 1
    if (!is.finite(fc_cap) || fc_cap <= 0) fc_cap <- 1
    dr$FC_plot <- dr$FC
    dr$FC_plot[is.infinite(dr$FC_plot) & dr$FC_plot > 0] <- fc_cap * 1.05
    dr$FC_plot[is.infinite(dr$FC_plot) & dr$FC_plot < 0] <- -fc_cap * 1.05
    dr$FC_plot[is.na(dr$FC_plot)] <- 0

    dr$significant <- !is.na(dr$p.adj) & dr$p.adj < fdr_cutoff
    dr$KO_gene <- gene
    dr$n_cells_input <- ncol(countMatrix)
    dr$n_genes_input <- nrow(countMatrix)

    # Rank by adjusted P first, then perturbation FC
    dr <- dr[order(dr$p.adj, -dr$FC, -dr$distance, na.last = TRUE), ]
    rownames(dr) <- NULL
    dr$rank <- seq_len(nrow(dr))

    sig <- dr[dr$significant, , drop = FALSE]
    top <- utils::head(dr, top_n)

    # Save tables/results
    saveRDS(
      result,
      file = file.path(run_dir, "05_scTenifoldKnk_result.rds")
    )

    write_xlsx_safe(
      dr,
      file.path(run_dir, "02_diffRegulation_all.xlsx"),
      sheet_name = "diffRegulation_all"
    )

    write_xlsx_safe(
      sig,
      file.path(run_dir, "03_diffRegulation_significant.xlsx"),
      sheet_name = "significant"
    )

    write_xlsx_safe(
      top,
      file.path(run_dir, "04_top_perturbed_genes.xlsx"),
      sheet_name = "top_perturbed"
    )

    # Parameters
    params <- c(
      paste0("KO_gene=", gene),
      paste0("subset=", subset_tag),
      paste0("assay=", assay),
      paste0("layer=", layer),
      paste0("n_cells=", ncol(countMatrix)),
      paste0("n_genes=", nrow(countMatrix)),
      paste0("KO_positive_cells=", ko_positive[[gene]]),
      paste0("qc=", qc),
      paste0("prefilter_genes=", prefilter_genes),
      paste0("min_cells_gene=", min_cells_gene),
      paste0("max_genes=", ifelse(is.null(max_genes), "NULL", max_genes)),
      paste0("nc_nNet=", nc_nNet),
      paste0("nc_nCells=", nc_nCells_effective),
      paste0("nc_nComp=", nc_nComp),
      paste0("nc_lambda=", nc_lambda),
      paste0("nc_q=", nc_q),
      paste0("td_K=", td_K),
      paste0("ma_nDim=", ma_nDim),
      paste0("fdr_cutoff=", fdr_cutoff),
      paste0("seed=", seed),
      paste0("nCores=", nCores),
      paste0("elapsed_minutes=", round(elapsed, 3)),
      paste0(
        "scTenifoldKnk_version=",
        as.character(utils::packageVersion("scTenifoldKnk"))
      )
    )
    writeLines(params, file.path(run_dir, "00_run_parameters.txt"))

    # ---------------------------
    # 7. Plot 1: perturbation significance bubble plot
    # ---------------------------
    label_dat <- utils::head(
      dr[order(dr$p.adj, -dr$FC, na.last = NA), , drop = FALSE],
      label_n
    )

    p1 <- ggplot2::ggplot(
      dr,
      ggplot2::aes(x = .data$FC_plot, y = .data$neg_log10_padj)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(size = .data$distance, color = .data$significant),
        alpha = 0.70
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(fdr_cutoff),
        linetype = "dashed",
        alpha = 0.6
      ) +
      ggplot2::labs(
        x = "scTenifoldKnk FC statistic",
        y = "-log10(FDR)",
        size = "Manifold distance",
        color = paste0("FDR < ", fdr_cutoff),
        title = paste0(gene, " virtual KO: perturbed genes"),
        subtitle = "FC is a perturbation statistic, not signed expression logFC"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 14
        ),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        legend.position = "right"
      )

    p1 <- add_labels(
      p = p1,
      dat = label_dat,
      xvar = "FC_plot",
      yvar = "neg_log10_padj"
    )

    save_gg(
      p1,
      file.path(fig_dir, "01_perturbation_significance_bubble"),
      width = 8.5, height = 6.5
    )

    # ---------------------------
    # 8. Plot 2: Top genes by perturbation FC
    # ---------------------------
    top_fc <- utils::head(
      dr[order(dr$FC, decreasing = TRUE, na.last = NA), , drop = FALSE],
      top_n
    )

    top_fc$gene_plot <- factor(
      top_fc$gene,
      levels = rev(top_fc$gene)
    )

    p2 <- ggplot2::ggplot(
      top_fc,
      ggplot2::aes(x = .data$FC_plot, y = .data$gene_plot)
    ) +
      ggplot2::geom_col(
        ggplot2::aes(fill = .data$neg_log10_padj),
        width = 0.72,
        alpha = 0.85
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", .data$FC)),
        hjust = -0.15,
        size = 3.2
      ) +
      ggplot2::scale_fill_gradient(
        low = "lightblue",
        high = "darkblue",
        name = "-log10(FDR)"
      ) +
      ggplot2::expand_limits(x = max(top_fc$FC_plot, na.rm = TRUE) * 1.18) +
      ggplot2::labs(
        x = "scTenifoldKnk FC statistic",
        y = "Gene",
        title = paste0(gene, " KO: Top ", nrow(top_fc), " perturbed genes")
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 14
        ),
        axis.text.y = ggplot2::element_text(face = "italic"),
        legend.position = "top",
        panel.grid.minor = ggplot2::element_blank()
      )

    save_gg(
      p2,
      file.path(fig_dir, "02_top_perturbed_genes_FC"),
      width = 8, height = max(6, 0.28 * nrow(top_fc) + 2)
    )

    # ---------------------------
    # 9. Plot 3: Top significant genes by FDR
    # ---------------------------
    sig_for_plot <- sig[sig$gene != gene, , drop = FALSE]

    if (nrow(sig_for_plot) > 0) {
      top_sig <- utils::head(
        sig_for_plot[
          order(sig_for_plot$p.adj, -sig_for_plot$FC, na.last = NA),
          ,
          drop = FALSE
        ],
        top_n
      )

      top_sig$gene_plot <- factor(
        top_sig$gene,
        levels = rev(top_sig$gene)
      )

      p3 <- ggplot2::ggplot(
        top_sig,
        ggplot2::aes(x = .data$neg_log10_padj, y = .data$gene_plot)
      ) +
        ggplot2::geom_col(
          ggplot2::aes(fill = .data$FC),
          width = 0.72,
          alpha = 0.85
        ) +
        ggplot2::labs(
          x = "-log10(FDR)",
          y = "Gene",
          fill = "FC statistic",
          title = paste0(
            gene, " KO: Top significant perturbed genes"
          )
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5, face = "bold", size = 14
          ),
          axis.text.y = ggplot2::element_text(face = "italic"),
          legend.position = "top",
          panel.grid.minor = ggplot2::element_blank()
        )

      save_gg(
        p3,
        file.path(fig_dir, "03_top_significant_genes_FDR"),
        width = 8,
        height = max(6, 0.28 * nrow(top_sig) + 2)
      )
    } else {
      p3 <- NULL
      message("[Plot] No genes pass FDR < ", fdr_cutoff,
              "; significant-gene bar plot skipped.")
    }

    # ---------------------------
    # 10. Plot 4: manifold displacement
    # ---------------------------
    manifold_out <- make_manifold_plot(
      result = result,
      dr = dr,
      fig_dir = fig_dir,
      gene = gene,
      top_n = top_n
    )

    # ---------------------------
    # 11. Plot 5: official KO-centered network
    # ---------------------------
    save_official_network_plot(
      result = result,
      gene = gene,
      fig_dir = fig_dir
    )

    # ---------------------------
    # 12. Optional enrichment
    # ---------------------------
    enrichment <- run_enrichr(
      genes = setdiff(sig$gene, gene),
      outdir = enr_dir
    )

    # ---------------------------
    # 13. Human-readable summary
    # ---------------------------
    summary_lines <- c(
      paste0("scTenifoldKnk virtual KO: ", gene),
      paste0("Subset: ", subset_tag),
      paste0("Cells used: ", ncol(countMatrix)),
      paste0("Genes used: ", nrow(countMatrix)),
      paste0("KO-positive cells: ", ko_positive[[gene]]),
      paste0("FDR cutoff: ", fdr_cutoff),
      paste0("Significantly perturbed genes: ", nrow(sig)),
      paste0("Elapsed time (min): ", round(elapsed, 2)),
      "",
      "Top perturbed genes:",
      paste(
        paste0(
          utils::head(dr$gene, min(20, nrow(dr))),
          " | FC=",
          signif(utils::head(dr$FC, min(20, nrow(dr))), 4),
          " | FDR=",
          signif(utils::head(dr$p.adj, min(20, nrow(dr))), 4)
        ),
        collapse = "\n"
      ),
      "",
      "Interpretation note:",
      "scTenifoldKnk FC is a non-directional perturbation statistic derived",
      "from WT-vs-KO manifold distance. It should not be interpreted as",
      "conventional signed gene-expression log2FC (up/down regulation)."
    )

    writeLines(
      summary_lines,
      file.path(run_dir, "analysis_summary.txt")
    )

    utils::capture.output(
      utils::sessionInfo(),
      file = file.path(run_dir, "sessionInfo.txt")
    )

    bundle <- list(
      gene = gene,
      run_dir = run_dir,
      input_summary = input_summary,
      result = result,
      diffRegulation = dr,
      significant = sig,
      top = top,
      plots = list(
        perturbation_bubble = p1,
        top_FC = p2,
        top_significant = p3,
        manifold = if (is.null(manifold_out)) NULL else manifold_out$plot
      ),
      manifold_data = if (is.null(manifold_out)) NULL else manifold_out$data,
      enrichment = enrichment,
      parameters = params
    )

    saveRDS(
      bundle,
      file = file.path(run_dir, "06_analysis_bundle.rds")
    )

    all_runs[[gene]] <- bundle

    message(
      "[Done] ", gene, " | significant genes: ", nrow(sig),
      " | output: ", normalizePath(run_dir, winslash = "/", mustWork = FALSE)
    )
  }

  invisible(all_runs)
}


# -------------------------------------------------------------------------
# EXAMPLES
# -------------------------------------------------------------------------
#
# 1) A merged Seurat object, analyze one cell type (recommended):
#
# res <- run_scTenifoldKnk_KO(
#   object = sce,
#   gKO = "CASP4",
#   subset_col = "celltype_major",
#   subset_values = "Macrophage",
#   output_dir = "scTenifoldKnk_CASP4",
#   assay = "RNA",
#   layer = "counts",
#   qc = FALSE,
#   max_genes = 5000,
#   nCores = 2
# )
#
# 2) Read directly from an RDS:
#
# res <- run_scTenifoldKnk_KO(
#   object = "/path/to/seurat_object.rds",
#   gKO = "CASP4",
#   output_dir = "KO_CASP4",
#   qc = FALSE,
#   nCores = 2
# )
#
# 3) Multiple virtual KOs (the wrapper loops over genes one-by-one):
#
# res <- run_scTenifoldKnk_KO(
#   object = sce,
#   gKO = c("CASP4", "NLRP3", "GSDMD"),
#   subset_col = "celltype_major",
#   subset_values = "Macrophage",
#   output_dir = "multi_KO",
#   nCores = 2
# )
#
# 4) Optional Enrichr enrichment:
#
# res <- run_scTenifoldKnk_KO(
#   object = sce,
#   gKO = "CASP4",
#   subset_col = "celltype_major",
#   subset_values = "Macrophage",
#   run_enrichment = TRUE,
#   nCores = 2
# )
#
# Notes:
# - Keep layer = "counts".
# - If your Seurat object has already passed QC, qc = FALSE is reasonable.
# - For a very large merged object, analyze within a biologically coherent
#   cell type/state rather than across all cell types.
# - max_genes = 5000 is a practical computational cap. Set max_genes = NULL
#   if you intentionally want to retain every gene and have enough memory.

# Adapt the contributed QC interface to scTenifoldKnk before/after version 1.1.
.biomed_knk_qc_arguments <- function(args, backend_formals, mt, min_library, min_cells) {
  if ("qc_minLibSize" %in% backend_formals) {
    args$qc_maxMTratio <- mt
    args$qc_minLibSize <- min_library
    args$qc_minPCT <- min(1, min_cells / ncol(args$countMatrix))
  } else {
    args$qc_mtThreshold <- mt
    args$qc_minLSize <- min_library
    args$qc_minCells <- min_cells
  }
  args
}
