#' Batch-plot spatial pathway scores from RDS files
#'
#' Read pathway score tables, align barcodes and generate continuous spatial maps across selected pathways, RDS files and images.
#'
#' @param object A Seurat object. For virtual knockout analysis, also accepts a genes-by-cells count matrix or an RDS/QS file containing one of these objects.
#' @param score_dir Directory containing pathway score RDS files.
#' @param selected_rds RDS paths or filenames within score_dir; NULL selects all RDS files there.
#' @param selected_pathways Pathway names or patterns; NULL selects all numeric pathways.
#' @param selected_images Spatial image names; NULL selects all images.
#' @param output_dir Output directory; defaults and NULL behavior are shown in Usage. NULL disables exports for gene-expression statistics and selects the configuration directory for DEG.
#' @param pathway_match Match pathway names exactly, by case-insensitive literal substring, or by regular expression.
#' @param image_scale Spatial coordinate scale: lowres or hires.
#' @param smooth_k Number of nearest spots in Gaussian smoothing.
#' @param sigma_factor Gaussian bandwidth in multiples of median spot spacing.
#' @param mask_radius_factor Maximum distance from a measured spot in multiples of spot spacing.
#' @param grid_n Number of grid positions along x; y adapts to the aspect ratio.
#' @param colors A color vector or named palette: viridis, magma, inferno, plasma, cividis, turbo, blue_red, blue_yellow_red, white_red or navy_yellow_red.
#' @param show_he Overlay the stored H&E image when available.
#' @param he_alpha H&E image opacity.
#' @param surface_alpha Interpolated surface opacity.
#' @param crop Crop to the observed spot extent.
#' @param crop_padding Relative padding around the cropped extent.
#' @param xlim Optional x limits, overriding automatic cropping.
#' @param ylim Optional y limits, overriding automatic cropping.
#' @param contour Overlay contour lines.
#' @param contour_bins Number of contour levels.
#' @param show_spots Overlay original spot centers.
#' @param transform Value transformation: none, log1p, sqrt or zscore.
#' @param cap_quantile Two ordered clipping probabilities, or NULL to disable clipping.
#' @param save_pdf Save each spatial map as PDF.
#' @param save_png Save each spatial map as PNG.
#' @param width Saved plot width in inches.
#' @param height Saved plot height in inches.
#' @param dpi Resolution for raster output.
#' @param keep_plots Retain ggplot objects in the returned list.
#' @param verbose Print progress messages.
#' @details Accepts RDS matrices, data frames, Seurat metadata, or lists containing a score table. Spot-by-pathway and pathway-by-spot orientations are recognized by barcode matches, and barcode columns are supported. Only numeric score columns with finite values are used. Each pathway needs at least ten valid spots. Temporary metadata is added to a local copy. PNG and PDF files are saved by default, with pathway_plot_summary.csv and pathway_plot_log.csv when results exist. Calls sce_plot_spatial_continuous directly; no manual source call is needed. Missing files and individual failed plots are reported and skipped; inspect the returned log.
#' @return Invisibly, a list with plots (retained only when keep_plots is TRUE), summary, log and images.
#' @export
sce_plot_pathway_rds_continuous <- function(
    object,
    score_dir = "gmt_pathway_scores",
    selected_rds = NULL,
    selected_pathways = NULL,
    selected_images = NULL,
    output_dir = "Spatial_pathway_continuous",
    pathway_match = c("exact", "contains", "regex"),
    image_scale = c("lowres", "hires"),
    smooth_k = 8,
    sigma_factor = 0.8,
    mask_radius_factor = 1.3,
    grid_n = 300,
    colors = "magma",
    show_he = TRUE,
    he_alpha = 0.40,
    surface_alpha = 0.82,
    crop = TRUE,
    crop_padding = 0.02,
    xlim = NULL,
    ylim = NULL,
    contour = FALSE,
    contour_bins = 8,
    show_spots = FALSE,
    transform = "none",
    cap_quantile = c(0.01, 0.99),
    save_pdf = TRUE,
    save_png = TRUE,
    width = 7,
    height = 6,
    dpi = 300,
    keep_plots = FALSE,
    verbose = TRUE
) {
  
  pathway_match <- match.arg(pathway_match)
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.")
  image_scale <- match.arg(image_scale)
  
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("\u8bf7\u5148\u5b89\u88c5 ggplot2\u3002")
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  object_cells <- colnames(object)
  if (is.null(object_cells) || length(object_cells) == 0) {
    stop("object \u4e2d\u6ca1\u6709 cell / spot names\u3002")
  }
  
  # ============================================================
  # spatial images
  # ============================================================
  
  all_images <- names(object@images)
  if (length(all_images) == 0) stop("Seurat object \u4e2d\u6ca1\u6709 spatial images\u3002")
  
  if (is.null(selected_images)) {
    image_names <- all_images
  } else {
    missing_images <- setdiff(selected_images, all_images)
    if (length(missing_images) > 0) {
      warning("\u4ee5\u4e0b selected_images \u4e0d\u5b58\u5728\uff0c\u5c06\u8df3\u8fc7\uff1a\n", paste(missing_images, collapse = "\n"))
    }
    image_names <- intersect(selected_images, all_images)
    if (length(image_names) == 0) stop("selected_images \u4e2d\u6ca1\u6709\u6709\u6548 spatial image\u3002")
  }
  
  if (verbose) message("Spatial images to plot: ", length(image_names))
  
  # ============================================================
  # RDS files
  # ============================================================
  
  if (is.null(selected_rds)) {
    selected_rds <- list.files(score_dir, pattern = "\\.rds$", full.names = FALSE)
    if (length(selected_rds) == 0) stop("\u5728\u76ee\u5f55\u4e2d\u6ca1\u6709\u627e\u5230 RDS \u6587\u4ef6\uff1a", score_dir)
  }
  
  # ============================================================
  # helper functions
  # ============================================================
  
  make_safe_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_.-]", "_", x)
    gsub("_+", "_", x)
  }
  
  extract_score_table <- function(x, rds_file) {
    
    if (is.matrix(x)) return(as.data.frame(x, check.names = FALSE))
    if (is.data.frame(x)) return(as.data.frame(x, check.names = FALSE))
    if (inherits(x, "Seurat")) return(as.data.frame(x@meta.data, check.names = FALSE))
    
    if (is.list(x)) {
      
      candidates <- which(vapply(
        x,
        function(z) is.matrix(z) || is.data.frame(z),
        logical(1)
      ))
      
      if (length(candidates) == 0) {
        stop("RDS \u4e2d\u627e\u4e0d\u5230 matrix/data.frame\uff1a", rds_file)
      }
      
      common_names <- c(
        "scores", "score", "meta", "meta.data", "metadata",
        "result", "results", "ucell", "UCell"
      )
      
      if (!is.null(names(x))) {
        hit <- intersect(
          which(tolower(names(x)) %in% tolower(common_names)),
          candidates
        )
        
        if (length(hit) > 0) {
          use_id <- hit[1]
          if (verbose) message("RDS list -> using element: ", names(x)[use_id])
          return(as.data.frame(x[[use_id]], check.names = FALSE))
        }
      }
      
      match_number <- vapply(candidates, function(i) {
        z <- x[[i]]
        r_match <- if (!is.null(rownames(z))) sum(rownames(z) %in% object_cells) else 0
        c_match <- if (!is.null(colnames(z))) sum(colnames(z) %in% object_cells) else 0
        max(r_match, c_match)
      }, numeric(1))
      
      use_id <- candidates[which.max(match_number)]
      
      if (verbose) {
        nm <- if (!is.null(names(x)) && length(names(x)) >= use_id && nzchar(names(x)[use_id])) {
          names(x)[use_id]
        } else {
          paste0("[[", use_id, "]]")
        }
        message("RDS list -> using element: ", nm)
      }
      
      return(as.data.frame(x[[use_id]], check.names = FALSE))
    }
    
    stop("\u6682\u4e0d\u652f\u6301 RDS \u7c7b\u578b\uff1a", paste(class(x), collapse = "/"))
  }
  
  standardize_score_table <- function(score_df, rds_file) {
    
    score_df <- as.data.frame(score_df, check.names = FALSE)
    
    barcode_candidates <- c(
      "cell", "cells", "barcode", "barcodes",
      "spot", "spots", "cell_id", "spot_id"
    )
    
    hit_barcode <- which(tolower(colnames(score_df)) %in% barcode_candidates)
    
    if (length(hit_barcode) > 0) {
      bc_col <- hit_barcode[1]
      possible_bc <- as.character(score_df[[bc_col]])
      if (sum(possible_bc %in% object_cells) > 0) {
        rownames(score_df) <- possible_bc
        score_df <- score_df[, -bc_col, drop = FALSE]
      }
    }
    
    row_match <- if (!is.null(rownames(score_df))) {
      sum(rownames(score_df) %in% object_cells)
    } else 0
    
    col_match <- if (!is.null(colnames(score_df))) {
      sum(colnames(score_df) %in% object_cells)
    } else 0
    
    if (verbose) {
      message("Barcode matching: rows = ", row_match, "; columns = ", col_match)
    }
    
    if (col_match > row_match) {
      if (verbose) message("Detected pathway \u00d7 spot matrix -> transpose")
      score_df <- as.data.frame(t(as.matrix(score_df)), check.names = FALSE)
      row_match <- sum(rownames(score_df) %in% object_cells)
    }
    
    if (row_match == 0) {
      stop(
        "RDS \u4e0e spatial object \u6ca1\u6709\u5339\u914d\u5230 spot barcode\uff1a", rds_file,
        "\nobject example:\n", paste(utils::head(object_cells), collapse = "\n"),
        "\nRDS example:\n", paste(utils::head(rownames(score_df)), collapse = "\n")
      )
    }
    
    numeric_cols <- vapply(score_df, is.numeric, logical(1))
    score_df <- score_df[, numeric_cols, drop = FALSE]
    
    if (ncol(score_df) == 0) {
      stop("RDS \u4e2d\u6ca1\u6709 numeric pathway scores\uff1a", rds_file)
    }
    
    keep <- vapply(score_df, function(z) any(is.finite(z)), logical(1))
    score_df <- score_df[, keep, drop = FALSE]
    
    if (ncol(score_df) == 0) {
      stop("RDS \u4e2d\u6240\u6709 pathway scores \u90fd\u662f NA\uff1a", rds_file)
    }
    
    score_df
  }
  
  select_pathways_fun <- function(pathway_names) {
    
    if (is.null(selected_pathways)) return(pathway_names)
    
    if (pathway_match == "exact") {
      return(pathway_names[pathway_names %in% selected_pathways])
    }
    
    if (pathway_match == "contains") {
      keep <- rep(FALSE, length(pathway_names))
      for (pat in selected_pathways) {
        keep <- keep | grepl(tolower(pat), tolower(pathway_names), fixed = TRUE)
      }
      return(pathway_names[keep])
    }
    
    if (pathway_match == "regex") {
      keep <- rep(FALSE, length(pathway_names))
      for (pat in selected_pathways) {
        keep <- keep | grepl(pat, pathway_names, ignore.case = TRUE)
      }
      return(pathway_names[keep])
    }
  }
  
  # ============================================================
  # containers
  # ============================================================
  
  plot_list <- list()
  summary_list <- list()
  plot_log_list <- list()
  
  total_pathways <- 0
  plots_generated <- 0
  
  # ============================================================
  # loop RDS
  # ============================================================
  
  for (rds_id in seq_along(selected_rds)) {
    
    rds_file <- selected_rds[rds_id]
    
    if (file.exists(rds_file)) {
      rds_path <- rds_file
    } else {
      rds_path <- file.path(score_dir, rds_file)
    }
    
    if (!file.exists(rds_path)) {
      warning("RDS \u6587\u4ef6\u4e0d\u5b58\u5728\uff0c\u8df3\u8fc7\uff1a", rds_path)
      next
    }
    
    if (verbose) {
      message(
        "\n============================================================",
        "\nRDS [", rds_id, "/", length(selected_rds), "]: ",
        basename(rds_path),
        "\n============================================================"
      )
    }
    
    x <- readRDS(rds_path)
    score_df <- extract_score_table(x, basename(rds_path))
    score_df <- standardize_score_table(score_df, basename(rds_path))
    
    all_pathways <- colnames(score_df)
    pathways_use <- select_pathways_fun(all_pathways)
    
    if (verbose) {
      message("Total pathways: ", length(all_pathways))
      message("Selected pathways: ", length(pathways_use))
    }
    
    if (length(pathways_use) == 0) {
      warning("\u5728 ", basename(rds_path), " \u4e2d\u6ca1\u6709\u5339\u914d\u5230 pathway\u3002")
      next
    }
    
    rds_prefix <- tools::file_path_sans_ext(basename(rds_path))
    subdir <- file.path(output_dir, make_safe_name(rds_prefix))
    dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
    
    # ==========================================================
    # loop pathways
    # ==========================================================
    
    for (pathway_id in seq_along(pathways_use)) {
      
      pathway <- pathways_use[pathway_id]
      total_pathways <- total_pathways + 1
      
      if (verbose) {
        message(
          "\nPathway [", pathway_id, "/", length(pathways_use), "]: ",
          pathway
        )
      }
      
      common_cells <- intersect(object_cells, rownames(score_df))
      
      score_vector <- rep(NA_real_, length(object_cells))
      names(score_vector) <- object_cells
      
      score_vector[common_cells] <- as.numeric(
        score_df[common_cells, pathway, drop = TRUE]
      )
      
      n_valid <- sum(is.finite(score_vector))
      
      if (n_valid < 10) {
        warning(pathway, ": \u6709\u6548 spots \u592a\u5c11\uff1a", n_valid)
        next
      }
      
      # ========================================================
      # temporary metadata
      # ========================================================
      
      temp_feature <- paste0(".TEMP_UCELL_", total_pathways)
      
      object_tmp <- object
      object_tmp@meta.data[[temp_feature]] <- score_vector[rownames(object_tmp@meta.data)]
      
      pathway_dir <- file.path(subdir, make_safe_name(pathway))
      dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)
      
      n_images_success <- 0
      
      # ========================================================
      # loop spatial images
      # ========================================================
      
      for (image_id in seq_along(image_names)) {
        
        image_name <- image_names[image_id]
        
        if (verbose) {
          message("  Image [", image_id, "/", length(image_names), "]: ", image_name)
        }
        
        p <- tryCatch(
          sce_plot_spatial_continuous(
            object = object_tmp,
            feature = temp_feature,
            image = image_name,
            image_scale = image_scale,
            smooth_k = smooth_k,
            sigma_factor = sigma_factor,
            mask_radius_factor = mask_radius_factor,
            grid_n = grid_n,
            transform = transform,
            cap_quantile = cap_quantile,
            colors = colors,
            show_he = show_he,
            he_alpha = he_alpha,
            surface_alpha = surface_alpha,
            crop = crop,
            crop_padding = crop_padding,
            xlim = xlim,
            ylim = ylim,
            contour = contour,
            contour_bins = contour_bins,
            show_spots = show_spots,
            title = paste0(pathway, "\n", image_name),
            legend_title = pathway
          ),
          error = function(e) {
            warning(
              "\nFailed:",
              "\n  pathway = ", pathway,
              "\n  image = ", image_name,
              "\n  error = ", conditionMessage(e)
            )
            NULL
          }
        )
        
        if (is.null(p)) {
          plot_log_list[[length(plot_log_list) + 1]] <- data.frame(
            rds_file = basename(rds_path),
            pathway = pathway,
            image = image_name,
            status = "failed",
            stringsAsFactors = FALSE
          )
          next
        }
        
        plots_generated <- plots_generated + 1
        n_images_success <- n_images_success + 1
        
        if (keep_plots) {
          list_name <- paste(rds_prefix, pathway, image_name, sep = "__")
          plot_list[[list_name]] <- p
        }
        
        safe_image <- make_safe_name(image_name)
        
        if (save_png) {
          ggplot2::ggsave(
            filename = file.path(pathway_dir, paste0(safe_image, "_continuous.png")),
            plot = p,
            width = width,
            height = height,
            dpi = dpi,
            bg = "white"
          )
        }
        
        if (save_pdf) {
          ggplot2::ggsave(
            filename = file.path(pathway_dir, paste0(safe_image, "_continuous.pdf")),
            plot = p,
            width = width,
            height = height,
            bg = "white"
          )
        }
        
        plot_log_list[[length(plot_log_list) + 1]] <- data.frame(
          rds_file = basename(rds_path),
          pathway = pathway,
          image = image_name,
          status = "success",
          stringsAsFactors = FALSE
        )
        
        if (!keep_plots) rm(p)
      }
      
      # ========================================================
      # pathway summary
      # ========================================================
      
      score_values <- score_vector[is.finite(score_vector)]
      
      summary_list[[length(summary_list) + 1]] <- data.frame(
        rds_file = basename(rds_path),
        pathway = pathway,
        matched_spots = length(common_cells),
        valid_spots = n_valid,
        images_requested = length(image_names),
        images_plotted = n_images_success,
        min = min(score_values, na.rm = TRUE),
        q01 = as.numeric(stats::quantile(score_values, 0.01, na.rm = TRUE)),
        median = stats::median(score_values, na.rm = TRUE),
        mean = mean(score_values, na.rm = TRUE),
        q99 = as.numeric(stats::quantile(score_values, 0.99, na.rm = TRUE)),
        max = max(score_values, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      
      rm(object_tmp)
      invisible(gc(verbose = FALSE))
    }
  }
  
  # ============================================================
  # summary
  # ============================================================
  
  if (length(summary_list) > 0) {
    plot_summary <- do.call(rbind, summary_list)
    utils::write.csv(
      plot_summary,
      file.path(output_dir, "pathway_plot_summary.csv"),
      row.names = FALSE
    )
  } else {
    plot_summary <- data.frame()
  }
  
  if (length(plot_log_list) > 0) {
    plot_log <- do.call(rbind, plot_log_list)
    utils::write.csv(
      plot_log,
      file.path(output_dir, "pathway_plot_log.csv"),
      row.names = FALSE
    )
  } else {
    plot_log <- data.frame()
  }
  
  # ============================================================
  # missing pathways
  # ============================================================
  
  if (!is.null(selected_pathways) && pathway_match == "exact") {
    
    found_pathways <- if (nrow(plot_summary) > 0) {
      unique(plot_summary$pathway)
    } else {
      character(0)
    }
    
    not_found <- setdiff(selected_pathways, found_pathways)
    
    if (length(not_found) > 0) {
      warning(
        "\u4ee5\u4e0b selected_pathways \u6ca1\u6709\u5728\u4efb\u4f55 RDS \u4e2d\u627e\u5230\uff1a\n",
        paste(not_found, collapse = "\n")
      )
    }
  }
  
  # ============================================================
  # finish
  # ============================================================
  
  if (verbose) {
    message(
      "\n============================================================",
      "\nFinished",
      "\nRDS files: ", length(selected_rds),
      "\nSpatial images: ", length(image_names),
      "\nPathways processed: ", nrow(plot_summary),
      "\nSpatial maps generated: ", plots_generated,
      "\nOutput: ", normalizePath(output_dir, mustWork = FALSE),
      "\n============================================================"
    )
  }
  
  invisible(list(
    plots = plot_list,
    summary = plot_summary,
    log = plot_log,
    images = image_names
  ))
}
