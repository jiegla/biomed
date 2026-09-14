
#' Plot a continuous spatial expression or metadata surface
#'
#' Align a gene or numeric metadata feature with spatial spots, apply Gaussian nearest-neighbor smoothing, mask distant regions and optionally overlay an H&E image.
#'
#' @param object A Seurat object. For virtual knockout analysis, also accepts a genes-by-cells count matrix or an RDS/QS file containing one of these objects.
#' @param feature One gene or numeric metadata column to visualize.
#' @param image Spatial image name; NULL selects the first image.
#' @param assay Assay name. For spatial plotting, NULL uses the default assay.
#' @param layer Expression layer or slot. Virtual knockout requires raw counts. Spatial plotting can choose an available layer when NULL.
#' @param image_scale Spatial coordinate scale: lowres or hires.
#' @param smooth_k Number of nearest spots in Gaussian smoothing.
#' @param sigma_factor Gaussian bandwidth in multiples of median spot spacing.
#' @param grid_n Number of grid positions along x; y adapts to the aspect ratio.
#' @param mask_radius_factor Maximum distance from a measured spot in multiples of spot spacing.
#' @param transform Value transformation: none, log1p, sqrt or zscore.
#' @param cap_quantile Two ordered clipping probabilities, or NULL to disable clipping.
#' @param colors A color vector or named palette: viridis, magma, inferno, plasma, cividis, turbo, blue_red, blue_yellow_red, white_red or navy_yellow_red.
#' @param n_colors Number of interpolated colors.
#' @param na_color Color for missing or masked surface values.
#' @param show_he Overlay the stored H&E image when available.
#' @param he_alpha H&E image opacity.
#' @param surface_alpha Interpolated surface opacity.
#' @param crop Crop to the observed spot extent.
#' @param crop_padding Relative padding around the cropped extent.
#' @param xlim Optional x limits, overriding automatic cropping.
#' @param ylim Optional y limits, overriding automatic cropping.
#' @param contour Overlay contour lines.
#' @param contour_bins Number of contour levels.
#' @param contour_color Contour line color.
#' @param contour_alpha Contour line opacity.
#' @param contour_size Contour line width.
#' @param show_spots Overlay original spot centers.
#' @param spot_size Spot point size.
#' @param spot_alpha Spot point opacity.
#' @param spot_color Spot point color.
#' @param reverse_y Reverse the y axis to match image coordinates.
#' @param title Optional plot title; defaults to the feature name.
#' @param legend_title Optional legend title; defaults to the feature name.
#' @param legend_position Legend position passed to ggplot2.
#' @param save_path Optional filename for saving the individual spatial plot.
#' @param width Saved plot width in inches.
#' @param height Saved plot height in inches.
#' @param dpi Resolution for raster output.
#' @param return_data Return the spot table, interpolated surface and parameters with the plot.
#' @details Metadata columns take precedence over gene names. At least ten matched, finite spots are required. Smoothing and masking are measured in multiples of median spot spacing. Transformations precede quantile clipping. The surface is a visualization of interpolated values, not additional measurements. The selected image scale must match the resolution of the image stored in the Seurat spatial object. Custom color vectors are accepted, including biomed_colors().
#' @return A ggplot, or when return_data is TRUE a list with plot, spot_data, surface_data and parameters.
#' @export
sce_plot_spatial_continuous <- function(
    object,
    feature,
    
    # -------------------------
    # Seurat / spatial settings
    # -------------------------
    image = NULL,
    assay = NULL,
    layer = NULL,
    image_scale = c("lowres", "hires"),
    
    # -------------------------
    # smoothing
    # -------------------------
    smooth_k = 12,
    sigma_factor = 1.0,
    grid_n = 300,
    
    # 防止跨越组织空洞进行插值
    mask_radius_factor = 1.5,
    
    # -------------------------
    # value processing
    # -------------------------
    transform = c("none", "log1p", "sqrt", "zscore"),
    cap_quantile = c(0.01, 0.99),
    
    # -------------------------
    # color
    # -------------------------
    colors = "magma",
    n_colors = 100,
    na_color = "transparent",
    
    # -------------------------
    # H&E
    # -------------------------
    show_he = TRUE,
    he_alpha = 0.55,
    
    # 连续信号层透明度
    surface_alpha = 0.80,
    
    # -------------------------
    # crop
    # -------------------------
    crop = TRUE,
    crop_padding = 0.03,
    
    # 手工裁剪，优先级高于 crop
    xlim = NULL,
    ylim = NULL,
    
    # -------------------------
    # contour
    # -------------------------
    contour = FALSE,
    contour_bins = 8,
    contour_color = "white",
    contour_alpha = 0.35,
    contour_size = 0.3,
    
    # -------------------------
    # 原始 spots
    # -------------------------
    show_spots = FALSE,
    spot_size = 0.15,
    spot_alpha = 0.35,
    spot_color = "black",
    
    # -------------------------
    # orientation
    # -------------------------
    reverse_y = TRUE,
    
    # -------------------------
    # labels/theme
    # -------------------------
    title = NULL,
    legend_title = NULL,
    legend_position = "right",
    
    # -------------------------
    # output
    # -------------------------
    save_path = NULL,
    width = 8,
    height = 7,
    dpi = 300,
    return_data = FALSE
) {
  
  # ============================================================
  # 0. packages
  # ============================================================
  
  pkgs <- c(
    "Seurat",
    "SeuratObject",
    "ggplot2",
    "FNN",
    "viridisLite"
  )
  
  miss <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(miss) > 0) {
    stop(
      "\u8bf7\u5148\u5b89\u88c5\u4ee5\u4e0b\u5305:\n",
      paste(miss, collapse = ", ")
    )
  }
  
  
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.")
  image_scale <- match.arg(image_scale)
  transform <- match.arg(transform)
  
  
  # ============================================================
  # 1. image
  # ============================================================
  
  if (!length(names(object@images))) stop("The Seurat object has no spatial images.")
  if (!is.character(feature) || length(feature) != 1L || is.na(feature) || !nzchar(feature))
    stop("feature must be one gene or numeric metadata column name.")
  for (value in list(sigma_factor, mask_radius_factor))
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      stop("sigma_factor and mask_radius_factor must be positive finite numbers.")
  for (value in list(smooth_k, grid_n, n_colors))
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value))
      stop("smooth_k, grid_n and n_colors must be positive integers.")
  if (grid_n < 2 || n_colors < 2) stop("grid_n and n_colors must be at least two.")
  if (!is.null(cap_quantile) && (!is.numeric(cap_quantile) || length(cap_quantile) != 2L ||
      any(!is.finite(cap_quantile)) || cap_quantile[1] < 0 || cap_quantile[2] > 1 ||
      cap_quantile[1] > cap_quantile[2])) stop("cap_quantile must be ordered probabilities.")
  if (is.null(image)) {
    image <- names(object@images)[1]
  }
  
  if (!is.character(image) || length(image) != 1L || is.na(image) || !image %in% names(object@images)) {
    stop(
      "\u627e\u4e0d\u5230 image: ", image,
      "\n\u53ef\u7528 images:\n",
      paste(names(object@images), collapse = "\n")
    )
  }
  
  img_obj <- object@images[[image]]
  
  
  # ============================================================
  # 2. coordinates
  # ============================================================
  
  coords <- tryCatch(
    {
      Seurat::GetTissueCoordinates(
        img_obj,
        scale = image_scale,
        cols = c("imagecol", "imagerow")
      )
    },
    error = function(e) {
      Seurat::GetTissueCoordinates(
        img_obj,
        scale = image_scale
      )
    }
  )
  
  coords <- as.data.frame(coords)
  
  
  # ------------------------------
  # determine barcode / cell names
  # ------------------------------
  
  if ("cell" %in% colnames(coords)) {
    spot_names <- as.character(coords$cell)
    
  } else if ("barcode" %in% colnames(coords)) {
    spot_names <- as.character(coords$barcode)
    
  } else {
    
    spot_names <- rownames(coords)
    
    if (
      is.null(spot_names) ||
      all(spot_names %in% as.character(seq_len(nrow(coords))))
    ) {
      
      tmp_cells <- tryCatch(
        SeuratObject::Cells(img_obj),
        error = function(e) NULL
      )
      
      if (
        !is.null(tmp_cells) &&
        length(tmp_cells) == nrow(coords)
      ) {
        spot_names <- tmp_cells
      }
    }
  }
  
  
  # ------------------------------
  # determine x/y columns
  # ------------------------------
  
  if (
    all(c("imagecol", "imagerow") %in% colnames(coords))
  ) {
    
    x <- coords$imagecol
    y <- coords$imagerow
    
  } else if (
    all(c("x", "y") %in% colnames(coords))
  ) {
    
    x <- coords$x
    y <- coords$y
    
  } else {
    
    numeric_cols <- which(
      vapply(coords, is.numeric, logical(1))
    )
    
    if (length(numeric_cols) < 2) {
      stop("\u65e0\u6cd5\u4ece spatial coordinates \u4e2d\u8bc6\u522b x/y \u5750\u6807\u3002")
    }
    
    x <- coords[[numeric_cols[1]]]
    y <- coords[[numeric_cols[2]]]
  }
  
  
  coord_df <- data.frame(
    spot = spot_names,
    x = as.numeric(x),
    y = as.numeric(y),
    stringsAsFactors = FALSE
  )
  
  
  # ============================================================
  # 3. extract feature
  # ============================================================
  
  meta <- object[[]]
  
  if (feature %in% colnames(meta)) {
    
    value <- meta[[feature]]
    if (!is.numeric(value)) stop("The metadata feature must be numeric.")
    names(value) <- rownames(meta)
    
    feature_source <- "metadata"
    
  } else {
    
    if (is.null(assay)) {
      assay <- SeuratObject::DefaultAssay(object)
    }
    
    if (!assay %in% names(object@assays)) {
      stop("\u627e\u4e0d\u5230 assay: ", assay)
    }
    
    if (!feature %in% rownames(object[[assay]])) {
      stop(
        "\u627e\u4e0d\u5230 feature: ", feature,
        "\n\u65e2\u4e0d\u5728 meta.data\uff0c\u4e5f\u4e0d\u5728 assay ", assay, " \u4e2d\u3002"
      )
    }
    
    
    available_layers <- tryCatch(
      SeuratObject::Layers(object[[assay]]),
      error = function(e) character(0)
    )
    
    
    # 自动选择 layer
    if (is.null(layer)) {
      
      priority <- c(
        "data",
        "scvi_normalized",
        "X",
        "counts"
      )
      
      candidate <- priority[
        priority %in% available_layers
      ]
      
      if (length(candidate) == 0) {
        stop(
          "\u65e0\u6cd5\u81ea\u52a8\u9009\u62e9 layer\u3002\n\u53ef\u7528 layers:\n",
          paste(available_layers, collapse = "\n")
        )
      }
      
      layer <- candidate[1]
      
      message(
        "Auto-selected layer = ",
        layer
      )
    }
    
    
    if (!layer %in% available_layers) {
      stop(
        "\u627e\u4e0d\u5230 layer = ", layer,
        "\n\u53ef\u7528 layers:\n",
        paste(available_layers, collapse = "\n")
      )
    }
    
    
    mat <- SeuratObject::LayerData(
      object[[assay]],
      layer = layer,
      features = feature
    )
    
    
    value <- as.numeric(mat[1, ])
    names(value) <- colnames(mat)
    
    feature_source <- paste0(
      assay,
      "/",
      layer
    )
  }
  
  
  # ============================================================
  # 4. combine coordinates + value
  # ============================================================
  
  coord_df$value <- value[
    match(coord_df$spot, names(value))
  ]
  
  coord_df <- coord_df[
    is.finite(coord_df$value) &
      is.finite(coord_df$x) &
      is.finite(coord_df$y),
    ,
    drop = FALSE
  ]
  
  
  if (nrow(coord_df) < 10) {
    stop(
      "\u6709\u6548 spatial spots \u592a\u5c11\uff0c\u4ec5\u6709 ",
      nrow(coord_df),
      " \u4e2a\u3002\u8bf7\u68c0\u67e5 barcode / coordinate \u662f\u5426\u5339\u914d\u3002"
    )
  }
  
  
  # ============================================================
  # 5. transform
  # ============================================================
  
  if (transform == "log1p") {
    
    if (any(coord_df$value < 0)) {
      stop("log1p transform \u4e0d\u80fd\u7528\u4e8e\u8d1f\u503c\u3002")
    }
    
    coord_df$value <- log1p(coord_df$value)
    
  } else if (transform == "sqrt") {
    
    if (any(coord_df$value < 0)) {
      stop("sqrt transform \u4e0d\u80fd\u7528\u4e8e\u8d1f\u503c\u3002")
    }
    
    coord_df$value <- sqrt(coord_df$value)
    
  } else if (transform == "zscore") {
    
    value_sd <- stats::sd(coord_df$value)
    if (!is.finite(value_sd) || value_sd == 0) stop("Cannot z-score a constant feature.")
    coord_df$value <- as.numeric(scale(coord_df$value))
  }
  
  
  # ============================================================
  # 6. Winsorization / outlier clipping
  # ============================================================
  
  if (!is.null(cap_quantile)) {
    
    if (
      length(cap_quantile) != 2 ||
      cap_quantile[1] < 0 ||
      cap_quantile[2] > 1
    ) {
      stop(
        "cap_quantile \u5e94\u4e3a\u7c7b\u4f3c c(0.01, 0.99)"
      )
    }
    
    q <- stats::quantile(
      coord_df$value,
      probs = cap_quantile,
      na.rm = TRUE
    )
    
    coord_df$value_plot <- pmax(
      q[1],
      pmin(q[2], coord_df$value)
    )
    
  } else {
    
    coord_df$value_plot <- coord_df$value
  }
  
  
  # ============================================================
  # 7. estimate spot spacing
  # ============================================================
  
  spot_xy <- as.matrix(
    coord_df[, c("x", "y")]
  )
  
  
  nn_spot <- FNN::get.knn(
    spot_xy,
    k = 1
  )
  
  spot_spacing <- stats::median(
    nn_spot$nn.dist[, 1],
    na.rm = TRUE
  )
  
  
  if (!is.finite(spot_spacing) ||
      spot_spacing <= 0) {
    
    stop("\u65e0\u6cd5\u4f30\u8ba1 spot spacing\u3002")
  }
  
  
  sigma <- sigma_factor * spot_spacing
  
  
  # ============================================================
  # 8. surface bounding box
  # ============================================================
  
  xr <- range(coord_df$x, na.rm = TRUE)
  yr <- range(coord_df$y, na.rm = TRUE)
  
  dx <- diff(xr)
  dy <- diff(yr)
  if (dx <= 0 || dy <= 0) stop("Spatial coordinates must span both axes.")
  
  
  surface_xr <- xr + c(
    -1,
    1
  ) * dx * 0.01
  
  surface_yr <- yr + c(
    -1,
    1
  ) * dy * 0.01
  
  
  # ============================================================
  # 9. construct grid
  # ============================================================
  
  nx <- grid_n
  
  # 保持真实长宽比
  ny <- max(
    50,
    round(
      grid_n *
        diff(surface_yr) /
        diff(surface_xr)
    )
  )
  
  
  gx <- seq(
    surface_xr[1],
    surface_xr[2],
    length.out = nx
  )
  
  gy <- seq(
    surface_yr[1],
    surface_yr[2],
    length.out = ny
  )
  
  
  grid_df <- expand.grid(
    x = gx,
    y = gy
  )
  
  
  # ============================================================
  # 10. kNN Gaussian smoothing
  # ============================================================
  
  smooth_k <- max(
    1,
    min(
      smooth_k,
      nrow(coord_df)
    )
  )
  
  
  knn <- FNN::get.knnx(
    data = spot_xy,
    query = as.matrix(
      grid_df[, c("x", "y")]
    ),
    k = smooth_k
  )
  
  
  idx <- knn$nn.index
  dst <- knn$nn.dist
  
  
  zmat <- matrix(
    coord_df$value_plot[idx],
    nrow = nrow(idx),
    ncol = ncol(idx)
  )
  
  
  # Gaussian weight
  weight <- exp(
    -(dst^2) /
      (2 * sigma^2)
  )
  
  
  z <- rowSums(
    weight * zmat,
    na.rm = TRUE
  ) /
    rowSums(
      weight,
      na.rm = TRUE
    )
  
  
  # ============================================================
  # 11. spatial mask
  # ============================================================
  
  # 如果 grid 点离任何真实 spot 太远，则不绘制。
  # 这样可以避免在组织外、孔洞、切片断裂区域造出假信号。
  
  nearest_dist <- dst[, 1]
  
  z[
    nearest_dist >
      mask_radius_factor * spot_spacing
  ] <- NA
  
  
  grid_df$z <- z
  
  
  # ============================================================
  # 12. color palette
  # ============================================================
  
  make_palette <- function(colors, n) {
    
    if (length(colors) > 1) {
      return(
        grDevices::colorRampPalette(colors)(n)
      )
    }
    
    
    pal_name <- tolower(colors)
    
    
    if (pal_name == "viridis") {
      
      return(
        viridisLite::viridis(n)
      )
      
    } else if (pal_name == "magma") {
      
      return(
        viridisLite::magma(n)
      )
      
    } else if (pal_name == "inferno") {
      
      return(
        viridisLite::inferno(n)
      )
      
    } else if (pal_name == "plasma") {
      
      return(
        viridisLite::plasma(n)
      )
      
    } else if (pal_name == "cividis") {
      
      return(
        viridisLite::cividis(n)
      )
      
    } else if (pal_name == "turbo") {
      
      return(
        viridisLite::turbo(n)
      )
      
    } else if (pal_name == "blue_red") {
      
      return(
        grDevices::colorRampPalette(
          c(
            "#2166AC",
            "#F7F7F7",
            "#B2182B"
          )
        )(n)
      )
      
    } else if (pal_name == "blue_yellow_red") {
      
      return(
        grDevices::colorRampPalette(
          c(
            "#313695",
            "#74ADD1",
            "#FFFFBF",
            "#F46D43",
            "#A50026"
          )
        )(n)
      )
      
    } else if (pal_name == "white_red") {
      
      return(
        grDevices::colorRampPalette(
          c(
            "#FFFFFF",
            "#FDBE85",
            "#EF6548",
            "#990000"
          )
        )(n)
      )
      
    } else if (pal_name == "navy_yellow_red") {
      
      return(
        grDevices::colorRampPalette(
          c(
            "#081D58",
            "#225EA8",
            "#41B6C4",
            "#FFFFCC",
            "#FD8D3C",
            "#BD0026"
          )
        )(n)
      )
      
    } else {
      
      stop(
        "\u672a\u77e5 colors = ", colors,
        "\n\u53ef\u9009\uff1aviridis / magma / inferno / plasma / ",
        "cividis / turbo / blue_red / blue_yellow_red / ",
        "white_red / navy_yellow_red",
        "\n\u6216\u8005\u76f4\u63a5\u4f20\u5165\u989c\u8272\u5411\u91cf\uff0c\u4f8b\u5982\uff1a",
        "\ncolors = c('navy','white','red')"
      )
    }
  }
  
  
  pal <- make_palette(
    colors,
    n_colors
  )
  
  
  # ============================================================
  # 13. H&E image
  # ============================================================
  
  he_raster <- NULL
  he_width <- NULL
  he_height <- NULL
  
  
  if (show_he) {
    
    he_array <- tryCatch(
      img_obj@image,
      error = function(e) NULL
    )
    
    
    if (is.null(he_array)) {
      
      warning(
        "\u65e0\u6cd5\u4ece image object \u4e2d\u63d0\u53d6 H&E image\uff1b",
        "\u5c06\u53ea\u7ed8\u5236 continuous surface\u3002"
      )
      
      show_he <- FALSE
      
    } else {
      
      he_height <- dim(he_array)[1]
      he_width <- dim(he_array)[2]
      
      
      # RGB
      if (
        length(dim(he_array)) == 3 &&
        dim(he_array)[3] >= 3
      ) {
        
        r <- he_array[, , 1]
        g <- he_array[, , 2]
        b <- he_array[, , 3]
        
        
        # 如果是 0~255
        if (max(he_array, na.rm = TRUE) > 1) {
          
          r <- r / 255
          g <- g / 255
          b <- b / 255
        }
        
        
        rgb_vector <- grDevices::rgb(
          r,
          g,
          b,
          alpha = he_alpha,
          maxColorValue = 1
        )
        
        
        he_raster <- grDevices::as.raster(
          matrix(
            rgb_vector,
            nrow = he_height,
            ncol = he_width
          )
        )
        
      } else {
        
        warning(
          "H&E image \u4e0d\u662f\u6807\u51c6 RGB array\uff1b",
          "\u8df3\u8fc7 H&E \u80cc\u666f\u3002"
        )
        
        show_he <- FALSE
      }
    }
  }
  
  
  # ============================================================
  # 14. plot crop region
  # ============================================================
  
  if (!is.null(xlim)) {
    
    plot_xr <- xlim
    
  } else if (crop) {
    
    plot_xr <- xr + c(
      -1,
      1
    ) * dx * crop_padding
    
  } else if (
    show_he &&
    !is.null(he_width)
  ) {
    
    plot_xr <- c(
      0,
      he_width
    )
    
  } else {
    
    plot_xr <- surface_xr
  }
  
  
  if (!is.null(ylim)) {
    
    plot_yr <- ylim
    
  } else if (crop) {
    
    plot_yr <- yr + c(
      -1,
      1
    ) * dy * crop_padding
    
  } else if (
    show_he &&
    !is.null(he_height)
  ) {
    
    plot_yr <- c(
      0,
      he_height
    )
    
  } else {
    
    plot_yr <- surface_yr
  }
  
  
  # ============================================================
  # 15. ggplot
  # ============================================================
  
  p <- ggplot2::ggplot()
  
  
  # ----------------------------
  # H&E
  # ----------------------------
  
  if (show_he) {
    
    p <- p +
      ggplot2::annotation_raster(
        raster = he_raster,
        xmin = 0,
        xmax = he_width,
        ymin = 0,
        ymax = he_height
      )
  }
  
  
  # ----------------------------
  # continuous surface
  # ----------------------------
  
  p <- p +
    ggplot2::geom_raster(
      data = grid_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        fill = .data$z
      ),
      interpolate = TRUE,
      alpha = surface_alpha
    )
  
  
  # ----------------------------
  # contour
  # ----------------------------
  
  if (contour) {
    
    p <- p +
      ggplot2::geom_contour(
        data = grid_df,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          z = .data$z
        ),
        bins = contour_bins,
        color = contour_color,
        alpha = contour_alpha,
        linewidth = contour_size
      )
  }
  
  
  # ----------------------------
  # original spot centers
  # ----------------------------
  
  if (show_spots) {
    
    p <- p +
      ggplot2::geom_point(
        data = coord_df,
        ggplot2::aes(
          x = .data$x,
          y = .data$y
        ),
        inherit.aes = FALSE,
        size = spot_size,
        alpha = spot_alpha,
        color = spot_color
      )
  }
  
  
  # ----------------------------
  # color
  # ----------------------------
  
  p <- p +
    ggplot2::scale_fill_gradientn(
      colours = pal,
      na.value = na_color
    )
  
  
  # ----------------------------
  # orientation
  # ----------------------------
  
  if (reverse_y) {
    
    p <- p +
      ggplot2::scale_y_reverse()
  }
  
  
  # ----------------------------
  # title
  # ----------------------------
  
  if (is.null(title)) {
    title <- feature
  }
  
  if (is.null(legend_title)) {
    legend_title <- feature
  }
  
  
  p <- p +
    ggplot2::coord_fixed(
      xlim = plot_xr,
      ylim = plot_yr,
      expand = FALSE
    ) +
    ggplot2::labs(
      title = title,
      fill = legend_title
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = legend_position,
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      )
    )
  
  
  # ============================================================
  # 16. save
  # ============================================================
  
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }
  
  
  # ============================================================
  # 17. message
  # ============================================================
  
  message(
    "\nFeature: ", feature,
    "\nSource: ", feature_source,
    "\nSpots: ", nrow(coord_df),
    "\nMedian spot spacing: ",
    round(spot_spacing, 3),
    "\nSmooth k: ", smooth_k,
    "\nSigma: ",
    round(sigma, 3),
    " (",
    sigma_factor,
    " \u00d7 spot spacing)",
    "\nMask radius: ",
    mask_radius_factor,
    " \u00d7 spot spacing"
  )
  
  
  # ============================================================
  # 18. return
  # ============================================================
  
  if (return_data) {
    
    return(
      list(
        plot = p,
        spot_data = coord_df,
        surface_data = grid_df,
        parameters = list(
          feature = feature,
          image = image,
          assay = assay,
          layer = layer,
          smooth_k = smooth_k,
          sigma_factor = sigma_factor,
          sigma = sigma,
          spot_spacing = spot_spacing,
          mask_radius_factor = mask_radius_factor,
          grid_n = grid_n
        )
      )
    )
    
  } else {
    
    return(p)
  }
}



