#' Compare gene expression across groups within cell types
#'
#' Summarize gene expression by sample or by cell, perform omnibus and pairwise tests, and optionally export reports and faceted plots.
#'
#' @param object A Seurat object. For virtual knockout analysis, also accepts a genes-by-cells count matrix or an RDS/QS file containing one of these objects.
#' @param genes Character vector of gene names in the selected assay.
#' @param sample_col Metadata column identifying biological samples.
#' @param celltype_col Metadata column identifying cell types. For DEG, NULL or 'all' combines all cells.
#' @param group_col Metadata column defining comparison groups.
#' @param assay Assay name. For spatial plotting, NULL uses the default assay.
#' @param layer Expression layer or slot. Virtual knockout requires raw counts. Spatial plotting can choose an available layer when NULL.
#' @param statistic_by Use sample_level summaries or individual cell_level observations.
#' @param sample_summary_method Sample expression summary: mean or median.
#' @param min_cells Minimum retained cells per sample, cell type, group and gene.
#' @param remove_zero Remove nonpositive expression before summarization and filtering.
#' @param p_adjust_method Multiple-testing method accepted by stats::p.adjust.
#' @param make_plot Whether to construct per-gene faceted plots.
#' @param output_dir Output directory; defaults and NULL behavior are shown in Usage. NULL disables exports for gene-expression statistics and selects the configuration directory for DEG.
#' @details Sample-level analysis (the default) tests one mean or median per sample and cell type. Cell-level analysis treats cells as observations and does not model their dependence within samples. Two groups use an unpaired Wilcoxon test; more groups use Kruskal-Wallis, with pairwise Wilcoxon tests. Adjusted P values are provided within gene and globally; pairwise results additionally adjust within gene and cell type. Differences are group2 minus group1. With remove_zero = TRUE, nonpositive cells are removed before summaries and min_cells filtering, so detection fractions then describe the retained positive cells. Exports include an Excel workbook when writexl is available (otherwise CSV tables), analysis-data and expression-long RDS files, and one PNG/PDF plot per gene.
#' @return A list containing parameters, genes_found, genes_not_found, expression_long, sample_summary, analysis_data, descriptive_statistics, overall_statistics, pairwise_statistics and plots.
#' @export
sce_gene_expression_statistic <- function(
    object,
    genes,
    sample_col = "sampleid",
    celltype_col = "main.celltype",
    group_col = "response_new",
    assay = "RNA",
    layer = "data",
    statistic_by = c("sample_level", "cell_level"),
    sample_summary_method = c("mean", "median"),
    min_cells = 20,
    remove_zero = FALSE,
    p_adjust_method = "BH",
    make_plot = TRUE,
    output_dir = NULL
) {
  
  #============================================================
  # 0. 参数设置
  #============================================================
  for (pkg in c("SeuratObject", "dplyr", "tidyr", "tidyselect", "purrr")) .biomed_require(pkg)
  if (!is.numeric(min_cells) || length(min_cells) != 1L ||
      !is.finite(min_cells) || min_cells < 1 || min_cells != floor(min_cells))
    stop("min_cells must be a positive integer.")
  if (!p_adjust_method %in% stats::p.adjust.methods) stop("Unknown p_adjust_method.")
  if (!is.character(genes) || !length(genes) || anyNA(genes) || any(!nzchar(genes)))
    stop("genes must be nonempty gene names.")
  # Local bindings for data-masked column symbols used by dplyr.
  cell <- sample_id <- celltype <- group <- gene <- expression <- n_groups <-
    n_cells <- mean_expression <- median_expression <- positive_cells <-
    expression_value <- p_value <- p_adjust_gene <- p_adjust_global <-
    p_adjust_within_celltype <- NULL
  statistic_by <- match.arg(statistic_by)
  sample_summary_method <- match.arg(sample_summary_method)
  
  genes <- unique(genes)
  
  if (!inherits(object, "Seurat")) {
    stop("object\u5fc5\u987b\u662fSeurat\u5bf9\u8c61\u3002")
  }
  
  if (!assay %in% names(object@assays)) {
    stop(
      "Seurat\u5bf9\u8c61\u4e2d\u4e0d\u5b58\u5728assay\uff1a", assay,
      "\n\u53ef\u7528assay\u5305\u62ec\uff1a",
      paste(names(object@assays), collapse = ", ")
    )
  }
  
  meta <- object@meta.data
  
  required_cols <- c(
    sample_col,
    celltype_col,
    group_col
  )
  
  missing_meta_cols <- setdiff(
    required_cols,
    colnames(meta)
  )
  
  if (length(missing_meta_cols) > 0) {
    stop(
      "metadata\u4e2d\u7f3a\u5c11\u4ee5\u4e0b\u5217\uff1a",
      paste(missing_meta_cols, collapse = ", ")
    )
  }
  
  #============================================================
  # 1. 检查基因
  #============================================================
  assay_genes <- rownames(object[[assay]])
  
  genes_found <- intersect(
    genes,
    assay_genes
  )
  
  genes_not_found <- setdiff(
    genes,
    genes_found
  )
  
  if (length(genes_not_found) > 0) {
    warning(
      "\u4ee5\u4e0b\u57fa\u56e0\u672a\u5728 ", assay,
      " assay\u4e2d\u627e\u5230\uff0c\u5c06\u88ab\u5ffd\u7565\uff1a",
      paste(genes_not_found, collapse = ", ")
    )
  }
  
  if (length(genes_found) == 0) {
    stop("gene\u5217\u8868\u4e2d\u7684\u57fa\u56e0\u5747\u672a\u5728\u6307\u5b9aassay\u4e2d\u627e\u5230\u3002")
  }
  
  message(
    "\u5171\u8f93\u5165 ", length(genes),
    " \u4e2a\u57fa\u56e0\uff0c\u5176\u4e2d ",
    length(genes_found),
    " \u4e2a\u57fa\u56e0\u53ef\u7528\u4e8e\u5206\u6790\u3002"
  )
  
  #============================================================
  # 2. 提取基因表达
  # 兼容Seurat v4和v5
  #============================================================
  expression_matrix <- .biomed_sce_expression(object, assay = assay, layer = layer)
  expression_wide <- as.data.frame(
    t(as.matrix(expression_matrix[genes_found, , drop = FALSE])),
    check.names = FALSE
  )

  #============================================================
  # 3. 整理metadata
  #============================================================
  reserved <- c("cell", "sample_id", "celltype", "group", "gene", "expression")
  if (any(genes_found %in% reserved))
    stop("Gene names conflict with report columns: ", paste(intersect(genes_found, reserved), collapse = ", "))
  metadata_use <- meta[
    rownames(expression_wide),
    required_cols,
    drop = FALSE
  ] |>
    tibble::rownames_to_column("cell") |>
    dplyr::transmute(
      cell = cell,
      sample_id = as.character(.data[[sample_col]]),
      celltype = as.character(.data[[celltype_col]]),
      group = .data[[group_col]]
    )
  
  expression_tbl <- expression_wide |>
    tibble::rownames_to_column("cell")
  
  expression_long <- expression_tbl |>
    dplyr::left_join(
      metadata_use,
      by = "cell"
    ) |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(genes_found),
      names_to = "gene",
      values_to = "expression"
    ) |>
    dplyr::filter(
      !is.na(sample_id),
      !is.na(celltype),
      !is.na(group),
      !is.na(expression)
    )
  
  # 根据参数决定是否仅分析阳性细胞
  if (remove_zero) {
    expression_long <- expression_long |>
      dplyr::filter(expression > 0)
  }
  
  if (nrow(expression_long) == 0) {
    stop("\u8fc7\u6ee4\u540e\u6ca1\u6709\u53ef\u7528\u4e8e\u5206\u6790\u7684\u8868\u8fbe\u6570\u636e\u3002")
  }
  
  #============================================================
  # 4. 检查一个sample是否属于多个group
  #============================================================
  sample_group_check <- expression_long |>
    dplyr::select(sample_id, group) |>
    dplyr::distinct() |>
    dplyr::count(sample_id, name = "n_groups") |>
    dplyr::filter(n_groups > 1)
  
  if (nrow(sample_group_check) > 0) {
    warning(
      nrow(sample_group_check),
      " \u4e2asample\u5bf9\u5e94\u591a\u4e2a\u5206\u7ec4\uff0c\u8bf7\u68c0\u67e5sample_col\u6216group_col\u3002",
      "\n\u5982\u679c\u662f\u914d\u5bf9\u7684\u6cbb\u7597\u524d\u540e\u6837\u672c\uff0c\u5e94\u4fdd\u8bc1\u6bcf\u4e2asample ID\u552f\u4e00\u3002"
    )
  }
  
  #============================================================
  # 5. 计算sample-level表达汇总
  #============================================================
  sample_summary <- expression_long |>
    dplyr::group_by(
      sample_id,
      celltype,
      group,
      gene
    ) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      
      mean_expression = mean(
        expression,
        na.rm = TRUE
      ),
      
      median_expression = stats::median(
        expression,
        na.rm = TRUE
      ),
      
      sd_expression = stats::sd(
        expression,
        na.rm = TRUE
      ),
      
      positive_cells = sum(
        expression > 0,
        na.rm = TRUE
      ),
      
      positive_ratio = positive_cells / n_cells,
      
      mean_positive_expression = ifelse(
        positive_cells > 0,
        mean(
          expression[expression > 0],
          na.rm = TRUE
        ),
        NA_real_
      ),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      expression_value = ifelse(
        sample_summary_method == "mean",
        mean_expression,
        median_expression
      )
    )
  
  #============================================================
  # 6. 根据statistic_by生成统计数据
  #============================================================
  if (statistic_by == "sample_level") {
    
    analysis_data <- sample_summary |>
      dplyr::filter(
        n_cells >= min_cells,
        !is.na(expression_value)
      ) |>
      dplyr::transmute(
        observation_id = sample_id,
        sample_id = sample_id,
        celltype = celltype,
        group = group,
        gene = gene,
        expression_value = expression_value,
        n_cells = n_cells
      )
    
  } else {
    
    # 首先过滤每个sample-celltype中细胞数不足的组合
    valid_sample_celltype <- expression_long |>
      dplyr::count(
        sample_id,
        celltype,
        group,
        gene,
        name = "n_cells"
      ) |>
      dplyr::filter(n_cells >= min_cells)
    
    analysis_data <- expression_long |>
      dplyr::inner_join(
        valid_sample_celltype,
        by = c(
          "sample_id",
          "celltype",
          "group",
          "gene"
        )
      ) |>
      dplyr::transmute(
        observation_id = cell,
        sample_id = sample_id,
        celltype = celltype,
        group = group,
        gene = gene,
        expression_value = expression,
        n_cells = n_cells
      )
  }
  
  if (nrow(analysis_data) == 0) {
    stop(
      "\u7ecf\u8fc7min_cells\u8fc7\u6ee4\u540e\uff0c\u6ca1\u6709\u6570\u636e\u53ef\u7528\u4e8e\u5206\u6790\u3002",
      "\n\u8bf7\u5c1d\u8bd5\u964d\u4f4emin_cells\u3002"
    )
  }
  
  #============================================================
  # 7. 描述性统计
  #============================================================
  descriptive_statistics <- analysis_data |>
    dplyr::group_by(
      gene,
      celltype,
      group
    ) |>
    dplyr::summarise(
      n_observations = dplyr::n(),
      
      n_samples = dplyr::n_distinct(
        sample_id
      ),
      
      mean_expression = mean(
        expression_value,
        na.rm = TRUE
      ),
      
      median_expression = stats::median(
        expression_value,
        na.rm = TRUE
      ),
      
      sd_expression = stats::sd(
        expression_value,
        na.rm = TRUE
      ),
      
      q1_expression = stats::quantile(
        expression_value,
        probs = 0.25,
        na.rm = TRUE
      ),
      
      q3_expression = stats::quantile(
        expression_value,
        probs = 0.75,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  #============================================================
  # 8. 总体统计检验
  # 两组：Wilcoxon
  # 多组：Kruskal-Wallis
  #============================================================
  overall_statistics <- analysis_data |>
    dplyr::group_by(
      gene,
      celltype
    ) |>
    dplyr::group_modify(
      ~ {
        
        dat <- .x |>
          dplyr::filter(
            !is.na(group),
            !is.na(expression_value)
          )
        
        group_values <- unique(
          as.character(dat$group)
        )
        
        group_values <- group_values[
          !is.na(group_values)
        ]
        
        n_groups <- length(group_values)
        
        if (n_groups < 2) {
          return(
            tibble::tibble(
              n_groups = n_groups,
              n_observations = nrow(dat),
              n_samples = dplyr::n_distinct(dat$sample_id),
              test = NA_character_,
              statistic = NA_real_,
              p_value = NA_real_
            )
          )
        }
        
        dat$group_test <- factor(
          as.character(dat$group),
          levels = group_values
        )
        
        if (n_groups == 2) {
          
          test_result <- tryCatch(
            suppressWarnings(
              stats::wilcox.test(
                expression_value ~ group_test,
                data = dat,
                exact = FALSE
              )
            ),
            error = function(e) NULL
          )
          
          tibble::tibble(
            n_groups = n_groups,
            n_observations = nrow(dat),
            n_samples = dplyr::n_distinct(dat$sample_id),
            test = "Wilcoxon rank-sum test",
            statistic = ifelse(
              is.null(test_result),
              NA_real_,
              unname(test_result$statistic)
            ),
            p_value = ifelse(
              is.null(test_result),
              NA_real_,
              test_result$p.value
            )
          )
          
        } else {
          
          test_result <- tryCatch(
            stats::kruskal.test(
              expression_value ~ group_test,
              data = dat
            ),
            error = function(e) NULL
          )
          
          tibble::tibble(
            n_groups = n_groups,
            n_observations = nrow(dat),
            n_samples = dplyr::n_distinct(dat$sample_id),
            test = "Kruskal-Wallis test",
            statistic = ifelse(
              is.null(test_result),
              NA_real_,
              unname(test_result$statistic)
            ),
            p_value = ifelse(
              is.null(test_result),
              NA_real_,
              test_result$p.value
            )
          )
        }
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(gene) |>
    dplyr::mutate(
      p_adjust_gene = stats::p.adjust(
        p_value,
        method = p_adjust_method
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      p_adjust_global = stats::p.adjust(
        p_value,
        method = p_adjust_method
      ),
      
      significance = dplyr::case_when(
        p_adjust_global < 0.0001 ~ "****",
        p_adjust_global < 0.001  ~ "***",
        p_adjust_global < 0.01   ~ "**",
        p_adjust_global < 0.05   ~ "*",
        TRUE ~ "ns"
      )
    ) |>
    dplyr::arrange(p_adjust_global)
  
  #============================================================
  # 9. 两两比较
  # 即使是多组，也进行两两Wilcoxon检验
  #============================================================
  pairwise_statistics <- analysis_data |>
    dplyr::group_by(
      gene,
      celltype
    ) |>
    dplyr::group_modify(
      ~ {
        
        dat <- .x |>
          dplyr::filter(
            !is.na(group),
            !is.na(expression_value)
          )
        
        group_values <- unique(
          as.character(dat$group)
        )
        
        group_values <- group_values[
          !is.na(group_values)
        ]
        
        if (length(group_values) < 2) {
          return(
            tibble::tibble(
              group1 = NA_character_,
              group2 = NA_character_,
              n_group1 = NA_integer_,
              n_group2 = NA_integer_,
              n_sample_group1 = NA_integer_,
              n_sample_group2 = NA_integer_,
              mean_group1 = NA_real_,
              mean_group2 = NA_real_,
              median_group1 = NA_real_,
              median_group2 = NA_real_,
              mean_difference = NA_real_,
              median_difference = NA_real_,
              statistic = NA_real_,
              p_value = NA_real_
            )
          )
        }
        
        comparison_list <- utils::combn(
          group_values,
          2,
          simplify = FALSE
        )
        
        purrr::map_dfr(
          comparison_list,
          function(comparison) {
            
            group1 <- comparison[1]
            group2 <- comparison[2]
            
            dat_group1 <- dat |>
              dplyr::filter(
                as.character(group) == group1
              )
            
            dat_group2 <- dat |>
              dplyr::filter(
                as.character(group) == group2
              )
            
            x <- dat_group1$expression_value
            y <- dat_group2$expression_value
            
            test_result <- tryCatch(
              suppressWarnings(
                stats::wilcox.test(
                  x,
                  y,
                  exact = FALSE
                )
              ),
              error = function(e) NULL
            )
            
            tibble::tibble(
              group1 = group1,
              group2 = group2,
              
              n_group1 = length(x),
              n_group2 = length(y),
              
              n_sample_group1 = dplyr::n_distinct(
                dat_group1$sample_id
              ),
              
              n_sample_group2 = dplyr::n_distinct(
                dat_group2$sample_id
              ),
              
              mean_group1 = mean(
                x,
                na.rm = TRUE
              ),
              
              mean_group2 = mean(
                y,
                na.rm = TRUE
              ),
              
              median_group1 = stats::median(
                x,
                na.rm = TRUE
              ),
              
              median_group2 = stats::median(
                y,
                na.rm = TRUE
              ),
              
              # group2 - group1
              mean_difference =
                mean(y, na.rm = TRUE) -
                mean(x, na.rm = TRUE),
              
              median_difference =
                stats::median(y, na.rm = TRUE) -
                stats::median(x, na.rm = TRUE),
              
              statistic = ifelse(
                is.null(test_result),
                NA_real_,
                unname(test_result$statistic)
              ),
              
              p_value = ifelse(
                is.null(test_result),
                NA_real_,
                test_result$p.value
              )
            )
          }
        )
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      gene,
      celltype
    ) |>
    dplyr::mutate(
      p_adjust_within_celltype = stats::p.adjust(
        p_value,
        method = p_adjust_method
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(gene) |>
    dplyr::mutate(
      p_adjust_gene = stats::p.adjust(
        p_value,
        method = p_adjust_method
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      p_adjust_global = stats::p.adjust(
        p_value,
        method = p_adjust_method
      ),
      
      significance = dplyr::case_when(
        p_adjust_global < 0.0001 ~ "****",
        p_adjust_global < 0.001  ~ "***",
        p_adjust_global < 0.01   ~ "**",
        p_adjust_global < 0.05   ~ "*",
        TRUE ~ "ns"
      )
    ) |>
    dplyr::arrange(p_adjust_global)
  
  #============================================================
  # 10. 绘图
  #============================================================
  plot_list <- list()
  
  if (make_plot) {
    
    for (current_gene in unique(analysis_data$gene)) {
      
      plot_data <- analysis_data |>
        dplyr::filter(gene == current_gene)
      
      if (statistic_by == "sample_level") {
        
        p <- ggplot2::ggplot(
          plot_data,
          ggplot2::aes(
            x = .data$group,
            y = .data$expression_value
          )
        ) +
          ggplot2::geom_boxplot(
            width = 0.65,
            outlier.shape = NA
          ) +
          ggplot2::geom_jitter(
            ggplot2::aes(shape = .data$group),
            width = 0.15,
            height = 0,
            size = 2.2
          )
        
      } else {
        
        p <- ggplot2::ggplot(
          plot_data,
          ggplot2::aes(
            x = .data$group,
            y = .data$expression_value
          )
        ) +
          ggplot2::geom_violin(
            scale = "width",
            trim = TRUE
          ) +
          ggplot2::geom_boxplot(
            width = 0.18,
            outlier.shape = NA
          )
      }
      
      p <- p +
        ggplot2::facet_wrap(
          ~ celltype,
          scales = "free_y",
          ncol = 4
        ) +
        ggplot2::labs(
          title = paste0(
            current_gene,
            " expression comparison"
          ),
          subtitle = paste0(
            "Statistical unit: ",
            statistic_by
          ),
          x = group_col,
          y = ifelse(
            statistic_by == "sample_level",
            paste0(
              sample_summary_method,
              " expression per sample"
            ),
            "Expression per cell"
          )
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          strip.text = ggplot2::element_text(
            face = "bold"
          ),
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1
          ),
          legend.position = "none"
        )
      
      plot_list[[current_gene]] <- p
    }
  }
  #============================================================
  # 11. 保存结果
  #============================================================
  if (!is.null(output_dir)) {
    
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    output_prefix <- paste0(
      "gene_expression_",
      statistic_by
    )
    
    #----------------------------------------------------------
    # 11.1 Excel只保存体积较小的汇总结果
    #----------------------------------------------------------
    if (requireNamespace("writexl", quietly = TRUE)) {
      
      writexl::write_xlsx(
        list(
          descriptive_statistics = descriptive_statistics,
          overall_statistics = overall_statistics,
          pairwise_statistics = pairwise_statistics,
          sample_summary = sample_summary
        ),
        path = file.path(
          output_dir,
          paste0(
            output_prefix,
            "_summary_results.xlsx"
          )
        )
      )
      
    } else {
      
      warning(
        "\u672a\u5b89\u88c5writexl\uff0c\u5c06\u7edf\u8ba1\u7ed3\u679c\u4fdd\u5b58\u4e3aCSV\u6587\u4ef6\u3002"
      )
      
      utils::write.csv(
        descriptive_statistics,
        file.path(
          output_dir,
          paste0(
            output_prefix,
            "_descriptive_statistics.csv"
          )
        ),
        row.names = FALSE
      )
      
      utils::write.csv(
        overall_statistics,
        file.path(
          output_dir,
          paste0(
            output_prefix,
            "_overall_statistics.csv"
          )
        ),
        row.names = FALSE
      )
      
      utils::write.csv(
        pairwise_statistics,
        file.path(
          output_dir,
          paste0(
            output_prefix,
            "_pairwise_statistics.csv"
          )
        ),
        row.names = FALSE
      )
      
      utils::write.csv(
        sample_summary,
        file.path(
          output_dir,
          paste0(
            output_prefix,
            "_sample_summary.csv"
          )
        ),
        row.names = FALSE
      )
    }
    
    #----------------------------------------------------------
    # 11.2 大型数据保存为RDS
    #----------------------------------------------------------
    saveRDS(
      analysis_data,
      file = file.path(
        output_dir,
        paste0(
          output_prefix,
          "_analysis_data.rds"
        )
      ),
      compress = TRUE
    )
    
    saveRDS(
      expression_long,
      file = file.path(
        output_dir,
        paste0(
          output_prefix,
          "_expression_long.rds"
        )
      ),
      compress = TRUE
    )
    
    #----------------------------------------------------------
    # 11.3 保存图片
    #----------------------------------------------------------
    if (make_plot && length(plot_list) > 0) {
      
      for (current_gene in names(plot_list)) {
        
        safe_gene_name <- gsub(
          "[^A-Za-z0-9_.-]",
          "_",
          current_gene
        )
        
        ggplot2::ggsave(
          filename = file.path(
            output_dir,
            paste0(
              safe_gene_name,
              "_",
              statistic_by,
              ".png"
            )
          ),
          plot = plot_list[[current_gene]],
          width = 13,
          height = 9,
          dpi = 300
        )
        
        ggplot2::ggsave(
          filename = file.path(
            output_dir,
            paste0(
              safe_gene_name,
              "_",
              statistic_by,
              ".pdf"
            )
          ),
          plot = plot_list[[current_gene]],
          width = 13,
          height = 9
        )
      }
    }
  }
  #============================================================
  # 12. 返回结果
  #============================================================
  result <- list(
    parameters = list(
      statistic_by = statistic_by,
      sample_summary_method = sample_summary_method,
      assay = assay,
      layer = layer,
      min_cells = min_cells,
      remove_zero = remove_zero
    ),
    
    genes_found = genes_found,
    genes_not_found = genes_not_found,
    
    expression_long = expression_long,
    sample_summary = sample_summary,
    analysis_data = analysis_data,
    
    descriptive_statistics = descriptive_statistics,
    overall_statistics = overall_statistics,
    pairwise_statistics = pairwise_statistics,
    
    plots = plot_list
  )
  
  return(result)
}

