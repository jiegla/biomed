#' Add sample metadata to a Seurat object
#'
#' Single-cell utility adapted from the contributed workflow, with explicit package dependencies and cell alignment.
#' @param sce A Seurat object; for SCENIC, also accepts the path to an RDS object.
#' @param metadata_file Path to CSV, XLSX, or XLS metadata. The first column holds sample IDs.
#' @param cols Tidyselect expression selecting incoming metadata columns; the first ID column is excluded.
#' @param sample_col Sample identifier column in Seurat metadata.
#' @param replace_sample Replace existing metadata columns when TRUE; otherwise create .new, .new2, etc.
#' @param sheet Excel sheet name or index.
#' @param file_encoding CSV encoding, or NULL to try detected encoding, UTF-8, GB18030 and GBK.
#' @param strict_match Require exactly the same sample ID sets in Seurat and the metadata file.
#' @param allow_missing Explicitly permit unmatched Seurat samples to receive NA. Default FALSE requires all Seurat samples to match.
#' @return The updated Seurat object, with cells in their original order.
#' @export
sce_add_seurat_metadata <- function(
    sce,
    metadata_file,
    cols = -1,
    sample_col = "orig.ident",
    replace_sample = FALSE,
    sheet = 1,
    file_encoding = NULL,
    strict_match = FALSE,
    allow_missing = FALSE
) {

  if (!inherits(sce, "Seurat")) {
    stop("`sce` \u5fc5\u987b\u662f\u4e00\u4e2a Seurat \u5bf9\u8c61\u3002")
  }

  if (!is.character(metadata_file) || length(metadata_file) != 1L) {
    stop("`metadata_file` \u5fc5\u987b\u662f\u5355\u4e2a\u6587\u4ef6\u8def\u5f84\u3002")
  }

  if (!file.exists(metadata_file)) {
    stop("\u627e\u4e0d\u5230 metadata \u6587\u4ef6\uff1a", metadata_file)
  }

  if (!is.character(sample_col) || length(sample_col) != 1L) {
    stop("`sample_col` \u5fc5\u987b\u662f\u5355\u4e2a\u5217\u540d\u3002")
  }

  seurat_meta <- sce[[]]

  if (!sample_col %in% colnames(seurat_meta)) {
    stop(
      "Seurat metadata \u4e2d\u4e0d\u5b58\u5728\u6837\u672c\u5217\uff1a", sample_col,
      "\n\u5f53\u524d\u53ef\u7528\u5217\uff1a\n",
      paste(colnames(seurat_meta), collapse = ", ")
    )
  }

  if (!is.logical(replace_sample) || length(replace_sample) != 1L || is.na(replace_sample)) {
    stop("`replace_sample` \u5fc5\u987b\u4e3a TRUE \u6216 FALSE\u3002")
  }

  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop("\u8bf7\u5148\u5b89\u88c5 rlang\uff1ainstall.packages('rlang')")
  }

  if (!requireNamespace("tidyselect", quietly = TRUE)) {
    stop("\u8bf7\u5148\u5b89\u88c5 tidyselect\uff1ainstall.packages('tidyselect')")
  }


  ext <- tolower(tools::file_ext(metadata_file))

  if (ext == "csv") {

    if (is.null(file_encoding)) {

      guessed_encoding <- NULL

      if (requireNamespace("readr", quietly = TRUE)) {
        guessed_encoding <- tryCatch(
          {
            enc_guess <- readr::guess_encoding(metadata_file)
            if (nrow(enc_guess) > 0L) {
              as.character(enc_guess$encoding[[1]])
            } else {
              NULL
            }
          },
          error = function(e) NULL
        )
      }

      encoding_candidates <- unique(
        c(
          guessed_encoding,
          "UTF-8-BOM",
          "UTF-8",
          "GB18030",
          "GBK"
        )
      )

    } else {

      if (!is.character(file_encoding) ||
          length(file_encoding) != 1L ||
          is.na(file_encoding) ||
          trimws(file_encoding) == "") {
        stop(
          "`file_encoding` \u5fc5\u987b\u4e3a NULL \u6216\u5355\u4e2a\u6709\u6548\u7f16\u7801\u540d\u79f0\uff0c",
          "\u4f8b\u5982 'UTF-8'\u3001'UTF-8-BOM'\u3001'GB18030'\u3002"
        )
      }

      encoding_candidates <- file_encoding
    }

    meta_new <- NULL
    encoding_used <- NULL

    for (enc in encoding_candidates) {

      tmp <- tryCatch({
        raw_text <- readBin(metadata_file, "raw", n = file.info(metadata_file)$size)
        text <- iconv(rawToChar(raw_text), from = sub("-BOM$", "", enc), to = "UTF-8")
        if (is.na(text)) stop("Invalid encoding")
        text <- sub("^\ufeff", "", text)
        dat <- utils::read.csv(text = text, stringsAsFactors = FALSE,
                               check.names = FALSE, colClasses = "character")
        if (ncol(dat) > 1L) {
          dat[-1L] <- lapply(dat[-1L], utils::type.convert, as.is = TRUE)
        }
        dat
      }, error = function(e) NULL)

      if (is.null(tmp)) {
        next
      }

      text_to_check <- colnames(tmp)

      char_cols <- vapply(
        tmp,
        is.character,
        logical(1)
      )

      if (any(char_cols)) {
        text_to_check <- c(
          text_to_check,
          unlist(tmp[char_cols], use.names = FALSE)
        )
      }

      valid_utf8 <- all(
        is.na(text_to_check) |
          !is.na(
            iconv(
              text_to_check,
              from = "",
              to = "UTF-8",
              sub = NA
            )
          )
      )

      if (valid_utf8) {
        meta_new <- tmp
        encoding_used <- enc
        break
      }
    }

    if (is.null(meta_new)) {
      stop(
        "\u65e0\u6cd5\u6b63\u786e\u8bfb\u53d6 CSV \u7f16\u7801\u3002\n",
        "\u8bf7\u5c1d\u8bd5\u624b\u52a8\u6307\u5b9a\uff0c\u4f8b\u5982\uff1a\n",
        "  file_encoding = 'GB18030'\n",
        "\u6216\uff1a\n",
        "  file_encoding = 'UTF-8'"
      )
    }

    colnames(meta_new) <- enc2utf8(colnames(meta_new))

    char_cols <- vapply(
      meta_new,
      is.character,
      logical(1)
    )

    if (any(char_cols)) {
      meta_new[char_cols] <- lapply(
        meta_new[char_cols],
        enc2utf8
      )
    }

    message(
      "CSV \u7f16\u7801\uff1a",
      encoding_used,
      if (is.null(file_encoding)) "\uff08\u81ea\u52a8\u8bc6\u522b/\u5c1d\u8bd5\uff09" else "\uff08\u7528\u6237\u6307\u5b9a\uff09"
    )

  } else if (ext %in% c("xlsx", "xls")) {

    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(
        "\u8bfb\u53d6 Excel \u6587\u4ef6\u9700\u8981 readxl \u5305\u3002\n",
        "\u8bf7\u5148\u8fd0\u884c\uff1ainstall.packages('readxl')"
      )
    }

    meta_new <- readxl::read_excel(
      path = metadata_file,
      sheet = sheet,
      .name_repair = "minimal"
    )

    meta_new <- as.data.frame(meta_new)

  } else {

    stop(
      "\u4e0d\u652f\u6301\u6587\u4ef6\u683c\u5f0f\uff1a.", ext,
      "\n\u76ee\u524d\u652f\u6301 csv\u3001xlsx \u548c xls\u3002"
    )
  }

  if (nrow(meta_new) == 0L) {
    stop("metadata \u6587\u4ef6\u4e3a\u7a7a\u3002")
  }

  if (ncol(meta_new) < 2L) {
    stop(
      "metadata \u6587\u4ef6\u81f3\u5c11\u9700\u8981\u4e24\u5217\uff1a",
      "\u7b2c\u4e00\u5217\u4e3a\u6837\u672c\u540d\u79f0\uff0c\u540e\u7eed\u5217\u4e3a\u5f85\u6dfb\u52a0\u4fe1\u606f\u3002"
    )
  }


  if (any(is.na(colnames(meta_new))) || any(trimws(colnames(meta_new)) == "")) {
    stop("metadata \u6587\u4ef6\u4e2d\u5b58\u5728\u7a7a\u5217\u540d\uff0c\u8bf7\u5148\u8865\u5168\u5217\u540d\u3002")
  }

  if (anyDuplicated(colnames(meta_new))) {
    dup_names <- unique(
      colnames(meta_new)[duplicated(colnames(meta_new))]
    )

    stop(
      "metadata \u6587\u4ef6\u4e2d\u5b58\u5728\u91cd\u590d\u5217\u540d\uff1a",
      paste(dup_names, collapse = ", "),
      "\n\u8bf7\u5148\u4fee\u6539\u91cd\u590d\u5217\u540d\u3002"
    )
  }


  file_sample_col <- colnames(meta_new)[1]

  file_sample <- trimws(as.character(meta_new[[1]]))
  seurat_sample <- trimws(as.character(seurat_meta[[sample_col]]))

  meta_new[[1]] <- file_sample


  invalid_file_sample <- is.na(file_sample) | file_sample == ""

  if (any(invalid_file_sample)) {
    stop(
      "\u5916\u90e8 metadata \u7b2c\u4e00\u5217 `",
      file_sample_col,
      "` \u4e2d\u5b58\u5728 NA \u6216\u7a7a\u767d\u6837\u672c\u540d\u3002"
    )
  }

  invalid_seurat_sample <- is.na(seurat_sample) | seurat_sample == ""

  if (any(invalid_seurat_sample)) {
    stop(
      "Seurat metadata \u5217 `",
      sample_col,
      "` \u4e2d\u5b58\u5728 NA \u6216\u7a7a\u767d\u6837\u672c\u540d\u3002"
    )
  }


  if (anyDuplicated(file_sample)) {
    duplicated_samples <- unique(
      file_sample[duplicated(file_sample)]
    )

    stop(
      "\u5916\u90e8 metadata \u7b2c\u4e00\u5217\u5b58\u5728\u91cd\u590d\u6837\u672c\u540d\u3002\n",
      "\u6bcf\u4e2a\u6837\u672c\u53ea\u80fd\u51fa\u73b0\u4e00\u6b21\u3002\n",
      "\u91cd\u590d\u6837\u672c\uff1a",
      paste(duplicated_samples, collapse = ", ")
    )
  }


  sample_seurat_unique <- unique(seurat_sample)
  sample_file_unique <- unique(file_sample)

  common_samples <- intersect(
    sample_seurat_unique,
    sample_file_unique
  )

  missing_in_file <- setdiff(
    sample_seurat_unique,
    sample_file_unique
  )

  extra_in_file <- setdiff(
    sample_file_unique,
    sample_seurat_unique
  )

  if (!is.logical(strict_match) || length(strict_match) != 1L || is.na(strict_match) ||
      !is.logical(allow_missing) || length(allow_missing) != 1L || is.na(allow_missing)) {
    stop("`strict_match` and `allow_missing` must be TRUE or FALSE.")
  }
  if (strict_match && (length(missing_in_file) || length(extra_in_file))) {
    stop("Sample IDs must match exactly when strict_match = TRUE.")
  }
  if (!allow_missing && length(missing_in_file)) {
    stop("Metadata missing for Seurat samples: ", paste(missing_in_file, collapse = ", "),
         ". Use allow_missing = TRUE to explicitly fill with NA.")
  }
  message("")
  message("========================================")
  message("\u6837\u672c\u5339\u914d\u68c0\u67e5")
  message("========================================")
  message("Seurat \u6837\u672c\u5217\uff1a", sample_col)
  message("\u5916\u90e8 metadata \u6837\u672c\u5217\uff1a", file_sample_col)
  message("Seurat \u6837\u672c\u6570\uff1a", length(sample_seurat_unique))
  message("\u5916\u90e8 metadata \u6837\u672c\u6570\uff1a", length(sample_file_unique))
  message("\u5171\u540c\u6837\u672c\u6570\uff1a", length(common_samples))

  if (length(extra_in_file) > 0L) {
    message("")
    message(
      "\u26a0 \u5916\u90e8 metadata \u4e2d\u6709 ",
      length(extra_in_file),
      " \u4e2a\u6837\u672c\u5728 Seurat \u4e2d\u4e0d\u5b58\u5728\u3002"
    )
    message("\u8fd9\u4e9b\u6837\u672c\u5c06\u88ab\u5ffd\u7565\uff1a")
    message("  ", paste(extra_in_file, collapse = ", "))
  }

  if (length(missing_in_file) > 0L) {
    message("")
    message(
      "\u26a0 Seurat \u4e2d\u6709 ",
      length(missing_in_file),
      " \u4e2a\u6837\u672c\u5728\u5916\u90e8 metadata \u4e2d\u4e0d\u5b58\u5728\u3002"
    )
    message("\u8fd9\u4e9b\u6837\u672c\u5bf9\u5e94\u7ec6\u80de\u7684\u65b0\u589e metadata \u5c06\u8bbe\u7f6e\u4e3a NA\uff1a")
    message("  ", paste(missing_in_file, collapse = ", "))
  }

  if (length(common_samples) == 0L) {
    stop(
      "\nSeurat \u4e0e\u5916\u90e8 metadata \u6ca1\u6709\u4efb\u4f55\u5171\u540c\u6837\u672c\u3002",
      "\n\u8bf7\u68c0\u67e5\u6837\u672c\u540d\u79f0\u662f\u5426\u4e00\u81f4\u3002"
    )
  }

  message("")
  message(
    "\u2713 \u5c06\u4f7f\u7528 intersect() \u5f97\u5230\u7684 ",
    length(common_samples),
    " \u4e2a\u5171\u540c\u6837\u672c\u8fdb\u884c metadata \u5339\u914d\u3002"
  )


  meta_use <- meta_new[
    meta_new[[1]] %in% common_samples,
    ,
    drop = FALSE
  ]

  meta_use_sample <- as.character(meta_use[[1]])


  cols_quo <- rlang::enquo(cols)

  selected_index <- tidyselect::eval_select(
    expr = cols_quo,
    data = meta_new
  )

  if (length(selected_index) == 0L) {
    stop("\u6ca1\u6709\u9009\u62e9\u4efb\u4f55\u9700\u8981\u6dfb\u52a0\u7684 metadata \u5217\u3002")
  }

  selected_cols <- names(selected_index)

  if (file_sample_col %in% selected_cols) {
    message("")
    message(
      "\u26a0 \u7b2c\u4e00\u5217 `",
      file_sample_col,
      "` \u662f\u6837\u672c\u5339\u914d\u5217\uff0c\u4e0d\u4f1a\u4f5c\u4e3a\u65b0\u589e metadata \u6dfb\u52a0\u3002"
    )

    selected_cols <- setdiff(
      selected_cols,
      file_sample_col
    )
  }

  if (length(selected_cols) == 0L) {
    stop("\u9664\u6837\u672c ID \u5217\u5916\uff0c\u6ca1\u6709\u5176\u4ed6\u9700\u8981\u6dfb\u52a0\u7684 metadata\u3002")
  }

  message("")
  message("========================================")
  message("\u51c6\u5907\u6dfb\u52a0 metadata")
  message("========================================")
  message("\u9009\u62e9\u7684\u5217\uff1a", paste(selected_cols, collapse = ", "))


  match_index <- match(
    seurat_sample,
    meta_use_sample
  )

  cell_names <- colnames(sce)

  if (length(match_index) != ncol(sce)) {
    stop("\u5185\u90e8\u68c0\u67e5\u5931\u8d25\uff1a\u6837\u672c\u5339\u914d\u5411\u91cf\u957f\u5ea6\u4e0e Seurat \u7ec6\u80de\u6570\u4e0d\u4e00\u81f4\u3002")
  }


  added_cols <- character()
  replaced_cols <- character()
  renamed_cols <- character()
  result_map <- character()

  for (src_col in selected_cols) {

    target_col <- src_col
    current_meta_cols <- colnames(sce[[]])

    if (target_col %in% current_meta_cols) {

      if (replace_sample) {

        message("")
        message(
          "\u26a0 \u68c0\u6d4b\u5230\u91cd\u590d\u5217\uff1a",
          src_col,
          "\uff1breplace_sample = TRUE\uff0c\u5c06\u66ff\u6362 Seurat \u4e2d\u539f\u6709\u5217\u3002"
        )

        replaced_cols <- c(
          replaced_cols,
          src_col
        )

      } else {

        candidate <- paste0(src_col, ".new")

        if (candidate %in% current_meta_cols) {
          i <- 2L

          while (paste0(src_col, ".new", i) %in% current_meta_cols) {
            i <- i + 1L
          }

          candidate <- paste0(src_col, ".new", i)
        }

        target_col <- candidate

        message("")
        message("\u26a0 \u68c0\u6d4b\u5230\u91cd\u590d\u5217\uff1a", src_col)
        message(
          "  replace_sample = FALSE\uff1b\u81ea\u52a8\u4fee\u6539\u4e3a\uff1a",
          target_col
        )

        renamed_cols <- c(
          renamed_cols,
          paste0(src_col, " -> ", target_col)
        )
      }
    }

    values <- meta_use[[src_col]][match_index]
    names(values) <- cell_names

    sce[[target_col]] <- values

    added_cols <- c(
      added_cols,
      target_col
    )

    result_map[src_col] <- target_col

    message("")
    message("----------------------------------------")

    if (src_col == target_col) {
      message("Metadata: ", target_col)
    } else {
      message("Metadata: ", src_col, " -> ", target_col)
    }

    message("----------------------------------------")

    tab <- table(
      sce[[]][[target_col]],
      useNA = "ifany"
    )

    message(
      paste(
        utils::capture.output(print(tab)),
        collapse = "\n"
      )
    )
  }


  message("")
  message("========================================")
  message("Metadata \u6dfb\u52a0\u5b8c\u6210")
  message("========================================")
  message("")
  message(
    "\u2713 \u6210\u529f\u6dfb\u52a0/\u66f4\u65b0 ",
    length(added_cols),
    " \u4e2a metadata\uff1a"
  )
  message("  ", paste(added_cols, collapse = ", "))

  if (length(renamed_cols) > 0L) {
    message("")
    message("\u26a0 \u4ee5\u4e0b metadata \u56e0\u5217\u540d\u91cd\u590d\u800c\u81ea\u52a8\u6539\u540d\uff1a")

    for (x in renamed_cols) {
      message("  ", x)
    }
  }

  if (length(replaced_cols) > 0L) {
    message("")
    message("\u26a0 \u4ee5\u4e0b Seurat \u539f\u6709 metadata \u5df2\u88ab\u66ff\u6362\uff1a")
    message(
      "  ",
      paste(unique(replaced_cols), collapse = ", ")
    )
  }

  if (length(extra_in_file) > 0L) {
    message("")
    message(
      "\u26a0 \u5916\u90e8 metadata \u4e2d\u4ee5\u4e0b\u6837\u672c\u672a\u5728 Seurat \u4e2d\u627e\u5230\uff0c\u5df2\u5ffd\u7565\uff1a"
    )
    message("  ", paste(extra_in_file, collapse = ", "))
  }

  if (length(missing_in_file) > 0L) {
    n_missing_cells <- sum(
      seurat_sample %in% missing_in_file
    )

    message("")
    message(
      "\u26a0 Seurat \u4e2d\u4ee5\u4e0b\u6837\u672c\u6ca1\u6709\u5bf9\u5e94\u7684\u5916\u90e8 metadata\uff1a"
    )
    message("  ", paste(missing_in_file, collapse = ", "))
    message(
      "\u5bf9\u5e94 ",
      n_missing_cells,
      " \u4e2a\u7ec6\u80de\u7684\u65b0\u589e metadata \u5df2\u8bbe\u7f6e\u4e3a NA\u3002"
    )
  }

  message("")
  message("\u5217\u540d\u5bf9\u5e94\u5173\u7cfb\uff1a")

  for (x in names(result_map)) {
    message(
      "  ",
      x,
      " -> ",
      result_map[[x]]
    )
  }

  message("")
  message(
    "Sample \u5339\u914d\uff1a",
    length(common_samples),
    " / ",
    length(sample_seurat_unique),
    " \u4e2a Seurat samples \u627e\u5230\u5916\u90e8 metadata\u3002"
  )
  message("========================================")

  return(sce)
}
