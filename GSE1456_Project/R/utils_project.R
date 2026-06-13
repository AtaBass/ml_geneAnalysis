# Shared utilities for GSE1456 project R scripts

get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(file.path(dirname(sub("--file=", "", file_arg[1])), ".."), mustWork = FALSE))
  }
  if (dir.exists("data")) return(normalizePath(getwd(), mustWork = FALSE))
  if (dir.exists("../data")) return(normalizePath("..", mustWork = FALSE))
  normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
}

check_packages <- function(pkgs, bioc = character()) {
  missing_cran <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_cran) > 0) {
    cat("Eksik CRAN paketleri:", paste(missing_cran, collapse = ", "), "\n")
    cat("Kurulum: install.packages(c(", paste0('"', missing_cran, '"', collapse = ", "), "))\n")
  }
  if (!requireNamespace("BiocManager", quietly = TRUE) && length(bioc) > 0) {
    cat("BiocManager gerekli. install.packages('BiocManager')\n")
    return(invisible(FALSE))
  }
  missing_bioc <- bioc[!bioc %in% rownames(installed.packages())]
  if (length(missing_bioc) > 0) {
    cat("Eksik Bioconductor paketleri:", paste(missing_bioc, collapse = ", "), "\n")
    cat("Kurulum: BiocManager::install(c(", paste0('"', missing_bioc, '"', collapse = ", "), "))\n")
  }
  invisible(length(c(missing_cran, missing_bioc)) == 0)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

parse_er_status <- function(x) {
  if (is.null(x) || all(is.na(x))) return(rep(NA_character_, length(x)))
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  pos <- grepl("^(ER\\+?|ER\\s*pos|positive|pos|1|yes|true)$", x, ignore.case = TRUE) |
    grepl("ER.*\\+", x, ignore.case = TRUE) |
    grepl("^\\+$", x)
  neg <- grepl("^(ER-?|ER\\s*neg|negative|neg|0|no|false)$", x, ignore.case = TRUE) |
    grepl("ER.*-", x, ignore.case = TRUE) & !grepl("\\+", x)
  out[pos] <- "ER+"
  out[neg] <- "ER-"
  unknown <- grepl("unknown|na|not available|missing|nd", x, ignore.case = TRUE)
  out[unknown] <- NA_character_
  out
}

infer_er_from_esr1 <- function(expr_mat, gene_col = rownames(expr_mat)) {
  if (!"ESR1" %in% gene_col) {
    stop("ESR1 geni ifade matrisinde bulunamadi; ER durumu cikarilamadi.")
  }
  esr1 <- expr_mat["ESR1", , drop = TRUE]
  med <- median(esr1, na.rm = TRUE)
  ifelse(esr1 >= med, "ER+", "ER-")
}

er_from_subtype <- function(subtype) {
  subtype <- trimws(as.character(subtype))
  out <- rep(NA_character_, length(subtype))
  out[subtype %in% c("Luminal A", "Luminal B")] <- "ER+"
  out[subtype %in% c("Basal", "ERBB2")] <- "ER-"
  out
}

parse_characteristic_value <- function(x, key) {
  x <- as.character(x)
  pattern <- paste0("^", key, ":\\s*")
  ifelse(grepl(pattern, x, ignore.case = TRUE),
         sub(pattern, "", x, ignore.case = TRUE),
         NA_character_)
}

extract_metadata_from_series_matrix <- function(series_file) {
  if (grepl("\\.gz$", series_file, ignore.case = TRUE)) {
    con <- gzfile(series_file, "r")
  } else {
    con <- file(series_file, "r")
  }
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  table_start <- grep("^!series_matrix_table_begin", lines, ignore.case = TRUE)[1]
  if (!is.na(table_start)) lines <- lines[seq_len(table_start - 1)]
  gsm_line <- grep("^!Sample_geo_accession", lines, value = TRUE)[1]
  if (is.na(gsm_line)) stop("Sample geo_accession satiri bulunamadi.")
  gsm_ids <- strsplit(sub("\t", "\t", gsm_line), "\t")[[1]][-1]
  gsm_ids <- gsub('"', "", gsm_ids)

  char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  char_mat <- do.call(rbind, lapply(char_lines, function(ln) {
    vals <- strsplit(ln, "\t")[[1]][-1]
    gsub('"', "", vals)
  }))

  meta <- data.frame(sample_id = gsm_ids, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(char_mat))) {
    first <- char_mat[i, 1]
    if (grepl(":", first)) {
      key <- sub(":.*", "", first)
      vals <- vapply(char_mat[i, ], parse_characteristic_value, character(1), key = key)
      col_name <- tolower(gsub("[^a-z0-9]+", "_", tolower(key)))
      col_name <- gsub("_+", "_", col_name)
      col_name <- gsub("^_|_$", "", col_name)
      if (nchar(col_name) == 0) next
      meta[[col_name]] <- vals
    }
  }
  meta
}

load_processed_data <- function(project_root) {
  proc_dir <- file.path(project_root, "data", "processed")
  expr <- readRDS(file.path(proc_dir, "expression_matrix.rds"))
  meta <- readRDS(file.path(proc_dir, "metadata.rds"))
  list(expression = expr, metadata = meta)
}

save_figure <- function(path, plot_obj, width = 8, height = 6) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, plot = plot_obj, width = width, height = height, dpi = 300)
}

save_base_figure <- function(path, plot_expr, width = 8, height = 6) {
  ensure_dir(dirname(path))
  png(path, width = width, height = height, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  plot_expr
}
