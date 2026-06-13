#!/usr/bin/env Rscript
# GSE1456 — 01 Preprocessing: CEL → RMA → gene symbols + metadata

set.seed(42)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("--file=", "", file_arg[1]), mustWork = FALSE))
} else {
  file.path(getwd(), "R")
}
source(file.path(script_dir, "utils_project.R"))

PROJECT_ROOT <- get_project_root()
cat("Proje dizini:", PROJECT_ROOT, "\n")

cran_pkgs <- c("ggplot2")
bioc_pkgs <- c("affy", "hgu133a.db", "AnnotationDbi")
if (!check_packages(cran_pkgs, bioc_pkgs)) {
  stop("Eksik paketleri kurun: Rscript 00_install_packages.R")
}

suppressPackageStartupMessages({
  library(affy)
  library(hgu133a.db)
  library(AnnotationDbi)
  library(ggplot2)
})

fig_dir <- file.path(PROJECT_ROOT, "results", "figures")
proc_dir <- file.path(PROJECT_ROOT, "data", "processed")
raw_dir <- file.path(PROJECT_ROOT, "data", "raw")
ensure_dir(fig_dir)
ensure_dir(proc_dir)

# --- 1. Extract RAW tar if needed & locate CEL files ---
tar_candidates <- c(
  file.path(raw_dir, "GSE1456_RAW.tar"),
  file.path(PROJECT_ROOT, "data", "GSE1456_RAW.tar"),
  file.path(PROJECT_ROOT, "..", "GSE1456_RAW.tar")
)
tar_file <- tar_candidates[file.exists(tar_candidates)][1]
extract_dir <- file.path(raw_dir, "GSE1456_RAW")
if (!is.na(tar_file) && !dir.exists(extract_dir)) {
  cat("GSE1456_RAW.tar cikariliyor ->", extract_dir, "\n")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  utils::untar(tar_file, exdir = extract_dir)
}

cel_dirs <- c(
  extract_dir,
  file.path(raw_dir, "GSE1456_RAW"),
  file.path(PROJECT_ROOT, "..", "GSE1456_RAW"),
  raw_dir
)
cel_dir <- cel_dirs[dir.exists(cel_dirs)][1]
if (is.na(cel_dir)) stop("CEL dosyalari bulunamadi. GSE1456_RAW.tar dosyasini data/raw/ altina cikarin.")

cel_files <- list.files(cel_dir, pattern = "\\.CEL(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
cat("Bulunan CEL dosyasi:", length(cel_files), "\n")
if (length(cel_files) == 0) stop("CEL dosyasi yok.")

# --- 2. Series matrix metadata (GEOquery + manual parse) ---
series_candidates <- c(
  file.path(PROJECT_ROOT, "data", "GSE1456-GPL96_series_matrix.txt"),
  file.path(PROJECT_ROOT, "data", "GSE1456-GPL96_series_matrix.txt.gz"),
  file.path(PROJECT_ROOT, "..", "GSE1456-GPL96_series_matrix.txt"),
  file.path(PROJECT_ROOT, "..", "GSE1456-GPL96_series_matrix.txt.gz")
)
series_file <- series_candidates[file.exists(series_candidates)][1]
if (is.na(series_file)) stop("series_matrix dosyasi bulunamadi.")

cat("Metadata okunuyor:", series_file, "\n")
meta <- extract_metadata_from_series_matrix(series_file)

# Standardize column names
if ("relapse" %in% names(meta)) meta$rfs_event <- as.integer(meta$relapse)
if ("surv_relapse" %in% names(meta)) meta$rfs_time <- as.numeric(meta$surv_relapse)
if ("elston" %in% names(meta)) {
  meta$grade <- suppressWarnings(as.numeric(meta$elston))
  meta$grade[meta$elston %in% c("NA", "unknown", "", NA) | is.na(meta$elston)] <- NA
}
if ("subtype" %in% names(meta)) meta$subtype <- meta$subtype

# Keep only U133A samples present in series matrix
gsm_in_matrix <- meta$sample_id

# Match CEL to GSM
cel_gsm <- sub("\\.CEL(\\.gz)?$", "", basename(cel_files), ignore.case = TRUE)
keep_cel <- cel_gsm %in% gsm_in_matrix
cel_files <- cel_files[keep_cel]
cat("Series matrix ile eslesen CEL:", length(cel_files), "\n")

log_step <- function(msg) {
  cat(format(Sys.time(), "%H:%M:%S"), "—", msg, "\n")
  flush.console()
}

# --- 3. Load & RMA ---
log_step(paste("CEL yukleniyor (", length(cel_files), " dosya, 5-15 dk)...", sep = ""))
affy_data <- ReadAffy(filenames = cel_files)
sampleNames(affy_data) <- sub("\\.CEL(\\.gz)?$", "", basename(sampleNames(affy_data)), ignore.case = TRUE)
rma_rds <- file.path(proc_dir, "expr_probe_rma.rds")
if (file.exists(rma_rds)) {
  log_step("Onceki RMA sonucu yukleniyor (expr_probe_rma.rds)...")
  expr_probe <- readRDS(rma_rds)
} else {
  log_step("RMA normalizasyonu (10-25 dk, en uzun adim)...")
  eset_rma <- rma(affy_data)
  expr_probe <- exprs(eset_rma)
  saveRDS(expr_probe, rma_rds)
  log_step("RMA tamamlandi ve ara kayit yapildi.")
}

# Pre-normalization QC — ozet istatistik (tum PM matrisini ggplot'a vermeyin; cok yavas)
log_step("Normalizasyon oncesi QC (boxplot ozet)...")
pm_pre <- pm(affy_data)
qc_pre <- data.frame(
  sample = colnames(pm_pre),
  median = apply(pm_pre, 2, median),
  q25 = apply(pm_pre, 2, quantile, 0.25),
  q75 = apply(pm_pre, 2, quantile, 0.75)
)
p_box_pre <- ggplot(qc_pre, aes(x = sample, y = median)) +
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.2, alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Pre-RMA PM intensity (median + IQR)", x = "Sample", y = "Intensity")
ggsave(file.path(fig_dir, "01_boxplot_pre_rma.png"), p_box_pre, width = 12, height = 5, dpi = 300)

# --- 4. Probe → SYMBOL ---
cat("Probeset anotasyonu (hgu133a.db)...\n")
probe_ids <- rownames(expr_probe)
symbols <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = "SYMBOL",
  keytype = "PROBEID"
)
symbols <- symbols[!is.na(symbols$SYMBOL) & symbols$SYMBOL != "", ]

expr_probe_ann <- expr_probe[symbols$PROBEID, , drop = FALSE]
probe_map <- symbols$PROBEID
gene_map <- symbols$SYMBOL

# Collapse: max mean expression per gene
cat("Gen sembollerine indirgeme (max mean probe)...\n")
mean_expr <- rowMeans(expr_probe_ann)
ord <- order(gene_map, -mean_expr)
expr_probe_ann <- expr_probe_ann[ord, , drop = FALSE]
gene_map <- gene_map[ord]
probe_map <- probe_map[ord]

keep_gene <- !duplicated(gene_map)
expr_gene <- expr_probe_ann[keep_gene, , drop = FALSE]
rownames(expr_gene) <- gene_map[keep_gene]
cat("Gen sayisi (anotasyonlu):", nrow(expr_gene), "\n")

# --- 5. Post-RMA QC (ornek basina ozet; milyon satirlik data.frame yok) ---
log_step("Normalizasyon sonrasi QC...")
qc_post <- data.frame(
  sample = colnames(expr_gene),
  median = apply(expr_gene, 2, median),
  q25 = apply(expr_gene, 2, quantile, 0.25),
  q75 = apply(expr_gene, 2, quantile, 0.75)
)
p_box_post <- ggplot(qc_post, aes(x = sample, y = median)) +
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.2, alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Post-RMA expression (median + IQR)", x = "Sample", y = "log2 expression")
ggsave(file.path(fig_dir, "01_boxplot_post_rma.png"), p_box_post, width = 12, height = 5, dpi = 300)

png(file.path(fig_dir, "01_density_post_rma.png"), width = 10, height = 6, units = "in", res = 300)
plot(density(expr_gene[, 1]), main = "Post-RMA density (all samples)", xlab = "log2 expression", col = 1)
if (ncol(expr_gene) > 1) {
  for (i in 2:ncol(expr_gene)) {
    lines(density(expr_gene[, i]), col = adjustcolor("steelblue", 0.15))
  }
}
dev.off()

# --- 6. ER status (once expr ile meta hizala) ---
common <- intersect(colnames(expr_gene), meta$sample_id)
expr_gene <- expr_gene[, common, drop = FALSE]
meta <- meta[match(common, meta$sample_id), , drop = FALSE]
rownames(meta) <- NULL

log_step("ER durumu (subtype + ESR1)...")
meta$er_status <- rep(NA_character_, nrow(meta))
meta$er_source <- rep(NA_character_, nrow(meta))

subtype_vec <- if ("subtype" %in% names(meta)) meta$subtype else rep(NA_character_, nrow(meta))
er_sub <- er_from_subtype(subtype_vec)
fill_sub <- !is.na(er_sub)
meta$er_status[fill_sub] <- er_sub[fill_sub]
meta$er_source[fill_sub] <- "subtype"

need_esr1 <- is.na(meta$er_status)
if (any(need_esr1) && "ESR1" %in% rownames(expr_gene)) {
  esr1_er <- infer_er_from_esr1(expr_gene)
  names(esr1_er) <- colnames(expr_gene)
  idx_esr1 <- which(need_esr1)
  meta$er_status[idx_esr1] <- unname(esr1_er[meta$sample_id[idx_esr1]])
  meta$er_source[idx_esr1] <- "ESR1"
}

cat("ER+:", sum(meta$er_status == "ER+", na.rm = TRUE),
    " ER-:", sum(meta$er_status == "ER-", na.rm = TRUE),
    " NA/unknown:", sum(is.na(meta$er_status) | meta$er_status == "unknown"), "\n")

# Filter unknown ER
keep_er <- !is.na(meta$er_status) & meta$er_status %in% c("ER+", "ER-")
cat("ER bilinmeyen ornekler cikariliyor:", sum(!keep_er), "\n")
expr_gene <- expr_gene[, keep_er, drop = FALSE]
meta <- meta[keep_er, ]
meta$er_status <- factor(meta$er_status, levels = c("ER+", "ER-"))

# --- HVG secimi (PDF: Top-3000 highly variable genes) ---
HVG_N <- 3000
log_step(paste("HVG secimi: top", HVG_N, "gen (varyansa gore)..."))
gene_var <- apply(expr_gene, 1, var)
top_hvg <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(HVG_N, length(gene_var)))]
expr_hvg <- expr_gene[top_hvg, , drop = FALSE]

hvg_df <- data.frame(
  gene = names(gene_var),
  variance = as.numeric(gene_var),
  is_hvg = names(gene_var) %in% top_hvg
)
hvg_df <- hvg_df[order(-hvg_df$variance), ]
write.csv(hvg_df, file.path(proc_dir, "hvg3000_gene_list.csv"), row.names = FALSE)

p_hvg <- ggplot(hvg_df, aes(x = rank(-variance), y = variance)) +
  geom_line() +
  geom_vline(xintercept = HVG_N, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    title = paste0("Highly variable genes (top ", HVG_N, " secildi)"),
    x = "Gen (varyansa gore siralama)",
    y = "Varyans"
  )
ggsave(file.path(fig_dir, "01_hvg3000_selection.png"), p_hvg, width = 8, height = 5, dpi = 300)
cat("HVG matris:", nrow(expr_hvg), "gen x", ncol(expr_hvg), "ornek\n")

log_step("PCA hesaplaniyor (tum genler)...")
pca <- prcomp(t(expr_gene), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  er_status = meta$er_status,
  sample_id = rownames(pca$x)
)
var_exp <- summary(pca)$importance[2, 1:2] * 100
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = er_status)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_bw() +
  labs(
    title = "PCA — ER status",
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    color = "ER status"
  )
ggsave(file.path(fig_dir, "01_pca_er_status.png"), p_pca, width = 8, height = 6, dpi = 300)

# --- 7. Save ---
log_step(paste("Kayit:", proc_dir))
saveRDS(expr_gene, file.path(proc_dir, "expression_matrix.rds"))
saveRDS(expr_hvg, file.path(proc_dir, "expression_matrix_hvg3000.rds"))
saveRDS(meta, file.path(proc_dir, "metadata.rds"))

# CSV for Python convenience
write.csv(expr_gene, file.path(proc_dir, "expression_matrix.csv"))
write.csv(expr_hvg, file.path(proc_dir, "expression_matrix_hvg3000.csv"))
write.csv(meta, file.path(proc_dir, "metadata.csv"), row.names = FALSE)

cat("01_preprocessing.R tamamlandi.\n")
cat("  Ornek:", ncol(expr_gene), " | Tum gen:", nrow(expr_gene),
    " | HVG3000:", nrow(expr_hvg), "\n")
