# Optional: run once to install dependencies

# CRAN mirror (required for non-interactive Rscript)
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran <- c("ggplot2", "ggpubr", "pheatmap", "igraph")
bioc <- c(
  "affy", "hgu133a.db", "limma", "AnnotationDbi",
  "survival", "survminer", "WGCNA", "clusterProfiler",
  "org.Hs.eg.db", "enrichplot"
)

to_install_cran <- setdiff(cran, rownames(installed.packages()))
if (length(to_install_cran) > 0) {
  cat("CRAN paketleri kuruluyor:", paste(to_install_cran, collapse = ", "), "\n")
  install.packages(to_install_cran, repos = getOption("repos"))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = getOption("repos"))
}

to_install_bioc <- setdiff(bioc, rownames(installed.packages()))
if (length(to_install_bioc) > 0) {
  cat("Bioconductor paketleri kuruluyor:", paste(to_install_bioc, collapse = ", "), "\n")
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
}

cat("Kurulum tamamlandi.\n")
