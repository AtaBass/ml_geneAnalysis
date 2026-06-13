#!/usr/bin/env Rscript
# GSE1456 — 02 Differential Gene Expression (limma)

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
if (!check_packages(c("ggplot2", "pheatmap"), c("limma"))) {
  stop("Eksik paketleri kurun: Rscript 00_install_packages.R")
}

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

fig_dir <- file.path(PROJECT_ROOT, "results", "figures")
tab_dir <- file.path(PROJECT_ROOT, "results", "tables")
ensure_dir(fig_dir)
ensure_dir(tab_dir)

cat("Islenmis veri yukleniyor...\n")
dat <- load_processed_data(PROJECT_ROOT)
expr <- dat$expression
meta <- dat$metadata

meta$er_status <- factor(meta$er_status, levels = c("ER+", "ER-"))
cat("Ornek sayisi:", ncol(expr), "| ER+:", sum(meta$er_status == "ER+"),
    " ER-:", sum(meta$er_status == "ER-"), "\n")

# Design: ER+ reference (R-safe column names: + / - not allowed in makeContrasts)
design <- model.matrix(~ 0 + meta$er_status)
colnames(design) <- c("ERpos", "ERneg")
cat("Design matrix — referans grup: ER+ (ERpos)\n")

contrast <- makeContrasts(ERneg - ERpos, levels = design)
fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

res <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "none")
res$gene <- rownames(res)
res <- res[order(-abs(res$logFC)), ]

sig <- res$adj.P.Val < 0.05
strong <- abs(res$logFC) > 1
res_sig <- res[sig & strong, ]
cat("DE gen (adj.P < 0.05 & |logFC| > 1):", nrow(res_sig), "\n")

up_genes <- res_sig$gene[res_sig$logFC > 1]
down_genes <- res_sig$gene[res_sig$logFC < -1]
cat("Up-regulated:", length(up_genes), " Down-regulated:", length(down_genes), "\n")

write.csv(res, file.path(tab_dir, "DGE_results.csv"), row.names = FALSE)

# Volcano
res$neg_log10_padj <- -log10(pmax(res$adj.P.Val, 1e-300))
res$regulation <- "NS"
res$regulation[sig & res$logFC > 1] <- "Up"
res$regulation[sig & res$logFC < -1] <- "Down"

p_volcano <- ggplot(res, aes(x = logFC, y = neg_log10_padj, color = regulation)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c(Up = "#E64B35", Down = "#4DBBD5", NS = "grey70")) +
  theme_bw() +
  labs(
    title = "Volcano plot — ER- vs ER+ (ref)",
    x = "log2 Fold Change",
    y = "-log10(adj. P-value)",
    color = "Regulation"
  )
ggsave(file.path(fig_dir, "02_volcano.png"), p_volcano, width = 8, height = 6, dpi = 300)

# Heatmap top 50 DE
top50 <- head(res_sig$gene, 50)
if (length(top50) >= 2) {
  mat_hm <- expr[top50, , drop = FALSE]
  mat_hm <- t(scale(t(mat_hm)))
  ann_col <- data.frame(ER = meta$er_status)
  rownames(ann_col) <- colnames(mat_hm)
  png(file.path(fig_dir, "02_heatmap_top50_DE.png"), width = 10, height = 8, units = "in", res = 300)
  pheatmap(
    mat_hm,
    annotation_col = ann_col,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Top 50 DE genes (z-score)"
  )
  dev.off()
}

# Top 5 up / down box & violin
plot_gene_comparison <- function(gene, direction) {
  df <- data.frame(
    expression = as.numeric(expr[gene, ]),
    er_status = meta$er_status
  )
  wt <- wilcox.test(expression ~ er_status, data = df)
  p_label <- paste0("Wilcoxon p = ", format(wt$p.value, digits = 3, scientific = TRUE))

  p_box <- ggplot(df, aes(x = er_status, y = expression, fill = er_status)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() +
    labs(
      title = paste0(gene, " (", direction, ")"),
      subtitle = p_label,
      y = "log2 expression",
      x = "ER status"
    ) +
    theme(legend.position = "none")

  p_violin <- ggplot(df, aes(x = er_status, y = expression, fill = er_status)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.size = 0.5) +
    theme_bw() +
    labs(title = paste0(gene, " — violin (", direction, ")"), y = "log2 expression", x = "ER status") +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, paste0("02_boxplot_", gene, ".png")), p_box, width = 6, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, paste0("02_violin_", gene, ".png")), p_violin, width = 6, height = 5, dpi = 300)
}

top5_up <- head(up_genes, 5)
top5_down <- head(down_genes, 5)
for (g in top5_up) plot_gene_comparison(g, "up")
for (g in top5_down) plot_gene_comparison(g, "down")

cat("02_DGE.R tamamlandi.\n")
