#!/usr/bin/env Rscript
# GSE1456 — 03 Survival analysis (KM + Cox)

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
if (!check_packages(c("ggplot2"), c("survival"))) {
  stop("Eksik paketleri kurun: Rscript 00_install_packages.R")
}

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

fig_dir <- file.path(PROJECT_ROOT, "results", "figures", "survival")
fig_main <- file.path(PROJECT_ROOT, "results", "figures")
tab_dir <- file.path(PROJECT_ROOT, "results", "tables")
ensure_dir(fig_dir)
ensure_dir(fig_main)

dat <- load_processed_data(PROJECT_ROOT)
expr <- dat$expression
meta <- dat$metadata

dge_file <- file.path(tab_dir, "DGE_results.csv")
if (!file.exists(dge_file)) stop("Once 02_DGE.R calistirin.")
dge <- read.csv(dge_file, stringsAsFactors = FALSE)

sig <- dge$adj.P.Val < 0.05 & abs(dge$logFC) > 1
up <- dge$gene[dge$logFC > 1 & sig]
down <- dge$gene[dge$logFC < -1 & sig]
up <- up[order(dge$logFC[match(up, dge$gene)], decreasing = TRUE)]
down <- down[order(dge$logFC[match(down, dge$gene)])]

genes_km <- unique(c(head(up, 5), head(down, 5)))
cat("Survival analizi genleri:", paste(genes_km, collapse = ", "), "\n")

if (!all(c("rfs_time", "rfs_event") %in% names(meta))) {
  stop("metadata'da rfs_time ve rfs_event bulunamadi.")
}

surv_df <- meta
surv_df$time <- as.numeric(surv_df$rfs_time)
surv_df$event <- as.integer(surv_df$rfs_event)
valid <- !is.na(surv_df$time) & !is.na(surv_df$event) & surv_df$time > 0
surv_df <- surv_df[valid, , drop = FALSE]
expr <- expr[, surv_df$sample_id, drop = FALSE]

plot_km <- function(gene) {
  gexpr <- as.numeric(expr[gene, surv_df$sample_id])
  surv_df$group <- factor(
    ifelse(gexpr >= median(gexpr, na.rm = TRUE), "high", "low"),
    levels = c("low", "high")
  )

  fit <- survfit(Surv(time, event) ~ group, data = surv_df)
  lr <- survdiff(Surv(time, event) ~ group, data = surv_df)
  pval <- 1 - pchisq(lr$chisq, length(lr$n) - 1)

  png(file.path(fig_dir, paste0("03_KM_", gene, ".png")), width = 9, height = 7, units = "in", res = 300)
  plot(
    fit,
    col = c("#4DBBD5", "#E64B35"),
    lwd = 2,
    xlab = "Time (years)",
    ylab = "Relapse-free survival probability",
    main = paste0(gene, " — relapse-free survival\nlog-rank p = ", signif(pval, 3))
  )
  legend(
    "topright",
    legend = c("Low expression", "High expression"),
    col = c("#4DBBD5", "#E64B35"),
    lwd = 2,
    bty = "n"
  )
  dev.off()
  cat("KM kaydedildi:", gene, " (p =", signif(pval, 3), ")\n")
}

for (gene in genes_km) {
  if (gene %in% rownames(expr)) plot_km(gene)
}

# --- Univariant Cox: tum genler ---
ensure_dir(tab_dir)
uni_rows <- list()
for (gene in genes_km) {
  if (!gene %in% rownames(expr)) next
  gexpr_tmp <- as.numeric(expr[gene, surv_df$sample_id])
  surv_df$gene_high_tmp <- as.integer(gexpr_tmp >= median(gexpr_tmp, na.rm = TRUE))
  cf <- coxph(Surv(time, event) ~ gene_high_tmp, data = surv_df)
  s  <- summary(cf)
  ci <- confint(cf)
  uni_rows[[gene]] <- data.frame(
    gene  = gene,
    hr    = exp(s$coefficients[, "coef"]),
    lower = exp(ci[, 1]),
    upper = exp(ci[, 2]),
    p     = s$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  cat("Univariant Cox —", gene, ": HR =", round(exp(s$coefficients[, "coef"]), 3),
      " p =", round(s$coefficients[, "Pr(>|z|)"], 4), "\n")
}
uni_df <- do.call(rbind, uni_rows)
rownames(uni_df) <- NULL
write.csv(uni_df, file.path(tab_dir, "cox_univariate_results.csv"), row.names = FALSE)
cat("cox_univariate_results.csv kaydedildi.\n")

# Birlesik univariant forest plot
uni_df$sig <- uni_df$p < 0.05
p_uni_forest <- ggplot(uni_df, aes(x = hr, y = reorder(gene, hr))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = sig), size = 3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, orientation = "y") +
  scale_x_log10() +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey60"), name = "p < 0.05") +
  theme_bw() +
  labs(
    title = "Univariant Cox — tum DE genler (relapse-free survival)",
    x = "Hazard Ratio (log scale)",
    y = ""
  )
ggsave(file.path(fig_main, "03_cox_forest_univariate.png"), p_uni_forest,
       width = 10, height = 6, dpi = 300)
cat("03_cox_forest_univariate.png kaydedildi.\n")

# --- Multivariable Cox: en dusuk p-degerli gen + klinik degiskenler ---
best_cox_gene <- uni_df$gene[which.min(uni_df$p)]
cat("\nMultivariable Cox geni (en dusuk univariant p):", best_cox_gene, "\n")

gexpr <- as.numeric(expr[best_cox_gene, surv_df$sample_id])
surv_df$gene_high <- as.integer(gexpr >= median(gexpr, na.rm = TRUE))
surv_df$er_neg <- as.integer(as.character(surv_df$er_status) == "ER-")

cox_vars <- c("gene_high", "er_neg")
if ("grade" %in% names(surv_df) && any(!is.na(surv_df$grade))) {
  surv_df$grade <- as.numeric(surv_df$grade)
  cox_vars <- c(cox_vars, "grade")
}

fml <- as.formula(paste("Surv(time, event) ~", paste(cox_vars, collapse = " + ")))
cox_fit <- coxph(fml, data = surv_df)
cat("\nCox multivariable model ozeti (", best_cox_gene, "):\n", sep = "")
print(summary(cox_fit))

# Multivariable forest plot
cox_tab <- summary(cox_fit)
ci_mv <- confint(cox_fit)
coef_df <- data.frame(
  term = rownames(cox_tab$coefficients),
  hr = exp(cox_tab$coefficients[, "coef"]),
  lower = exp(ci_mv[, 1]),
  upper = exp(ci_mv[, 2]),
  p = cox_tab$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)
coef_df$term <- gsub("_", " ", coef_df$term)

p_forest <- ggplot(coef_df, aes(x = hr, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, color = "#E64B35") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, orientation = "y") +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = paste0("Multivariable Cox — ", best_cox_gene, " + klinik degiskenler"),
    x = "Hazard ratio (log scale)",
    y = ""
  )
ggsave(file.path(fig_main, "03_cox_forest.png"), p_forest, width = 10, height = 6, dpi = 300)

cat("03_survival.R tamamlandi.\n")
