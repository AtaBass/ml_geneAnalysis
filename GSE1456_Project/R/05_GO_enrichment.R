#!/usr/bin/env Rscript
# GSE1456 — 05 GO enrichment (clusterProfiler)

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
if (!check_packages(c("ggplot2"), c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))) {
  stop("Eksik paketleri kurun: Rscript 00_install_packages.R")
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

fig_dir <- file.path(PROJECT_ROOT, "results", "figures")
tab_dir <- file.path(PROJECT_ROOT, "results", "tables")
proc_dir <- file.path(PROJECT_ROOT, "data", "processed")
ensure_dir(fig_dir)
ensure_dir(tab_dir)

wgcna_mod_file <- file.path(proc_dir, "wgcna_gene_modules.rds")
wgcna_sum_file <- file.path(proc_dir, "wgcna_summary.rds")
if (!file.exists(wgcna_mod_file)) stop("Once 04_WGCNA.R calistirin (wgcna_gene_modules.rds).")

gene_modules <- readRDS(wgcna_mod_file)
wgcna_sum <- if (file.exists(wgcna_sum_file)) readRDS(wgcna_sum_file) else list(best_module = NULL)
best_color <- wgcna_sum$best_module
if (is.null(best_color)) {
  tab <- table(gene_modules$module)
  tab <- tab[names(tab) != "grey"]
  best_color <- names(which.max(tab))
}
cat("GO enrichment modulu:", best_color, "\n")

module_genes <- gene_modules$gene[gene_modules$module == best_color]
cat("Modul gen sayisi:", length(module_genes), "\n")

dat <- load_processed_data(PROJECT_ROOT)
expr <- dat$expression

gene_entrez <- bitr(module_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe <- bitr(rownames(expr), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("Entrez eslesen modul geni:", nrow(gene_entrez), "/", length(module_genes), "\n")

run_enrichGO <- function(p_cut, q_cut, use_universe = TRUE) {
  args <- list(
    gene = gene_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = p_cut,
    qvalueCutoff = q_cut,
    minGSSize = 10,
    readable = TRUE
  )
  if (use_universe) args$universe <- universe$ENTREZID
  do.call(enrichGO, args)
}

enrich_method <- "universe_p005"
ego <- run_enrichGO(0.05, 0.2, use_universe = TRUE)
ego_df <- as.data.frame(ego)

if (nrow(ego_df) == 0) {
  cat("Universe ile terim yok; modul gen listesi uzerinde (universe olmadan) deneniyor...\n")
  enrich_method <- "module_only_p005"
  ego <- run_enrichGO(0.05, 0.2, use_universe = FALSE)
  ego_df <- as.data.frame(ego)
}

if (nrow(ego_df) == 0) {
  cat("GO enrichment sonucu bos.\n")
  write.csv(
    data.frame(note = "No significant GO BP terms for this module"),
    file.path(tab_dir, "GO_enrichment_results.csv"),
    row.names = FALSE
  )
} else {
  ego_df$enrichment_method <- enrich_method
  write.csv(ego_df, file.path(tab_dir, "GO_enrichment_results.csv"), row.names = FALSE)
  cat("Kaydedilen GO terimi:", nrow(ego_df), "\n")

  p_bar <- barplot(ego, showCategory = min(20, nrow(ego_df))) +
    ggtitle(paste0("GO Biological Process — WGCNA module ", best_color))
  ggsave(file.path(fig_dir, "05_GO_barplot.png"), p_bar, width = 10, height = 8, dpi = 300)

  p_dot <- dotplot(ego, showCategory = min(20, nrow(ego_df))) +
    ggtitle(paste0("GO BP dot plot — module ", best_color))
  ggsave(file.path(fig_dir, "05_GO_dotplot.png"), p_dot, width = 10, height = 8, dpi = 300)

  top5 <- head(ego_df[order(ego_df$p.adjust), ], 5)
  cat("\n--- En anlamli 5 biyolojik surec (rapor yorumu) ---\n")
  for (i in seq_len(nrow(top5))) {
    cat(sprintf("%d) %s (adj.p=%.2e, gene ratio=%s)\n",
                i, top5$Description[i], top5$p.adjust[i], top5$GeneRatio[i]))
  }
}

# Rapor icin kisa yorum cercevesi:
# - Modul ER ile iliskiliyse estrojen/ luminal veya proliferasyon terimleri beklenebilir.
# - Immune/ECM terimleri tumor mikrocevresi veya invazyon ile ilgili olabilir.
# - Sonuclar veri seti ve modul buyuklugune bagli olarak az sayida terim dondurmus olabilir.

cat("05_GO_enrichment.R tamamlandi.\n")
