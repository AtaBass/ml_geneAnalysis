#!/usr/bin/env Rscript
# GSE1456 — 04 WGCNA co-expression modules

set.seed(42)
tryCatch(enableWGCNAThreads(nThreads = 2), error = function(e) NULL)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("--file=", "", file_arg[1]), mustWork = FALSE))
} else {
  file.path(getwd(), "R")
}
source(file.path(script_dir, "utils_project.R"))

PROJECT_ROOT <- get_project_root()
check_packages(c("ggplot2", "pheatmap"), c("WGCNA"))

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(pheatmap)
})

fig_dir <- file.path(PROJECT_ROOT, "results", "figures")
tab_dir <- file.path(PROJECT_ROOT, "results", "tables")
ensure_dir(fig_dir)
ensure_dir(tab_dir)

dat <- load_processed_data(PROJECT_ROOT)
expr <- dat$expression
meta <- dat$metadata

# Transpose: samples x genes
datExpr0 <- t(expr)
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]

# Top 5000 HVG
vars <- apply(datExpr0, 2, var)
top_genes <- names(sort(vars, decreasing = TRUE))[seq_len(min(5000, ncol(datExpr0)))]
datExpr <- datExpr0[, top_genes]
cat("WGCNA — HVG gen sayisi:", ncol(datExpr), " Ornek:", nrow(datExpr), "\n")

proc_dir <- file.path(PROJECT_ROOT, "data", "processed")
wgcna_cache <- file.path(proc_dir, "wgcna_blockwise.rds")

# Hub/GO adimlarini tekrar denemek icin onbellek (tam WGCNA'yi atlar)
if (file.exists(wgcna_cache) && Sys.getenv("WGCNA_FORCE", "") != "1") {
  cat("Onbellek yukleniyor (hub/GO adimi):", wgcna_cache, "\n")
  cached <- readRDS(wgcna_cache)
  datExpr <- cached$datExpr
  module_colors <- cached$module_colors
  soft_power <- cached$soft_power
  skip_network <- TRUE
} else {
  skip_network <- FALSE
}

if (!skip_network) {
# Soft threshold
powers <- c(1:30)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "unsigned")
fit <- sft$fitIndices
r2 <- fit[, "SFT.R.sq"]
powers_ok <- fit[, "Power"][r2 > 0.85]
soft_power <- if (length(powers_ok) > 0) min(powers_ok) else fit[which.max(r2), "Power"]
cat("Secilen soft power:", soft_power, "(R^2 > 0.85 min veya max R^2)\n")

png(file.path(fig_dir, "04_soft_threshold.png"), width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(fit[, "Power"], r2, xlab = "Soft Threshold", ylab = "Scale Free Topology R^2",
     main = "Scale independence", type = "b")
abline(h = 0.85, col = "red", lty = 2)
plot(fit[, "Power"], fit[, "mean.k."], xlab = "Soft Threshold", ylab = "Mean connectivity",
     main = "Mean connectivity", type = "b")
dev.off()

# Modules
net <- blockwiseModules(
  datExpr,
  power = soft_power,
  TOMType = "unsigned",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  maxBlockSize = 5000,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)
module_labels <- net$colors
module_colors <- labels2colors(module_labels)
cat("Modul sayisi:", length(unique(module_colors)), "\n")

saveRDS(
  list(
    datExpr = datExpr,
    module_colors = module_colors,
    soft_power = soft_power
  ),
  wgcna_cache
)
cat("WGCNA ara sonuc kaydedildi.\n")
} # skip_network

# Traits
traits <- data.frame(
  ER_neg = as.numeric(meta$er_status == "ER-"),
  RFS_event = as.numeric(meta$rfs_event)
)
rownames(traits) <- meta$sample_id
traits <- traits[rownames(datExpr), , drop = FALSE]

MEs0 <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

png(file.path(fig_dir, "04_module_trait_heatmap.png"), width = 8, height = 6, units = "in", res = 300)
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitP, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()

# Best module for ER (grey = atanmamis genler, haric tut)
er_cors <- moduleTraitCor[, "ER_neg", drop = TRUE]
names(er_cors) <- sub("^ME", "", names(er_cors))
non_grey <- names(er_cors) != "grey"
best_color <- names(which.max(abs(er_cors[non_grey])))
cat("ER ile en yuksek |korelasyon| modul:", best_color,
    "(r =", round(er_cors[best_color], 3), ")\n")

# Hub genes per module
gene_module_membership <- as.data.frame(cor(datExpr, MEs, use = "p"))

hub_list <- list()
for (mod in unique(module_colors)) {
  if (mod == "grey") next
  genes_mod <- colnames(datExpr)[module_colors == mod]
  if (length(genes_mod) == 0) next
  me_col <- paste0("ME", mod)
  if (!me_col %in% colnames(gene_module_membership)) next
  kme <- gene_module_membership[genes_mod, me_col, drop = TRUE]
  names(kme) <- genes_mod
  if (length(kme) == 0) next
  top5 <- names(sort(kme, decreasing = TRUE))[seq_len(min(5, length(kme)))]
  if (length(top5) == 0) next
  hub_list[[mod]] <- data.frame(
    module = mod,
    gene = top5,
    kME = as.numeric(kme[top5]),
    stringsAsFactors = FALSE
  )
}
hub_df <- do.call(rbind, hub_list)
write.csv(hub_df, file.path(tab_dir, "WGCNA_hub_genes.csv"), row.names = FALSE)

# Save module assignment for downstream GO script
gene_modules <- data.frame(
  gene = colnames(datExpr),
  module = module_colors,
  stringsAsFactors = FALSE
)
saveRDS(gene_modules, file.path(proc_dir, "wgcna_gene_modules.rds"))
saveRDS(list(best_module = best_color, soft_power = soft_power), file.path(proc_dir, "wgcna_summary.rds"))

# Cytoscape export — best module top 50 edges
best_genes <- colnames(datExpr)[module_colors == best_color]
if (length(best_genes) >= 2) {
  datExpr_mod <- datExpr[, best_genes, drop = FALSE]
  tom <- TOMsimilarityFromExpr(datExpr_mod, power = soft_power, TOMType = "unsigned")
  adj <- tom
  adj[lower.tri(adj, diag = TRUE)] <- 0
  edge_idx <- which(adj > 0, arr.ind = TRUE)
  edges <- data.frame(
    from = best_genes[edge_idx[, 1]],
    to = best_genes[edge_idx[, 2]],
    weight = adj[edge_idx],
    stringsAsFactors = FALSE
  )
  edges <- edges[order(-edges$weight), ]
  edges <- head(edges, 50)
  nodes <- data.frame(
    id = unique(c(edges$from, edges$to)),
    module = best_color,
    stringsAsFactors = FALSE
  )
  write.csv(edges, file.path(tab_dir, "WGCNA_cytoscape_edges.csv"), row.names = FALSE)
  write.csv(nodes, file.path(tab_dir, "WGCNA_cytoscape_nodes.csv"), row.names = FALSE)
  cat("Cytoscape edge/node dosyalari kaydedildi.\n")

  # Network grafigi (rapor PNG; CSV dosyalari Cytoscape'e import edilebilir)
  if (requireNamespace("igraph", quietly = TRUE)) {
    suppressPackageStartupMessages(library(igraph))
    g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
    w <- edges$weight
    E(g)$width <- pmax(0.5, (w - min(w)) / max(max(w) - min(w), 1e-6) * 4)
    V(g)$size <- pmax(6, degree(g) * 1.5)
    png(file.path(fig_dir, "04_cytoscape_network.png"), width = 10, height = 8, units = "in", res = 300)
    set.seed(42)
    plot(
      g,
      vertex.label = V(g)$name,
      vertex.label.cex = 0.55,
      vertex.color = adjustcolor("tomato", 0.75),
      edge.color = "grey40",
      main = paste0(
        "Co-expression network — module ", best_color,
        "\n(top ", nrow(edges), " TOM edges; import: WGCNA_cytoscape_*.csv"
      )
    )
    dev.off()
    cat("Network PNG kaydedildi: 04_cytoscape_network.png\n")
  } else {
    cat("igraph yok — network PNG atlandi. Kurulum: install.packages('igraph')\n")
  }
}

cat("04_WGCNA.R tamamlandi.\n")
