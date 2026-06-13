#!/usr/bin/env Rscript
# Cytoscape CSV -> network PNG (04 sonrasi veya bagimsiz)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("--file=", "", file_arg[1]), mustWork = FALSE))
} else {
  file.path(getwd(), "R")
}
source(file.path(script_dir, "utils_project.R"))

PROJECT_ROOT <- get_project_root()
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("install.packages('igraph')")
}

tab_dir <- file.path(PROJECT_ROOT, "results", "tables")
fig_dir <- file.path(PROJECT_ROOT, "results", "figures")
edges <- read.csv(file.path(tab_dir, "WGCNA_cytoscape_edges.csv"), stringsAsFactors = FALSE)
nodes <- read.csv(file.path(tab_dir, "WGCNA_cytoscape_nodes.csv"), stringsAsFactors = FALSE)
mod_name <- unique(nodes$module)[1]

suppressPackageStartupMessages(library(igraph))
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
w <- edges$weight
E(g)$width <- pmax(0.5, (w - min(w)) / max(max(w) - min(w), 1e-6) * 4)
V(g)$size <- pmax(6, degree(g) * 1.5)

png(file.path(fig_dir, "04_cytoscape_network.png"), width = 10, height = 8, units = "in", res = 300)
set.seed(42)
plot(g, vertex.label = V(g)$name, vertex.label.cex = 0.55,
     vertex.color = adjustcolor("tomato", 0.75), edge.color = "grey40",
     main = paste0("Co-expression network — module ", mod_name))
dev.off()
cat("04_cytoscape_network.png kaydedildi.\n")
