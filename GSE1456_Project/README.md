# GSE1456 Biyoenformatik Projesi (BLM3810)

Meme kanseri **GSE1456** veri seti üzerinde ER+ / ER- sınıflandırması, diferansiyel ifade, survival, WGCNA, GO zenginleştirme ve makine öğrenmesi analizleri.

## Klasör yapısı

```
GSE1456_Project/
├── data/
│   ├── raw/GSE1456_RAW/     # CEL dosyaları (.CEL.gz)
│   ├── GSE1456-GPL96_series_matrix.txt
│   └── processed/           # RDS + CSV çıktıları
├── R/01_preprocessing.R … 05_GO_enrichment.R
├── Python/06_ML_classification.py
├── results/figures/ | tables/
└── report/
```

## Gereksinimler

### R (RStudio önerilir)

```r
install.packages(c("ggplot2", "ggpubr", "pheatmap", "GEOquery", "BiocManager"))
BiocManager::install(c(
  "affy", "hgu133a.db", "limma", "AnnotationDbi",
  "survival", "survminer", "WGCNA", "clusterProfiler",
  "org.Hs.eg.db", "enrichplot"
))
```

### Python

```bash
pip install -r requirements.txt
```

## Çalıştırma sırası

Proje kökünden (`GSE1456_Project/`):

```bash
cd R
Rscript 01_preprocessing.R   # ~10–30 dk (159 CEL, RMA)
Rscript 02_DGE.R
Rscript 03_survival.R
Rscript 04_WGCNA.R           # bellek yoğun; maxBlockSize=5000
Rscript 05_GO_enrichment.R

cd ../Python
python3 06_ML_classification.py
```

RStudio’da her scripti `setwd(".../GSE1456_Project/R")` ile açıp **Source** edebilirsiniz.

## ER durumu notu

GEO series matrix’te doğrudan ER alanı yoktur. `01_preprocessing.R` sırasıyla:

1. GEOquery `pData` içinde ER alanı aranır  
2. Luminal A/B → ER+, Basal/ERBB2 → ER- (subtype)  
3. Kalan örnekler için **ESR1** medyan eşiği  

Bilinmeyen ER örnekleri filtrelenir.

## Çıktılar

| Script | Tablolar | Figürler |
|--------|----------|----------|
| 01 | `expression_matrix.rds/csv`, `expression_matrix_hvg3000.*`, `hvg3000_gene_list.csv`, `metadata` | boxplot, density, PCA, **01_hvg3000_selection.png** |
| 02 | `DGE_results.csv` | volcano, heatmap, box/violin |
| 03 | — | `survival/03_KM_*.png`, Cox forest |
| 04 | hub genes, Cytoscape CSV | soft threshold, module-trait, **04_cytoscape_network.png** |
| 05 | `GO_enrichment_results.csv` | **05_GO_barplot.png**, **05_GO_dotplot.png** |
| 06 | `ML_results.csv` | `ML/` confusion + ROC |

Tüm figürler **PNG, 300 dpi**.

## Seed

- R: `set.seed(42)`  
- Python: `random_state=42`, `np.random.seed(42)`
