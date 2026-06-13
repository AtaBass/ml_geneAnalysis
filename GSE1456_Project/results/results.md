# GSE1456 Projesi — Sonuç Dosyaları Rehberi

Bu belge, `results/` klasöründeki **tüm tabloları** ve **tüm görselleri** tek tek açıklar. Amaç: rapor yazarken veya sunum hazırlarken her dosyanın **ne işe yaradığını**, **nasıl okunacağını** ve **biyolojik/istatistiksel olarak ne anlama geldiğini** netleştirmektir.

**Proje bağlamı (kısa):**
- Veri seti: **GSE1456** (meme kanseri, Affymetrix HG-U133A)
- **159 tümör örneği**, sınıflandırma hedefi: **ER+** (n=92) vs **ER−** (n=67)
- DGE kontrastı: **ER− − ER+** (ER+ referans grup) → pozitif logFC = gen ER−’de daha yüksek; negatif logFC = ER+’de daha yüksek

---

## İçindekiler

1. [Tablolar (`results/tables/`)](#1-tablolar-resultstables)
2. [Genel görseller (`results/figures/`)](#2-genel-görseller-resultsfigures)
   - [01 — Ön işleme](#21--ön-işleme-01_preprocessingr)
   - [02 — Diferansiyel gen ifadesi](#22--diferansiyel-gen-ifadesi-02_dger)
   - [03 — Cox modeli](#23--cox-modeli-03_survivalr)
   - [04 — WGCNA](#24--wgna-04_wgcnar)
   - [05 — GO enrichment](#25--go-enrichment-05_go_enrichmentr)
3. [Survival görselleri (`results/figures/survival/`)](#3-survival-görselleri-resultsfiguressurvival)
4. [Makine öğrenmesi görselleri (`results/figures/ML/`)](#4-makine-öğrenmesi-görselleri-resultsfiguresml)
5. [Dosyaları raporda nasıl eşleştirirsiniz?](#5-dosyaları-raporda-nasıl-eşleştirirsiniz)

---

## 1. Tablolar (`results/tables/`)

### `DGE_results.csv`

**Üreten script:** `R/02_DGE.R`  
**Ne işe yarar:** ER+ ve ER− grupları arasındaki **tüm genlerin** limma diferansiyel ifade sonuçlarının tam listesi. Volcano plot, heatmap ve “kaç gen DE?” sorusunun kaynağıdır.

| Sütun | Anlamı |
|--------|--------|
| `gene` | Gen sembolü (HG-U133A anotasyonu sonrası) |
| `logFC` | Log2 fold change (**ER− − ER+**). Pozitif → ER−’de daha yüksek ifade; negatif → ER+’de daha yüksek |
| `AveExpr` | Ortalama ifade düzeyi (log2 ölçekte) |
| `t` | limma t istatistiği |
| `P.Value` | Ham p-değeri |
| `adj.P.Val` | BH ile çoklu test düzeltmesi sonrası adjusted p |
| `B` | limma B istatistiği (efekt büyüklüğü / güven ile ilgili) |

**Projede kullanılan filtre:** `adj.P.Val < 0.05` ve `|logFC| > 1` → **46 diferansiyel gen** (14 up / ER− yüksek, 32 down / ER+ yüksek).  
**Örnek yorum:** `ESR1` negatif ve büyük logFC ile ER+ tümörlerde beklenen şekilde daha yüksektir; `S100A8` pozitif logFC ile ER− / basal benzeri profille uyumludur.

**Rapor için:** Ana DGE tablosu veya ekte; volcano/heatmap bu dosyadan türetilir.

---

### `WGCNA_hub_genes.csv`

**Üreten script:** `R/04_WGCNA.R`  
**Ne işe yarar:** Her WGCNA modülünden seçilen **hub genler** (modülün “merkez” genleri). Setup-1 makine öğrenmesinde biyolojik özellik listesinin bir parçasıdır.

| Sütun | Anlamı |
|--------|--------|
| `module` | Modül rengi (ör. `turquoise`, `red`, `blue`) |
| `gene` | Hub gen sembolü |
| `kME` | Module Membership (kME): gen ifadesi ile modül eigengene’i arasındaki korelasyon; **yüksek kME = o modülün tipik geni** |

**Satır sayısı:** Modül başına en fazla 5 gen → toplam **85 satır** (17 renkli modül × 5, grey hariç).  
**Örnek:** `turquoise` modülünde `BUB1B`, `CDC20` → hücre döngüsü / proliferasyon ile ilişkili modül adayı.

---

### `WGCNA_cytoscape_edges.csv` ve `WGCNA_cytoscape_nodes.csv`

**Üreten script:** `R/04_WGCNA.R`  
**Ne işe yarar:** **Cytoscape** (veya benzeri ağ yazılımı) ile modül alt ağını çizmek için kenar ve düğüm listesi. ER ile en güçlü ilişkili modül (**`red`**) için en yüksek **TOM** (Topological Overlap Matrix) ağırlıklı **50 kenar** seçilmiştir.

**`WGCNA_cytoscape_edges.csv`:**

| Sütun | Anlamı |
|--------|--------|
| `from` | Kenarın başlangıç geni |
| `to` | Kenarın bitiş geni |
| `weight` | TOM benzerlik ağırlığı (yüksek = güçlü ko-ifade bağı) |

**`WGCNA_cytoscape_nodes.csv`:**

| Sütun | Anlamı |
|--------|--------|
| `id` | Düğüm = gen adı |
| `module` | Hangi modüle ait (`red`) |

**Görsel karşılığı:** `figures/04_cytoscape_network.png` (aynı kenarların R/igraph özeti).

---

### `GO_enrichment_results.csv`

**Üreten script:** `R/05_GO_enrichment.R`  
**Ne işe yarar:** Seçilen WGCNA modülü (**`red`**, 185 gen) için **Gene Ontology — Biological Process (BP)** zenginleştirme sonuçları.

| Sütun | Anlamı |
|--------|--------|
| `ID` | GO terim kimliği |
| `Description` | Biyolojik sürecin adı (okunabilir) |
| `GeneRatio` | Terimle ilişkili modül geni sayısı / test edilen gen sayısı |
| `BgRatio` | Arka plandaki oran |
| `pvalue` / `p.adjust` / `qvalue` | Anlamlılık (raporda genelde **p.adjust** kullanın) |
| `geneID` | Bu terime katkı veren genler |
| `enrichment_method` | `module_only_p005`: tüm gen universe’i sonuç vermediği için modül gen seti üzerinde çalıştırıldı |

**Sizin koşunuzda 8 anlamlı terim** (ör. fenol içeren bileşik metabolizması, miyelinasyon). Bar ve dot plot bu tablodan üretilir.

---

### `ML_results.csv`

**Üreten script:** `Python/06_ML_classification.py`  
**Ne işe yarar:** ER+ / ER− sınıflandırmasında **5 katlı stratified cross-validation** ile elde edilen performans özeti.

| Sütun | Anlamı |
|--------|--------|
| `setup` | `Setup1_Biological` (DGE + WGCNA hub genleri) veya `Setup2_Statistical` (varyans filtresi + SelectKBest 100 gen) |
| `model` | LogisticRegression, RandomForest, SVC |
| `metric` | accuracy, precision, recall, f1 |
| `mean` | 5 fold ortalaması |
| `std` | 5 fold standart sapması |

**Metrik kısa açıklama:**
- **Accuracy:** Doğru tahmin oranı
- **Precision:** “ER+ dediğimiz” örneklerin ne kadarı gerçekten ER+
- **Recall:** Gerçek ER+ örneklerin ne kadarını yakaladık
- **F1:** Precision ve recall harmonik ortalaması (dengesiz sınıflarda önemli)

**Görsel karşılıklar:** `figures/ML/` altındaki ROC ve confusion matrix dosyaları.

---

## 2. Genel görseller (`results/figures/`)

### 2.1 — Ön işleme (`01_preprocessing.R`)

#### `01_boxplot_pre_rma.png`

**Ne gösterir:** RMA **öncesi** her örneğin PM (Perfect Match) prob intensitelerinin özeti (medyan + IQR).  
**Ne işe yarar:** Normalizasyon gerekli mi? Örnekler arası sistematik fark (batch, RNA kalitesi) var mı?  
**Nasıl okunur:** Örneklerin medyanları birbirinden çok kopuksa batch/teknik etki olabilir. RMA bunu düzeltmeyi hedefler.

#### `01_boxplot_post_rma.png`

**Ne gösterir:** RMA **sonrası** log2 ifade düzeylerinin örnek bazlı özeti.  
**Ne işe yarar:** Normalizasyonun başarılı olduğunu göstermek — örnek medyanları birbirine yaklaşmalı.  
**Rapor cümlesi örneği:** “RMA sonrası örnekler arası dağılım hizalandı.”

#### `01_density_post_rma.png`

**Ne gösterir:** Tüm örneklerin ifade yoğunluk eğrileri üst üste.  
**Ne işe yarar:** Normalizasyon sonrası dağılım şekillerinin benzerleştiğini görsel kanıt.  
**Dikkat:** Çok sayıda eğri olduğu için tek tek örnek ayırt etmek zordur; genel örtüşme önemlidir.

#### `01_hvg3000_selection.png`

**Ne gösterir:** Genlerin varyansa göre sıralanmış eğrisi; **kırmızı kesik çizgi = top 3000 HVG** eşiği.  
**Ne işe yarar:** PDF’de istenen “highly variable gene” seçim adımının belgelenmesi.  
**Anlam:** En değişken 3000 gen, sonraki ağ analizleri / özellik seçimi için biyolojik sinyal taşıyan alt küme olarak işaretlenir (tam matris DGE ve ML için saklanır).

#### `01_pca_er_status.png`

**Ne gösterir:** Tüm genler (veya işlenmiş matris) üzerinden **PCA**; noktalar **ER+ / ER−** ile renkli.  
**Ne işe yarar:** Sınıflar ifade uzayında ayrışıyor mu? ML’nin neden yüksek accuracy verebileceğini destekler.  
**Nasıl okunur:** ER+ ve ER− kümeleri ayrılıyorsa moleküler fark güçlüdür. Eksenlerde % varyans (PC1, PC2) yazılıdır.

---

### 2.2 — Diferansiyel gen ifadesi (`02_DGE.R`)

#### `02_volcano.png`

**Ne gösterir:** Her gen için **x = logFC**, **y = −log10(adj.P.Val)**.  
**Renkler (tipik):** Up (ER− yüksek), Down (ER+ yüksek), gri = anlamsız veya eşik dışı.  
**Ne işe yarar:** Binlerce geni tek grafikte özetlemek; en ayırt edici genleri göstermek.  
**Rapor:** DGE bölümünün ana figürü.

#### `02_heatmap_top50_DE.png`

**Ne gösterir:** En güçlü **50 DE gen** × tüm örnekler; hücre rengi = z-score (satır bazlı).  
**Ne işe yarar:** Hangi örneklerin hangi gen grubuna benzediğini görmek; ER gruplarının gen ifadesi blokları oluşuyor mu?  
**Nasıl okunur:** Aynı renk bloklarına sahip örnek kümeleri benzer biyolojik profildir; sütun annotasyonunda ER durumu beklenir.

#### `02_boxplot_<GEN>.png` ve `02_violin_<GEN>.png` (10 gen × 2 = 20 dosya)

**Genler:** `ESR1`, `NAT1`, `SCUBE2`, `TFF1`, `CA12` (ER+’de relatif yüksek / “down” listesi) ve `S100A8`, `S100A9`, `GABRP`, `PROM1`, `RARRES1` (ER−’de relatif yüksek / “up” listesi).

| Dosya | Ne işe yarar |
|--------|----------------|
| **Boxplot** | ER+ vs ER− medyan ve çeyreklikler; alt başlıkta **Wilcoxon p-değeri** |
| **Violin** | Aynı karşılaştırmanın dağılım şekli (tek mod, çok modlu fark) |

**Ne anlama gelir:** Bu genler DGE listesinde en uçlarda olduğu için **biyolojik olarak en güçlü ayırt edicilerden** seçilmiştir; raporda “örnek DE gen” olarak anlatılır.

- **ESR1, GATA3 yolu (TFF1, NAT1):** Luminal / ER+ programı  
- **S100A8/A9, GABRP:** Basal / agresif alt tip sinyali  
- **CA12, SCUBE2:** Ek luminal veya diferansiyasyon ile ilişkili adaylar  

---

### 2.3 — Cox modeli (`03_survival.R`)

#### `03_cox_forest.png`

**Ne gösterir:** Cox regresyon modelindeki değişkenlerin **hazard ratio (HR)** ve güven aralıkları.  
**Modelde tipik değişkenler:** Seçilen bir DE genin yüksek/düşük ifadesi (`gene_high`), **ER−** göstergesi (`er_neg`), **grade (Elston)**.  
**Nasıl okunur:**
- **HR = 1** (dikey kesik çizgi): etki yok  
- **HR > 1:** relaps riski artışı (o değişken arttıkça)  
- **HR < 1:** relaps riski azalışı  

**Ne işe yarar:** Gen ifadesinin survival’a etkisini **klinik değişkenlerden bağımsız** (kısmen) değerlendirmek.  
**Not:** GSE1456’da yaş/cinsiyet metadata’da yok; grade ve ER kullanılmıştır.

---

### 2.4 — WGCNA (`04_WGCNA.R`)

#### `04_soft_threshold.png`

**Ne gösterir:** İki panel — (1) soft threshold **power** vs scale-free topology **R²**, (2) power vs ortalama bağlantı.  
**Ne işe yarar:** WGCNA’da ağ inşası için **soft power** seçiminin gerekçesi (sizde **power = 4**, R² > 0.85).  
**Rapor:** “Scale-free ağ varsayımına uygun minimum power seçildi.”

#### `04_module_trait_heatmap.png`

**Ne gösterir:** Her modül eigengene’inin **ER_neg** ve **RFS_event** ile korelasyonu; hücrede r ve parantez içinde p.  
**Ne işe yarar:** Hangi modülün ER durumu veya relaps olayı ile en çok ilişkili olduğunu göstermek (**red** modülü ER için öne çıkar).  
**Nasıl okunur:** Kırmızı = pozitif korelasyon, mavi = negatif; mutlak değer büyük + düşük p = güçlü ilişki.

#### `04_cytoscape_network.png`

**Ne gösterir:** **`red` modülünden** en güçlü 50 TOM kenarının ağ grafiği (gen = düğüm, kenar kalınlığı = benzerlik).  
**Ne işe yarar:** PDF’de istenen “modül alt ağı” görseli; hub genlerin birbirine nasıl bağlandığını gösterir.  
**İlişkili tablolar:** `WGCNA_cytoscape_edges.csv`, `WGCNA_cytoscape_nodes.csv` (Cytoscape’te düzenlenebilir sürüm).

---

### 2.5 — GO enrichment (`05_GO_enrichment.R`)

#### `05_GO_barplot.png`

**Ne gösterir:** En anlamlı GO BP terimlerinin **gen sayısı / zenginleşme** çubuk grafiği (clusterProfiler `barplot`).  
**Ne işe yarar:** Modülün hangi biyolojik süreçlere yoğunlaştığını tek bakışta özetlemek.

#### `05_GO_dotplot.png`

**Ne gösterir:** Aynı terimler; nokta büyüklüğü = gen oranı, renk = p.adjust (tipik).  
**Ne işe yarar:** Bar plot’a alternatif, raporda ikinci GO figürü; **GeneRatio** ve anlamlılığı birlikte okumak kolay.

**Birlikte yorum (sizin sonuçlarınız):** `red` modülü fenol/melanin metabolizması ve miyelinasyon/akson çevreleme terimleri ile zengin — raporda modülün biyolojik temasını anlatmak için kullanın.

---

## 3. Survival görselleri (`results/figures/survival/`)

Her dosya: **`03_KM_<GEN>.png`** — Kaplan–Meier **relapse-free survival** eğrisi.

**Ortak yöntem:**
1. DGE’den seçilen gen (top-5 up + top-5 down)
2. Gen ifadesi **medyan**e göre **high** / **low** grup
3. **Log-rank test** p-değeri başlıkta
4. Y ekseni: relaps **olmama** olasılığı; X ekseni: süre (yıl)

| Dosya | Gen | Tipik biyolojik bağlam | Grafikte ne aranır? |
|--------|-----|------------------------|---------------------|
| `03_KM_ESR1.png` | ESR1 | Östrojen reseptörü | ER+ biyolojisi; high/low grupların relaps eğrileri ayrılıyor mu? |
| `03_KM_NAT1.png` | NAT1 | ER ilişkili | Luminal program; prognostik fark var mı? |
| `03_KM_TFF1.png` | TFF1 | Luminal | Benzer |
| `03_KM_SCUBE2.png` | SCUBE2 | Diferansiyasyon | Süreç-spesifik prognostik sinyal |
| `03_KM_CA12.png` | CA12 | Karbonik anhidraz | pH / metabolizma ile ilişkili aday |
| `03_KM_S100A8.png` | S100A8 | Basal / inflamasyon | ER− profiline yakın örneklerde yüksek olabilir |
| `03_KM_S100A9.png` | S100A9 | S100A8 ile ilişkili | Aynı |
| `03_KM_GABRP.png` | GABRP | Basal benzeri | Agresif alt tip işareti |
| `03_KM_PROM1.png` | PROM1 | Kök hücre / progenitor | Prognostik heterojenlik |
| `03_KM_RARRES1.png` | RARRES1 | Retinoik asit yanıtı | Farklılaşma / büyüme kontrolü |

**Nasıl yorumlanır:**
- İki eğri **belirgin ayrılıyorsa** ve **p < 0.05** ise: gen ifadesi relapse-free survival ile ilişkili olabilir (tek gen, çok değişkenli analiz değil — dikkatli dil kullanın).
- Eğriler iç içe ve p büyükse: bu kohortta tek başına güçlü prognostik ayrım yok.

**Cox forest ile ilişki:** KM tek gen; `03_cox_forest.png` gen + ER + grade birlikte modelledi.

---

## 4. Makine öğrenmesi görselleri (`results/figures/ML/`)

Dosya adı şablonu:

```text
06_{roc|cm}_{Setup}_{Model}.png
```

- **Setup1_Biological:** DGE (up+down) ∪ WGCNA hub genleri (~131 özellik)  
- **Setup2_Statistical:** VarianceThreshold + SelectKBest (k=100)  
- **Model:** LogisticRegression, RandomForest, SVC  

### ROC eğrileri (`06_roc_*.png`)

**Ne gösterir:** False Positive Rate vs True Positive Rate; eğri altı alan = **AUC** (başlıkta veya lejantta).  
**Ne işe yarar:** ER+ sınıfını ne kadar iyi ayırdığımızı **olasılık eşiği** üzerinden göstermek.  
**Nasıl okunur:**
- **AUC ≈ 1:** Mükemmel ayırım  
- **AUC ≈ 0.5:** Rastgele tahmin  
- **Setup2** genelde **Setup1**’e yakın veya biraz daha iyi (tabloda LR ~0.89 accuracy)

**Üretim:** 5-fold CV ile out-of-fold olasılık tahminleri (`cross_val_predict`).

### Confusion matrix (`06_cm_*.png`)

**Ne gösterir:** **Son CV fold** için gerçek vs tahmin sayıları (ısı haritası).  
**Eksenler:** Gerçek (True) vs Tahmin (Predicted); sınıflar **ER−** (0) ve **ER+** (1).  
**Ne işe yarar:** Hangi sınıfta daha çok hata yapıldığını göstermek (ör. ER−’yi ER+ sanma).  

**Dikkat:** Tek fold’un anlık görüntüsüdür; `ML_results.csv` tüm fold’ların ortalamasıdır — raporda **tablo = özet**, **CM = örnek görsel** diye belirtin.

**Örnek okuma:** Köşegen (ER−→ER−, ER+→ER+) yüksekse model doğru; off-diagonal hatalar klinik karışıklık riskini gösterir.

---

## 5. Dosyaları raporda nasıl eşleştirirsiniz?

| Rapor bölümü | Önerilen dosyalar |
|--------------|-------------------|
| Veri ve ön işleme | `01_boxplot_pre/post`, `01_density`, `01_hvg3000`, `01_pca` |
| DGE | `02_volcano`, `02_heatmap`, 2–4 adet `02_boxplot` + `02_violin`, `DGE_results.csv` (ek) |
| Survival | 2–3 adet `03_KM_*`, `03_cox_forest` |
| WGCNA | `04_soft_threshold`, `04_module_trait`, `04_cytoscape_network`, `WGCNA_hub_genes.csv` |
| GO | `05_GO_barplot`, `05_GO_dotplot`, `GO_enrichment_results.csv` |
| Makine öğrenmesi | `ML_results.csv`, 1 ROC + 1 CM (en iyi model), kısa setup karşılaştırması |

---

## Hızlı dosya sayımı

| Klasör | Dosya sayısı |
|--------|----------------|
| `tables/` | 6 CSV |
| `figures/` (kök) | 32 PNG |
| `figures/survival/` | 10 PNG |
| `figures/ML/` | 12 PNG |
| **Toplam** | **6 tablo + 54 görsel** |

---

*Bu rehber, GSE1456_Project pipeline çıktılarına (Mayıs 2026) göre yazılmıştır. Script güncellenirse dosya listesi değişebilir; güncel liste için: `find results -type f | sort`.*
