============================================================
  #  PIPELINE: Microglia NSG vs BL6 Control (PSEUDOBULK)
  #  Estrategia correcta cuando BL6 = 8 columnas (muestras bulk)
  #
  #  LÓGICA:
  #  - NSG  → single-cell Seurat → extraer microglia → 
  #            agregar por muestra (pseudobulk) → DESeq2
  #  - BL6  → Excel ya es pseudobulk → DESeq2 directamente
  #  - Comparación: DESeq2 NSG_microglia vs BL6_control
  # ============================================================

# ============================================================
#  BLOQUE A — NSG SINGLE-CELL
#  (igual que ANGEL: cargar, integrar, extraer microglia)
# ============================================================

# Enter commands in R (or R studio, if installed)
install.packages('Seurat')

setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

library(Seurat)
library (dplyr)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Input files will be filtered_feature_bc_matrix.h5 (H5 files) from 10X Genomics Cell Ranger pipeline output.
ctx1.data <- Read10X_h5("C:/Users/mdmsa/Desktop/Transcriptomics/Transcriptomic for moseq/filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)

ctx11.data <- Read10X_h5("C:/Users/mdmsa/Desktop/Transcriptomics/Transcriptomic for moseq/filtered_feature_bc_matrix_11.h5", use.names = TRUE, unique.features = TRUE)

wm2.data <- Read10X_h5("C:/Users/mdmsa/Desktop/Transcriptomics/Transcriptomic for moseq/filtered_feature_bc_matrix_2.h5", use.names = TRUE, unique.features = TRUE)

wm12.data <- Read10X_h5("C:/Users/mdmsa/Desktop/Transcriptomics/Transcriptomic for moseq/filtered_feature_bc_matrix_12.h5", use.names = TRUE, unique.features = TRUE)

#Create seurat objects for each sample
ctx1 <- CreateSeuratObject(counts = ctx1.data, project = "ctx_7dpi", min.cells = 3, min.features = 200)

ctx11 <- CreateSeuratObject(counts = ctx11.data, project = "ctx_30dpi", min.cells = 3, min.features = 200)

wm2 <- CreateSeuratObject(counts = wm2.data, project = "wm_7dpi", min.cells = 3, min.features = 200)

wm12 <- CreateSeuratObject(counts = wm12.data, project = "wm_30dpi", min.cells = 3, min.features = 200)

ctx1[["percent.mt"]] <- PercentageFeatureSet(ctx1, pattern = "^mt-")

ctx11[["percent.mt"]] <- PercentageFeatureSet(ctx11, pattern = "^mt-")

wm2[["percent.mt"]] <- PercentageFeatureSet(wm2, pattern = "^mt-")

wm12[["percent.mt"]] <- PercentageFeatureSet(wm12, pattern = "^mt-")


#QC for no. of genes, molecules and mito percent
ctx1_f <- subset(ctx1, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & nCount_RNA < 12000 & percent.mt < 2)
ctx11_f <- subset(ctx11, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & nCount_RNA < 12000 & percent.mt < 2)
wm2_f <- subset(wm2, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & nCount_RNA < 12000 & percent.mt < 2)
wm12_f <- subset(wm12, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & nCount_RNA < 12000 & percent.mt < 2)



#Normalization with SCT Method. 

library(sctransform)

ctx1_f <- PercentageFeatureSet(ctx1_f, pattern = "^mt-", col.name = "percent.mt")

ctx11_f <- PercentageFeatureSet(ctx11_f, pattern = "^mt-", col.name = "percent.mt")

wm2_f <- PercentageFeatureSet(wm2_f, pattern = "^mt-", col.name = "percent.mt")

wm12_f <- PercentageFeatureSet(wm12_f, pattern = "^mt-", col.name = "percent.mt")


ctx1_f <- SCTransform(ctx1_f, vars.to.regress = "percent.mt", verbose = FALSE)

ctx11_f <- SCTransform(ctx11_f, vars.to.regress = "percent.mt", verbose = FALSE)

wm2_f <- SCTransform(wm2_f, vars.to.regress = "percent.mt", verbose = FALSE)

wm12_f <- SCTransform(wm12_f, vars.to.regress = "percent.mt", verbose = FALSE)


#Select Integration Features: Select the features that will be used for integration. These are usually the most variable features across the datasets.
features <- SelectIntegrationFeatures(object.list = list(ctx1_f, ctx11_f, wm2_f, wm12_f), nfeatures = 3000)

#Prepare for Integration: Prepare the Seurat objects for integration by running SCT normalization again on the selected features. This step harmonizes the scale of the data between the datasets.
PrepSCTIntegration(object = list(ctx1_f, ctx11_f, wm2_f, wm12_f), anchor.features = features)

#Find Integration Anchors: Identify anchors between the datasets, which are pairs of cells that are similar across the datasets. This is a critical step where the integration actually takes place.
anchors <- FindIntegrationAnchors(object.list = list(ctx1_f, ctx11_f, wm2_f, wm12_f), anchor.features = features)

#Integrate data
integrated_all <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Run downstream analyses
integrated_all <- RunPCA(integrated_all)
ElbowPlot(integrated_all)
integrated_all <- FindNeighbors(object = integrated_all, dims = 1:10)
integrated_all <- FindClusters(object = integrated_all)
integrated_all <- RunUMAP(object = integrated_all, reduction = "pca", dims = 1:10)
DimPlot(integrated_all, reduction = "umap")

VlnPlot(integrated_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(integrated_all, features = c("nFeature_RNA"))
VlnPlot(integrated_all, features = c("nCount_RNA"))
VlnPlot(integrated_all, features = c("percent.mt"))

plot1 <- FeatureScatter(integrated_all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(integrated_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

VizDimLoadings(integrated_all, dims = 1:2, reduction = "pca")

DimPlot(integrated_all, reduction = "pca")

Idents(integrated_all) = "integrated_snn_res.0.8"

VlnPlot(integrated_all, features = c("Cx3cr1", "Hexb", "C1qa", "C1qb" , "C1qc" , "P2ry12"), ncol = 3) #Mic/Mac

#Paso 1 — Extraer solo microglia
microglia <- subset(integrated_all, idents = c(13, 18))
cat("Células microglia:", ncol(microglia), "\n")
#Coge solo las células de los clusters 12 y 22 y crea un objeto Seurat nuevo. El cat te dice cuántas células tienes — espera entre 500 y 3000 células típicamente.

# SC TRANSFORM
microglia <- SCTransform(microglia, 
                         vars.to.regress = "percent.mt",
                         verbose = FALSE)
#Paso 3 — PCA
microglia <- RunPCA(microglia, npcs = 30, verbose = FALSE)
ElbowPlot(microglia, ndims = 30)
#Reduce la dimensionalidad. El ElbowPlot te muestra dónde cae la curva — busca el "codo". Para microglia suele ser entre 10 y 20 PCs.

microglia <- FindNeighbors(microglia,
                           reduction = "pca",
                           dims = 1:7)

microglia <- FindClusters(microglia, resolution = 0.2)

microglia <- RunUMAP(microglia,
                     reduction = "pca",
                     dims = 1:7,
                     min.dist = 0.3,
                     n.neighbors = 20)

DimPlot(microglia, label = TRUE, pt.size = 1.5) +
  ggtitle("Subpoblaciones microglia NSG") + NoLegend()

# Ver con números de cluster
DimPlot(microglia, label = TRUE, pt.size = 1.5) + 
  ggtitle("Subclusters microglia") + NoLegend()

# Ver por región
DimPlot(microglia, group.by = "region", pt.size = 1.5,
        cols = c("CTX" = "#534AB7", "WM" = "#E74C3C")) +
  ggtitle("CTX vs WM")

# Ver marcadores de cada subcluster
DotPlot(microglia,
        features = c("P2ry12", "Tmem119", "Cx3cr1",   # homeostáticos
                     "Trem2", "Apoe", "Lpl", "Spp1",   # DAM/WAM
                     "Mki67", "Top2a",                  # proliferación
                     "Il1b", "Tnf", "Nos2",             # pro-inflamatorio
                     "Ifit1", "Irf7"),                  # interferón
        dot.scale = 6,
        cols = c("lightgrey", "#534AB7")) +
  RotatedAxis() +
  ggtitle("Identidad subclusters microglia NSG")

microglia <- PrepSCTFindMarkers(microglia)

FindMarkers(microglia, ident.1 = 2,
            min.pct = 0.25,
            logfc.threshold = 0.5) %>% head(20)

# Eliminar cluster 2 contaminante
microglia <- subset(microglia, idents = c(0, 1))
cat("Células microglia limpias:", ncol(microglia), "\n")

# Ver distribución final
table(microglia$orig.ident)
table(microglia$region)

# Guardar
saveRDS(microglia, "microglia_NSG_limpia.rds")

# UMAP final limpio
DimPlot(microglia, label = TRUE, pt.size = 1.5) +
  ggtitle("Microglia NSG — cluster 0 y 1") + NoLegend()

DimPlot(microglia, group.by = "region", pt.size = 1.5,
        cols = c("CTX" = "#534AB7", "WM" = "#E74C3C")) +
  ggtitle("CTX vs WM")

DimPlot(microglia, group.by = "orig.ident", pt.size = 1.5) +
  ggtitle("Por muestra")

# Aislar solo microglia WM de tus datos NSG
microglia_wm <- subset(microglia, subset = region == "WM")
cat("Células microglia WM NSG:", ncol(microglia_wm), "\n")

# ── B1. Pseudobulk NSG: sumar counts por muestra ─────────────
#  AggregateExpression suma los raw counts de todas las células
#  de cada muestra (orig.ident) → resultado: genes x muestras
pseudo_nsg <- AggregateExpression(
  microglia_nsg,
  assays       = "RNA",
  group.by     = "orig.ident",   # agrupa por muestra (ctx_7dpi, wm_7dpi, etc.)
  return.seurat = FALSE
)$RNA

# Verificar estructura
cat("Dimensiones pseudobulk NSG:", dim(pseudo_nsg), "\n")
# Esperado: ~15000-20000 genes x 4 muestras
print(head(pseudo_nsg))
print(colnames(pseudo_nsg))   # ctx_7dpi, ctx_30dpi, wm_7dpi, wm_30dpi

# ── B2. Leer Excel BL6 ───────────────────────────────────────
BL6data <- read_excel("C:/Users/mdmsa/Desktop/Transcriptomics/Transcriptomic for moseq/GSE298799_ProcessedData_microglia.xlsx")


# Explorar estructura
# Ver dimensiones
dim(BL6data)

# Ver primeras filas y columnas
head(BL6data)

# Ver todos los nombres de columnas
colnames(BL6data)

# Nombre de columna de genes — ajusta si no es "GeneID"
col_genes <- "GeneID"   # <── cambia si necesario: "gene", "Gene", "gene_name"
genes_bl6 <- datos_bl6[[col_genes]]

# Seleccionar las 8 columnas control
# Ajusta el patrón según cómo se llamen: "^Con_", "^Ctrl_", "^WT_", "^Control"
patron_control <- "^Con_"   # <── cambia si necesario
cols_control   <- grep(patron_control, colnames(datos_bl6), value = TRUE)
cat("Muestras control BL6:", length(cols_control), "\n")
print(cols_control)

# Si no detecta nada, selección manual:
# cols_control <- c("nombre_col1", "nombre_col2", ...)

stopifnot("No encontré 8 muestras control" = length(cols_control) == 8)

# Crear matriz BL6
pseudo_bl6 <- as.matrix(datos_bl6[, cols_control])
rownames(pseudo_bl6) <- genes_bl6
storage.mode(pseudo_bl6) <- "integer"
pseudo_bl6 <- pseudo_bl6[!is.na(rownames(pseudo_bl6)), ]

cat("Dimensiones pseudobulk BL6:", dim(pseudo_bl6), "\n")

# ── B3. Homogeneizar nombres de genes ────────────────────────
#  NSG: genes en formato mouse (Cx3cr1, Trem2...)
#  BL6: verificar si está igual o en mayúsculas (CX3CR1, TREM2)
cat("\nGenes NSG (5 ejemplos):", head(rownames(pseudo_nsg), 5), "\n")
cat("Genes BL6 (5 ejemplos):", head(rownames(pseudo_bl6), 5), "\n")

#  Si BL6 tiene genes en MAYÚSCULAS y NSG en Title Case → convertir:
#  rownames(pseudo_bl6) <- str_to_title(rownames(pseudo_bl6))
#  (descomentar la línea de arriba si es necesario)

# Genes comunes
genes_comunes <- intersect(rownames(pseudo_nsg), rownames(pseudo_bl6))
cat("Genes en común:", length(genes_comunes), "\n")

pseudo_nsg_c <- pseudo_nsg[genes_comunes, ]
pseudo_bl6_c <- pseudo_bl6[genes_comunes, ]

# ── B4. Combinar matrices NSG + BL6 ──────────────────────────
counts_combinados <- cbind(pseudo_nsg_c, pseudo_bl6_c)
cat("Matriz combinada:", dim(counts_combinados), "\n")
# Esperado: genes x 12 (4 NSG + 8 BL6)


# ── B5. Metadata para DESeq2 ─────────────────────────────────
col_data <- data.frame(
  sample    = colnames(counts_combinados),
  genotype  = c(rep("NSG", ncol(pseudo_nsg_c)),
                rep("BL6", ncol(pseudo_bl6_c))),
  condition = c(rep("experimental", ncol(pseudo_nsg_c)),
                rep("control",      ncol(pseudo_bl6_c))),
  row.names = colnames(counts_combinados)
)

# Añadir región y timepoint para NSG (opcional, para análisis estratificado)
col_data$region    <- "unknown"
col_data$timepoint <- "unknown"
col_data[grep("ctx",  col_data$sample, ignore.case=TRUE), "region"]    <- "CTX"
col_data[grep("wm",   col_data$sample, ignore.case=TRUE), "region"]    <- "WM"
col_data[grep("7dpi", col_data$sample, ignore.case=TRUE), "timepoint"] <- "7dpi"
col_data[grep("30dpi",col_data$sample, ignore.case=TRUE), "timepoint"] <- "30dpi"
col_data[col_data$genotype == "BL6", "region"]    <- "WM"     # ajusta si BL6 es CTX o ambos
col_data[col_data$genotype == "BL6", "timepoint"] <- "control"

# Genotipo como factor, BL6 = referencia
col_data$genotype <- factor(col_data$genotype, levels = c("BL6", "NSG"))

cat("\nMetadata DESeq2:\n")
print(col_data)


# ============================================================
#  BLOQUE C — DESeq2: NSG microglia vs BL6 control
# ============================================================

# ── C1. Crear objeto DESeq2 ───────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_combinados,
  colData   = col_data,
  design    = ~ genotype     # comparación principal: NSG vs BL6
)

# Filtro mínimo: quitar genes con muy pocas lecturas
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genes tras filtrado:", nrow(dds), "\n")

# ── C2. Normalización y análisis DESeq2 ──────────────────────
dds <- DESeq(dds)
cat("\nResultNames:\n")
print(resultsNames(dds))

# ── C3. Extraer resultados NSG vs BL6 ────────────────────────
res <- results(
  dds,
  contrast  = c("genotype", "NSG", "BL6"),  # NSG vs BL6 (BL6=referencia)
  alpha     = 0.05
)

# Resumen
cat("\n=== RESUMEN DESeq2: NSG microglia vs BL6 control ===\n")
summary(res)

# Convertir a dataframe y ordenar
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  filter(!is.na(padj))

# Top genes up en NSG
cat("\nTop 20 genes upregulados en NSG (vs BL6):\n")
print(res_df %>% filter(log2FoldChange > 0) %>% head(20))

# Top genes down en NSG (es decir, más altos en BL6)
cat("\nTop 20 genes downregulados en NSG (más altos en BL6):\n")
print(res_df %>% filter(log2FoldChange < 0) %>% head(20))

# Guardar tabla completa
write.csv(res_df, "DEG_microglia_NSG_vs_BL6_DESeq2.csv", row.names = FALSE, quote = FALSE)
cat("\nTabla guardada: DEG_microglia_NSG_vs_BL6_DESeq2.csv\n")


# ============================================================
#  BLOQUE D — VISUALIZACIONES
# ============================================================

# ── D1. Volcano plot ─────────────────────────────────────────
res_df$sig <- "NS"
res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "Up en NSG"
res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Up en BL6"

# Genes a etiquetar: marcadores conocidos de microglia
genes_label <- c("P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Siglech",   # homeostáticos
                 "Trem2", "Apoe", "Lpl", "Spp1", "Cd9",               # DAM/WAM
                 "Il1b", "Tnf", "Nos2", "Ifit1", "Irf7",              # inflamación/IFN
                 "Mki67", "Top2a")                                     # proliferación

res_label <- res_df %>%
  filter(gene %in% genes_label | sig != "NS") %>%
  slice_head(n = 40)

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c(
    "NS"       = "grey70",
    "Up en NSG" = "#E74C3C",
    "Up en BL6" = "#2980B9"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data = res_label,
    aes(label = gene),
    size = 3, max.overlaps = 25, color = "black"
  ) +
  labs(
    title    = "Microglia NSG vs BL6 control",
    subtitle = "Positivo = up en NSG  |  Negativo = up en BL6",
    x        = "log2 Fold Change",
    y        = "-log10 (adj. p-value)",
    color    = ""
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(volcano)
ggsave("volcano_NSG_vs_BL6.pdf", volcano, width = 8, height = 7)

# ── D2. PCA de muestras (todos los 12) ───────────────────────
vst_data <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vst_data, intgroup = c("genotype", "condition"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                 color = genotype, shape = genotype, label = name)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text(vjust = -0.8, size = 3, color = "black") +
  scale_color_manual(values = c("NSG" = "#E74C3C", "BL6" = "#2980B9")) +
  labs(
    title = "PCA — microglia NSG vs BL6 (VST normalizado)",
    x     = paste0("PC1: ", pct_var[1], "% varianza"),
    y     = paste0("PC2: ", pct_var[2], "% varianza")
  ) +
  theme_classic(base_size = 12)

print(pca_plot)
ggsave("PCA_NSG_vs_BL6.pdf", pca_plot, width = 7, height = 6)

# ── D3. Heatmap: microglia homeostática vs DAM ───────────────
genes_heatmap <- c(
  # Homeostáticos (deben ser más altos en BL6)
  "P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Siglech", "Fcrls",
  # DAM / activados (pueden ser más altos en NSG)
  "Trem2", "Apoe", "Lpl", "Spp1", "Cd9", "Lgals3",
  # Inflamación
  "Il1b", "Tnf", "Il6",
  # Interferón
  "Ifit1", "Ifit3", "Irf7",
  # Proliferación
  "Mki67", "Top2a"
)

# Filtrar genes presentes en el dataset
genes_heatmap_ok <- genes_heatmap[genes_heatmap %in% rownames(vst_data)]
cat("Genes de heatmap presentes:", length(genes_heatmap_ok), "/", length(genes_heatmap), "\n")

if (length(genes_heatmap_ok) >= 5) {
  mat_heatmap  <- assay(vst_data)[genes_heatmap_ok, ]
  mat_scaled   <- t(scale(t(mat_heatmap)))   # z-score por gen
  
  ann_col <- data.frame(
    Genotipo = col_data$genotype,
    row.names = rownames(col_data)
  )
  
  pheatmap(
    mat_scaled,
    annotation_col  = ann_col,
    cluster_rows     = FALSE,
    cluster_cols     = TRUE,
    color            = colorRampPalette(c("#2980B9", "white", "#E74C3C"))(100),
    fontsize_row     = 9,
    fontsize_col     = 8,
    main             = "Microglia NSG vs BL6 — genes clave (z-score)",
    annotation_colors = list(Genotipo = c("NSG" = "#E74C3C", "BL6" = "#2980B9")),
    filename         = "heatmap_microglia_NSG_vs_BL6.pdf",
    width = 8, height = 8
  )
  cat("Heatmap guardado.\n")
}

# ── D4. Score de homeostasis vs DAM por muestra ──────────────
#  Útil para ver si NSG pierde el perfil homeostático
genes_homeo <- c("P2ry12","Tmem119","Cx3cr1","Hexb","Siglech")
genes_dam   <- c("Trem2","Apoe","Lpl","Spp1","Cd9")

genes_homeo_ok <- genes_homeo[genes_homeo %in% rownames(vst_data)]
genes_dam_ok   <- genes_dam[genes_dam %in% rownames(vst_data)]

if (length(genes_homeo_ok) >= 2 & length(genes_dam_ok) >= 2) {
  vst_mat <- assay(vst_data)
  
  score_df <- data.frame(
    sample      = colnames(vst_mat),
    genotype    = col_data$genotype,
    score_homeo = colMeans(vst_mat[genes_homeo_ok, ]),
    score_DAM   = colMeans(vst_mat[genes_dam_ok, ])
  )
  
  p_homeo <- ggplot(score_df, aes(x = genotype, y = score_homeo, color = genotype)) +
    geom_jitter(width = 0.15, size = 4) +
    scale_color_manual(values = c("NSG" = "#E74C3C", "BL6" = "#2980B9")) +
    labs(title = "Score homeostasis", y = "Expresión media VST", x = "") +
    theme_classic() + theme(legend.position = "none")
  
  p_dam <- ggplot(score_df, aes(x = genotype, y = score_DAM, color = genotype)) +
    geom_jitter(width = 0.15, size = 4) +
    scale_color_manual(values = c("NSG" = "#E74C3C", "BL6" = "#2980B9")) +
    labs(title = "Score DAM", y = "Expresión media VST", x = "") +
    theme_classic() + theme(legend.position = "none")
  
  scores_plot <- grid.arrange(p_homeo, p_dam, ncol = 2,
                              top = "Microglia NSG vs BL6 — estado funcional")
  ggsave("scores_homeo_DAM_NSG_vs_BL6.pdf", scores_plot, width = 8, height = 5)
  cat("Score plots guardados.\n")
}


# ============================================================
#  BLOQUE E — ANÁLISIS ESTRATIFICADO (opcional)
#  NSG WM vs BL6  |  NSG CTX vs BL6  |  7dpi vs 30dpi
# ============================================================

# ── E1. Solo WM: NSG_WM vs BL6 ───────────────────────────────
idx_wm_nsg <- which(col_data$genotype == "NSG" & col_data$region == "WM")
idx_bl6    <- which(col_data$genotype == "BL6")

counts_wm <- counts_combinados[, c(idx_wm_nsg, idx_bl6)]
meta_wm   <- col_data[c(idx_wm_nsg, idx_bl6), ]
meta_wm$genotype <- droplevels(meta_wm$genotype)

dds_wm <- DESeqDataSetFromMatrix(
  countData = counts_wm,
  colData   = meta_wm,
  design    = ~ genotype
)
keep_wm <- rowSums(counts(dds_wm) >= 5) >= 2
dds_wm  <- dds_wm[keep_wm, ]
dds_wm  <- DESeq(dds_wm, quiet = TRUE)

res_wm <- results(dds_wm, contrast = c("genotype", "NSG", "BL6"), alpha = 0.05)
cat("\n=== NSG WM vs BL6 ===\n")
summary(res_wm)

res_wm_df <- as.data.frame(res_wm) %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  filter(!is.na(padj))
write.csv(res_wm_df, "DEG_NSG_WM_vs_BL6.csv", row.names = FALSE)

# ── E2. Solo CTX: NSG_CTX vs BL6 ─────────────────────────────
idx_ctx_nsg <- which(col_data$genotype == "NSG" & col_data$region == "CTX")
counts_ctx  <- counts_combinados[, c(idx_ctx_nsg, idx_bl6)]
meta_ctx    <- col_data[c(idx_ctx_nsg, idx_bl6), ]
meta_ctx$genotype <- droplevels(meta_ctx$genotype)

dds_ctx <- DESeqDataSetFromMatrix(
  countData = counts_ctx,
  colData   = meta_ctx,
  design    = ~ genotype
)
keep_ctx <- rowSums(counts(dds_ctx) >= 5) >= 2
dds_ctx  <- dds_ctx[keep_ctx, ]
dds_ctx  <- DESeq(dds_ctx, quiet = TRUE)

res_ctx <- results(dds_ctx, contrast = c("genotype", "NSG", "BL6"), alpha = 0.05)
cat("\n=== NSG CTX vs BL6 ===\n")
summary(res_ctx)

res_ctx_df <- as.data.frame(res_ctx) %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  filter(!is.na(padj))
write.csv(res_ctx_df, "DEG_NSG_CTX_vs_BL6.csv", row.names = FALSE)


# ── E3. 7dpi vs 30dpi dentro de NSG (contexto temporal) ──────
idx_7dpi  <- which(col_data$genotype == "NSG" & col_data$timepoint == "7dpi")
idx_30dpi <- which(col_data$genotype == "NSG" & col_data$timepoint == "30dpi")

counts_time <- counts_combinados[, c(idx_7dpi, idx_30dpi)]
meta_time   <- col_data[c(idx_7dpi, idx_30dpi), ]
meta_time$timepoint <- factor(meta_time$timepoint, levels = c("7dpi", "30dpi"))

dds_time <- DESeqDataSetFromMatrix(
  countData = counts_time,
  colData   = meta_time,
  design    = ~ timepoint
)
keep_time <- rowSums(counts(dds_time) >= 5) >= 2
dds_time  <- dds_time[keep_time, ]
dds_time  <- DESeq(dds_time, quiet = TRUE)

res_time <- results(dds_time, contrast = c("timepoint", "30dpi", "7dpi"), alpha = 0.05)
cat("\n=== NSG 30dpi vs 7dpi ===\n")
summary(res_time)

res_time_df <- as.data.frame(res_time) %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  filter(!is.na(padj))
write.csv(res_time_df, "DEG_NSG_30dpi_vs_7dpi.csv", row.names = FALSE)


# ============================================================
#  RESUMEN FINAL
# ============================================================
cat("\n")
cat("=============================================\n")
cat(" PIPELINE COMPLETADO\n")
cat("=============================================\n")
cat(" Archivos generados:\n")
cat("   microglia_NSG_limpia.rds\n")
cat("   DEG_microglia_NSG_vs_BL6_DESeq2.csv\n")
cat("   DEG_NSG_WM_vs_BL6.csv\n")
cat("   DEG_NSG_CTX_vs_BL6.csv\n")
cat("   DEG_NSG_30dpi_vs_7dpi.csv\n")
cat("   volcano_NSG_vs_BL6.pdf\n")
cat("   PCA_NSG_vs_BL6.pdf\n")
cat("   heatmap_microglia_NSG_vs_BL6.pdf\n")
cat("   scores_homeo_DAM_NSG_vs_BL6.pdf\n")
cat("=============================================\n")
