############################################################
# TITLE: scRNA-seq Analysis of Colorectal Cancer (Seurat)
# AUTHOR: Jibin Varghese
# DESCRIPTION: Preprocessing, clustering, annotation,
#              and export for AI-based biomarker analysis
############################################################

############################################################
# STEP 1 — LOAD LIBRARIES
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

theme_set(theme_classic(base_size = 14))

############################################################
# STEP 2 — LOAD DATA
############################################################

data_dir <- "C:/Users/Jibin Varghese/Downloads/Parent_Visium_Human_ColorectalCancer_filtered_feature_bc_matrix/filtered_feature_bc_matrix"

data <- Read10X(data.dir = data_dir)

seurat_obj <- CreateSeuratObject(counts = data)

cat("Initial dimensions:\n")
print(dim(seurat_obj))

############################################################
# STEP 3 — QUALITY CONTROL
############################################################

# Mitochondrial gene percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# QC visualization
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

############################################################
# STEP 4 — FILTER CELLS
############################################################

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 5
)

cat("After QC filtering:\n")
print(dim(seurat_obj))

############################################################
# STEP 5 — REMOVE MITOCHONDRIAL GENES
############################################################

mt_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
cat("MT genes removed:", length(mt_genes), "\n")

seurat_obj <- seurat_obj[!rownames(seurat_obj) %in% mt_genes, ]

############################################################
# STEP 6 — NORMALIZATION
############################################################

seurat_obj <- NormalizeData(seurat_obj)

############################################################
# STEP 7 — VARIABLE FEATURES
############################################################

seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)

top10 <- head(VariableFeatures(seurat_obj), 10)

plot1 <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

############################################################
# STEP 8 — SCALING
############################################################

seurat_obj <- ScaleData(seurat_obj)

############################################################
# STEP 9 — PCA
############################################################

seurat_obj <- RunPCA(seurat_obj)

DimPlot(seurat_obj, reduction = "pca") +
  ggtitle("PCA Projection")

ElbowPlot(seurat_obj)

############################################################
# STEP 10 — CLUSTERING
############################################################

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:12)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

############################################################
# STEP 11 — UMAP
############################################################

seurat_obj <- RunUMAP(seurat_obj, dims = 1:12)

DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
  ggtitle("UMAP Clusters")

############################################################
# STEP 12 — MARKER GENES
############################################################

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

############################################################
# STEP 13 — CELL TYPE ANNOTATION
############################################################

tumor_markers <- c("EPCAM","KRT8","KRT18","KRT19")
immune_markers <- c("PTPRC","CD68","LYZ")
fibro_markers <- c("COL1A1","COL1A2","DCN","LUM")

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(tumor_markers, immune_markers, fibro_markers),
  name = c("Tumor","Immune","Fibro")
)

scores <- FetchData(seurat_obj, vars = c("Tumor1","Immune2","Fibro3"))

seurat_obj$cell_type <- apply(scores, 1, function(x) {
  c("Tumor","Immune","Fibro")[which.max(x)]
})

DimPlot(seurat_obj, group.by = "cell_type", label = TRUE) +
  ggtitle("Cell Type Annotation")

############################################################
# STEP 14 — DATA ANALYSIS & RESULTS
############################################################

# Expression matrix
expr_data <- GetAssayData(seurat_obj, layer = "data")

# Summary statistics
summary_table <- data.frame(
  Mean = rowMeans(expr_data),
  Variance = apply(expr_data, 1, var)
)

write.csv(summary_table, "gene_summary_statistics.csv")

# Histogram
hist(as.vector(expr_data),
     breaks = 50,
     col = "skyblue",
     main = "Gene Expression Distribution",
     xlab = "Expression")

# Cluster distribution
cluster_counts <- table(Idents(seurat_obj))
write.csv(cluster_counts, "cluster_distribution.csv")

barplot(cluster_counts,
        col = "steelblue",
        main = "Cells per Cluster",
        ylab = "Cell Count")

# Heatmap of top markers
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 3)

DoHeatmap(seurat_obj, features = top_markers$gene) +
  ggtitle("Top Marker Genes Heatmap")

# Save markers
write.csv(markers, "marker_genes.csv")

############################################################
# STEP 15 — EXPORT FOR PYTHON (ML/DL)
############################################################

expr_matrix <- as.matrix(expr_data)
expr_matrix <- t(expr_matrix)
expr_matrix <- expr_matrix[, VariableFeatures(seurat_obj)]

cat("Final matrix shape:\n")
print(dim(expr_matrix))

write.csv(expr_matrix, "expr_matrix.csv")

clusters <- Idents(seurat_obj)
write.csv(clusters, "clusters.csv")

############################################################
# END OF PIPELINE
############################################################