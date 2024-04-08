# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 2k Sorted Cells from Human Glioblastoma Multiforme, 3’ v3.1
# data source: https://www.10xgenomics.com/datasets/2-k-sorted-cells-from-human-glioblastoma-multiforme-3-v-3-1-3-1-standard-6-0-0         

# setwd("~/Documents/Rtutorial/Brain_Tumor_3p_raw_feature_bc_matrix")


# import the libraries
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
Brain_Tumor.data <- Read10X_h5(filename = 'C:/Users/Dr Medite/Documents/Rtutorial/Brain_Tumor_3p_raw_feature_bc_matrix.h5')

# Initialize the Seurat object with the raw (non-normalized data).
Brain_Tumor <- CreateSeuratObject(counts = Brain_Tumor.data, project = "Brain_Tumor", min.cells = 3, min.features = 200)
Brain_Tumor

# QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Brain_Tumor[["percent.mt"]] <- PercentageFeatureSet(Brain_Tumor, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Brain_Tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(Brain_Tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Brain_Tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Normalizing the data
Brain_Tumor <- NormalizeData(Brain_Tumor)

# Identification of highly variable features (feature selection)
Brain_Tumor <- FindVariableFeatures(Brain_Tumor, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Brain_Tumor), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Brain_Tumor)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#Scaling the data
all.genes <- rownames(Brain_Tumor)
Brain_Tumor <- ScaleData(Brain_Tumor, features = all.genes)

#Perform linear dimensional reduction
Brain_Tumor <- RunPCA(Brain_Tumor, features = VariableFeatures(object = Brain_Tumor))

# Examine and visualize PCA results a few different ways
print(Brain_Tumor[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Brain_Tumor, dims = 1:2, reduction = "pca")

DimPlot(Brain_Tumor, reduction = "pca") + NoLegend()

DimHeatmap(Brain_Tumor, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(Brain_Tumor, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset
ElbowPlot(Brain_Tumor)

# Cluster the cells
Brain_Tumor <- FindNeighbors(Brain_Tumor, dims = 1:10)
Brain_Tumor <- FindClusters(Brain_Tumor, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(Brain_Tumor), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
Brain_Tumor <- RunUMAP(Brain_Tumor, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Brain_Tumor, reduction = "umap", label = TRUE)

saveRDS(Brain_Tumor, file = "../Brain_Tumor_tutorial.rds")

# Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(Brain_Tumor, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(Brain_Tumor, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Brain_Tumor.markers <- FindAllMarkers(Brain_Tumor, only.pos = TRUE)
Brain_Tumor.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# ROC test returns the ‘classification power’ for any individual marker
cluster0.markers <- FindMarkers(Brain_Tumor, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(Brain_Tumor, features = c("MTRNR2L12", "MT-ATP6"))

# you can plot raw counts as well
VlnPlot(Brain_Tumor, features = c("MT-CO1", "MT-CO2"), slot = "counts", log = TRUE)

# you can plot raw counts as well
VlnPlot(Brain_Tumor, features = c("MTRNR2L12", "MT-ATP6"), slot = "counts", log = TRUE)

FeaturePlot(Brain_Tumor, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "MT-CO1",
                               "CD8A"))

Brain_Tumor.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(Brain_Tumor, features = top10$gene) + NoLegend()


