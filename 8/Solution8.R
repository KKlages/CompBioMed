install.packages("dplyr")
install.packages("patchwork")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Seurat")

library(dplyr)
library(Seurat)
library(patchwork)

pbmc.data <- Read10X(data.dir = "data_Seurat/")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

dense.size <- object.size(as.matrix(pbmc.data))
dense.size

sparse.size <- object.size(pbmc.data)
sparse.size
dense.size/sparse.size

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. In Seurat v5, Normalized values are stored in pbmc[["RNA"]]$data.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "pbmc_tutorial.rds")

cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
saveRDS(pbmc, file = "pbmc3k_final.rds")

##############
install.packages("hdf5r")
library(hdf5r)

D1_data <- Read10X_h5("data_new/GSM5022599_D1_filtered_feature_bc_matrix.h5")
D5_data <- Read10X_h5("data_new/GSM5022604_D6-D10_Pool_CMG_filtered_feature_bc_matrix.h5")
D1 <- CreateSeuratObject(counts = D1_data, project = "D1")
D5 <- CreateSeuratObject(counts = D5_data, project = "D5")

D1[["percent.mt"]] <- PercentageFeatureSet(D1, pattern = "^MT-")
head(D1@meta.data, 5)
D5[["percent.mt"]] <- PercentageFeatureSet(D5, pattern = "^MT-")
head(D5@meta.data, 5)
D1 <- subset(D1, subset =  percent.mt < 10)
D5 <- subset(D5, subset =  percent.mt < 10)

D1 <- NormalizeData(D1)
D5 <- NormalizeData(D5)
datasets <- list(D1, D5)

anchors <- FindIntegrationAnchors(object.list = datasets, normalization.method = c("LogNormalize"))

integrated_data <- IntegrateData(anchorset = anchors)
D1D5 = integrated_data
VlnPlot(D1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(D5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(D1D5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


D1D5 <- FindVariableFeatures(D1D5, selection.method = "vst", nfeatures = 2000)
top10D1D5 <- head(VariableFeatures(D1D5), 10)

plot1 <- VariableFeaturePlot(D1D5)
plot2 <- LabelPoints(plot = plot1, points = top10D1D5, repel = TRUE)
plot1 + plot2


all.genes <- rownames(D1D5)
D1D5 <- ScaleData(D1D5, features = all.genes)

D1D5 <- RunPCA(D1D5, features = VariableFeatures(object = D1D5))
print(D1D5[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(D1D5, dims = 1:2, reduction = "pca")
DimPlot(D1D5, reduction = "pca") + NoLegend()
DimHeatmap(D1D5, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(D1D5, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(D1D5)


D1D5 <- FindNeighbors(D1D5, dims = 1:10)
D1D5 <- FindClusters(D1D5, resolution = 0.5)
head(Idents(D1D5), 5)

D1D5 <- RunUMAP(D1D5, dims = 1:10)
DimPlot(D1D5, reduction = "umap")
saveRDS(D1D5, file = "D1D5_tutorial.rds")


D1D5.markers <- FindAllMarkers(D1D5)
D1D5.markers %>%
  group_by(cluster) %>%R
dplyr::filter(avg_log2FC > 1)

# Get the unique cluster IDs
cluster_ids <- unique(Idents(D1D5))

# Create new, more descriptive cluster names
new_cluster_names <- paste0("Cluster ", 1:length(cluster_ids))

# Rename the cluster identities
Idents(D1D5) <- plyr::mapvalues(Idents(D1D5), 
                                from = cluster_ids, 
                                to = new_cluster_names)

# Plot the UMAP visualization with the new cluster names
DimPlot(D1D5, reduction = "umap", label = TRUE, label.size = 6) + 
  ggtitle("Integrated D1D5 Clustering") +
  theme(plot.title = element_text(size = 16, face = "bold"))

# Visualize cell markers
FeaturePlot(D1D5, features = c("MS4A1", "PECAM1", "EPCAM", "PDGFRB", "CD68", "JCHAIN", "MKI67", "CD3D"), 
            cols = c("lightgrey", "blue"), 
            pt.size = 1.5)


############
library(Seurat)
library(ggplot2)

spatial_data <- readRDS("data_new_2.rds")

str(spatial_data)
typeof(spatial_data)

VlnPlot(spatial_data, features = c("nCount_Spatial", "nFeature_Spatial"))


summary(spatial_data$nCount_Spatial)
summary(spatial_data$nFeature_Spatial)


length(Cells(spatial_data))


length(Features(spatial_data))

SpatialFeaturePlot(spatial_data, features = spatial_data@assays$SCT@var.features[1:3])