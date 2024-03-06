#scRNA seq analysis + Seurat
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)

#Load data
pbmc.data <- Read10X(data.dir = "C:/Users/ignac/Desktop/R/Seurat")

#Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200, project = "pbmc3k")
#min cell es la cantidad minima de cells que tiene >0 de ese gen y min features es la cantidad minima de genes por cell

pbmc

#QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#Violin Plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Scatter Plot
FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

#Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")

#Find Highly Variable Genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
top10

#Plot genes
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#Scaling data -> Put the mean at 0
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#PCA viz -> Different ways
print(pbmc[["pca"]], dims = 1:3, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, nfeatures = 15, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1, cells = 800, balanced = TRUE)
DimHeatmap(pbmc, dims = 2, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:2, cells = 500, balanced = TRUE)

#Elbow Plot -> Choose the ammount of PC to take into consideration according the variance
ElbowPlot(pbmc)
#In this case, I chose 10

#Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

#Find markers for cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
VlnPlot(pbmc, features = c(row.names(cluster1.markers[1]), row.names(cluster1.markers[2])))

#Find markers for cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
VlnPlot(pbmc, features = c(row.names(cluster2.markers[1]), row.names(cluster2.markers[2])))

#Find all markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Visulization
x <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
FeaturePlot(pbmc, features = x$gene[1:4])
FeaturePlot(pbmc, features = x$gene[5:8])

#Another way of visualizing 
p <- FeaturePlot(pbmc, features = c("CCR7", "FOLR3", "CD40LG", "LINC00926", "GZMK", "CKB", "AKR1C3", "SERPINF1"), combine = FALSE)
p <- lapply(X = p, FUN = function(x) x + 
              theme(plot.title = element_text(size = 8)) +
              theme(axis.title.y = element_text(size = 5)) +
              theme(axis.title.x = element_text(size = 5)) +
              theme(axis.text.y = element_text(size = 5)) +
              theme(axis.text.x = element_text(size = 5)) +
              theme(legend.position = "none"))
CombinePlots(plots = p)

#Heatmap
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10
p2 <- DoHeatmap(pbmc, features = top10$gene, group.bar.height = 0.01,size=3,combine = FALSE) 
p2 <- lapply(X = p2, FUN = function(x) x + 
               theme(plot.title = element_text(size = 8)) +
               theme(axis.title.y = element_text(size = 5)) +
               theme(axis.title.x = element_text(size = 5)) +
               theme(axis.text.y = element_text(size = 3)) +
               theme(legend.position = "none")  )
CombinePlots(plots = p2)

#Assign cell to cluster
new.cluster.id <- c("Naive CD4 T", "Memory CD4 T", "CD14 + MONO", "B", "CD8 T", "FCGR3A + MONO", "NK", "DC", "Platelet")
names(new.cluster.id) <- levels(pbmc)

pbmc <- RenameIdents(pbmc, new.cluster.id)
DimPlot(pbmc, reduction = "pca", label = TRUE, pt.size = 0.5)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)