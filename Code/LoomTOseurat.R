
rm(list = ls())

library(loomR)
library(Seurat)

# Load the processed seurat object from Anxo
Seurat_Intestine <- readRDS("./Data/RDS_umap_Epithelium_AMO1")

## Load loom data
ldat <- connect("~/Documents/Research/Projects/Velocity/Data/BigLoomFiltered.loom", mode = "r+", skip.validate = T)

#ldat$row.attrs

Seurat <- as.Seurat(ldat, cells = "obs_names", features = "gene_name")

#Seurat[["RNA"]] <- Seurat[["spliced"]]

VlnPlot(Seurat, features = c("Sparc", "Ftl1", "Junb", "Ccl4"), ncol = 2, pt.size = 0.1)

ldat$close_all()


Seurat[["percent.mt"]] <- PercentageFeatureSet(Seurat, pattern = "^mt-")

head(Seurat@meta.data, 5)

VlnPlot(Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Seurat <- subset(Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Seurat <- NormalizeData(Seurat, normalization.method = "LogNormalize")

#Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(Seurat), 10)


# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(Seurat)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2


#all.genes <- rownames(Seurat)

#Seurat <- ScaleData(Seurat, features = all.genes)

## Instead of NormalizeData , findVariableFeatures, and ScaleData : we will use SCTransform
Seurat <- SCTransform(object = Seurat, assay = "RNA", vars.to.regress = "percent.mt", verbose = T)


## Run PCA
Seurat <- RunPCA (Seurat, verbose = T)

print(Seurat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Seurat, dims = 1:2, reduction = "pca")


DimPlot(Seurat, reduction = "pca")


DimHeatmap(Seurat, dims = 1:15, cells = 500, balanced = TRUE)


#Seurat <- JackStraw(Seurat, num.replicate = 100)

#Seurat <- ScoreJackStraw(Seurat, dims = 1:20)

#JackStrawPlot(Seurat, dims = 1:15)

ElbowPlot(Seurat)

Seurat <- FindNeighbors(Seurat, dims = 1:15)

Seurat <- FindClusters(Seurat, resolution = 0.3)

head(Idents(Seurat), 5)

Seurat <- RunUMAP(Seurat, dims = 1:15)

DimPlot(Seurat, reduction = "umap")


#Clusters <- IntestineSeurat$seurat_clusters
#Cells <- intersect(colnames(IntestineSeurat), colnames(Seurat))
#Seurat <- Seurat[, Cells]
  
#Seurat$seurat_clusters <- Clusters  
  
cluster1.markers <- FindMarkers(Seurat, ident.1 = 0, min.pct = 0.25)
cluster2.markers <- FindMarkers(Seurat, ident.1 = 1, min.pct = 0.25)

cluster1.markers_x <- FindMarkers(Seurat_Intestine, ident.1 = "Tumor_dkot", min.pct = 0.25)

  
##############################################3
## YEAH BABY !!!!
  

