

# we will export all of the necessary meta-data from our Seurat object that is needed for our RNA Velocity object. This includes:
  
  
rm(list = ls())

library(Seurat)

IntestineSeurat <- readRDS("~/Desktop/RDS_umap_Epithelium_AMO1")


write.csv(Cells(IntestineSeurat), file = "cellID_obs.csv", row.names = FALSE)

CellID <- Cells(IntestineSeurat)
CellID <- gsub("\\-.+", "", CellID)

write.csv(Cells(IntestineSeurat), file = "cellID_obs.csv", row.names = FALSE)

write.csv(Embeddings(IntestineSeurat, reduction = "umap"), file = "cell_embeddings.csv")

Clusters <- data.frame(IntestineSeurat$sample.cell_type)
Clusters$CellID <- rownames(Clusters)
write.csv(Clusters, file = "clusters.csv", row.names = F)

data.frame(IntestineSeurat$sample.cell_type)


Idents(IntestineSeurat)
