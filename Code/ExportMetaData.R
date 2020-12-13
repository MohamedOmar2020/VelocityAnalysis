

# we will export all of the necessary meta-data from our Seurat object that is needed for our RNA Velocity object. This includes:
  
  
rm(list = ls())

library(Seurat)


# Load the seurat objec
IntestineSeurat <- readRDS("~/Desktop/RDS_umap_Epithelium_AMO1")



# Save cell ID
write.csv(Cells(IntestineSeurat), file = "cellID_obs.csv", row.names = FALSE)

# Save UMAP coords
write.csv(Embeddings(IntestineSeurat, reduction = "umap"), file = "cell_embeddings.csv")

# Save Cluster annotations
Clusters <- data.frame(Seurat_Intestine$seurat_clusters)
Clusters$CellID <- rownames(Clusters)
write.csv(Clusters, file = "clusters.csv", row.names = F)




#############################3
## to extract the colors >> not working yet
p <- Seurat::DimPlot(IntestineSeurat) # Generate the tSNE plot, but save it as an object
pbuild <- ggplot2::ggplot_build(p) # Use ggplot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # Pull the data used for the plot
#The colors used, in hexadecimal, are in the colour column of pdata, the groups in the group column.

#If you want a vector of the colors used you can do the following:
  
pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
ucols <- unique(pdata$colour) # Get a vector of unique colors
names(ucols) <- unique(pdata$group) # Add the groups to the vector of colors as names
names(ucols) <- unique(IntestineSeurat$sample.cell_type)

write.csv(ucols, file = "cluster_colors.csv")


