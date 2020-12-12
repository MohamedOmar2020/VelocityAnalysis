
library(loomR)


#looms <- list("~/Desktop/Looms/MAD1.loom", "~/Desktop/Looms/MAD2.loom", "~/Desktop/Looms/MAD3.loom", "~/Desktop/Looms/MAD4.loom", "~/Desktop/Looms/MAD5.loom", "~/Desktop/Looms/MAD6.loom")

#combine(looms = list(lfile1, lfile2), "LoomCombined.loom", chunk.size = 1000, skip.validate = TRUE)


#lfile1 <- connect(filename = "~/Desktop/Looms/MAD1.loom", mode = "r+", skip.validate = T)
#lfile2 <- connect(filename = "~/Desktop/Looms/MAD2.loom", mode = "r+", skip.validate = T)

lfile1$row.attrs


#Access all gene names
gene.names <- lfile[["row_attrs/gene_name"]][]
head(x = gene.names)

clusters = lfile1[["col_attrs/cluster_seurat"]]

lfile$close_all()
lfile1$close_all()
lfile2$close_all()




