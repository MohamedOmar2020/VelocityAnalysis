

library(devtools)
install_github("velocyto-team/velocyto.R")
library(velocyto.R)
library(loomR)
library(SCopeLoomR)

## Load loom data
ldat <- connect("~/Documents/Research/Projects/Velocity/Data/BigLoomFiltered.loom", mode = "r+", skip.validate = T)

fit.quantile <- 0.1

## pull out spliced and unspliced matrices from AnnData
emat <- as.matrix(t(ldat$layers['spliced']))
nmat <- as.matrix(t(adata$layers['unspliced']))
cells <- adata$obs_names$values
genes <- adata$var_names$values
colnames(emat) <- colnames(nmat) <- cells
rownames(emat) <- rownames(nmat) <- genes

## pull out PCA 
pcs <- adata$obsm['X_pca']
rownames(pcs) <- cells
cell.dist <- as.dist(1-cor(t(pcs))) ## cell distance in PC space

## filter genes
gexp1 <- log10(rowSums(emat)+1)
gexp2 <- log10(rowSums(nmat)+1)
#plot(gexp1, gexp2)
good.genes <- genes[gexp1 > 2 & gexp2 > 1]



## Close loom connection
ldat$close_all()