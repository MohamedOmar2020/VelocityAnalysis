

library(devtools)
install_github("velocyto-team/velocyto.R")
library(velocyto.R)
library(loomR)
library(SCopeLoomR)

## Load loom data
ldat <- open_loom("~/Documents/Research/Projects/Velocity/Loom_Combined.loom", mode = "r+")

##Normalize and cluster cells using pagoda2
# Using spliced expression matrix as input to pagoda2.
emat <- ldat$spliced
# this is where one woudl do some filtering
emat <- emat[,colSums(emat) >= 1000]

ldat$close_all()

gr_all <- merge(IntestineSeurat, as.Seurat(ldat)[, colnames(IntestineSeurat)])



