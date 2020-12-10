#Install devtools if it's not already installed

if (!require(devtools)) {
  install.packages("devtools")
}

if (!require(BiocManager)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DropletUtils", "BSgenome.Mmusculus.UCSC.mm10", 
                       "AnnotationHub", "SingleR"))

# Install from GitHub
devtools::install_github("BUStools/BUSpaRse")
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github("velocyto-team/velocyto.R")


library(BUSpaRse)
library(Seurat)
library(SeuratWrappers)
library(BSgenome.Mmusculus.UCSC.mm10)
library(AnnotationHub)
library(zeallot) # For %<-% that unpacks lists in the Python manner
library(DropletUtils)
library(tidyverse)
library(GGally) # For ggpairs
library(velocyto.R)
library(SingleR)
library(scales)
library(plotly)
theme_set(theme_bw())



#setwd("~/Desktop")

d <- "~/Desktop/MAD2_Velo/counts_unfiltered"


c(spliced, unspliced) %<-% read_velocity_output(spliced_dir = d,
                                                spliced_name = "spliced",
                                                unspliced_dir = d, 
                                                unspliced_name = "unspliced")

sum(unspliced@x) / (sum(unspliced@x) + sum(spliced@x))
dim(spliced)
dim(unspliced)

tot_count <- Matrix::colSums(spliced)
summary(tot_count)



bc_rank <- barcodeRanks(spliced)
bc_uns <- barcodeRanks(unspliced)

#Knee plot for filtering empty droplets
#' 
#' Visualizes the inflection point to filter empty droplets. This function plots 
#' different datasets with a different color. Facets can be added after calling
#' this function with `facet_*` functions.
#' 
#' @param bc_ranks A named list of output from `DropletUtil::barcodeRanks`.
#' @return A ggplot2 object.
#' @importFrom tibble tibble
#' @importFrom purrr map map_dbl
#' @importFrom dplyr distinct
#' @importFrom ggplot2 geom_line geom_hline geom_vline scale_x_log10 scale_y_log10
#' @importFrom tidyr unnest
#' @export
knee_plot <- function(bc_ranks) {
  # purrr pluck shorthand doesn't work on S4Vector DataFrame
  knee_plt <- tibble(rank = map(bc_ranks, ~ .x[["rank"]]), 
                     total = map(bc_ranks, ~ .x[["total"]]),
                     dataset = names(bc_ranks)) %>% 
    unnest(cols = c(rank, total)) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = map_dbl(bc_ranks, ~ metadata(.x)[["inflection"]]),
                  rank_cutoff = map_dbl(bc_ranks, 
                                        ~ max(.x$rank[.x$total >
                                                        metadata(.x)[["inflection"]]])),
                  dataset = names(bc_ranks))
  p <- ggplot(knee_plt, aes(rank, total, color = dataset)) +
    geom_line() +
    geom_hline(aes(yintercept = inflection, color = dataset), 
               data = annot, linetype = 2) +
    geom_vline(aes(xintercept = rank_cutoff, color = dataset),
               data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank", y = "Total UMIs")
  return(p)
}


knee_plot(list(spliced = bc_rank, unspliced = bc_uns)) +
  coord_flip()

# Can only plot barcodes with both spliced and unspliced counts
bcs_inter <- intersect(colnames(spliced), colnames(unspliced))
s <- colSums(spliced[,bcs_inter])
u <- colSums(unspliced[,bcs_inter])
# Grid points
sr <- sort(unique(exp(round(log(s)*100)/100)))
ur <- sort(unique(exp(round(log(u)*100)/100)))
# Run naive approach
bc2 <- bc_ranks2(s, u, sr, ur)

# can't turn color to lot scale unless log values are plotted
z_use <- log10(bc2)
z_use[is.infinite(z_use)] <- NA
plot_ly(x = sr, y = ur, z = z_use) %>% add_surface() %>% 
  layout(scene = list(xaxis = list(title = "Total spliced UMIs", type = "log"),
                      yaxis = list(title = "Total unspliced UMIs", type = "log"),
                      zaxis = list(title = "Rank (log10)")))

bcs_use <- colnames(spliced)[tot_count > metadata(bc_rank)$inflection]
# Remove genes that aren't detected
tot_genes <- Matrix::rowSums(spliced)
genes_use <- rownames(spliced)[tot_genes > 0]
sf <- spliced[genes_use, bcs_use]
uf <- unspliced[genes_use, bcs_use]


dim(sf)

rownames(sf) <- str_remove(rownames(sf), "\\.\\d+")
rownames(uf) <- str_remove(rownames(uf), "\\.\\d+")


seu <- CreateSeuratObject(sf, assay = "sf") %>% 
  SCTransform(assay = "sf", new.assay.name = "spliced")



mouse.rnaseq <- MouseRNAseqData(ensembl = TRUE)



annots <- SingleR(GetAssayData(seu, assay = "spliced", slot = "data"), 
                  ref = mouse.rnaseq, labels = colData(mouse.rnaseq)$label.fine,
                  de.method = "wilcox", method = "single", BPPARAM = MulticoreParam(4))


seu$cell_type <- annots$pruned.labels


seu[["uf"]] <- CreateAssayObject(uf)
seu <- SCTransform(seu, assay = "uf", new.assay.name = "unspliced")

cols_use <- c("nCount_sf", "nFeature_sf", "nCount_uf", "nFeature_uf")
VlnPlot(seu, cols_use, pt.size = 0.1, ncol = 1, group.by = "cell_type")


# Helper functions for ggpairs
log10_diagonal <- function(data, mapping, ...) {
  ggally_densityDiag(data, mapping, ...) + scale_x_log10()
}
log10_points <- function(data, mapping, ...) {
  ggally_points(data, mapping, ...) + scale_x_log10() + scale_y_log10()
}

ggpairs(seu@meta.data, columns = cols_use,
        upper = list(continuous = "cor"),
        diag = list(continuous = log10_diagonal),
        lower = list(continuous = wrap(log10_points, alpha = 0.1, size=0.3)),
        progress = FALSE)



DefaultAssay(seu) <- "spliced"
seu <- RunPCA(seu, verbose = FALSE, npcs = 70)
ElbowPlot(seu, ndims = 70)


DimPlot(seu, reduction = "pca",
        group.by = "cell_type", pt.size = 0.5, label = TRUE, repel = TRUE) +
  scale_color_brewer(type = "qual", palette = "Set2")




























