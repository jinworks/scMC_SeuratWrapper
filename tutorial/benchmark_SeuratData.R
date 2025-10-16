#
rm(list=ls())
library(scMC)
library(Seurat)
library(patchwork)
library(dplyr)

setwd("/Users/suoqinjin/Documents/scMC")

################## Systematic comparative analysis of human PBMC #################
# import data
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
InstallData("pbmcsca")
data("pbmcsca")

######### Part I: Setup the a list of Seurat objects, one per dataset ##############
pbmcsca.list <- SplitObject(pbmcsca, split.by = "Method")
# step1. pre-processing each dataset
## Normalizing, scaling the data and feature selection
pbmcsca.list <- lapply(X = pbmcsca.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # perform scaling on the previously identified variable features
  x <- ScaleData(x, verbose = FALSE)
})

########### Part II: Perform an integrated analysis using scMC ###########
future::plan("multiprocess", workers = 4)
ptm = Sys.time()
combined <- RunscMC(object.list = pbmcsca.list)
execution.time = Sys.time() - ptm
Misc(combined, slot = 'execution.time') <- execution.time


########### Part III: Run the standard workflow for visualization and clustering ###########
nPC = 40
combined <- FindNeighbors(combined, reduction = "scMC", dims = 1:nPC)
combined <- FindClusters(combined, algorithm = 4, resolution = 0.8)
combined <- RunUMAP(combined, reduction='scMC', dims = 1:nPC)

## Visualization
DimPlot(combined, group.by = c("Method", "ident", "CellType"), ncol = 3)

save.image(combined, file = "scMC_pbmcsca.RData")


################## Interferon-stimulated and control PBMC #################
# import data
library(SeuratData)
InstallData("ifnb")
data("ifnb")

######### Part I: Setup the a list of Seurat objects, one per dataset ##############
ifnb.list <- SplitObject(ifnb, split.by = "stim")
# step1. pre-processing each dataset
## Normalizing, scaling the data and feature selection
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # perform scaling on the previously identified variable features
  x <- ScaleData(x, verbose = FALSE)
})

########### Part II: Perform an integrated analysis using scMC ###########
future::plan("multiprocess", workers = 4)
ptm = Sys.time()
combined <- RunscMC(object.list = ifnb.list)
execution.time = Sys.time() - ptm
Misc(combined, slot = 'execution.time') <- execution.time

########### Part III: Run the standard workflow for visualization and clustering ###########
nPC = 40
combined <- FindNeighbors(combined, reduction = "scMC", dims = 1:nPC)
combined <- FindClusters(combined, algorithm = 4, resolution = 0.05)
combined <- RunUMAP(combined, reduction='scMC', dims = 1:nPC)

## Visualization
DimPlot(combined, group.by = c("stim", "ident", "seurat_annotations"), ncol = 3)

save.image(combined, file = "scMC_ifnb.RData")



################## Eight human pancreatic islet datasets #################
# import data
library(SeuratData)
InstallData("panc8")
data("panc8")

######### Part I: Setup the a list of Seurat objects, one per dataset ##############
panc8.list <- SplitObject(panc8, split.by = "replicate")
# step1. pre-processing each dataset
## Normalizing, scaling the data and feature selection
panc8.list <- lapply(X = panc8.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # perform scaling on the previously identified variable features
  x <- ScaleData(x, verbose = FALSE)
})

########### Part II: Perform an integrated analysis using scMC ###########
future::plan("multiprocess", workers = 4)
ptm = Sys.time()
combined <- RunscMC(object.list = panc8.list)
execution.time = Sys.time() - ptm
Misc(combined, slot = 'execution.time') <- execution.time

########### Part III: Run the standard workflow for visualization and clustering ###########
nPC = 40
combined <- FindNeighbors(combined, reduction = "scMC", dims = 1:nPC)
combined <- FindClusters(combined, algorithm = 4, resolution = 0.05)
combined <- RunUMAP(combined, reduction='scMC', dims = 1:nPC)

## Visualization
DimPlot(combined, group.by = c("replicate", "ident", "celltype"), ncol = 3)

save.image(combined, file = "scMC_panc8.RData")




