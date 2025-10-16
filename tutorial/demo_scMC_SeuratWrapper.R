#
rm(list=ls())
library(scMC)
library(Seurat)
library(patchwork)
library(dplyr)
setwd("/Users/suoqinjin/Documents/scMC")
# import data
load("/Users/suoqinjin/Documents/scMC/tutorial/data_dermis.rda")
data.input <- data_dermis
sample.name <- names(data.input)

######### Part I: Setup the a list of Seurat objects, one per dataset ##############
future::plan("multiprocess", workers = 4)
object.list <- vector("list", length(sample.name))
names(object.list) <- sample.name
for (i in 1:length(object.list)) {
  # Initialize the Seurat object with the raw (non-normalized data) for each dataset
  object.list0 <- CreateSeuratObject(counts = data.input[[i]], min.cells = 3, min.features = 200)
  # calculate mitochondrial QC metrics
  object.list0[["percent.mt"]] <- PercentageFeatureSet(object.list0, pattern = "^mt-")
  object.list0$sample.name <- sample.name[i]
  object.list[[i]] <-  object.list0
  rm(object.list0)
}

# step1. pre-processing each dataset
##  QC and selecting cells for further analysis
for (i in 1:length(object.list)) {
  # Visualize QC metrics as a violin plot
  gg <- VlnPlot(object.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001,cols=c("#a6cee3"))
  print(gg)
  Sys.sleep(2)
  # plot1 <- FeatureScatter(object.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1,cols=c("black"))
  # plot2 <- FeatureScatter(object.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1,cols=c("black"))
  # plot1 + plot2
  object.list[[i]] <- subset(object.list[[i]], subset = nFeature_RNA < 7000 & nCount_RNA < 40000 & percent.mt < 10)
}

## Normalizing, scaling the data and feature selection
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # perform scaling on the previously identified variable features
  x <- ScaleData(x, verbose = FALSE)
})

########### Part II: Perform an integrated analysis using scMC ###########
integrated.scMC <- RunscMC(object.list)

########### Part III: Run the standard workflow for visualization and clustering ###########
nPC = 40
integrated.scMC <- FindNeighbors(integrated.scMC, reduction = "scMC", dims = 1:nPC)
integrated.scMC <- FindClusters(integrated.scMC, algorithm = 4, resolution = 0.05)
levels(Idents(integrated.scMC))
integrated.scMC <- BuildClusterTree(integrated.scMC, reorder = T, reorder.numeric = T, verbose = F)
integrated.scMC <- RunUMAP(integrated.scMC, reduction='scMC', dims = 1:nPC)

integrated.scMC <- ScaleData(integrated.scMC, feature = rownames(object), verbose = FALSE)
markers <- FindAllMarkers(integrated.scMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(integrated.scMC, features = top10$gene) + NoLegend()

## Visualization
DimPlot(integrated.scMC, reduction = "umap", label = F)
p1 <- DimPlot(integrated.scMC, reduction = "umap", group.by = "sample.name")
p2 <- DimPlot(integrated.scMC, reduction = "umap", label = F)
patchwork::wrap_plots(p1, p2, ncol = 2)
DimPlot(integrated.scMC, group.by = c("sample.name", "ident"), ncol = 2)

features <- c('Pdgfra','Lox','Ptch1','Gli1')
FeaturePlot(integrated.scMC, features)

save.image(integrated.scMC, file = "scMC_dermis.RData")


