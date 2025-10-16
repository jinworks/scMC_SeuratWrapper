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
for (i in 1:length(object.list)) {
  # Initialize the Seurat object with the raw (non-normalized data) for each dataset
  x <- CreateSeuratObject(counts = data.input[[i]], min.cells = 3, min.features = 200, project = sample.name[i])
  # calculate mitochondrial QC metrics
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  x$sample.name <- sample.name[i]
  x <- RenameCells(x, new.names = paste0(Cells(x), "_", x$sample.name))
  object.list[[i]] <-  x
  rm(x)
}
lapply(object.list, function(x) dim(x@assays$RNA@data))

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
  x <- FindVariableFeatures(x, verbose = FALSE, nfeatures = 3000)
  #x <- FindVariableFeatures(x, selection.method = "mean.var.plot", verbose = FALSE, mean.cutoff = c(0.01, 5), dispersion.cutoff = c(0.25, Inf))
  # perform scaling on the previously identified variable features
  x <- ScaleData(x, verbose = FALSE)
})


########### Part II: Perform an integrated analysis using scMC ###########
future::plan("multiprocess", workers = 4)
options(future.rng.onMisuse="ignore")

# step2. identify clusters with different resolution for each condition
# compute SNN
object.list <- identifyNeighbors(object.list)
# identify clusters
object.list <- identifyClusters(object.list)

## step3. detect cluster-specific cells with high confident
features.integration = identifyIntegrationFeatures(object.list)
object.list <- identifyConfidentCells(object.list, features.integration)

## step4. Identify marker genes associated with the putative cell clusters in each dataset
object.list <- identifyMarkers(object.list)

## step 5. Learn technical variation between any two datasets
structured_mat <- learnTechnicalVariation(object.list, features.integration)

## step 6. Learn a shared embedding of cells across all datasets after removing technical variation
combined <- merge(x = object.list[[1]],y = object.list[2:length(x = object.list)])
combined@meta.data <- combined@meta.data %>% select(-starts_with("RNA_snn_res"))
combined$sample.name <- factor(combined$sample.name, levels = sample.name)
VariableFeatures(combined) <- features.integration
combined <- integrateData(combined, structured_mat)

########### Part III: Run the standard workflow for visualization and clustering ###########
nPC = 40
combined <- FindNeighbors(combined, reduction = "scMC", dims = 1:nPC)
combined <- FindClusters(combined, algorithm = 4, resolution = 0.05)
levels(Idents(combined))
combined <- BuildClusterTree(combined, reorder = T, reorder.numeric = T, verbose = F)
combined <- RunUMAP(combined, reduction='scMC', dims = 1:nPC)

combined <- ScaleData(combined, feature = rownames(object), verbose = FALSE)
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(combined, features = top10$gene) + NoLegend()

## Visualization
DimPlot(combined, reduction = "umap", label = F)
DimPlot(combined, group.by = c("sample.name", "ident"), ncol = 2)

features <- c('Pdgfra','Lox','Ptch1','Gli1','Ptprc')
FeaturePlot(combined, features)

save.image(combined, file = "scMC_dermis.RData")


