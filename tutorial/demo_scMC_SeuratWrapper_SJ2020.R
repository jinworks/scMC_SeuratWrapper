#
rm(list=ls())
reticulate::use_python("/Users/suoqinjin/opt/miniconda3/bin/python", required=T)
library(scMC)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
setwd("/Users/suoqinjin/Documents/scMC_SeuratWrapper")

# import data
load("/Users/suoqinjin/Documents/scMC_SeuratWrapper/tutorial/data_dermis.rda")
data.input <- data_dermis
sample.name <- names(data.input)

######### Part I: Setup the a list of Seurat objects, one per dataset ##############
object.list <- vector("list", length(sample.name))
names(object.list) <- sample.name
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
nFeature_RNA1 = c(7000, 7000); nCount_RNA1 = c(40000, 40000); percent.mt1 = c(10, 10)
for (i in 1:length(object.list)) {
  # Visualize QC metrics as a violin plot
  gg <- VlnPlot(object.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001,cols=c("#a6cee3")) + geom_hline(yintercept = percent.mt1[i], linetype = 2)
  print(gg)
  cowplot::save_plot(filename=paste0("QC1_", sample.name[i],"_.pdf"), plot=gg, base_width = 7, base_height = 4.5)
  # VlnPlot(object.list[[i]], features = c("percent.mt"))+ geom_hline(yintercept = 15, linetype = 2)
  Sys.sleep(2)
  plot1 <- FeatureScatter(object.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1,cols=c("black")) + geom_hline(yintercept = percent.mt1[i], linetype = 2)
  plot2 <- FeatureScatter(object.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1,cols=c("black")) + geom_hline(yintercept = nFeature_RNA1[i], linetype = 2)
  gg <- wrap_plots(plots = plot1, plot2, ncol = 2)
  print(gg)
  cowplot::save_plot(filename=paste0("QC2_", sample.name[i],"_.pdf"), plot=gg, base_width = 10, base_height = 3.5)
  object.list[[i]] <- subset(object.list[[i]], subset = (nFeature_RNA < nFeature_RNA1[i]) & (nCount_RNA < nCount_RNA1[i]) & (percent.mt < percent.mt1[i]))
}
lapply(object.list, function(x) dim(x@assays$RNA@data))

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
dimPlot(combined, reduction = "umap", group.by = c("sample.name","ident"))

combined <- ScaleData(combined, feature = rownames(combined), verbose = FALSE)
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Misc(combined, slot = 'markers') <- markers
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf('dermis_heatmap112.pdf')
#DoHeatmap(combined, features = top10$gene) + scale_fill_gradientn(colors = color.heatmap.use) + NoLegend()
#doHeatmap(combined, features = top10$gene)
doHeatmap(combined, features=top10$gene, group.by='ident', additional.group.by = c('sample.name'))+
  theme(axis.text.y = element_text(size = 6))
dev.off()

## Visualization and annotation
### Feature plot of known marker genes
features = c('Lox','Ptch1','Cd68','Pecam1','Myh11','Plp1')
gg <- featurePlot(combined, features = features, show.legend = F, show.axes = F)
gg <- patchwork::wrap_plots(plots = gg, ncol = 3)
gg
cowplot::save_plot(filename=paste0("overlayKnownMarkers_integration_dermis", "_umap.pdf"), plot=gg, base_width = 5.5, base_height = 6)


## Annotation
## Labeling the clusters by cell type ##
new.cluster.ids <- c("Immune", "Hh-inactive Fib", "Hh-active Fib",  "Schwann", "Endotheial","Muscle")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
new.order <- c("Hh-inactive Fib", "Hh-active Fib","Immune","Endotheial", "Muscle", "Schwann")
combined@active.ident <- factor(combined@active.ident, levels = new.order)


## Visualize cells onto the low-dimensional space
p1 <- dimPlot(combined, reduction = "umap", group.by = "sample.name", colors.ggplot = T)
p2 <- dimPlot(combined, reduction = "umap", label = F)
gg <- patchwork::wrap_plots(p1, p2, ncol = 2)
gg
cowplot::save_plot(filename=paste0("integration_dermis", "_umap.pdf"), plot=gg, base_width = 8, base_height = 3)

# Split the plot into each dataset
dimPlot(combined, reduction = "umap", split.by = "sample.name",  combine = T)


### Violin plot
#### Stacked violin plot
features = c('Lox','Ptch1','Cd68','Pecam1','Myh11','Plp1')
gg <- StackedVlnPlot(combined, features = features)
gg
cowplot::save_plot(paste0("violin_markers", "_integration_dermis", ".pdf"), gg, base_height = 3.5, base_width = 2)

#### Splitted violin plot
features = c('Lox','Ptch1','Cd68','Pecam1','Myh11','Plp1')
gg <- StackedVlnPlot(combined, features = features, split.by = "sample.name")
gg

#### Violin plot with statistical test
features = c('Pdgfra','Lox','Ptch1','Gli1')
gg <- vlnPlot(combined, features = features, stat.add = T, comparisons = list(c("Hh-inactive Fib", "Hh-active Fib")))
patchwork::wrap_plots(plots = gg, ncol = 4)

### Dot plot
dotPlot(combined, features =c('Lox','Ptch1','Cd68','Pecam1','Myh11','Plp1'))

#### Splitted Dot plot
dotPlot(combined, features = c('Pdgfra','Lox','Ptch1','Gli1'), split.by = "sample.name", idents = c("Hh-inactive Fib", "Hh-active Fib"))


# Part IV: Downstream analysis
combined$clusters.final <- Idents(combined)

## Compute the proportion
gg1 <- computeProportion(combined, x = "clusters.final", fill = "sample.name")
gg2 <- computeProportion(combined, x = "sample.name", fill = "clusters.final", colors.use = scPalette(6))
gg1 + gg2

## DEG analysis
markers.sample <- identifyDEG(combined, group.by = "sample.name")

## GO entichment analysis
res.go <- getEnrichedGO(markers, category = "BP", do.simplify = F, idents = c("Hh-inactive Fib", "Hh-active Fib"))


# save data
combined <- ScaleData(combined)
combined$clusters.final <- Idents(combined)
save(combined, file = "scMC_dermis_SeuratV4.RData")
# combined.loom <- as.loom(combined, filename = "scMC_dermis_SeuratV4.loom", verbose = FALSE)
# combined.loom$close_all()
write.table(GetAssayData(combined),file = "preprocessedData_integration_dermis.txt",sep = '\t')
write.table(combined@meta.data,file = "metaData_integration_dermis.txt",sep = '\t')
write.table(Idents(combined),file = "identity_integration_dermis.txt",sep = '\t')
write.table(Embeddings(combined, "scMC"),file = "integratedSpace_scMC_integration.txt",sep = '\t')
write.table(combined[["umap"]]@cell.embeddings,file = "projectedData_umap_integration_dermis.txt",sep = '\t')
write.table(markers,file = "markers_integration_dermis.txt",sep = '\t')
write.table(top10,file = "markersTop10_integration_dermis.txt",sep = '\t')

# combined.loom <- loomR::connect(filename = "scMC_dermis_SeuratV4.loom", mode = "r")
# combined <- as.Seurat(combined.loom)
load("scMC_dermis_SeuratV4.RData")
       
       
 ############# Downstream analysis ##################
## Compute the proportion
gg1 <- computeProportion(combined, x = "clusters.final", fill = "conditions", ylabel = "Cellular composition (%)")
gg2 <- computeProportion(combined, x = "conditions", fill = "clusters.final", colors.use = scPalette(nlevels(combined$clusters.final)),ylabel = "Cellular composition (%)")
gg <- patchwork::wrap_plots(gg1, gg2, ncol = 2,widths = c(3, 1))
gg
cowplot::save_plot(filename=paste0("percent_integration_iPSC_D12_LacZ&KD", "_scMC_K9.pdf"), plot=gg, base_width = 6, base_height = 4)


## DEG analysis
markers.conditions <- identifyDEG(combined, group.by = "conditions")
top10.conditions <- markers.conditions %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
pdf('heatmap_DEG__conditions_iPSC_D12_LacZ_vs_KD.pdf', width = 4, height = 5)
doHeatmap(combined, features=top10.conditions$gene, group.by='conditions')
dev.off()

unique(markers$cluster)
res.KEGG <- getEnrichedGO(markers.conditions, GO.old = res.KEGG, category = "KEGG", do.simplify = F, OrgDb= "org.Hs.eg.db", organism = "hsa",
                          file.name = "KEGG_comparison_DEG_conditions_iPSC_D12_LacZ_vs_KD", width = 6, height = 5)


# gene scoring analysis
library(tidyverse)
library(readxl)

s.genes <- Seurat::cc.genes$s.genes
s.genes <- s.genes[s.genes %in% rownames(combined)] # genes in dataset
g2m.genes <- Seurat::cc.genes$g2m.genes
g2m.genes <- g2m.genes[g2m.genes %in% rownames(combined)] # genes in dataset
combined <- CellCycleScoring(object = combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
combined$Phase <- factor(combined$Phase, levels = c("G1","S","G2M"))
colors.phase <- c("#bdbdbd", "#e7298a","#78c679")
gg <- dimPlot(combined, reduction = "umap", group.by = "Phase",split.by = "conditions", colors.use = colors.phase, combine = T)
cowplot::save_plot(filename=paste0("integration_iPSC_D12_LacZ&KD", "_umap_cellcycle_scMC.pdf"), plot=gg, base_width = 6, base_height = 3)


geneSets <- excel_sheets("/Users/suoqinjin/Google Drive/projects/project_Tim/scRNA-seq/signature_genes_iPSC_reprogramming.xlsx")
geneSets <- geneSets[15:18]
for (i in 1:length(geneSets)) {
  genes <- read_excel("/Users/suoqinjin/Google Drive/projects/project_Tim/scRNA-seq/signature_genes_iPSC_reprogramming.xlsx", sheet =  geneSets[i])
  genes.use <- toupper(genes[,1][[1]])
  genes.use <- genes.use[genes.use %in% rownames(combined)]
  combined <- AddModuleScore(combined, features = list(genes.use), name = geneSets[i],assay = "RNA")
}
names(combined@meta.data) <- gsub("1","",names(combined@meta.data))
names(combined@meta.data) <- gsub("signaling1","signaling",names(combined@meta.data))

gg <- featurePlot(combined, features = geneSets, show.legend = T, show.axes = F)
gg <- patchwork::wrap_plots(plots = gg, ncol = 4)
gg
cowplot::save_plot(filename=paste0("overlaySignalingScores_integration_iPSC_D12_LacZ&KD", "_umap_scMC.png"), plot=gg, base_width = 9, base_height = 2.5)

gg <- vlnPlot(combined, group.by = 'conditions',features = geneSets, show.median = T ,colors.ggplot = T )
gg <- patchwork::wrap_plots(plots = gg, ncol = 4)
cowplot::save_plot(filename=paste0("violinSignalingScores_iPSC_D12_LacZ_vs_KD", ".png"), plot=gg, base_width = 6, base_height = 3.5)





