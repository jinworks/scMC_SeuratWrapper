#
rm(list=ls())
library(scMC)
setwd("/Users/suoqinjin/Documents/scMC")
# import data
load("/Users/suoqinjin/Documents/scMC/tutorial/data_dermis.rda")

scMC <- create_scMC(raw.data = data_dermis)

## step1. process data
scMC <- preprocessing(scMC, add.names = c("control","mutant"))

## step2. identify clusters with different resolution for each condition
# compute SNN
scMC <- identifyNeighbors(scMC, mode = "separate")
# identify clusters
scMC <- identifyClusters(scMC, mode = "separate")

## step3. detect cluster-specific cells with high confident
scMC <- identifyConfidentCells(scMC)

## step4. Identify marker genes associated with the putative cell clusters in each dataset
scMC <- identifyMarkers(scMC, mode = "separate")

## step 5. Learn technical variation between any two datasets
structured_mat <- learnTechnicalVariation(scMC)

## step 7. Learn a shared embedding of cells across all datasets after removing technical variation
scMC <- integrateData(scMC, structured_mat)

## step 8. Perform an integrated analysis
# Run the standard workflow for visualization and clustering
scMC <- runUMAP(scMC, reducedDims = "integrated")
scMC <- identifyNeighbors(scMC, mode = "integrated")
scMC <- identifyClusters(scMC, mode = "integrated", resolution = 0.04)

## step 9. Visualization
gg1 <- cellVisualization(scMC, color.by = "sample.id", colors.ggplot = T)
gg2 <- cellVisualization(scMC, color.by = "Integrated clusters")
cowplot::plot_grid(gg1, gg2)

save.image(scMC, file = "scMC_dermis.RData")
