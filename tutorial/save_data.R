rm(list = ls())
setwd("/Users/suoqinjin/Downloads/simulation_data/data2")
library(data.table)
# X1 <- read.table("Kidney_input_RNA.txt",sep = '\t',row.names=1,header=T)
X1 <- fread("X1.txt")
X1 <- as.matrix(X1, rownames = TRUE)
X2 <- fread("X2.txt")
X2 <- as.matrix(X2, rownames = TRUE)
X1 <- as(X1, "dgCMatrix")
X2 <- as(X2, "dgCMatrix")
data <- list("dataset1" = X1, "dataset2" = X2)

labels1 = read.table(file = "label1.txt", header = F, sep = "\t")
row.names(labels1) <- colnames(X1); colnames(labels1) <- "labels"
labels2 = read.table(file = "label2.txt", header = F, sep = "\t")
row.names(labels2) <- colnames(X2); colnames(labels2) <- "labels"
labels <- list("dataset1" = labels1, "dataset2" = labels2)

data_simulation <- list("data" = data, "labels" = labels)

setwd("/Users/suoqinjin/Documents/scMC/data")
save(data_simulation, file = "data_simulation.rda")

## Dermis
setwd("/Users/suoqinjin/Downloads")
load("GSE112671_dermis_control_rawData.RData")
data1 <- w10x.data
X1 <- as(as.matrix(data1), "dgCMatrix")
load("GSE112671_dermis_mutant_rawData.RData")
data2 <- w10x.data
X2 <- as(as.matrix(data2), "dgCMatrix")
data_dermis <- list("Control" = X1, "Hh activation" = X2)
setwd("/Users/suoqinjin/Documents/scMC/data")
save(data_dermis, file = "data_dermis.rda")


