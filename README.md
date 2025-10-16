# scMC: Integrating and comparing multiple single cell genomic datasets

## Capabilities
- scMC is an R toolkit for integrating and comparing multiple single cell genomic datasets from single cell RNA-seq and ATAC-seq experiments across different conditions, time points and tissues. 
- scMC exhibits superior performance in detecting context-shared and -specific biological signals, particularly noticeable for the datasets with imbalanced cell population compositions across interrelated biological conditions. 
- scMC learns a shared reduced dimensional embedding of cells that retains the biological variation while removing the technical variation. This shared embedding can enhance a variety of single cell analysis tasks, such as low-dimensional visualization, cell clustering and pseudotemporal trajectory inference. 

## Installation
scMC R package can be easily installed from Github using devtools:  

```
devtools::install_github("jinworks/scMC_SeuratWrapper")
```
### Installation of other dependencies
- Install Leiden python pacakge for identifying cell clusters: ```pip install leidenalg```. Please check [here](https://github.com/vtraag/leidenalg) if you encounter any issue.


## Tutorials
Please check the tutorial directory of the repo.

- [R Markdown: Demo of scMC Seurat Wrapper (in-house, SJ2020)](https://htmlpreview.github.io/?https://github.com/sqjin/scMC_SeuratWrapper/blob/main/tutorial/demo_scMC_SeuratWrapper_SJ2020.html)

- [R code: Demo of scMC Seurat Wrapper (in-house, SJ2020)](https://github.com/sqjin/scMC_SeuratWrapper/blob/main/tutorial/demo_scMC_SeuratWrapper_SJ2020.R)


## Other great tools for single-cell data processing and downstream analysis
### R package
1. Seurat (https://satijalab.org/seurat/articles/get_started_v5_new)
2. SCP (https://github.com/zhanghao-njmu/SCP) ( Very nice visualization; Pipelines embedded with multiple integration methods for scRNA-seq or scATAC-seq data, including Uncorrected, Seurat, scVI, MNN, fastMNN, Harmony, Scanorama, BBKNN, CSS, LIGER, Conos, ComBat.; Multiple single-cell downstream analyses such as identification of differential features, enrichment analysis, GSEA analysis, identification of dynamic features, PAGA, RNA velocity, Palantir, Monocle2, Monocle3, etc.)


