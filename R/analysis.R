#' Stacked bar plot showing the proportion of cells across certrain cell groups
#'
#' @param object seurat object
#' @param x Name of one metadata column to show on the x-axis
#' @param fill Name of one metadata column to compare the proportion of cells
#' @param facet Name of one metadata column defining faceting groups
#' @param colors.use defining the color of stacked bar plot; either a char vector defining a color for each cell group or a palette name from brewer.pal
#' @param n.colors Number of colors when setting colors.use to be a palette name from brewer.pal
#' @param n.row Number of rows in facet_grid()
#' @param title.name Name of the main title
#' @param legend.title Name of legend
#' @param xlabel Name of x label
#' @param ylabel Name of y label
#' @param width bar width
#' @param show.legend Whether show the legend
#' @param x.lab.rot Whether rorate the xtick labels
#' @param text.size font size
#' @param flip Whether flip the cartesian coordinates so that horizontal becomes vertical
#' @return ggplot2 object
#' @export
#'
#' @import ggplot2
#' @importFrom plyr ddply as.quoted
computeProportion <- function(object, x = "cellType", fill = "condition", facet = NULL, colors.use = NULL, n.colors = 8, n.row = 1, title.name = NULL, legend.title = NULL,
                              xlabel = NULL, ylabel = "Cellular composition (%)", width = 0.8,
                              show.legend = TRUE, x.lab.rot = TRUE, text.size = 10, flip = FALSE) {

  df <- plyr::ddply(object@meta.data, plyr::as.quoted(c(x,fill)), nrow)

  gg <- ggplot(df, aes_string(x = x, y = "V1", fill = fill)) +
    geom_bar(stat="identity", position="fill", width = width) +
    scale_y_continuous(name = ylabel, labels = c(0,25,50,75,100))
  if (!is.null(facet)) {
    gg <- gg + facet_wrap(facet, nrow = n.row)
  }

  gg <- gg + theme_classic() + ylab(ylabel) + xlab(xlabel) +
    labs(title = title.name) +  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = text.size), axis.text = element_text(colour="black"))

  if (!is.null(colors.use)) {
    if (length(colors.use) == 1) {
      colors.use <- RColorBrewer::brewer.pal(n.colors, colors.use)[1:length(unique(df[, fill]))]
    }
    gg <- gg + scale_fill_manual(values = alpha(colors.use, alpha = 1), drop = FALSE)
    #   gg <- gg + scale_color_manual(values = alpha(color.use, alpha = 1), drop = FALSE) + guides(colour = FALSE)
  }
  if (is.null(legend.title)) {
    gg <- gg + theme(legend.title = element_blank())
  } else {
    gg <- gg + guides(fill=guide_legend(legend.title))
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    gg <- gg + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=text.size))
  }
  if (flip) {
    gg <- gg + coord_flip()
  }
  gg
  return(gg)
}

#' Identify differential expressed genes
#'
#' @param object Seurat object
#' @param ident.1 Identity class to define markers for. If NULL (default) - use FindAllMarkers.
#' @param ident.2 A second identity class for comparison. If NULL (default) - use all other cells for comparison.
#' @param group.by A vector of variables to group cells by; pass NULL to group by cell identity classes
#' @param grouping.var grouping variable when finding the markers that are conserved between the groups, e.g, grouping.var = "batch"
#' @param test.use which test to use
#' @param only.pos Only return positive markers
#' @param min.pct Threshold of the percent of cells enriched in one cluster
#' @param logfc.threshold Threshold of Log Fold Change
#' @param do.plot whether plot the heatmap
#' @param gg.return whether return the ggplot object
#' @param file.name if not NULL, it will automatically save a PDF figure of the heatmap
#' @param width,height set the figure size when saving figures
#' @param n.top the number of genes per cell group shown in the heatmap
#' @param additional.group.by A vector of variables to group cells by
#' @param colors.use A list of colors to use for the color bar. e.g., colors.use = list(ident = c("red","blue"))
#' @param ... other parameters in FindAllMarkers, FindMarkers, FindConservedMarkers or doHeatmap
#'
#' @import Seurat
#' @import dplyr
#' @export
#'
identifyDEG <- function(object, ident.1 = NULL, ident.2 = NULL, group.by = NULL, grouping.var = NULL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox",
                        do.plot = TRUE, gg.return = FALSE, file.name = "heatmap_DEG", width = NA, height = NA,
                        n.top = 10, additional.group.by = NULL, colors.use = NULL, ...) {
  if (is.null(group.by)) {
    message("The identity classes (i.e., active.ident) is used as default for testing")
  } else if (group.by %in% colnames(object@meta.data)){
    Idents(object) <- group.by
  } else {
    stop("Please provide a valid char name of one metadata column!")
  }

  if (is.null(grouping.var)) {
    if (is.null(ident.1)) {
      message("Identify DEG using FindAllMarkers")
      markers <- tryCatch({
        FindAllMarkers(object, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use, ...)
      }, error = function(e) {
        FindAllMarkers(object, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use)
      })
    } else {
      message("Identify DEG using FindMarkers")
      markers <- tryCatch({
        FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use, ...)
      }, error = function(e) {
        FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use)
      })
    }
  } else {
    message("Identify DEG using FindConservedMarkers")
    markers <- tryCatch({
      FindConservedMarkers(object, ident.1 = ident.1, ident.2 = ident.2, grouping.var = grouping.var, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use, ...)
    }, error = function(e) {
      FindConservedMarkers(object, ident.1 = ident.1, ident.2 = ident.2, grouping.var = grouping.var, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use)
    })
  }
  if (do.plot) {
    top.markers <- markers %>% group_by(cluster) %>% top_n(n = n.top, wt = avg_log2FC)
    gg <- tryCatch({
      doHeatmap(object, features = top.markers$gene, group.by = group.by, additional.group.by = additional.group.by, colors.use = colors.use, ...)
    }, error = function(e) {
      doHeatmap(object, features = top.markers$gene, group.by = group.by, additional.group.by = additional.group.by, colors.use = colors.use)
    })
    if (!is.null(file.name)) {
      ggsave(filename = paste0(file.name,"_", group.by, ".pdf"), plot=gg,  units = 'in', height = height, width = width, device = "pdf")
    }
    print(gg)
  }
  if (gg.return) {
    return(markers = markers, gg.obj = gg)
  } else {
    return(markers)
  }

}



#' Perform GO or KEGG enrichment analysis
#'
#' @param df A data frame containing one column named 'cluster' and one column named 'gene'
#' @param OrgDb Genome wide annotation for organism when performing GO analysis
#' @param category Type of GO or KEGG enrichment analysis
#' @param organism organism when performing KEGG analysis
#' @param GO.old Use previously calculated enrichResult object to directly plot the results
#' @param idents A vector of cell clusters from df$cluster to perform analysis; default will perform GO comparison analsyis for all the cell clusters
#' @param invert invert the cell clusters defined by idents
#' @param color.heatmap A character string or vector indicating the colormap option to use. It can be the avaibale color palette in viridis_pal() or brewer.pal()
#' @param direction Sets the order of colors in the scale. If 1, the default colors are used. If -1, the order of colors is reversed.
#' @param n.colors number of basic colors to generate from color palette
#' @param dot.size the dot size
#' @param title.name title name
#' @param legend.remove whether remove the figure legend
#' @param file.name if not NULL, it will automatically save a PDF figure and .csv file of the results
#' @param width,height set the figure size when saving figures
#' @param gg.return whether return the ggplot object
#' @param showCategory Number of GO terms to keep
#' @param pvalueCutoff pvalue cutoff on enrichment tests to report
#' @param qvalueCutoff value cutoff on enrichment tests to report as significant.
#' @param do.simplify whether simply the GO terms by removing redundant terms.
#' @param simplifyCutoff cutoff to use for simplifying the GO terms
#' @param label_format a numeric value sets wrap length, alternatively a custom function to format axis labels
#' @param p_val cutoff associated with each gene when performing the DEG analysis
#' @param text.size Font size of the figure
#'
#' @import clusterProfiler
#' @export
#'
getEnrichedGO <- function(df, OrgDb = c("org.Mm.eg.db", "org.Hs.eg.db"), category = c("BP", "MF", "KEGG"), organism = c("mouse", "hsa"), GO.old = NULL, idents = NULL, invert = FALSE,
                          color.heatmap = "Reds", n.colors = 8, direction = 1, dot.size = 3,
                          title.name = " enrichment analysis", legend.remove = FALSE,
                          file.name = "GO_comparison_DEG_clusters", width = NA, height = NA, gg.return = FALSE,
                          showCategory = 10, pvalueCutoff  = 0.01, qvalueCutoff = 0.05,
                          do.simplify = TRUE, simplifyCutoff = 0.7, label_format = 1000,
                          p_val = 0.01, text.size = 10) {
  category <- match.arg(category)
  OrgDb <- match.arg(OrgDb)
  organism <- match.arg(organism)
  if (!is.null(title.name)) {
    title.name <- paste0(category,title.name)
  }
  if (!("p_val" %in% colnames(df))) {
    warning("p_val is not a colname in df! Set it to be zero!!")
    df$p_val <- 0
  }
  if (is.null(idents) | length(idents) >= 2) {
    clusters <- levels(df$cluster)
    if (!is.null(idents)) {
      if (invert) {
        clusters <- clusters[!(clusters %in% idents)]
      } else {
        clusters <- clusters[clusters %in% idents]
      }
    }
    # Create a list with genes from each cell group
    
    if (is.null(GO.old)) {
      gene.sets <- list()
      for (i in 1:length(clusters)) {
        name <- clusters[i]
        genes <- df$gene[df$cluster == clusters[i] & df$p_val < p_val]
        if (category == "KEGG") {
          ID <- bitr(as.character(genes), fromType="SYMBOL", toType="ENTREZID", OrgDb = OrgDb)
          gene.sets[[name]] <- ID$ENTREZID
        } else {
          gene.sets[[name]] <- as.character(genes)
        }
        
      }
      if (category == "KEGG") {
        compGO <- compareCluster(gene.sets, fun = "enrichKEGG", keyType = "ncbi-geneid", organism = organism, pvalueCutoff  = pvalueCutoff, qvalueCutoff = qvalueCutoff)
      } else {
        compGO <- compareCluster(gene.sets, fun = "enrichGO", keyType = "SYMBOL", OrgDb = OrgDb, ont = category, pvalueCutoff  = pvalueCutoff, qvalueCutoff = qvalueCutoff)
        if (do.simplify) {
          compGO <- simplify(compGO, cutoff = simplifyCutoff, by="p.adjust", select_fun=min)
        }
      }
      
    } else {
      compGO <- GO.old
    }
  } else {
    if (is.null(GO.old)) {
      genes <- as.character(df$gene[df$cluster == idents & df$p_val < p_val])
      gene.sets <- list()
      if (category == "KEGG") {
        ID <- bitr(as.character(genes), fromType="SYMBOL", toType="ENTREZID", OrgDb = OrgDb)
        gene.sets[[idents]] <- ID$ENTREZID
      } else {
        gene.sets[[idents]] <- as.character(genes)
      }
      if (category == "KEGG") {
        compGO <- enrichKEGG(gene.sets[[idents]], keyType = "ncbi-geneid", organism = organism, pvalueCutoff  = pvalueCutoff, qvalueCutoff = qvalueCutoff)
      } else {
        # compGO <- enrichGO(gene.sets[[idents]], OrgDb = OrgDb,  keyType = "ENTREZID", ont = category, pvalueCutoff  = pvalueCutoff, qvalueCutoff = qvalueCutoff)
        compGO <- enrichGO(gene.sets[[idents]], OrgDb = OrgDb,  keyType = "SYMBOL", ont = category, pvalueCutoff  = pvalueCutoff, qvalueCutoff = qvalueCutoff)
      }
      
    } else {
      compGO <- GO.old
    }
  }
  
  gg <- dotplot(compGO, showCategory = showCategory, title = title.name, label_format = label_format) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1,size = rel(1))) +
    theme(axis.text.x = element_text(size = text.size), axis.text.y = element_text(size = text.size)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  if (is.null(idents) | length(idents) >= 2) {
    gg <- gg + scale_radius(range = c(1, dot.size),  name = "Gene ratio")
  } else {
    gg <- gg + scale_radius(range = c(1, dot.size),  name = "Gene count")
  }
  
  
  if (!is.null(color.heatmap)) {
    if (length(color.heatmap) == 1) {
      color.heatmap.use <- tryCatch({
        RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
      }, error = function(e) {
        scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
      })
      color.heatmap.use <- color.heatmap.use[3:(length(color.heatmap.use)-1)]
    } else if (length(color.heatmap) > 1) {
      color.heatmap.use <- color.heatmap
    }
    if (direction == -1) {
      color.heatmap.use <- rev(color.heatmap.use)
    }
    color.heatmap.use <- colorRampPalette(color.heatmap.use)(99)
    gg <- gg + scale_colour_gradientn(colours = color.heatmap.use, trans = 'reverse',
                                      guide = guide_colorbar(title = "Adjusted p-value", ticks = T, label = T, barwidth = 0.5), na.value = "white")
  }
  if (legend.remove) {
    gg <- gg + theme(legend.position = "none")
  }
  print(gg)
  if (!is.null(file.name)) {
    ggsave(filename = paste0(file.name, ".pdf"), plot=gg,  units = 'in', height = height, width = width, device = "pdf")
    # Output results from GO analysis to a table
    GO_summary <- data.frame(compGO)
    utils::write.csv(GO_summary, file = paste0(file.name, ".csv"))
  }
  if (gg.return) {
    return(list(res.GO = compGO, gg.obj = gg))
  } else {
    return(compGO)
  }
}


#' Classify New Data
#'
#' Classify new data based on the cluster information of the provided object.
#' Random Forests are used as the basis of the classification.
#'
#' @param object Seurat object on which to train the classifier
#' @param classifier Random Forest classifier from BuildRFClassifier. If not provided,
#' it will be built from the training data provided.
#' @param training.genes Vector of genes to build the classifier on
#' @param training.classes Vector of classes to build the classifier on
#' @param new.data New data to classify
#' @param assay Specific assay to get data from the training data; default = "RNA"
#' @param ... additional parameters passed to ranger
#'
#' @return Vector of cluster ids
#'
#' @import Matrix
#' @importFrom stats predict
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # take the first 10 cells as test data and train on the remaining 70 cells
#' test.pbmc <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[1:10])
#' train.pbmc <- SubsetData(object = pbmc_small, cells.use = pbmc_small@cell.names[11:80])
#' predicted.classes <- ClassifyCells(
#'   object = train.pbmc,
#'   training.classes = Idents(train.pbmc),
#'   new.data = test.pbmc[["RNA"]]@data
#' )
#'
classifyCells <- function(
  object,
  classifier,
  training.genes = NULL,
  training.classes = NULL,
  new.data = NULL,
  assay = "RNA",
  ...
) {
  PackageCheck('ranger')
  # build the classifier
  if (missing(classifier)){
    classifier <- buildRFClassifier(
      object = object,
      training.genes = training.genes,
      training.classes = training.classes,
      assay = "RNA",
      ...
    )
  }
  # run the classifier on the new data
  features <- classifier$forest$independent.variable.names
  genes.to.add <- setdiff(x = features, y = rownames(x = new.data))
  data.to.add <- matrix(
    data = 0,
    nrow = length(x = genes.to.add),
    ncol = ncol(x = new.data)
  )
  rownames(x = data.to.add) <- genes.to.add
  new.data <- rbind(new.data, data.to.add)
  new.data <- new.data[features, ]
  new.data <- as.matrix(x = t(x = new.data))
  message("Running Classifier ...")
  prediction <- predict(classifier, new.data)
  new.classes <- prediction$predictions
  return(new.classes)
}

#' Build Random Forest Classifier
#'
#' Train the random forest classifier
#'
#'
#' @param object Seurat object on which to train the classifier
#' @param training.genes Vector of genes to build the classifier on
#' @param training.classes Vector of classes to build the classifier on
#' @param assay Specific assay to get data from; default = "RNA"
#' @param verbose Additional progress print statements
#' @param ... additional parameters passed to ranger
#'
#' @return Returns the random forest classifier
#'
#' @import Matrix
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Builds the random forest classifier to be used with ClassifyCells
#' # Useful if you want to use the same classifier with several sets of new data
#' classifier <- BuildRFClassifier(pbmc_small, training.classes = pbmc_small@ident)
#'
buildRFClassifier <- function(
  object,
  training.genes = NULL,
  training.classes = NULL,
  assay = "RNA",
  verbose = TRUE,
  ...
) {
  PackageCheck('ranger')
  training.classes <- as.vector(x = training.classes)
  data.training = GetAssayData(object, slot = "data", assay = assay)
  training.genes <- SetIfNull(
    x = training.genes,
    default = rownames(x = data.training)
  )
  training.data <- as.data.frame(
    x = as.matrix(
      x = t(
        x = data.training[training.genes, ]
      )
    )
  )
  training.data$class <- factor(x = training.classes)
  if (verbose) {
    message("Training Classifier ...")
  }
  classifier <- ranger::ranger(
    data = training.data,
    dependent.variable.name = "class",
    classification = TRUE,
    write.forest = TRUE,
    ...
  )
  return(classifier)
}

# gg <- EnhancedVolcano::EnhancedVolcano(markers.sample,
#                                        lab = markers.sample$gene,
#                                        selectLab = NULL,
#                                        x = 'avg_log2FC',
#                                        y = 'p_val_adj',
#                                        pCutoff = 10e-10,
#                                        FCcutoff = 1,
#                                        labSize = 3,
#                                        title = 'DEG analysis',
#                                        col=c('grey50', 'grey50', 'grey50', '#e41a1c'),
#                                        colAlpha = 0.5,
#                                        subtitle = NULL,
#                                        legend = NULL,
#                                        axisLabSize = 10,
#                                        titleLabSize = 12,
#                                        captionLabSize = 10,
#                                        boxedLabels = F,drawConnectors = F) + theme(legend.position = "none")


#' Convert ExpressionSet object or data matrix to Seutat object
#'
#' @param eset an ExpressionSet object
#' @param data.mat a data matrix
#' @param gene.names a char vector giving the gene names
#' @param labels a char vector giving the cell labels of each cell
#' @param min.cells Include features detected in at least this many cells.
#'
#' @return
#' @export
ExpressionSet2seurat <- function(eset = NULL, data.mat = NULL, gene.names = NULL, labels = NULL, min.cells = 10) {
 if (!is.null(eset)) {
   mat <- exprs(eset)
   feature <- fData(eset)
   meta <- pData(eset)
 } else if (!is.null(data.mat)) {
   mat <- data.mat
   feature <- data.frame(gene_short_name = gene.names, gene_id = rownames(data.mat))
   meta <- data.frame(labels = labels, row.names = colnames(data.mat))
 }

  idx <- which(rowSums(mat > 0) <= min.cells)
  mat <- mat[-idx, ]
  feature <- feature[-idx,]
  
  identical(feature$gene_id, rownames(mat))
  # length(unique(feature$gene_short_name))
  gene.rep = feature$gene_short_name[duplicated(feature$gene_short_name)]
  
  if (length(gene.rep) > 0) {
    mat1 <- mat[feature$gene_id[feature$gene_short_name %in% gene.rep == FALSE], ]
    rownames(mat1) <- feature$gene_short_name[match(rownames(mat1), feature$gene_id)]
    mat2 <- c()
    for (i in 1:length(gene.rep)) {
      mat2 <- rbind(mat2, Matrix::colMeans(mat[feature$gene_id[feature$gene_short_name == gene.rep[i]], ]))
    }
    rownames(mat2) <- gene.rep
    mat <- rbind(mat1, mat2)
    rownames(mat) <- c(rownames(mat1), rownames(mat2))
  }
  
  seu = CreateSeuratObject(mat, meta.data = meta, min.cells = min.cells, min.features = 0)
  print(seu)
  return(seu)
}



