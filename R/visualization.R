#' ggplot theme in scMC
#'
#' @export
#'
#' @examples
#' @importFrom ggplot2 theme_classic element_rect theme element_blank element_line element_text
scMC_theme_opts <- function() {
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme_classic() +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(color = "black")) +
    theme(axis.line.y = element_line(color = "black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank()) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}


#' Generate ggplot2 colors
#'
#' @param n number of colors to generate
#' @importFrom grDevices hcl
#' @export
#'
ggPalette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Generate colors from a customed color palette
#'
#' @param n number of colors
#'
#' @return A color palette for plotting
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}


#' visualize cells in 2D-dimensional space
#'
#' @param object seurat object
#' @param colors.use defining the color for each cell group
#' @param colors.ggplot whether use ggplot color scheme; default: colors.ggplot = FALSE
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' @param shape.by If NULL, all points are circles (default).
#' @param label whether label the clusters in the 2D-space
#' @param shuffle Whether to randomly shuffle the order of points. This can be useful for crowded plots if points of interest are being buried.
#' @param seed Sets the seed if randomly shuffling the order of points.
#' @param pt.size size of the dots
#' @param combine Combine plots into a single patchworked ggplot object. If FALSE, return a list of ggplot objects
#' @param ... Extra parameters passed to DimPlot from Seurat package
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom Seurat DimPlot
dimPlot <- function(object, colors.use = NULL, colors.ggplot = FALSE, reduction = NULL,
                    group.by = NULL, split.by = NULL, shape.by = NULL, label = FALSE, shuffle = TRUE, seed = 1,
                    pt.size = 0.0001, combine = TRUE,...) {
  if (is.null(colors.use)) {
    numCluster <- length(levels(Idents(object)))
    if (colors.ggplot) {
      colors.use <- NULL
    } else {
      colors.use <- scPalette(numCluster)
    }
  }
  if (!is.null(split.by)) {shuffle = FALSE}
  if (combine) {
    gg <- DimPlot(object, reduction = reduction, group.by = group.by, split.by = split.by, shape.by = shape.by, label = label,
                  shuffle = shuffle, seed = seed, pt.size = pt.size, cols = colors.use, combine = combine,...) +
      theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
      theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.line = element_line(colour = 'black')) +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  } else {
    gg <- DimPlot(object, reduction = reduction, group.by = group.by, split.by = split.by, shape.by = shape.by, label = label, shuffle = shuffle, seed = seed,
                  pt.size = pt.size, cols = colors.use, combine = combine, ...)
    for(i in 1:length(gg)) {
      gg[[i]] <- gg[[i]] +
        theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
        theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.line = element_line(colour = 'black')) +
        theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
      xlims <- ggplot_build(gg[[i]])$layout$panel_scales_x[[1]]$range$range
      xlims <- c(xlims[1]+max(xlims[1]*0.1, -1), xlims[2]+min(xlims[2]*0.1, 1))
      # xlims[1] <- max(xlims[1]*1.1, floor(xlims[1]))
      # xlims[2] <- min(xlims[2]*1.1, ceiling(xlims[2]))
      ylims <- ggplot_build(gg[[i]])$layout$panel_scales_y[[1]]$range$range
      ylims <- c(ylims[1]+max(ylims[1]*0.1, -1), ylims[2]+min(ylims[2]*0.1, 1))
      gg[[i]] <- gg[[i]] + xlim(xlims) + ylim(ylims)
    }
  }
  return(gg)
}


#' visualize 'features' of single cells in 2D-dimensional space
#'
#' @param object seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param colormap RColorbrewer palette to use (check available palette using RColorBrewer::display.brewer.all()). default will use customed color palette
#' @param n.colors Number of colors when setting a color name based on RColorBrewer::brewer.all()
#' @param color.direction Sets the order of colours in the scale. If 1, the default, colours are as output by RColorBrewer::brewer.pal(). If -1, the order of colours is reversed.
#' @param show.text.axis whther show axis text
#' @param show.legend whether show individual legend
#' @param show.legend.combined  whether just show one legend
#' @param show.axes whether show the axes
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' @param shape.by If NULL, all points are circles (default).
#' @param pt.size size of the dots
#' @param combine Combine plots into a single patchworked ggplot object. If FALSE, return a list of ggplot objects
#' @param ... Extra parameters passed to DimPlot from Seurat package
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom Seurat FeaturePlot
featurePlot <- function(object, features, colormap = NULL, n.colors = 9, color.direction = 1,show.text.axis = FALSE, show.legend = T, show.legend.combined = F, show.axes = T,
                        reduction = NULL, group.by = NULL, split.by = NULL, shape.by = NULL,
                        pt.size = 0.1, combine = FALSE,...) {

  if (is.null(colormap)) {
    colormap.use <- grDevices::colorRampPalette(c("#FFFFEF", "#FFFF00", "#FF0000", "#0A0000"))(64)
    colormap.use[1] <- "#E5E5E5"
  } else {
    colormap.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = colormap)
    }, error = function(e) {
      scales::viridis_pal(option = colormap, direction = -1)(n.colors)
    })
    if (color.direction == -1) {
      colormap.use <- rev(colormap.use)
    }
  }

  if (combine) {
    gg <- FeaturePlot(object, features,
                      reduction = reduction, split.by = split.by, shape.by = shape.by,
                      pt.size = pt.size, combine = combine, ...)
    gg <- gg + theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in")) +
      theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
      theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.line = element_line(colour = 'black')) +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
    gg <- gg + scale_colour_gradientn(colours = colormap.use, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.3), na.value = "lightgrey")
    # if (is.null(colormap)) {
    #   gg <- gg + scale_colour_gradientn(colours = colormap.use, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5), na.value = "lightgrey")
    # } else {
    #   gg <- gg + scale_color_distiller(palette = colormap, direction = color.direction, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5), na.value = "lightgrey")
    # }
    if (!show.text.axis) {
      gg <- gg + theme(axis.ticks=element_blank(), axis.text=element_blank())
    }
    if (!show.axes) {
      gg <- gg + theme_void()
    }
    if (!show.legend) {
      gg <- gg + theme(legend.position = "none")
    }

    if (show.legend.combined & i == length(gg)) {
      gg <- gg + theme(legend.position = "right", legend.key.height = grid::unit(0.15, "in"), legend.key.width = grid::unit(0.3, "in"), legend.title = NULL)
    }

  } else {
    gg <- FeaturePlot(object, features,
                      reduction = reduction, split.by = split.by, shape.by = shape.by,
                      pt.size = pt.size, combine = combine, ...)
    for (i in 1:length(gg)) {
      gg[[i]] <- gg[[i]] + theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in")) +
        theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
        theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.line = element_line(colour = 'black')) +
        theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
      gg[[i]] <- gg[[i]] + scale_colour_gradientn(colours = colormap.use, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.3), na.value = "lightgrey")
      # if (is.null(colormap)) {
      #   gg[[i]] <- gg[[i]] + scale_colour_gradientn(colours = colormap.use, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.3), na.value = "lightgrey")
      # } else {
      #   gg[[i]] <- gg[[i]] + scale_color_distiller(palette = colormap, direction = color.direction, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5), na.value = "lightgrey")
      # }
      xlims <- ggplot_build(gg[[i]])$layout$panel_scales_x[[1]]$range$range
      xlims <- c(xlims[1]+max(xlims[1]*0.1, -1), xlims[2]+min(xlims[2]*0.1, 1))
      # xlims[1] <- max(xlims[1]*1.1, floor(xlims[1]))
      # xlims[2] <- min(xlims[2]*1.1, ceiling(xlims[2]))
      ylims <- ggplot_build(gg[[i]])$layout$panel_scales_y[[1]]$range$range
      ylims <- c(ylims[1]+max(ylims[1]*0.1, -1), ylims[2]+min(ylims[2]*0.1, 1))
      gg[[i]] <- gg[[i]] + xlim(xlims) + ylim(ylims)

      if (!show.text.axis) {
        gg[[i]] <- gg[[i]] + theme(axis.ticks=element_blank(), axis.text=element_blank())
      }
      if (!show.axes) {
        gg[[i]] <- gg[[i]] + theme_void()
      }
      if (!show.legend) {
        gg[[i]] <- gg[[i]] + theme(legend.position = "none")
      }

      if (show.legend.combined & i == length(gg)) {
        gg[[i]] <- gg[[i]] + theme(legend.position = "right", legend.key.height = grid::unit(0.15, "in"), legend.key.width = grid::unit(0.3, "in"), legend.title = NULL)
      }
    }
  }
  return(gg)
}



#' Dot plot
#'
#'The size of the dot encodes the percentage of cells within a class, while the color encodes the AverageExpression level across all cells within a class
#'
#' @param object seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param rotation whether rotate the plot
#' @param colors.use defining the color for each condition when split.by is not NULL
#' @param colormap RColorbrewer palette to use (check available palette using RColorBrewer::display.brewer.all()). default will use customed color palette
#' @param color.direction Sets the order of colours in the scale. If 1, the default, colours are as output by RColorBrewer::brewer.pal(). If -1, the order of colours is reversed.
#' @param idents Which classes to include in the plot (default is all)
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' @param legend.width legend width
#' @param legend.title legend title 
#' @param scale whther show x-axis text
#' @param col.min Minimum scaled average expression threshold (everything smaller will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger will be set to this)
#' @param dot.scale Scale the size of the points, similar to cex
#' @param assay Name of assay to use, defaults to the active assay
#' @param angle.x angle for x-axis text rotation
#' @param hjust.x adjust x axis text
#' @param angle.y angle for y-axis text rotation
#' @param hjust.y adjust y axis text
#' @param show.legend whether show the legend
#' @param ... Extra parameters passed to DotPlot from Seurat package
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom Seurat DotPlot
dotPlot <- function(object, features, rotation = TRUE, colors.use = NULL,colormap = "OrRd", color.direction = 1, scale = TRUE, col.min = -2.5, col.max = 2.5, dot.scale = 6, assay = NULL,
                    idents = NULL, group.by = NULL, split.by = NULL, legend.width = 0.5, legend.title = "Scaled expression",
                    angle.x = 45, hjust.x = 1, angle.y = 0, hjust.y = 0.5, show.legend = TRUE, ...) {
  if (is.null(colors.use)) {
    colors.use = c(ggPalette(2)[1], ggPalette(2)[2])
  }
  gg <- DotPlot(object, features = features, assay = assay, cols = colors.use,
                scale = scale, col.min = col.min, col.max = col.max, dot.scale = dot.scale,
                idents = idents, group.by = group.by, split.by = split.by,...)
  gg <- gg + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.line = element_line(colour = 'black')) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))+
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x), axis.text.y = element_text(angle = angle.y, hjust = hjust.y))
  
  gg <- gg + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))
  if (is.null(split.by)) {
    gg <- gg + guides(color = guide_colorbar(barwidth = legend.width, title = legend.title),size = guide_legend(title = 'Percent expressed'))
  }
  
  if (rotation) {
    gg <- gg + coord_flip()
  }
  if (!is.null(colormap)) {
    if (is.null(split.by)) {
      gg <- gg + scale_color_distiller(palette = colormap, direction = color.direction, guide = guide_colorbar(title = legend.title, ticks = T, label = T, barwidth = legend.width), na.value = "lightgrey")
    }
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  return(gg)
}




#' Violin plot
#'
#' @param object seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param colors.use defining the color for each cell group
#' @param colors.ggplot whether use ggplot color scheme; default: colors.ggplot = FALSE
#' @param idents Which classes to include in the plot (default is all)
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' @param fill.by Color violins based on either 'feature' or 'ident'
#' @param show.text.x whther show x-axis text
#' @param angle.x angle for x-axis text rotation
#' @param hjust.x adjust x axis text
#' @param show.median whether show the median value
#' @param median.size the shape size of the median
#' @param show.boxplot whether show box plot
#' @param pt.size size of the dots
#' @param stat.add whether add statistical test
#' @param comparisons A list of length-2 vectors. The entries in the vector are either the names of 2 values on the x-axis or the 2 integers that correspond to the index of the groups of interest, to be compared.
#' @param stat.method a character string indicating which method to be used for comparing
#' @param label.format character string specifying label type. Allowed values include "p.signif" (shows the significance levels), "p.format" (shows the formatted p value).
#' @param label.size label size
#' @param label.x numeric Coordinates (in data units) to be used for absolute positioning of the label. If too short they will be recycled.
#' @param y.adjust adjust the y.max when stat.add is TRUE
#' @param y.adjust.scale scale factor to adjust the y.max when stat.add is TRUE and y.adjust is NULL
#' @param combine Combine plots into a single patchworked ggplot object. If FALSE, return a list of ggplot objects
#' @param show.legend whether show the legend
#' @param ... Extra parameters passed to VlnPlot from Seurat package
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#' @importFrom stats median
#' @importFrom purrr map2 map_dbl
#' @importFrom ggpubr stat_compare_means
vlnPlot <- function(object, features, colors.use = NULL, colors.ggplot = FALSE,
                    idents = NULL, group.by = NULL, split.by = NULL, fill.by = "ident",
                    show.text.x = TRUE, angle.x = 45, hjust.x = 1,
                    show.median = FALSE, median.size = 1, show.boxplot = FALSE, pt.size = 0,
                    stat.add = FALSE, comparisons = list(c("1","2")), stat.method = "wilcox.test", label.format = "p.format", label.size = 3, label.x = NULL, y.adjust = NULL, y.adjust.scale = 1/5, 
                    combine = FALSE, show.legend = FALSE, ...) {
  if (is.null(colors.use)) {
    if (is.null(group.by)) {
      cell.levels <- levels(Idents(object))
    } else {
      cell.levels <- levels(object@meta.data[,group.by])
    }
    numCluster <- length(cell.levels)
    if (colors.ggplot) {
      colors.use <- c(ggPalette(2)[1], ggPalette(2)[2])
    } else {
      colors.use <- scPalette(numCluster)
    }
  }
  for (ii in 1:length(comparisons)) {
    comparisons[[ii]] <- plyr::mapvalues(comparisons[[ii]], from = comparisons[[ii]],to = cell.levels[as.numeric(comparisons[[ii]])])
  }
  if (combine) {
    gg <- VlnPlot(object, features, idents = idents, group.by = group.by, split.by = split.by,
                  pt.size = pt.size, cols = colors.use, combine = combine,...)
    # gg <- gg + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    #   theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.line = element_line(colour = 'black')) +
    #   theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))+ theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x))
    if (show.median) {
      gg <- gg + stat_summary(fun.y=median, geom="point", shape=3, size=median.size)
    }
    if (show.boxplot) {
      gg <- gg + geom_boxplot(width=0.1, fill="white", outlier.shape = NA)
    }
    if (!show.legend) {
      gg <- gg + theme(legend.position = "none")
    }
    if (!show.text.x) {
      gg <- gg + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
    }
  } else {
    gg <- VlnPlot(object, features, idents = idents, group.by = group.by, split.by = split.by,
                  pt.size = pt.size, cols = colors.use, combine = combine)
    if (stat.add) {
      # change the y-axis tick to only max value
      ymaxs<- purrr::map_dbl(gg, extract_max)
      if (is.null(y.adjust)) {
        y.adjust <- signif(ymaxs*y.adjust.scale,2)
      }
      gg<- purrr::map2(gg, ymaxs + y.adjust, function(x,y) x +
                         scale_y_continuous(breaks = c(y)) +
                         expand_limits(y = y))
    }
    for(i in 1:length(gg)) {
      gg[[i]] <- gg[[i]] + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(text = element_text(size = 10)) +
        theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.line = element_line(colour = 'black')) +
        theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x))
      if (show.median) {
        gg[[i]] <- gg[[i]] + stat_summary(fun.y=median, geom="point", shape=3, size=median.size)
      }
      if (show.boxplot) {
        gg[[i]] <- gg[[i]] + geom_boxplot(width=0.1, fill="white", outlier.shape = NA)
      }
      if (!show.legend) {
        gg[[i]] <- gg[[i]] + theme(legend.position = "none")
      }
      if (!show.text.x) {
        gg[[i]] <- gg[[i]] + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
      }
      if (stat.add) {
        gg[[i]] <- gg[[i]] + ggpubr::stat_compare_means(comparisons = comparisons, method = stat.method, label.x = label.x,
                                                        label = label.format, size = label.size, label.y = ymaxs[[i]])
      }
    }
  }
  return(gg)
}




#' Stacked Violin plot
#'
#' @param object seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param colors.use defining the color for each cell group
#' @param colors.ggplot whether use ggplot color scheme; default: colors.ggplot = FALSE
#' @param split.by Name of a metadata column to split plot by;
#' @param idents Which classes to include in the plot (default is all)
#' @param show.text.y whther show y-axis text
#' @param line.size line width in the violin plot
#' @param pt.size size of the dots
#' @param plot.margin adjust the white space between each plot
#' @param angle.x angle for x-axis text rotation
#' @param vjust.x adjust x axis text
#' @param hjust.x adjust x axis text
#' @param ... Extra parameters passed to VlnPlot from Seurat package
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom  patchwork wrap_plots
#' @importFrom purrr map map2 map_dbl
#' @importFrom Seurat VlnPlot
StackedVlnPlot<- function(object, features, idents = NULL, split.by = NULL,
                          colors.use = NULL, colors.ggplot = FALSE,
                          angle.x = 90, vjust.x = NULL, hjust.x = NULL, show.text.y = TRUE, line.size = NULL,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  options(warn=-1)
  if (is.null(colors.use)) {
    numCluster <- length(levels(Idents(object)))
    if (colors.ggplot) {
      colors.use <- c(ggPalette(2)[1], ggPalette(2)[2])
    } else {
      colors.use <- scPalette(numCluster)
    }
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 1, 1)
    vjust=c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(object = object, features = x, idents = idents, split.by = split.by, cols = colors.use, pt.size = pt.size,
                                                              show.text.y = show.text.y, line.size = line.size, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x)) +
    theme(axis.text.x = element_text(size = 10))
  
  # # change the y-axis tick to only max value
  # ymaxs<- purrr::map_dbl(plot_list, extract_max)
  # print(ymaxs)
  # plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
  #                           scale_y_continuous(breaks = c(y)) +
  #                           expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#' modified vlnplot
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param split.by Name of a metadata column to split plot by;
#' @param idents Which classes to include in the plot (default is all)
#' @param cols defining the color for each cell group
#' @param show.text.y whther show y-axis text
#' @param line.size line width in the violin plot
#' @param pt.size size of the dots
#' @param plot.margin adjust the white space between each plot
#' @param ... pass any arguments to VlnPlot in Seurat
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#'
modify_vlnplot<- function(object,
                          features,
                          idents = NULL,
                          split.by = NULL,
                          cols = NULL,
                          show.text.y = TRUE,
                          line.size = NULL,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  options(warn=-1)
  p<- VlnPlot(object, features = features, cols = cols, pt.size = pt.size, idents = idents, split.by = split.by,  ... )  +
    xlab("") + ylab(features) + ggtitle("")
  p <- p + theme(text = element_text(size = 10)) + theme(axis.line = element_line(size=line.size)) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.line.x = element_line(colour = 'black', size=line.size),axis.line.y = element_line(colour = 'black', size= line.size))
  # theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  p <- p + theme(legend.position = "none",
                 plot.title= element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(size = rel(1), angle = 0),
                 axis.text.y = element_text(size = rel(1)),
                 plot.margin = plot.margin ) +
    theme(axis.text.y = element_text(size = 8))
 
  p <- p + scale_y_continuous(labels = function(x) {
    idx0 = which(x == 0)
    if (idx0 > 1) {
      c(rep(x = "", times = idx0-1), "0",rep(x = "", times = length(x) -2-idx0), x[length(x) - 1], "")
    } else {
      c("0", rep(x = "", times = length(x)-3), x[length(x) - 1], "")
    }
  })
   #  #c(rep(x = "", times = length(x)-2), x[length(x) - 1], ""))
    
  p <- p + theme(element_line(size=line.size))
  
  if (!show.text.y) {
    p <- p + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())
  }
  return(p)
}

#' extract the max value of the y axis
#' @param p ggplot object
#' @importFrom  ggplot2 ggplot_build
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(signif(ymax,2))
}


#' Feature expression heatmap
#'
#' Draws a heatmap of single cell feature expression.
#'
#' @param object Seurat object
#' @param features A vector of features to plot, defaults to \code{VariableFeatures(object = object)}
#' @param cells A vector of cells to plot
#' @param group.bar Add a color bar showing group status for cells
#' @param group.by A vector of variables to group cells by; pass 'ident' to group by cell identity classes
#' @param additional.group.by A vector of variables to group cells by.
#' @param additional.group.sort.by which variable from additional.group.by to use for sorting cells
#' @param colors.use A list of colors to use for the color bar. e.g., colors.use = list(ident = c("red","blue"))
#' @param colors.ggplot whether use ggplot color scheme; default: colors.ggplot = FALSE
#' @param color.heatmap A character string or vector indicating the colormap option to use. It can be the avaibale color palette in viridis_pal() or brewer.pal()
#' @param direction Sets the order of colors in the scale. If 1, the default colors are used. If -1, the order of colors is reversed.
#' @param n.colors number of basic colors to generate from color palette
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped); defaults to 2.5
#' if \code{slot} is 'scale.data', 6 otherwise
#' @param assay Assay to pull from
# @param check.plot Check that plotting will finish in a reasonable amount of time
#' @param text.size text size of gene names
#' @param label Label the cell identies above the color bar
#' @param size Size of text above color bar
#' @param hjust Horizontal justification of text above color bar
#' @param angle Angle of text above color bar
#' @param raster If true, plot with geom_raster, else use geom_tile. geom_raster may look blurry on
#' some viewing applications such as Preview due to how the raster is interpolated. Set this to FALSE
#' if you are encountering that issue (note that plots may take longer to produce/render).
#' @param draw.lines Include white lines to separate the groups
#' @param lines.width Integer number to adjust the width of the separating white lines.
#' Corresponds to the number of "cells" between each group.
#' @param group.bar.height Scale the height of the color bar
#' @param legend.remove whether remove legend
#' @param legend.title legend title
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom stats median
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian scale_color_manual
#' ggplot_build aes_string
#' @importFrom patchwork wrap_plots
#' @importFrom grid textGrob gpar
#' @importFrom RColorBrewer brewer.pal
#' @import rlang
#' @import Seurat
#' @export
doHeatmap <- function (object,
                       features = NULL,
                       cells = NULL,
                       group.by = "ident",
                       colors.use = NULL,
                       colors.ggplot = FALSE,
                       group.bar = TRUE,
                       additional.group.by = NULL,
                       additional.group.sort.by = NULL,
                       color.heatmap = "RdBu",
                       n.colors = 11,
                       direction = -1,
                       disp.min = -2.1,
                       disp.max = 2.1,
                       slot = "scale.data",
                       assay = NULL,
                       label = TRUE,
                       text.size = 8,
                       size = 4,
                       hjust = 0,
                       angle = 45,
                       raster = TRUE,
                       draw.lines = TRUE,
                       lines.width = NULL,
                       group.bar.height = 0.02,
                       legend.remove = FALSE,
                       legend.title = NULL,
                       combine = TRUE)
{
  # if (is.null(color.heatmap)) {
  #   color.heatmap.use <- PurpleAndYellow()
  # } else

  group.colors <- colors.use
  if (!is.null(additional.group.by)) {
    if (is.null(additional.group.sort.by)) {
      additional.group.sort.by <- additional.group.by[1]
    }
  }

  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object,
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot,
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ",
            slot, " slot for the ", assay, " assay: ", paste(bad.features,
                                                             collapse = ", "))
  }

  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ",
                paste(bad.sorts, collapse = ", "))
      }
    }
  }

  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object,
                                                             slot = slot)[features, cells, drop = FALSE])))

  object <- suppressMessages(expr = StashIdent(object = object,
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"

  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }

    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]

    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }
    }

    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) *
                                           lines.width), FUN = function(x) {
                                             return(RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }

      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells

      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells

      group.use <- rbind(group.use, placeholder.groups)

      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }

      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }

    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])

    plot <- SingleRasterMap(data = data.group, raster = raster,
                            disp.min = disp.min, disp.max = disp.max, feature.order = features,
                            cell.order = rownames(x = group.use), group.by = group.use[[i]], legend.title = legend.title)

    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }

        # Default
        if (colors.ggplot) {
          cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))
        } else {
          cols[[colname]] <- scPalette(length(x = levels(x = group.use[[colname]])))
          if (colname %in% additional.group.by) {
            cols[[colname]] <- ggPalette(length(x = levels(x = group.use[[colname]])))
          }
        }


        #Overwrite if better value is provided
        if (!is_null(group.colors[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(group.colors[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(group.colors[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(group.colors[[colname]]))) {
                cols[[colname]] <- as.vector(group.colors[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(group.colors[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }

        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])

        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)

        plot <- suppressMessages(plot +
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) +
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = grid::gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off"))

        #  same time run  or not   ggplot version    ?
        if ((colname == i ) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          # x.divs <- pbuild$layout$panel_params[[1]]$x.major
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% attr(x = pbuild$layout$panel_params[[1]]$x$get_breaks(), which = "pos")
          # group.use$x <- x.divs
          # label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
          #                       FUN = median) * x.max
          # label.x.pos <- data.frame(group = names(x = label.x.pos),
          #                           label.x.pos)

          x <- data.frame(group = sort(x = group.use[[colname]]), x = x.divs)
          label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = function(y) {
            if (isTRUE(x = draw.lines)) {
              mean(x = y[-length(x = y)])
            } else {
              mean(x = y)
            }
          })
          label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)

          plot <- plot + geom_text(stat = "identity",
                                   data = label.x.pos, aes_string(label = "group",
                                                                  x = "label.x.pos"), y = y.max + y.max *
                                     0.03 * 0.5, angle = angle, hjust = hjust,
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0,
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) *
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank()) +
      theme(axis.text.y = element_text(size = text.size))
    if (!is.null(color.heatmap)) {
      if (length(color.heatmap) == 1) {
        color.heatmap.use <- tryCatch({
          RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
        }, error = function(e) {
          scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
        })
      } else if (length(color.heatmap) > 1) {
        color.heatmap.use <- color.heatmap
      }
      if (direction == -1) {
        color.heatmap.use <- rev(color.heatmap.use)
      }
      color.heatmap.use <- colorRampPalette(color.heatmap.use)(99)
      plot <- plot + scale_fill_gradientn(colors = color.heatmap.use)
    }

    if (legend.remove) {
      plot <-  plot + NoLegend()
    }
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}

# A single heatmap from ggplot2 using geom_raster
#
# @param data A matrix or data frame with data to plot
# @param raster switch between geom_raster and geom_tile
# @param cell.order ...
# @param feature.order ...
# @param cols A vector of colors to use
# @param disp.min Minimum display value (all values below are clipped)
# @param disp.max Maximum display value (all values above are clipped)
# @param limits A two-length numeric vector with the limits for colors on the plot
# @param group.by A vector to group cells by, should be one grouping identity per cell
# @param legend.title legend title
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradient
#' scale_fill_gradientn theme element_blank labs geom_point guides guide_legend geom_tile
#
SingleRasterMap <- function(
  data,
  raster = TRUE,
  cell.order = NULL,
  feature.order = NULL,
  colors = PurpleAndYellow(),
  disp.min = -2.5,
  disp.max = 2.5,
  limits = NULL,
  group.by = NULL,
  legend.title = NULL
) {
  data <- MinMax(data = data, min = disp.min, max = disp.max)
  data <- Melt(x = t(x = data))
  colnames(x = data) <- c('Feature', 'Cell', 'Expression')
  if (!is.null(x = feature.order)) {
    data$Feature <- factor(x = data$Feature, levels = unique(x = feature.order))
  }
  if (!is.null(x = cell.order)) {
    data$Cell <- factor(x = data$Cell, levels = unique(x = cell.order))
  }
  if (!is.null(x = group.by)) {
    data$Identity <- group.by[data$Cell]
  }
  limits <- limits %||% c(min(data$Expression), max(data$Expression))
  if (length(x = limits) != 2 || !is.numeric(x = limits)) {
    stop("limits' must be a two-length numeric vector")
  }
  my_geom <- ifelse(test = raster, yes = geom_raster, no = geom_tile)
  plot <- ggplot(data = data) +
    my_geom(mapping = aes_string(x = 'Cell', y = 'Feature', fill = 'Expression')) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_gradientn(limits = limits, colors = colors, na.value = "white") +
    labs(x = NULL, y = NULL, fill = group.by %iff% 'Expression') +
    WhiteBackground() + NoAxes(keep.text = TRUE)
  plot <- plot + guides(fill = guide_colourbar(barwidth = 0.5, title = legend.title))
  # if (!is.null(x = group.by)) {
  #   plot <- plot + geom_point(
  #     mapping = aes_string(x = 'Cell', y = 'Feature', color = 'Identity'),
  #     alpha = 0
  #   ) +
  #     guides(color = guide_legend(override.aes = list(alpha = 1)))
  # }
  return(plot)
}
# Melt a data frame
#
# @param x A data frame
#
# @return A molten data frame
#
Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}
# Generate a random name
#
# Make a name from randomly sampled lowercase letters,
# pasted together with no spaces or other characters
#
# @param length How long should the name be
# @param ... Extra parameters passed to sample
#
# @return A character with nchar == length of randomly sampled letters
#
# @seealso \code{\link{sample}}
#
RandomName <- function(length = 5L, ...) {
  return(paste(sample(x = letters, size = length, ...), collapse = ''))
}
# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

#' Run PHATE from Seurat object
#'
#' PHATE is a data reduction method specifically designed for visualizing
#' **high** dimensional data in **low** dimensional spaces.
#' To run, you must first install the `phate` python
#' package (e.g. via pip install phate). Details on this package can be
#' found here: \url{https://github.com/KrishnaswamyLab/PHATE}. For a more in depth
#' discussion of the mathematics underlying PHATE, see the Nature Biotechnology paper here:
#' \url{https://www.nature.com/articles/s41587-019-0336-3}.
#'
#' @param object Seurat object
#' @param dims Which dimensions to use as input features, used only if
#' \code{features} is NULL
#' @param reduction Which dimensional reduction (PCA or ICA) to use for the
#' PHATE input. Default is PCA
#' @param features If set, run PHATE on this subset of features (instead of running on a
#' set of reduced dimensions). Not set (NULL) by default; \code{dims} must be NULL to run
#' on features
#' @param assay Assay to pull data for when using \code{features}
#' @param n.components Total number of dimensions to embed in PHATE.
#' @param knn int, optional, default: 5
#' number of nearest neighbors on which to build kernel
#' @param decay int, optional, default: 40
#' sets decay rate of kernel tails.
#' If NA, alpha decaying kernel is not used
#' @param n.landmark int, optional, default: 2000
#' number of landmarks to use in fast PHATE
#' @param gamma float, optional, default: 1
#' Informational distance constant between -1 and 1.
#' `gamma=1` gives the PHATE log potential, `gamma=0` gives
#' a square root potential.
#' @param t int, optional, default: 'auto'
#' power to which the diffusion operator is powered
#' sets the level of diffusion
#' @param mds.solver {'sgd', 'smacof'}, optional, default: 'sgd'
#' which solver to use for metric MDS. SGD is substantially faster,
#' but produces slightly less optimal results. Note that SMACOF was used
#' for all figures in the PHATE paper.
#' @param knn.dist.method string, optional, default: 'euclidean'.
#' recommended values: 'euclidean', 'cosine', 'precomputed'
#' Any metric from `scipy.spatial.distance` can be used
#' distance metric for building kNN graph. If 'precomputed',
#' `data` should be an n_samples x n_samples distance or
#' affinity matrix. Distance matrices are assumed to have zeros
#' down the diagonal, while affinity matrices are assumed to have
#' non-zero values down the diagonal. This is detected automatically using
#' `data[0,0]`. You can override this detection with
#' `knn.dist.method='precomputed_distance'` or
#' `knn.dist.method='precomputed_affinity'`.
#' @param mds.method string, optional, default: 'metric'
#' choose from 'classic', 'metric', and 'nonmetric'
#' which MDS algorithm is used for dimensionality reduction
#' @param mds.dist.method string, optional, default: 'euclidean'
#' recommended values: 'euclidean' and 'cosine'
#' @param t.max int, optional, default: 100.
#' Maximum value of t to test for automatic t selection.
#' @param npca int, optional, default: 100
#' Number of principal components to use for calculating
#' neighborhoods. For extremely large datasets, using
#' n_pca < 20 allows neighborhoods to be calculated in
#' log(n_samples) time.
#' @param plot.optimal.t boolean, optional, default: FALSE
#' If TRUE, produce a plot showing the Von Neumann Entropy
#' curve for automatic t selection.
#' @param verbose `int` or `boolean`, optional (default : 1)
#' If `TRUE` or `> 0`, print verbose updates.
#' @param n.jobs `int`, optional (default: 1)
#' The number of jobs to use for the computation.
#' If -1 all CPUs are used. If 1 is given, no parallel computing code is
#' used at all, which is useful for debugging.
#' For n_jobs below -1, (n.cpus + 1 + n.jobs) are used. Thus for
#' n_jobs = -2, all CPUs but one are used
#' @param seed.use int or `NA`, random state (default: `42`)
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. phate by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. PHATE_ by default
#' @param k Deprecated. Use `knn`.
#' @param alpha Deprecated. Use `decay`.
#' @param ... other parameters to phateR::phate
#'
#' @return Returns a Seurat object containing a PHATE representation
#'
#' @references Moon K, van Dijk D, Wang Z, Gigante S,
#' Burkhardt D, Chen W, van den Elzen A,
#' Hirn M, Coifman R, Ivanova N, Wolf G and Krishnaswamy S (2017).
#' "Visualizing Transitions and Structure for High Dimensional Data
#' Exploration." _bioRxiv_, pp. 120378. doi: 10.1101/120378
#' (URL: http://doi.org/10.1101/120378),
#' <URL: https://www.biorxiv.org/content/early/2017/12/01/120378>.
#' @examples
#' if (reticulate::py_module_available("phate")) {
#'
#' # Load data
#' pbmc_small
#'
#' # Run PHATE with default parameters
#' pbmc_small <- RunPHATE(object = pbmc_small)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction = 'phate')
#'
#' # Try smaller `k` for a small dataset, and larger `t` for a noisy embedding
#' pbmc_small <- RunPHATE(object = pbmc_small, k = 4, t = 12)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction = 'phate')
#'
#' # For increased emphasis on local structure, use sqrt potential
#' pbmc_small <- RunPHATE(object = pbmc_small, gamma = 0)
#' # Plot results
#' DimPlot(object = pbmc_small, reduction = 'phate')
#'
#' }
#'
#' @export
#'
RunPHATE <- function(
  object,
  dims = NULL,
  reduction = 'pca',
  features = NULL,
  assay = 'RNA',
  n.components = 2L,
  knn = 20L,
  decay = 40L,
  n.landmark=2000L,
  gamma = 1,
  t = "auto",
  mds.solver = 'sgd',
  knn.dist.method = "euclidean",
  mds.method = "metric",
  mds.dist.method = "euclidean",
  t.max = 100,
  npca = 40,
  plot.optimal.t = FALSE,
  verbose = 1,
  n.jobs = 1,
  seed.use = 42,
  reduction.name = "phate",
  reduction.key = "PHATE_",
  k = NULL,
  alpha = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = dims) && !is.null(x = features)) {
    stop("Please specify only one of the following arguments: dims or features")
  }
  if (!is.null(x = features)) {
    data.use <- t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, ])
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- assay %||% DefaultAssay(object = object[[reduction]])
  } else {
    data.use <- Embeddings(object[[reduction]])
    assay <- assay %||% DefaultAssay(object = object[[reduction]])
  }
  object[[reduction.name]] <- RunPHATE_internal(
    object = data.use,
    assay = assay,
    n.components = n.components,
    knn = knn,
    decay = decay,
    n.landmark = n.landmark,
    gamma = gamma,
    t = t,
    knn.dist.method = knn.dist.method,
    mds.solver = mds.solver,
    mds.method = mds.method,
    mds.dist.method = mds.dist.method,
    t.max = t.max,
    npca = npca,
    plot.optimal.t = plot.optimal.t,
    verbose = verbose,
    n.jobs = n.jobs,
    seed.use = seed.use,
    reduction.key = reduction.key,
    k = k,
    alpha = alpha,
    ...
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}


#' Run PHATE from a data matrix
#'
#' @param object matrix (n_samples, n_dimensions) 2 dimensional input data array with n_samples samples and n_dimensions dimensions.
#' If knn.dist.method is 'precomputed', data is treated as a (n_samples, n_samples) distance or affinity matrix
#' @param assay Assay to pull data for when using \code{features}
#' @param n.components Total number of dimensions to embed in PHATE.
#' @param knn int, optional, default: 5
#' number of nearest neighbors on which to build kernel
#' @param decay int, optional, default: 40
#' sets decay rate of kernel tails.
#' If NA, alpha decaying kernel is not used
#' @param n.landmark int, optional, default: 2000
#' number of landmarks to use in fast PHATE
#' @param gamma float, optional, default: 1
#' Informational distance constant between -1 and 1.
#' `gamma=1` gives the PHATE log potential, `gamma=0` gives
#' a square root potential.
#' @param t int, optional, default: 'auto'
#' power to which the diffusion operator is powered
#' sets the level of diffusion
#' @param mds.solver {'sgd', 'smacof'}, optional, default: 'sgd'
#' which solver to use for metric MDS. SGD is substantially faster,
#' but produces slightly less optimal results. Note that SMACOF was used
#' for all figures in the PHATE paper.
#' @param knn.dist.method string, optional, default: 'euclidean'.
#' recommended values: 'euclidean', 'cosine', 'precomputed'
#' Any metric from `scipy.spatial.distance` can be used
#' distance metric for building kNN graph. If 'precomputed',
#' `data` should be an n_samples x n_samples distance or
#' affinity matrix. Distance matrices are assumed to have zeros
#' down the diagonal, while affinity matrices are assumed to have
#' non-zero values down the diagonal. This is detected automatically using
#' `data[0,0]`. You can override this detection with
#' `knn.dist.method='precomputed_distance'` or
#' `knn.dist.method='precomputed_affinity'`.
#' @param mds.method string, optional, default: 'metric'
#' choose from 'classic', 'metric', and 'nonmetric'
#' which MDS algorithm is used for dimensionality reduction
#' @param mds.dist.method string, optional, default: 'euclidean'
#' recommended values: 'euclidean' and 'cosine'
#' @param t.max int, optional, default: 100.
#' Maximum value of t to test for automatic t selection.
#' @param npca int, optional, default: 100
#' Number of principal components to use for calculating
#' neighborhoods. For extremely large datasets, using
#' n_pca < 20 allows neighborhoods to be calculated in
#' log(n_samples) time.
#' @param plot.optimal.t boolean, optional, default: FALSE
#' If TRUE, produce a plot showing the Von Neumann Entropy
#' curve for automatic t selection.
#' @param verbose `int` or `boolean`, optional (default : 1)
#' If `TRUE` or `> 0`, print verbose updates.
#' @param n.jobs `int`, optional (default: 1)
#' The number of jobs to use for the computation.
#' If -1 all CPUs are used. If 1 is given, no parallel computing code is
#' used at all, which is useful for debugging.
#' For n_jobs below -1, (n.cpus + 1 + n.jobs) are used. Thus for
#' n_jobs = -2, all CPUs but one are used
#' @param seed.use int or `NA`, random state (default: `42`)
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. PHATE_ by default
#' @param k Deprecated. Use `knn`.
#' @param alpha Deprecated. Use `decay`.
#' @param ... other parameters to phateR::phate

#' @export
#'
RunPHATE_internal <- function(
  object,
  assay = NULL,
  n.components = 2L,
  knn = 20L,
  decay = 40L,
  n.landmark=2000L,
  gamma = 1,
  t = "auto",
  mds.solver = "sgd",
  knn.dist.method = "euclidean",
  mds.method = "metric",
  mds.dist.method = "euclidean",
  t.max=100,
  npca = 40,
  plot.optimal.t = FALSE,
  verbose = 1,
  n.jobs = 1,
  seed.use = 42,
  reduction.key = "PHATE_",
  k = NULL,
  alpha = NULL,
  ...
) {
  # if (!PackageCheck('phateR', error = FALSE)) {
  #   stop("Please install phateR - learn more at https://github.com/KrishnaswamyLab/phateR")
  # }
  if (!is.null(k)) {
    message("Argument k is deprecated. Using knn instead.")
    knn <- k
  }
  if (!is.null(alpha)) {
    message("Argument alpha is deprecated. Using decay instead.")
    decay <- alpha
  }
  #  CheckDots(...)
  phate_output <- phateR::phate(
    object,
    ndim = n.components,
    knn = knn,
    decay = decay,
    n.landmark = n.landmark,
    gamma = gamma,
    t = t,
    knn.dist.method = knn.dist.method,
    init = NULL,
    mds.solver = mds.solver,
    mds.method = mds.method,
    mds.dist.method = mds.dist.method,
    t.max = t.max,
    npca = npca,
    plot.optimal.t = plot.optimal.t,
    verbose = verbose,
    n.jobs = n.jobs,
    seed = seed.use,
    ...
  )
  phate_output <- as.matrix(phate_output)
  colnames(x = phate_output) <- paste0(reduction.key, 1:ncol(x = phate_output))
  rownames(x = phate_output) <- rownames(object)
  phate.reduction <- CreateDimReducObject(
    embeddings = phate_output,
    key = reduction.key,
    assay = assay
  )
  return(phate.reduction)
}
                         
 
compareDistribution <- function(df, color.use = NULL, title = NULL, group.labels = c("Brachial", "Lumbar"), alpha = 0.6) {
  gg <- ggplot(df, aes(x=data, fill=group)) +
    geom_density(alpha=alpha) +
    labs(title = title, x = "Expression levels", y = "Density") + theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = 10)) + theme(legend.title = element_blank())
  if (is.null(color.use)) {
    gg <- gg + scale_color_manual(name = "", labels = group.labels)
  } else {
    gg <- gg + scale_color_manual(values = color.use, name = "", labels = group.labels)
  }
  return(gg)
}
                         
                         


