library(Seurat)
library(readxl)
library(patchwork)

output_file <- read_excel("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF7 and T47D KLF4 Data/Analysis of KLF4 ChIP-Seq Data/KLF4_intersect.anno.xlsx")

colnames(output_file)

output_file <- output_file[output_file$distance < 20000, ]

output_file = output_file[order(output_file[,'ensembl_gene_id'],output_file[,'distance']),]
output_file = output_file[!duplicated(output_file$ensembl_gene_id),]

nearby_KLF4_genes <- output_file$external_gene_name

MCF7SO <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF7 and T47D KLF4 Data/MCF7SO.rds")
MCF7SO.updated <- UpdateSeuratObject(object = MCF7SO)
saveRDS(MCF7SO.updated, "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF7 and T47D KLF4 Data/MCF7SO.updated")

MCF7SO.updated <- readRDS(file = "C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF7 and T47D KLF4 Data//MCF7SO.updated")
MCF7SO.updated@assays$RNA@var.features <- nearby_KLF4_genes


T47DSO <- readRDS("C:/Users/jainar/Documents/Gap Year/CEDAR Student Worker Position/Analysis of MCF7 and T47D KLF4 Data/T47DSO.rds")

T47DSO@assays$RNA@var.features <- nearby_KLF4_genes

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

SingleCorPlot <- function(
  data,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  smooth = FALSE,
  rows.highlight = NULL,
  legend.title = NULL,
  na.value = 'grey50',
  span = NULL
) {
  pt.size <- pt.size <- pt.size %||% AutoPointSize(data = data)
  orig.names <- colnames(x = data)
  names.plot <- colnames(x = data) <- gsub(
    pattern = '-',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  names.plot <- colnames(x = data) <- gsub(
    pattern = ':',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  if (ncol(x = data) < 2) {
    msg <- "Too few variables passed"
    if (ncol(x = data) == 1) {
      msg <- paste0(msg, ', only have ', colnames(x = data)[1])
    }
    stop(msg, call. = FALSE)
  }
  plot.cor <- round(x = cor(x = data[, 1], y = data[, 2]), digits = 2)
  if (!is.null(x = rows.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = rows.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size,
      cols.highlight = 'red',
      col.base = 'black',
      pt.size = pt.size
    )
    cols <- highlight.info$color
    col.by <- factor(
      x = highlight.info$highlight,
      levels = rev(x = highlight.info$plot.order)
    )
    plot.order <- order(col.by)
    data <- data[plot.order, ]
    col.by <- col.by[plot.order]
  }
  if (!is.null(x = col.by)) {
    data$colors <- col.by
  }
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = names.plot[1], y = names.plot[2])
  ) +
    labs(
      x = orig.names[1],
      y = orig.names[2],
      title = plot.cor,
      color = legend.title
    )
  if (smooth) {
    # density <- kde2d(x = data[, names.plot[1]], y = data[, names.plot[2]], h = Bandwidth(data = data[, names.plot]), n = 200)
    # density <- data.frame(
    #   expand.grid(
    #     x = density$x,
    #     y = density$y
    #   ),
    #   density = as.vector(x = density$z)
    # )
    plot <- plot + stat_density2d(
      mapping = aes(fill = ..density.. ^ 0.25),
      geom = 'tile',
      contour = FALSE,
      n = 200,
      h = Bandwidth(data = data[, names.plot])
    ) +
      # geom_tile(
      #   mapping = aes_string(
      #     x = 'x',
      #     y = 'y',
      #     fill = 'density'
      #   ),
      #   data = density
      # ) +
      scale_fill_continuous(low = 'white', high = 'dodgerblue4') +
      guides(fill = FALSE)
  }
  if (!is.null(x = col.by)) {
    plot <- plot + geom_point(
      mapping = aes_string(color = 'colors'),
      position = 'jitter',
      size = pt.size
    )
  } else {
    plot <- plot + geom_point(position = 'jitter', size = pt.size)
  }
  if (!is.null(x = cols)) {
    cols.scale <- if (length(x = cols) == 1 && cols %in% rownames(x = brewer.pal.info)) {
      scale_color_brewer(palette = cols)
    } else {
      scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + cols.scale
    if (!is.null(x = rows.highlight)) {
      plot <- plot + guides(color = FALSE)
    }
  }
  plot <- plot + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(x = span)) {
    plot <- plot + geom_smooth(
      mapping = aes_string(x = names.plot[1], y = names.plot[2]),
      method = 'loess',
      span = span
    )
  }
  return(plot)
}



VarFeaturePlot <- function(object, cols = c("black", "red"), pt.size = 1, log = NULL, selection.method = NULL, assay = NULL, gene_list) {
  
  if (length(x = cols) != 2) {
    stop("'cols' must be of length 2")
  }
  
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  
  hvf.info <- HVFInfo(object = MCF7SO.updated)
  hvf.info <- hvf.info[, c(1, 3)]
  hvf.info_genes <- rownames(hvf.info)
  is_nearby_gene <- ifelse(hvf.info_genes %in% nearby_KLF4_genes == TRUE, "nearby_genes", "other_genes")
  
  axis.labels <- switch(EXPR = colnames(x = hvf.info)[2], variance.standardized = c("Average Expression", "Standardized Variance"), dispersion.scaled = c("Average Expression", "Dispersion"), residual_variance = c("Geometric Mean of Expression", "Residual Variance"))
  log <- log %||% (any(c("variance.standardized", "residual_variance") %in% colnames(x = hvf.info)))
  
  plot <- SingleCorPlot(data = hvf.info, col.by = is_nearby_gene, pt.size = pt.size)
  
  plot <- plot + labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) 
  
  if (log) {
    plot <- plot + scale_x_log10()
  }
  return(plot)
}

VarFeaturePlot(MCF7SO.updated, gene_list = nearby_KLF4_genes)
VarFeaturePlot(T47DSO, gene_list = nearby_KLF4_genes)

