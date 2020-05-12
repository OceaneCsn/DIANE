library(pheatmap)
library(stringr)
library(ggplot2)
library(reshape2)
library(limma)

#' Draw heatmap
#'
#' @param data data to plot
#' @param subset subset of rows to display
#' @param show_rownames show rownames or not
#' @param title plot title
#' @param log Show log(expression+1) in the heatmap if TRUE, expression if FALSE
#' @param profiles Show expression/mean(expression) for each gene if TRUE, expression if FALSE
#' @param conditions if NULL, shows all the conditions, else if character vector, shows only the required ones
#'
#'
#' @importFrom pheatmap pheatmap
#' @importFrom stringr str_split_fixed
#' @export
#' @return plot of the heatmap

draw_heatmap <-
  function(data,
           subset = NULL,
           show_rownames = FALSE,
           title = "Random preview heatmap",
           log = TRUE,
           profiles = FALSE,
           conditions = NULL) {
    if (is.null(subset)) {
      sample_subset <- sample(rownames(data), size = 100)
    }
    else
      sample_subset <- subset
    
    if (is.null(conditions))
      conds <- colnames(data)
    else
      conds <- unique(grep(paste(conditions, collapse = "|"),
                           colnames(data), value = TRUE))
    
    if (log)
      data <- log(data + 1)
    if (profiles)
      data <- data / rowMeans(data)
    
    
    mat <- data[sample_subset, conds]
    
    sample <- stringr::str_split_fixed(colnames(mat), '_', 2) [, 1]
    samples <- data.frame(sample, row.names = colnames(mat))
    
    pheatmap::pheatmap(
      mat,
      color = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(n = 7, name = "YlGnBu"))(100),
      annotation_col = samples,
      show_rownames = show_rownames,
      main = title,
      fontsize = 17
    )
  }



#' draw_distributions of count data
#'
#' @param data count data
#' @param boxplot if TRUE, plot each sample as a boxplot, else, as a violin plot
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#'
#' @return plot
#' @export
#'
draw_distributions <- function(data, boxplot = TRUE) {
  d <-
    reshape2::melt(log(data[sample(rownames(data),
                                   replace = F,
                                   size = round(dim(data)[1] / 4, 0)),] + 1))
  
  colnames(d)[c(length(colnames(d)) - 1, length(colnames(d)))] <-
    c("sample", "logCount")
  
  d$condition <- str_split_fixed(d$sample, "_", 2)[, 1]
  g <- ggplot2::ggplot(data = d, aes(x = sample, y = logCount))
  
  if (boxplot) {
    g <-
      g + ggplot2::geom_boxplot(
        alpha = 0.5,
        lwd = 1,
        ggplot2::aes(fill = condition),
        outlier.color = "black",
        outlier.alpha = 0.1
      )
  } else{
    g <-
      g + ggplot2::geom_violin(alpha = 0.5, lwd = 1, aes(fill = condition))
  }
  
  g <-
    g + ggplot2::theme(
      plot.title = ggplot2::element_text(size = 22, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 20),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 20, face = "bold"),
      legend.text = ggplot2::element_text(size = 22, angle = 0),
      axis.text.y = ggplot2::element_text(size = 18, angle = 30),
      axis.text.x = ggplot2::element_text(
        size = 15,
        angle = -50,
        hjust = 0,
        colour = "grey50"
      ),
      legend.text.align = 1,
      axis.title = element_text(size = 24)
    )
  g
}

#' draw_MDS
#'
#' @param normalized.count data to plot for MDS
#' @importFrom limma plotMDS
#' @return MDS plot
#' @export
#'
draw_MDS <- function(normalized.count) {
  mds <- limma::plotMDS(normalized.count, plot = FALSE)
  d <-
    data.frame(
      dim1 = mds$x,
      dim2 = mds$y,
      sample = names(mds$x),
      condition = str_split_fixed(names(mds$x), '_', 2)[, 1]
    )
  g <-
    ggplot2::ggplot(data = d, aes(x = dim1, y = dim2, color = condition)) + 
    ggplot2::geom_point(size = 6)
  
  g <- g + ggplot2::ggtitle("Multi Dimensional Scaling plot")
  
  g + ggplot2::theme(
    plot.title = ggplot2::element_text(size = 22, face = "bold"),
    strip.text.x = ggplot2::element_text(size = 20),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 20, face = "bold"),
    legend.text = ggplot2::element_text(size = 20, angle = 0),
    axis.text.y = ggplot2::element_text(size = 20, angle = 0),
    axis.text.x = ggplot2::element_text(
      size = 20,
      angle = 0,
      hjust = 0
    ),
    legend.text.align = 1,
    axis.title = ggplot2::element_blank()
  )
}
