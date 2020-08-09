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
#' @examples
#' data("abiotic_stresses")
#' DIANE::draw_heatmap(abiotic_stresses$normalized_counts, subset = abiotic_stresses$heat_DEGs,
#' title = "Log expression for DE genes under heat stress")
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
    
    if (sum(stringr::str_detect(rownames(data), paste0(sample_subset, collapse = '|'))) == 0) {
      stop("The required subset of genes was not found in expression data rownames")
    }
    
    if (is.null(conditions))
      conds <- colnames(data)
    else
      conds <-
        colnames(data)[str_split_fixed(colnames(data), '_', 2)[, 1] %in% conditions]
    
    if (log)
      data <- log(data + 1)
    if (profiles)
      data <- data / rowMeans(data)
    
    if (length(conds) == 0) {
      stop("The required conditions were not found in the expression data")
    }
    
    mat <- data[sample_subset, conds]
    
    sample <- stringr::str_split_fixed(colnames(mat), '_', 2) [, 1]
    samples <- data.frame(sample, row.names = colnames(mat))
    
    pheatmap::pheatmap(
      mat,
      color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlGnBu"))(100),
      annotation_col = samples,
      show_rownames = show_rownames,
      main = title,
      fontsize = 17
    )
  }



#' draw_distributions of count data
#'
#' @param data count data, as a numeric dataframe with condition names as columns and genes as rows
#' @param boxplot if TRUE, plot each sample as a boxplot, else, it is shown as a violin plot
#' @export
#' @examples
#' data("abiotic_stresses")
#' DIANE::draw_distributions(abiotic_stresses$normalized_counts, boxplot = FALSE)
draw_distributions <- function(data, boxplot = TRUE) {
  d <-
    reshape2::melt(log(data[sample(rownames(data),
                                   replace = FALSE,
                                   size = round(dim(data)[1] / 4, 0)), ] + 1))
  
  colnames(d)[c(length(colnames(d)) - 1, length(colnames(d)))] <-
    c("sample", "logCount")
  
  d$condition <- stringr::str_split_fixed(d$sample, "_", 2)[, 1]
  
    
  
  if (boxplot) {
    g <- ggplot2::ggplot(data = d, ggplot2::aes(x = sample, y = logCount))
      g <- g + ggplot2::geom_boxplot(
        alpha = 0.5,
        lwd = 1,
        ggplot2::aes(fill = condition),
        outlier.color = "black",
        outlier.alpha = 0.1
      )
  } else{
    g <- ggplot2::ggplot(data = d, ggplot2::aes(y = sample, x = logCount, color = condition)) +
      ggridges::geom_density_ridges(size = 2, fill = "#d6dbdf")
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
      axis.title = ggplot2::element_text(size = 24)
    )
  g 
}

#' Multi-dimensional scaling plot
#'
#' @param normalized.count data to plot for MDS
#' @importFrom limma plotMDS
#' @export
#' @examples
#' data("abiotic_stresses")
#' DIANE::draw_MDS(abiotic_stresses$normalized_counts)

draw_MDS <- function(normalized.count) {
  normalized.count <- normalized.count / rowMeans(normalized.count)
  mds <- limma::plotMDS(normalized.count, plot = FALSE)
  d <-
    data.frame(
      dim1 = mds$x,
      dim2 = mds$y,
      sample = names(mds$x),
      condition = str_split_fixed(names(mds$x), '_', 2)[, 1]
    )
  g <-
    ggplot2::ggplot(data = d, ggplot2::aes(x = dim1, y = dim2, color = condition, label = sample)) +
    ggplot2::geom_point(size = 6) + geom_text(
      color = "black",
      size = 6,
      alpha = 0.5,
      nudge_x = 0.07,
      nudge_y = 0.07
    )
  
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



#' Draw variables along principal PCA components
#'
#' @param data normalized expression data
#'
#' @export
#' @import ggplot2
#'
#' @examples
#' data("abiotic_stresses")
#' draw_PCA(abiotic_stresses$normalized_counts)
draw_PCA <- function(data) {
  # PCA computation
  data <- log(data + 2)
  data <- data / rowMeans(data)
  acp <-
    ade4::dudi.pca(
      data,
      center = TRUE,
      scale = TRUE,
      scannf = FALSE,
      nf = 4
    )
  
  acp$co$condition = stringr::str_split_fixed(rownames(acp$co), '_', 2)[, 1]
  
  scree <-
    data.frame(
      component = seq(1:length(acp$eig)),
      eigen.values = acp$eig,
      explained.variance = round(acp$eig / sum(acp$eig) *
                                   100, 2)
    )
  scree <- scree[1:min(nrow(scree), 4),]
  
  # Plots
  g1_2 <-
    ggplot(data = acp$co,
           aes(
             x = Comp1,
             y = Comp2,
             color = condition,
             label = condition
           )) + geom_text(
             color = "black",
             size = 6,
             alpha = 0.5,
             nudge_x = 0.07,
             nudge_y = 0.07
           ) +
    geom_point(size = 6, alpha = 0.7) + xlim(-1, 1) +
    ylim(-1, 1) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    theme(legend.position = "none") +
    ggtitle("Principal components 1 and 2") +
    xlab(paste("x-axis : Comp1 ", scree[1, "explained.variance"], "%")) +
    ylab(paste("y-axis : Comp2 ", scree[2, "explained.variance"], "%"))
  
  g2_3 <-
    ggplot(data = acp$co,
           aes(
             x = Comp2,
             y = Comp3,
             color = condition,
             label = condition
           )) + geom_text(
             color = "black",
             size = 6,
             alpha = 0.5,
             nudge_x = 0.07,
             nudge_y = 0.07
           ) +
    geom_point(size = 6, alpha = 0.7) + xlim(-1, 1) +
    ylim(-1, 1) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    theme(legend.position = "none") +
    ggtitle("Principal components 2 and 3") +
    xlab(paste("x-axis : Comp2 ", scree[2, "explained.variance"], "%")) +
    ylab(paste("y-axis : Comp3 ", scree[3, "explained.variance"], "%"))
  
  g3_4 <-
    ggplot(data = acp$co,
           aes(
             x = Comp3,
             y = Comp4,
             color = condition,
             label = condition
           )) + geom_text(
             color = "black",
             size = 6,
             alpha = 0.5,
             nudge_x = 0.07,
             nudge_y = 0.07
           ) +
    geom_point(size = 6, alpha = 0.7) + xlim(-1, 1) +
    ylim(-1, 1) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 18),
      legend.text.align = 1
    ) +
    ggtitle("Principal components 3 and 4") +
    xlab(paste("x-axis : Comp3 ", scree[3, "explained.variance"], "%")) +
    ylab(paste("y-axis : Comp4 ", scree[4, "explained.variance"], "%"))
  
  screeplot <- ggplot(scree,
                      aes(
                        y = explained.variance,
                        x = component,
                        fill = component,
                        label = paste(round(explained.variance, 1), '%')
                      )) +
    geom_bar(stat = "identity") + geom_text(size = 6,
                                            vjust = 1.6,
                                            color = "white") +
    ggtitle("PCA Screeplot") + theme(legend.position = "none")
  
  
  gridExtra::grid.arrange(g1_2, g2_3, g3_4, screeplot, ncol = 2)
}


#' Draw gene expression levels
#'
#' The normalized counts of the required genes in the required conditions are
#' shown. Please limit the number of input genes for readability reasons (up to 10 genes).
#'
#' @param data normalized expression dataframe, with genes as rownames and
#' conditions as colnames.
#' @param genes character vector of genes to be plotted (must be contained in
#' the rownames of data, entierly or as a substring)
#' @param conds conditions to be shown on expression levels (must be contained in
#' the column names of data before the _rep suffix). Default : all conditions.
#' @param gene.name.size size of the facet plot title font for each gene. Default : 12 
#' 
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' genes <- sample(abiotic_stresses$heat_DEGs, 4)
#' conditions <- c("C", "M", 'H', 'SH')
#' DIANE::draw_expression_levels(abiotic_stresses$normalized_counts,
#' genes = genes, conds = conditions)
draw_expression_levels <-
  function(data,
           genes,
           conds = unique(stringr::str_split_fixed(colnames(data), '_', 2)[,1]),
           gene.name.size = 12) {
    
    
    if (sum(stringr::str_detect(rownames(data), paste0(genes, collapse = '|'))) == 0) {
      stop("The required genes were not found in expression data rownames")
    }
    
    if (sum(stringr::str_detect(rownames(data), paste0(genes, collapse = '|'))) > 10) {
      stop("Please specify less than 10 genes, for readability reasons.")
    }
    
    conditions <-
      colnames(data)[stringr::str_split_fixed(colnames(data), '_', 2)[, 1] %in% conds]
    if (length(conditions) == 0) {
      stop("The required conditions were not found in the expression data")
    }
    
    data <- data.frame(data)
    data$gene <- rownames(data)
    
    d <-
      reshape2::melt(as.data.frame(data[stringr::str_detect(
        rownames(data), paste0(genes, collapse = '|')), c(conditions, 'gene')]))
    d$condition <- stringr::str_split_fixed(d$variable, '_', 2)[, 1]
    d$replicate <- stringr::str_split_fixed(d$variable, '_', 2)[, 2]
    ggplot(d,
           aes(
             x = condition,
             y = value,
             color = replicate
           )) +
      geom_point(size = 4, alpha = 0.8) + facet_wrap( ~ gene, scales = "free") +
      ggtitle("Normalized expression levels") +
      theme(
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = gene.name.size),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.text.y = element_text(size = 22, angle = 320),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 20)
      ) + xlab("") + ylab("Normalized counts")
  }
