library(pheatmap)
library(stringr)




#' Draw heatmap
#'
#' @param normalized.count 
#' @param subset 
#' @param show_rownames 
#' @param title 
#' @importFrom pheatmap pheatmap
#' @importFrom stringr str_split_fixed
#'
#' @return plot of the heatmap

draw_heatmap <- function(normalized.count, subset = "random", show_rownames = F, 
                         title = "Random preview heatmap",
                         log = TRUE){
  
  if (subset == "random") sample_subset <- sample(rownames(normalized.count), size = 100)
  else sample_subset <- subset
  
  if(log) mat <- log(normalized.count[sample_subset,] + 1)
  else mat <- normalized.count[sample_subset,]
  
  sample <- stringr::str_split_fixed(colnames(mat), '_',2) [,1]
  samples <- data.frame(sample, row.names = colnames(mat))
  pheatmap::pheatmap(mat, annotation_col = samples, show_rownames = show_rownames, 
           main = title, fontsize = 17)
}
