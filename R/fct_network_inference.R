
#' Genie3 regulatory importances estimation
#'
#' @param normalized.count normalized expression matrix containing the regressors and target genes
#' @param regressors genes to be taken as regressors during the inference procedures (regulator genes)
#' @param targets genes to be included in the inferred network. Regressors can also be in the targets
#' @param nTrees Number of trees by Random Forest
#' @param nCores Number of CPU cores to use during the procedure
#'
#' @return matrix object

network_inference <- function(normalized.count, regressors, targets, nTrees=5000, nCores=1){
  regressors <- intersect(rownames(normalized.count), regressors)
  mat <- GENIE3::GENIE3(as.matrix(normalized.count), regulators = regressors, targets = targets, 
                treeMethod = "RF", K = "sqrt", nTrees = nTrees, nCores = nCores, verbose = T)
  
  return(mat)
}


#' Thresholds a non sparse adjascency matrix to return the strongest weights only
#'
#' @param mat matrix containing the importance values for each target and regulator
#' @param n_edges number of edges to keep in the final network
#'
#' @return igraph object

network_thresholding <- function(mat, n_edges){
  links <- GENIE3::getLinkList(mat, reportMax = n_edges)
  g <- igraph::graph_from_data_frame(links, directed = T)
  return(g)
}
