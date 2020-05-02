
library(TCC)
library(HTSFilter)

#' normalize raw count data
#'
#' @param data raw counts to be normalized (data frame or matrix)
#' @param conditions condition of each column of the data argument
#' @param norm_method method used for normalisation, between tmm or deseq2
#' @param deg_method method used for deg detection, between edgeR ou deseq2
#' @param fdr pvalue threshold for adjusted pvalues for degs
#' @param iteration weather or not to perform a prior removal of DEGs
#' 
#' @import TCC
#' @import HTSFilter
#' @return the normalized data
#' @export

normalize <- function(data, conditions, norm_method = "tmm", deg_method = "edgeR", fdr = 0.01,
                      iteration = TRUE){
  tcc <- new("TCC", data, conditions)
  tcc <- TCC::calcNormFactors(tcc, norm.method = norm_method, test.method = deg_method, 
                              iteration = iteration, FDR = 0.01, floorPDEG = 0.05)
  norm_data <- TCC::getNormalizedData(tcc)
  return(list("normalized.counts" = norm_data, "norm_factors" = tcc$norm.factors))
}




#' filter_sum
#'
#' @param data data to be filtered to remove low count genes
#' 
#' @param thr the sum of counts across all samples to be exceeded for a gene
#' @export
#' @return the filtered data

filter_sum <- function(data, thr){
  return(data[rowSums(data) > thr,])
}


#' filter_hts
#'
#' @param data data to be filtered to remove low count genes
#' 
#' @param conditions the sum of counts across all samples to be exceeded for a gene
#' @export
#' @return the filtered data

filter_hts <- function(data, conditions){
  filter <- HTSFilter::HTSFilter(data, conditions, normalization = 'none')
  return(filter$filteredData)
}


