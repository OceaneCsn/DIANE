#' normalize raw count data, using the TCC package
#'
#' @param data raw counts to be normalized (data frame or matrix)
#' @param conditions condition of each column of the data argument
#' @param norm_method method used for normalisation, between tmm or deseq2
#' @param deg_method method used for deg detection, between edgeR ou deseq2
#' @param fdr pvalue threshold for adjusted pvalues for degs
#' @param iteration weather or not to perform a prior removal of DEGs (TRUE or FALSE)
#' 
#' @import TCC
#' @return a TCC-Class object containing the normalized data (filtering is highly recommended after normalization)
#' @export
#' @examples
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)


normalize <- function(data, conditions, norm_method = "tmm", deg_method = "edgeR", fdr = 0.01,
                      iteration = TRUE){
  tcc <- new("TCC", data, conditions)
  tcc <- TCC::calcNormFactors(tcc, norm.method = norm_method, test.method = deg_method, 
                              iteration = iteration, FDR = 0.01, floorPDEG = 0.05)
  return(tcc)
}




#' filter_low_counts
#'
#' @param tcc data to be filtered to remove low count genes
#' 
#' @param thr the sum of counts across all samples to be exceeded for a gene
#' @export
#' @return a TCC-Class object containing the filtered data
#' @examples
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)

filter_low_counts <- function(tcc, thr){
  tcc <- TCC::filterLowCountGenes(tcc, low.count = thr)
  return(tcc)
}

