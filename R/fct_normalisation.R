#' Normalize raw count data
#' 
#' @description The function corrects for different sequencing depths bewteen samples.
#' It relies on the TCC package, to build a TCC-class object containing the raw counts 
#' and the conditions for each sample. The function calcNormFactors is then applied, 
#' and uses the method chosen by the user. It can, weather or not, proceed to a first step or 
#' removing potentially differentially expressed genes to have less biased normalisation
#' factors in the second normalization step. It returns a TCC object,
#' with an element norm_factors containing the computed normalization factors. 
#'
#' @param data raw counts to be normalized (data frame or matrix)
#' @param conditions condition of each column of the data argument
#' @param norm_method method used for normalization, between tmm or deseq2
#' @param deg_method method used for deg detection, between edgeR ou deseq2
#' @param fdr pvalue threshold for adjusted pvalues for degs
#' @param iteration weather or not to perform a prior removal of DEGs (TRUE or FALSE)
#' @details # Warning
#' Filtering is highly recommended after normalization, consider using the DIANE::filter_low_counts
#' function just after.
#' @details # Note
#' You can get the normalized expression matrix with TCC::getNormalizedData(tcc),
#' tcc being the result of DIANE::normalize or DIANE::filter_low_counts
#' @importFrom TCC calcNormFactors TCC
#' @return a TCC-Class object
#' @export
#' @examples
#' data("abiotic_stresses")
#' tcc_object <- DIANE::normalize(abiotic_stresses$raw_counts, 
#' abiotic_stresses$conditions, iteration = FALSE)
normalize <- function(data, conditions, norm_method = "tmm", deg_method = "edgeR", fdr = 0.01,
                      iteration = TRUE){
  tcc <- TCC::TCC(count =  data, group = conditions)
  tcc <- TCC::calcNormFactors(tcc, norm.method = norm_method, test.method = deg_method, 
                              iteration = iteration, FDR = 0.01, floorPDEG = 0.05)
  return(tcc)
}

#' Remove low count genes
#' 
#' @description Removes genes having a sum of counts accross all sample
#' inferior to the specified threshold. It returns un aupdated TCC object,
#' which count element contains the filtered expression matrix.
#'
#' @param tcc data to be filtered to remove low count genes
#' @importFrom TCC filterLowCountGenes
#' @param thr the sum of counts across all samples to be exceeded for a gene
#' @export
#' @return a TCC-Class object
#' @examples
#' data("abiotic_stresses")
#' tcc_object <- DIANE::normalize(abiotic_stresses$raw_counts, 
#' abiotic_stresses$conditions, iteration = FALSE)
#' threshold = 10*length(abiotic_stresses$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
filter_low_counts <- function(tcc, thr){
  tcc <- TCC::filterLowCountGenes(tcc, low.count = thr)
  return(tcc)
}

