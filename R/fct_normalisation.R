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


#' Detects weather or not the data contains splice variants 
#' instead of unique locus. Returns TRUE if all genes ids match 
#' the regex for splice variants
#'
#' @param gene_ids gene vector to be tested
#'
#' @return boolean
are_splice_variants <- function(gene_ids){
  return(sum(stringr::str_detect( gene_ids, pattern = "\\.[[:digit:]]+$")) == length(gene_ids))
}

#' Merge all splice variants of an expression dataset into 
#' unique locus, unaware of alternative splicing, by summing
#' all variants for the same gene
#'
#' @param data dataframe : expression data with splice variants as rownames
#'
#' @return dataframe with agregated rows
#' @export
#' @examples 
#' data("abiotic_stresses")
#' d <- aggregate_splice_variants(abiotic_stresses$normalized_counts)
aggregate_splice_variants <- function(data){
  if(are_splice_variants(rownames(data))){
    
    data <- data.frame(data)
    
    locus <- stringr::str_replace_all(rownames(data), 
                                      pattern = "\\.[[:digit:]]+$", "")
    data$locus <- locus
    data_locus <- aggregate(. ~ locus, data, sum)
    rownames(data_locus) <- unique(locus)
    return(data_locus[, colnames(data_locus) != "locus"])
  }
  else print("The input data did not contain splice variants.
             Or at least not enterely.")
}


#' Get the locus ids from splice variants ids
#'
#' @param gene_ids list of gene ids with splice variants information
#' @param unique boolean, weather or not to return unique locus vector
#' @return character vector
#' @export
#' @examples 
#' splice_variants <- rownames(abiotic_stresses$normalized_counts)[1:20]
#' get_locus(splice_variants)
get_locus <- function(gene_ids, unique = TRUE){
  if(unique){
    return(unique(
      stringr::str_replace_all(gene_ids, 
                               pattern = "\\.[[:digit:]]+$", "")))
  }
  else{
    return(
      stringr::str_replace_all(gene_ids, 
                               pattern = "\\.[[:digit:]]+$", ""))
  }
}

