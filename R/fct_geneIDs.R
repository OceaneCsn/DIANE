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



#' Check compatibility between gene IDs and organism
#'
#' @param ids character vector of gene identifiers to be tested
#' @param organism organism, should be betwwen Arabidopsis thaliana,
#' Homo sapiens and Mus musculus.
#'
#' @return boolean, true if all of the gene IDs match the expected regex.
#' @export
#'
#' @examples 
#' data("abiotic_stresses")
#' check_IDs(rownames(abiotic_stresses$raw_counts),
#' organism = "Arabidopsis thaliana")
#' check_IDs(rownames(abiotic_stresses$raw_counts),
#' organism = "Homo sapiens")
check_IDs <- function(ids, organism){
  if(organism == "Arabidopsis thaliana")
    pattern = "AT[[:alnum:]]G[[:digit:]]{5}"
  
  if(organism == "Homo sapiens")
    pattern = "ENSG[0-9]{11}"
  
  if(organism == "Mus musculus")
    pattern = "ENSMUSG[0-9]{11}"
  
  matched <- sum(stringr::str_detect(ids, pattern = pattern))
  if( matched == length(ids))
    return(TRUE)
  else{
    if(matched > 0 & matched < length(ids)){
      warning("Some of the gene IDs do not match the expected regex")
      return(FALSE)
    }
    else{
      warning("None of the gene IDs match the expected regex")
      return(FALSE)
    }
  }
}