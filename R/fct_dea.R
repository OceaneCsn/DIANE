

#' estimateDispersion
#'
#' The function computes the common, trend, and tag-wise dispersions of the dataset.
#' It then fits a negative binomial model for each gene using default design : ~group+0 and GLMs.
#' 
#' @param tcc object countaining count, and norm factors
#' @importFrom stringr str_split_fixed
#' @importFrom edgeR DGEList estimateDisp glmFit
#' @return fit object, containing teh fitted models
#' @export
#'
estimateDispersion <- function(tcc){
  
  groups <- str_split_fixed(colnames(tcc$count), '_', 2)[,1]
  design = model.matrix(~groups+0)
  rownames(design) <- colnames(tcc$count)
  colnames(design) <- str_split_fixed(colnames(design), 'groups', 2)[,2]
  
  dge <- edgeR::DGEList(counts = tcc$count,
                 lib.size = tcc$norm.factors,
                 norm.factors = colSums(tcc$count), 
                 group = groups, genes = rownames(tcc$count))
  
  y <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmFit(y, design)
  return(fit)
}
  

#' estimateDEGs
#' 
#' Given a glmFit, performs a log ratio test and
#' returns the genes that are differentially expressed 
#' between the reference and perturbation conditions.
#' log ratios will be expressed as log(perturbation/reference)
#' 
#' @param fit edgeR glmFit
#' @param reference condition being considered as "normal"
#' @param perturbation condition we can to compare to the reference
#' @param fdr adjusted pvalue for DEGs detection during glmLRT
#'
#' @return top tags dataframe
#' @export

estimateDEGs <- function(fit, reference, perturbation, fdr = 0.01){
  contrast <- ifelse(colnames(fit$design) == reference, -1, ifelse(colnames(fit$design) == perturbation, 1, 0))
  print(contrast)
  print(colnames(fit$design))
  lrt <- edgeR::glmLRT(fit, contrast = contrast)
  print(fdr)
  top <- edgeR::topTags(lrt, p.value = fdr, n = Inf)
  return(top)
}
  
  
  
