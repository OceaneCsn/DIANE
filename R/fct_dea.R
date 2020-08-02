#' estimateDispersion
#'
#' @description This function computes the common, trended, and tag-wise dispersions of 
#' the gene expression dataset. Once, the dispersion is estimated,
#' it fits a negative binomial model for each gene using default design : ~group+0.
#' In this configuration, the log average expression of each gene is approximated
#' by a linear combination of each of the conditions.
#' It returns a the value of a glmFit object from the package edgeR, containing the attributes
#' of the tcc object given, plus fitted values and GLM coefficients, among other indicators of the 
#' fitting procedure.
#'
#' @param tcc TCC-class object countaining counts, and norm factors, such as obtained after 
#' DIANE::normalize() and DIANE::filter_low_counts()
#' @param conditions if NULL (default), takes the split before the character '_' as condition names.
#' else, the conditions specified by the user, as a character vector. Its order should match the 
#' columns names of the expression matrix used to build the tcc object.
#' @importFrom stringr str_split_fixed
#' @importFrom edgeR DGEList estimateDisp glmFit
#' @importFrom stats model.matrix
#' @return glmFit
#' @export
#' 
#' @section # Note:
#' The return value is meant to be used directly as the parameter of
#' DIANE::estimateDEGs()
#' 
#' @examples
#' data("abiotic_stresses")
#' tcc_object <- DIANE::normalize(abiotic_stresses$raw_counts, 
#' abiotic_stresses$conditions, iteration = FALSE)
#' threshold = 10*length(abiotic_stresses$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
#' 
estimateDispersion <- function(tcc, conditions = NULL) {
  if (is.null(conditions)) {
    groups <- stringr::str_split_fixed(colnames(tcc$count), '_', 2)[, 1]
  }
  else{
    groups <- conditions
  }
  design = stats::model.matrix( ~ groups + 0)
  rownames(design) <- colnames(tcc$count)
  colnames(design) <-
    str_split_fixed(colnames(design), 'groups', 2)[, 2]


  dge <- edgeR::DGEList(
    counts = tcc$count,
    lib.size = tcc$norm.factors,
    norm.factors = colSums(tcc$count),
    group = groups,
    genes = rownames(tcc$count)
  )
  
  y <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmFit(y, design)
  return(fit)
}

#' estimateDEGs
#'
#' Given a glmFit, performs a log ratio test using edgeR function glmLRT and
#' returns the genes that are differentially expressed
#' between the reference and perturbation conditions.
#' log ratios are expressed as log(perturbation/reference)
#'
#' @param fit edgeR glmFit
#' @param reference condition being considered as teh reference for differential analysis
#' @param perturbation condition we compared to the reference for differential analysis
#' @param p.value numeric cutoff value for adjusted p-values. Only tags with adjusted p-values equal or 
#' lower than specified are returned
#' @param lfc minimal absolute log fold change required for a gene to be considered as 
#' differnetially expressed.
#' @return topTags object, which table element contains DEGs dataframe.
#' @export
#' @examples
#' data("abiotic_stresses")
#' tcc_object <- DIANE::normalize(abiotic_stresses$raw_counts, 
#' abiotic_stresses$conditions, iteration = FALSE)
#' threshold = 10*length(abiotic_stresses$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, 
#' conditions = abiotic_stresses$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 0.01)
#' DEGs <- topTags$table
estimateDEGs <- function(fit, reference, perturbation, p.value = 1, lfc = 0) {
  contrast <-
    ifelse(colnames(fit$design) == reference,
           -1,
           ifelse(colnames(fit$design) == perturbation, 1, 0))
  lrt <- edgeR::glmLRT(fit, contrast = contrast)
  top <- edgeR::topTags(lrt, p.value = p.value, n = Inf)
  top$table <- top$table[abs(top$table$logFC) > lfc,]
  return(top)
}


#' MA or volcano plot for Differential expression
#'
#' @param tags tags returned bu estimateDEGs, function, that is to say topTags from edgeR
#' @param fdr pvalue for DEGs
#' @param MA TRUE : MA plot (LogFC depending on average log expression), or else "Volcano" for FDR depending on logFC.
#' @param lfc absolute log fold change threshold for differentially expressed genes, default is 0.
#' @export
#' @examples
#' data("abiotic_stresses")
#' tcc_object <- DIANE::normalize(abiotic_stresses$raw_counts, abiotic_stresses$conditions, 
#' iteration = FALSE)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, 10*length(abiotic_stresses$conditions))
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
#' tags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 1)
#' DIANE::draw_DEGs(tags, fdr = 0.01, lfc = 1)
draw_DEGs <- function(tags,
                     fdr = 0.01,
                     lfc = 0,
                     MA = TRUE) {
  top <- tags$table
  top$isDE <-
    ifelse(top$FDR < fdr & abs(top$logFC) > lfc, TRUE, FALSE)
  if (MA)
    g <-
    ggplot2::ggplot(data = top, ggplot2::aes(x = logCPM, y = logFC, color = isDE)) + 
    ggplot2::ggtitle("M-A plot")
  else
    g <-
    ggplot2::ggplot(data = top, ggplot2::aes(
      y = -log10(FDR),
      x = logFC,
      color = isDE
    )) + ggplot2::ggtitle("Volcano plot")
  g <-
    g + ggplot2::geom_point(size = 1.5, alpha = 0.75) + 
    ggplot2::scale_color_manual(name = "DEG", values = c("#999999", "#25BA40"))
  g + ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 18),
    axis.text.y = ggplot2::element_text(size = 18),
    title = ggplot2::element_text(size = 20, face = "bold")
  )
}
