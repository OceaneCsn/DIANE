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
  # TODO : add a genes dataframe containing the annotations here, they will be kept furing all dea process in edgeR
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
#' 
#'
#' @return tags dataframe
#' @export

estimateDEGs <- function(fit, reference, perturbation){
  contrast <- ifelse(colnames(fit$design) == reference, -1, ifelse(colnames(fit$design) == perturbation, 1, 0))
  lrt <- edgeR::glmLRT(fit, contrast = contrast)
  top <- edgeR::topTags(lrt, p.value = 1, n = Inf)
  return(top)
}
  
  
#' Title
#'
#' @param tags tags returned bu estimateDEGs, function, that is to say topTags from edgeR
#' @param fdr pvalue for DEGs
#' @param MA TRUE : MA plot (LogFC depending on average log expression), or else "Vulcano" for
#' @param lfc absolute logFC threshold for DEGs
#' FDR depending on logFC.
#' @return plot object
#' @export
#'
plotDEGs <- function(tags, fdr = 0.01, lfc = 0, MA = TRUE){
  top <- tags$table
  top$isDE <- ifelse(top$FDR < fdr & abs(top$logFC) > lfc, TRUE, FALSE)
  if (MA) g <- ggplot2::ggplot(data = top, aes(x = logCPM, y = logFC, color = isDE)) + ggtitle("M-A plot")
  else g <- ggplot2::ggplot(data = top, aes(y = -log10(FDR), x = logFC, color = isDE)) + ggplot2::ggtitle("Vulcano plot")
  g <- g + ggplot2::geom_point(size = 1.2) + ggplot2::scale_color_manual(values = c("#999999", "#25BA40"))
  g + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18),
            title = element_text(size = 20, face = "bold"))
}
