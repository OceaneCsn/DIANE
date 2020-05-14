#' estimateDispersion
#'
#' The function computes the common, trend, and tag-wise dispersions of the dataset.
#' It then fits a negative binomial model for each gene using default design : ~group+0 and GLMs.
#'
#' @param tcc object countaining count, and norm factors
#' @param conditions if NULL (default, takes the split before the character '_' as condition names)
#' else, the conditions regardless of biological replicates, as a character vector. Its order should match the 
#' columns names of the expression matrix used to build the tcc object.
#' @importFrom stringr str_split_fixed
#' @importFrom edgeR DGEList estimateDisp glmFit
#' @importFrom stats model.matrix
#' @return fit object, containing teh fitted models
#' @export
#' 
#' @examples
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' 
estimateDispersion <- function(tcc, conditions = NULL) {
  if (is.null(conditions)) {
    groups <- stringr::str_split_fixed(colnames(tcc$count), '_', 2)[, 1]
  }
  else{
    groups <- conditions
  }
  print(groups)
  design = stats::model.matrix( ~ groups + 0)
  rownames(design) <- colnames(tcc$count)
  colnames(design) <-
    str_split_fixed(colnames(design), 'groups', 2)[, 2]
  print(design)
  # TODO : add a genes dataframe containing the annotations here, they will 
  # be kept furing all dea process in edgeR
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
#' @param p.value numeric cutoff value for adjusted p-values. Only tags with adjusted p-values equal or lower than specified are returned
#'
#' @return topTags object, which table element contains DEGs dataframe
#' @export
#' @examples
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' DEGs <- topTags$table

estimateDEGs <- function(fit, reference, perturbation, p.value = 1) {
  contrast <-
    ifelse(colnames(fit$design) == reference,
           -1,
           ifelse(colnames(fit$design) == perturbation, 1, 0))
  lrt <- edgeR::glmLRT(fit, contrast = contrast)
  top <- edgeR::topTags(lrt, p.value = p.value, n = Inf)
  return(top)
}


#' plotDEGs : MA or volcano plot for DEGs
#'
#' @param tags tags returned bu estimateDEGs, function, that is to say topTags from edgeR
#' @param fdr pvalue for DEGs
#' @param MA TRUE : MA plot (LogFC depending on average log expression), or else "Volcano" for FDR depending on logFC.
#' @param lfc absolute log fold change threshold for differentially expressed genes, default is 0.
#' @return plot object
#' @export
#' @examples
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, 10*length(demo_data_At$conditions))
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' tags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 1)
#' DIANE::plotDEGs(tags, fdr = 0.01, lfc = 1)

plotDEGs <- function(tags,
                     fdr = 0.01,
                     lfc = 0,
                     MA = TRUE) {
  top <- tags$table
  top$isDE <-
    ifelse(top$FDR < fdr & abs(top$logFC) > lfc, TRUE, FALSE)
  if (MA)
    g <-
    ggplot2::ggplot(data = top, aes(x = logCPM, y = logFC, color = isDE)) + 
    ggplot2::ggtitle("M-A plot")
  else
    g <-
    ggplot2::ggplot(data = top, aes(
      y = -log10(FDR),
      x = logFC,
      color = isDE
    )) + ggplot2::ggtitle("Vulcano plot")
  g <-
    g + ggplot2::geom_point(size = 1.2) + 
    ggplot2::scale_color_manual(values = c("#999999", "#25BA40"))
  g + theme(
    axis.text.x = ggplot2::element_text(size = 18),
    axis.text.y = ggplot2::element_text(size = 18),
    title = ggplot2::element_text(size = 20, face = "bold")
  )
}
