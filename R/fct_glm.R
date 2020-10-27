
is.unique <- function(vector){
  if (length(unique(vector)) == 1) return(TRUE)
  return(FALSE)
}

#' Returns the factors that 
#' don't stay at the same level between 
#' between the specified conditions
#'
#' @param conditions character vetcor of conditions (must be in design rownames)
#' @param design design dataframe
#'
#' @return character vector of "active", or "perturbated" factors
get_factors_from_conditions <- function(conditions, design){
  if (is.null(design)) return(conditions)
  
  factors <- colnames(design)
  
  uniques <- mapply(design[conditions, ], FUN = is.unique)
  
  undesired <- names(uniques[uniques == TRUE])
  
  return(factors[!factors %in% undesired])
}

#' Fits a Poisson generalized linear model to a set of genes
#'
#' @description The idea is to extract the importance and effect of each factor.
#' To do so, the expression of each gene is modeled as a Poisson distribution.
#' The log of its parameter (the expected value) is approximated by a linear
#' combination of the factors in the experiment. The coefficients associated to 
#' each factors are estimated to fit gene expression, and can be insightful to
#' characterize genes behavior in a particular cluster.
#' The model with interactions is considered. It your design in not a
#' complete crossed design, the interaction term will be null.
#'
#' @param normalized_counts normalized counts
#' @param genes genes belonging to a specific expression-based clusters
#' @param design experimental design as a dataframe
#' @param factors factors to use for the fit (defalut is
#' all the factors of the design)
#' @section Note :
#' Note that we can only apply a glm fit to a set of genes that have very close expression 
#' profiles accros conditions, else we would have to introduce a new variable related to the genes
#' themselves.
#' @return glm object
#' @export
#' @examples 
#' data("abiotic_stresses")
#' genes_cluster <- DIANE::get_genes_in_cluster(
#' abiotic_stresses$heat_DEGs_coseq_membership, cluster = 3)
#' glm <- DIANE::fit_glm(abiotic_stresses$normalized_counts, genes_cluster, 
#' abiotic_stresses$design)
#' summary(glm)
fit_glm <-
  function(normalized_counts,
           genes,
           design,
           factors = colnames(design)) {
    normalized_counts <- as.matrix(normalized_counts)
    glmData <- reshape2::melt(round(normalized_counts[genes, ], 0))
    glmData <- glmData[, 2:length(colnames(glmData))]
    colnames(glmData) <- c("Sample", "Counts")
    glmData$condition <- stringr::str_split_fixed(glmData$Sample, '_', 2)[, 1]
    for (factor in factors) {
      glmData[, factor] <- design[glmData$condition, factor]
    }
    formula <- paste("Counts ~ ", paste(factors, collapse = '*'))
    glm <-
      glm(formula , data = glmData, family = poisson(link = "log"))
    return(glm)
  }


#' Plots the coefficients value of a Poisson generalized linear model 
#'
#' @param glm glm object returned by \code{DIANE::fit_glm()}
#' @export
#' @examples 
#' data("abiotic_stresses")
#' genes_cluster <- DIANE::get_genes_in_cluster(
#' abiotic_stresses$heat_DEGs_coseq_membership, cluster = 3)
#' glm <- DIANE::fit_glm(abiotic_stresses$normalized_counts, genes_cluster, 
#' abiotic_stresses$design)
#' draw_glm(glm)
draw_glm <- function(glm) {
  #remove the intercept
  coefs <- glm$coefficients[2:length(glm$coefficients)]
  
  #removes the perturbation level from the factor names
  d <-
    data.frame(Coefficient = stringr::str_remove_all(as.character(names(coefs)), "1"),
               Values = coefs)
  ggplot2::ggplot(d, ggplot2::aes(x = Coefficient, y = Values, fill = Coefficient)) + 
    ggplot2::geom_bar(color = 'black',
             stat = "identity",
             alpha = 0.4) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 25, angle = 320),
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 22, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 20)
    ) +
    ggplot2::ggtitle("Generalized linear model's coefficients values")
}