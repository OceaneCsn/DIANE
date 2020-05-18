


#' Fits a Poisson generalized linear model to a set of genes
#'
#' @description The idea is to extract the importance and effect of each factor.
#' To do so, the expression of each gene is modelled as a Poisson distribution.
#' The log of its parameter (the expected value) is approximated by a linear
#' combination of the factors in the experiment. The coefficients associated to 
#' each factors are estimated to fit gene expression, and can be insightful to
#' characterize genes behavior in a particular cluster.
#' The model with interactions is considered. It your design in not a
#' complete corssed design, the interaction term will be null.
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
#' library(DIANE)
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' genes <- topTags$table$genes
#' clustering <- DIANE::run_coseq(conds = unique(demo_data_At$conditions), data = normalized_counts, genes = genes, K = 6:9)
#' genes_cluster <- DIANE::get_genes_in_cluster(clustering$membership, cluster = 3)
#' glm <- DIANE::fit_glm(normalized_counts, genes_cluster, demo_data_At$design)
#' summary(glm)

fit_glm <-
  function(normalized_counts,
           genes,
           design,
           factors = colnames(design)) {
    
    glmData <- melt(round(normalized_counts[genes, ], 0))
    glmData <- glmData[, 2:length(colnames(glmData))]
    colnames(glmData) <- c("Sample", "Counts")
    glmData$condition <- str_split_fixed(glmData$Sample, '_', 2)[, 1]
    for (factor in factors) {
      glmData[, factor] <- design[glmData$condition, factor]
    }
    formula <- paste("Counts ~ ", paste(factors, collapse = '*'))
    glm <-
      glm(formula , data = glmData, family = poisson(link = "log"))
    return(glm)
  }


#' Plots the value of a Poisson generalized linear model
#'
#' @param glm glm object returned by ```DIANE::fit_glm()```
#' @export
#' @examples 
#' library(DIANE)
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' genes <- topTags$table$genes
#' clustering <- DIANE::run_coseq(conds = unique(demo_data_At$conditions), data = normalized_counts, genes = genes, K = 6:9)
#' genes_cluster <- DIANE::get_genes_in_cluster(clustering$membership, cluster = 3)
#' glm <- DIANE::fit_glm(normalized_counts, genes_cluster, demo_data_At$design)
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