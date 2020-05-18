


#' Fits a Poisson generalized linear model to a set of genes
#'
#' @description more to come
#'
#' @param normalized_counts normalized counts
#' @param genes genes belonging to a specific expression-based clusters
#' @param design experimental design as a dataframe
#' @param factors factors to use for the fit (defalut is
#' all the factors of the design)
#'
#' @return glm object
#' @export
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
    summary(glm)
    return(glm)
  }


#' Plots the value of a Poisson generalized linear model
#'
#' @param glm glm object returned by ```DIANE::fit_glm()```
#' @export
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