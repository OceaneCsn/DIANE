#' Estimates gene dispersions prior to differential expression testing
#'
#' @description This function computes the common, trended, and tag-wise dispersions of 
#' the genes in an expression dataset. Once the dispersion is estimated,
#' it fits a negative binomial model for each gene using default design : ~group+0.
#' In this configuration, the log average expression of each gene is approximated
#' by a linear combination of each of the conditions.
#' It returns a the value of a glmFit object from the package edgeR, containing the attributes
#' of the tcc object given, plus fitted values and GLM coefficients, among other indicators of the 
#' fitting procedure.
#'
#' @param tcc TCC-class object containing counts, and norm factors, such as obtained after 
#' \code{DIANE::normalize()} and \code{DIANE::filter_low_counts()} functions
#' @param conditions if NULL (default), takes the split before the character '_' in 
#' your expression data column names as condition names.
#' (This is how condition names should be formatted to be used in DIANE.)
#' Else, the conditions are specified by the user, as a character vector. Its order should match the 
#' columns names of the expression matrix used to build the tcc object.
#' @importFrom stringr str_split_fixed
#' @importFrom edgeR DGEList estimateDisp glmFit
#' @importFrom stats model.matrix
#' @return glmFit object from edgeR, which is the resulting model estimation
#' @export
#' 
#' @section Note:
#' The return value is meant to be used directly as the parameter of
#' \code{DIANE::estimateDEGs()}
#' 
#' @examples
#' data("abiotic_stresses")
#' tcc_object <- DIANE::normalize(abiotic_stresses$raw_counts, 
#' abiotic_stresses$conditions, iteration = FALSE)
#' threshold = 10*length(abiotic_stresses$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
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
    stringr::str_split_fixed(colnames(design), 'groups', 2)[, 2]


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

#' Estimates Differentially Expression Genes. Legacy version.
#'
#' @description Given a glmFit model, performs log ratio tests using edgeR 
#' function glmLRT and
#' returns the genes that are differentially expressed
#' between the reference and perturbation conditions.
#' log ratios are expressed as log(perturbation/reference)
#'
#' @param fit edgeR glmFit
#' @param reference condition being considered as the reference for differential analysis.
#' It should corresponds to a condition name, e.g. the string before the underscore and 
#' replicate number in your sample names.
#' @param perturbation condition we compared to the reference for differential analysis.
#' It should corresponds to a condition name, e.g. the string before the underscore and 
#' replicate number in your sample names.
#' @param p.value numeric cutoff value for adjusted p-values. Only tags with adjusted p-values equal or 
#' lower than specified are returned
#' @param lfc minimal absolute log fold change required for a gene to be considered as 
#' differentially expressed.
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
#' head(DEGs)
estimateDEGs_legacy <- function(fit, reference, perturbation, p.value = 1, lfc = 0) {
  contrast <-
    ifelse(colnames(fit$design) == reference,
           -1,
           ifelse(colnames(fit$design) == perturbation, 1, 0))
  lrt <- edgeR::glmLRT(fit, contrast = contrast)
  top <- edgeR::topTags(lrt, p.value = p.value, n = Inf)
  top$table <- top$table[abs(top$table$logFC) > lfc,]
  return(top)
}

#' Estimates Differentially Expression Genes
#'
#' @description Given a glmFit model, performs log ratio tests using edgeR 
#' function glmLRT and
#' returns the genes that are differentially expressed
#' between the reference and perturbation conditions.
#' log ratios are expressed as log(perturbation/reference)
#'
#' @param fit edgeR glmFit
#' @param reference condition (or a vector of conditions) being considered as the reference for differential analysis.
#' It should corresponds to a condition name, e.g. the string before the underscore and 
#' replicate number in your sample names.
#' @param perturbation condition (or a vector of conditions) we compared to the reference for differential analysis.
#' It should corresponds to a condition name, e.g. the string before the underscore and 
#' replicate number in your sample names.
#' @param p.value numeric cutoff value for adjusted p-values. Only tags with adjusted p-values equal or 
#' lower than specified are returned
#' @param lfc minimal absolute log fold change required for a gene to be considered as 
#' differentially expressed.
#' @param systematic_orientation If multiple conditions are indicated in reference or perturbation, compute all
#' comparisons and return only genes that have the same orientation in all comparisons. If False (default), the mean
#' of conditions is used.
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
#' head(DEGs)
#' 
#' # For multiple conditions
#' topTags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = c("H","M","S"), p.value = 0.01)
#' DEGs <- topTags$table
#' head(DEGs)
estimateDEGs <- function(fit, reference, perturbation, p.value = 1, lfc = 0, systematic_orientation = FALSE, diagnostic_plot = FALSE) {
  contrast <-  create_versus_design(colnames(fit$design), reference, perturbation)
  print(contrast)
  lrt <- edgeR::glmLRT(fit, contrast = contrast)
  top <- edgeR::topTags(lrt, p.value = 1, n = Inf)
  
  top$table <- top$table[abs(top$table$logFC) > lfc & top$table$FDR <= p.value,]
  
  ###Systematic orientation case
  if (systematic_orientation == TRUE &&
      length(reference) + length(perturbation) > 2) { ###Only usefull for multiple condition comparison.
    all_vs_all_deg <- list()
    for (ref in colnames(fit$design)[contrast > 0]) {
      ###Reference
      for (pert in colnames(fit$design)[contrast < 0]) {
        ###vs perturbation
        print(paste0(ref," vs ",pert))
        targeted_contrast <-
          ifelse(colnames(fit$design) == ref,
                 -1,
                 ifelse(colnames(fit$design) == pert, 1, 0))
        targeted_lrt <-
          edgeR::glmLRT(fit, contrast = targeted_contrast)
        targeted_top <-
          edgeR::topTags(targeted_lrt, p.value = p.value, n = Inf) ###TODO : change pvalue cutoff for lest restictive one.
        all_vs_all_deg[[paste0(ref, "_vs_", pert)]] <-
          targeted_top$table[abs(targeted_top$table$logFC) > lfc, ]
      }
    }
    
    ###Remove genes that are not DE in each condition.
    DE_genes <- top$table$genes
    upreg = DE_genes
    downreg = DE_genes
    for (comp in names(all_vs_all_deg)) {
      print(comp)
      upreg <-
        upreg[upreg %in% rownames(all_vs_all_deg[[comp]][which(all_vs_all_deg[[comp]][2] >= 1), ])]
      downreg  <-
        downreg[downreg %in% rownames(all_vs_all_deg[[comp]][which(all_vs_all_deg[[comp]][2] <= -1), ])]
    }
    top$table <- top$table[c(upreg, downreg), ]
    
  }
  
  return(top)
}

#' Draw raw pvalue histogram
#'
#' @param tags returned by estimateDEGs, function, that is to say topTags from edgeR, 
#' used with \code{p.value = 1}
#' @param bins number of bar to display
#'
#' @import ggplot2
#' 
#' @return ggplot object corresponding to raw pvalue histogram.
#' @examples
#' data("abiotic_stresses")
#' tcc_object <- DIANE::normalize(abiotic_stresses$raw_counts, abiotic_stresses$conditions, 
#' iteration = FALSE)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, 10*length(abiotic_stresses$conditions))
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
#' tags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 1)
#' DIANE::draw_raw_pvalue_histogram(tags)
draw_raw_pvalue_histogram <- function(tags, bins=100){
  
  diag_plot <- ggplot2::ggplot(data=tags$table, ggplot2::aes(x=tags$table$PValue)) + ggplot2::theme_classic() +
    ggplot2::geom_histogram(breaks=seq(0, 1, by=1/bins), 
                            fill="#92D9A2", alpha = 1) +
    ggplot2::labs(title="Histogram of raw pvalues", x="Pvalue", y="Count") + 
    ggplot2::xlim(c(0,1))
  return(diag_plot)
}

#' MA or volcano plot for differential expression results
#'
#' @param tags tags returned by estimateDEGs, function, that is to say topTags from edgeR, 
#' used with \code{p.value = 1}
#' @param fdr pvalue for DEGs detection
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


#' Venn diagram
#' 
#' @description Draws a Venn diagram between lists of DEGs.
#' Possible numbers of lists handled by the function are 2, 3 or 4.
#'
#' @param gene_list named list of genes
#'
#' @export
#' @examples 
#' data("abiotic_stresses")
#' tcc_object <- DIANE::normalize(abiotic_stresses$raw_counts, abiotic_stresses$conditions, 
#' iteration = FALSE)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, 10*length(abiotic_stresses$conditions))
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
#' genes_heat <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H",
#'  p.value = 0.05)$table$genes
#' genes_heat_mannitol <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "MH", 
#' p.value = 0.05)$table$genes
#' draw_venn(list("C - H" = genes_heat, "C - MH" = genes_heat_mannitol))
draw_venn <- function(gene_list){
  if (length(gene_list) < 2 | length(gene_list)>4)
    stop("The number of gene lists must be between 2 and 4 to be shown in Venn diagram")
  venn <- ggVennDiagram::ggVennDiagram(gene_list, color = "#888888")
  suppressMessages(venn + ggplot2::scale_fill_gradient(low = "#EEEEEE", high = "#61BF45"))
}

#' create_versus_design
#' 
#' @description Given a vector containing an experiment design and two other vectors
#' containg reference and comparison condition, return a design vector for 
#' equilibrate comparison in edgeR
#'
#'
#' @param experiment_design A vector containing all the conditions of a studied dataset.
#' @param reference_point All the points in experiment_design to consider as reference
#' @param comparison_point All the points in experiment_design to consider as comparison points.
#'
#' @noRd
#' @return
#'
#' @examples create_versus_design(c("reference","condition1","condition2","condition3"), c("condition1", "condition2", "condition3"), "reference")
create_versus_design <- function(experiment_design, reference_point, comparison_point){
  experiment_design = unique(experiment_design) ###We just need the name of the different expriment design. We dont care about the replicates or so.
  reference_point = unique(reference_point)
  comparison_point = unique(comparison_point)
  
  ###Check that all the specified points are in the experiment design
  if (!all(c(reference_point, comparison_point) %in% experiment_design)) {
    stop(
      paste0(
        c(reference_point, comparison_point)[!c(reference_point, comparison_point) %in% experiment_design],
        "not in design (experiment_design). Please check your input condition list."
      )
    )
  }
  ###Check that there is no common reference and perturbation.
  if (any(reference_point %in% comparison_point) | any(comparison_point %in% reference_point)){
    stop("A condition cannot be used as reference and comparison point. Please check your comparison vector.")
  }
  
  comparison_matrix=rep(0, length(experiment_design))
  
  if(length(reference_point) > length(comparison_point)){
    weight_comparison=1
    weight_reference=length(comparison_point)/length(reference_point)
  } else if (length(reference_point) < length(comparison_point)) {
    weight_reference=1
    weight_comparison=length(reference_point)/length(comparison_point)
  } else {
    weight_comparison=1
    weight_reference=1
    
  }
  
  comparison_matrix[which(experiment_design %in% reference_point)] <- -weight_reference
  comparison_matrix[which(experiment_design %in% comparison_point)] <- weight_comparison
  
  return(comparison_matrix)
}
