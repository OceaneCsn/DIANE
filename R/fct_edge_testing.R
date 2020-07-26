#' test_edges
#'
#' @description
#' We designed a method to perform statistical testing on TF-target gene pairs
#' from the GENIE3 approach.
#' The idea is to build a first biologically relevant network with the
#' strongest importances diven by a prior GENIE3 run, that would then be
#' refined by statistical testing.
#' Those tests are performed by the rfPermute package, providing empirical
#' pvalues on observed importance values from response variable permutations.
#'
#' @param mat matrix containing the importance values for each target and regulator
#' (preferably computed with GENIE3 and the OOB importance metric)
#' @param normalized_counts normalized expression data containing the genes present in
#' mat argument, and such as used for the first network inference step.
#' @param nGenes number of total genes in the network, union of thetarget genes,
#' and regulators
#' @param nRegulators number of regulators used for the network inference step
#' @param density approximate desired density, that will be used to build
#' a first network, which edges are the one to be statistically tested. Default
#' is 0.02. Biological networks are known to have densities (ratio of edges over
#' total possible edges in the graph) between 0.1 and 0.001.
#' The number of genes and regulators are needed to compute the density.
#' @param nTrees number of trees used for random forest importance computations
#' @param nShuffle number of times the response variable (target gene expression)
#' is randomized in order to estimate the null distribution of the
#' predictive variables (regulators) importances.
#' @param nCores Number of CPU cores to use during the procedure.
#' Default is the detected number of cores minus one.
#' @param verbose If set to TRUE, a feedback on the progress of the calculations
#' is given. Default: TRUE
#'
#' @return named list containing the edges pvalues, as well as graphics
#' intended to guide the choice of a pvalue threshold for the final network:
#' \itemize{
#' \item links: a dataframe containing the links of the network before testing,
#' as built from the user defined prior density. All edges are associated to their
#' pvalue and fdr-adjusted pvalue.
#' \item fdr_nEdges_curve : relation between the fdr threshold, and the final
#' number of edges in the final network
#' }
#' @export
#'
#' @examples
test_edges <-
  function(mat,
           normalized_counts,
           nGenes,
           nRegulators,
           density = 0.02,
           nTrees = 1000,
           nShuffle = 1000,
           nCores = ifelse(is.na(parallel::detectCores()),
                           1,
                           max(parallel::detectCores() - 1, 1)),
           verbose = TRUE) {
    if (nGenes < 0)
      stop("nGenes must be a positive integer")
    if (nRegulators < 0)
      stop("nRegulators must be a positive integer")
    if (density <= 0 | density > 1)
      stop("density must be strictly positive between 0 and 1, preferably less than 0.1")
    
    nEdges = round(density * (nGenes - 1) * nRegulators, 0)
    
    # ou, avec cette formule trouvee dans un papier nature com des GRN sur les bacteries
    #nEdges = 0.4*nGenes**(-0.78)*(nGenes -1)*nRegulators
    
    if (verbose)
      message(
        paste(
          "First network thresholding : Selecting",
          nEdges,
          "edges for a gloabl density of",
          density,
          "."
        )
      )
    
    links <- GENIE3::getLinkList(mat, reportMax = nEdges)
    
    # assign to each gene its regulators
    targets <- as.vector(unique(links$targetGene))

    target_to_TF <- list()
    for (t in targets) {
      target_to_TF[[t]] <- as.vector(links[links$targetGene == t,
                                           "regulatoryGene"])
    }
    # estimate pvalues
    dataT <- t(normalized_counts)
    
    doParallel::registerDoParallel(cores = nCores)
    if (verbose)
      message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
    "%dopar%" <- foreach::"%dopar%"
    suppressPackageStartupMessages(result.reg <-
                                     doRNG::"%dorng%"(foreach::foreach(
                                       target = targets, .combine = rbind
                                     ),
                                     {
                                       # remove target gene from input genes
                                       theseRegulatorNames <-
                                         target_to_TF[[target]]
                                       numRegulators <-
                                         length(theseRegulatorNames)
                                       
                                       y <-
                                         dataT[, target]
                                       
                                       x <-
                                         data.frame(dataT[, target_to_TF[[target]]])
                                       colnames(x) <-
                                         target_to_TF[[target]]
                                       
                                       y <-
                                         (y - mean(y)) / sd(y)
                                       rf <-
                                         rfPermute::rfPermute(
                                           x,
                                           y,
                                           ntree = nTrees,
                                           replace = FALSE,
                                           nodesize = 1,
                                           nrep = nShuffle,
                                           num.cores = 1
                                         )
                                       if (length(target_to_TF[[target]]) == 1) {
                                         pval <- rf$pval[, , "scaled"]
                                         names(pval) <-
                                           paste(names(pval), ".pval", sep = "")
                                         
                                         pvals <-
                                           setNames(pval["%IncMSE.pval"], target_to_TF[[target]])
                                       }
                                       else{
                                         pvals <- rfPermute::rp.importance(rf)[, "%IncMSE.pval"]
                                       }
                                       
                                       res <-
                                         data.frame("regulatoryGene" = theseRegulatorNames,
                                                    "targetGene" = rep(target, numRegulators))
                                       res$pval <-
                                         pvals[match(res$regulatoryGene, names(pvals))]
                                       res
                                     }))
    attr(result.reg, "rng") <- NULL
    # It contains the whole sequence of RNG seeds
    links <- result.reg
    links$fdr <- stats::p.adjust(links$pval, method = "fdr")
    
    d <- reshape2::melt(links[, c("pval", "fdr")])
    distr <-
      ggplot2::ggplot(d, ggplot2::aes(x = value, fill = variable)) +
      ggplot2::geom_density(alpha = 0.5) + ggplot2::xlim(0, 0.5)
    
    nedges <-
      sapply(
        X = seq(0.0, 0.1, by = 0.001),
        FUN = function(fdr) {
          return(sum(links$fdr < fdr))
        }
      )
    d <-
      data.frame(FDRs = seq(0.0, 0.1, by = 0.001), n_edges = nedges)
    curve <-
      ggplot2::ggplot(d, ggplot2::aes(x = FDRs, y = n_edges)) +
      ggplot2::geom_point(size = 3, color = "darkgreen") +
      ggplot2::ggtitle("Number of edges depending on the FDR threshold") +
      ggplot2::xlab("fdr")
    
    return(list(
      links = links,
      fdr_nEdges_curve = curve,
      pvalues_distributions =
        distr
    ))
  }


#' Create network from edges statistical tests
#'
#' @param links dataframe of the network edges and associated pvalues,
#' as in the links attribute of the \code{test_edges} method.
#' @param fdr threshold value such as all edges wit hadjusted pvalues lesser than
#' the argument will be discarded for the final network construction
#'
#' @return an oriented weighted network as an igraph object
#' @export
#'
#' @examples
network_from_tests <- function(links, fdr) {
  if (fdr <= 0 | fdr > 1)
    stop(
      "fdr threshold must be strictly positive between 0 and 1,
         standard values being 0.01, 0.05, or O.1"
    )
  
  message(paste((sum(links$fdr < fdr)), "edges kept in final network"))
  
  links <- links[links$fdr < fdr,]
  net <- igraph::graph_from_data_frame(links, directed = TRUE)
  return(net)
}