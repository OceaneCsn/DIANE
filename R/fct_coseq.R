

#' Mixture Model clustering
#' 
#' @description Performs expression based clustering on genes.
#' Uses the coseq package to fit either Poisson or Gaussian 
#' mixtures to genes clusters, estimating their multivariate distribution
#' parameters via an EM algorithm.
#' Several numbers of clusters can be tested, and evaluated in terms of
#' of likelihood to help the user in the decision
#'
#' @param conds Condition names to be used for clustering.
#' Must be a unique vector containing the conditions you want to consider
#' for gene clustering, without the replicate information 
#' (string before the underscore in sample names)
#' @param genes Genes used as an input for the clustering. They must be present 
#' in the row names of data.
#' @param data normalized counts with genes as rownames and samples as columns
#' @param K range of number of clusters to test.
#' @param model Model to use for mixture models : to choose between
#' Poisson or Normal.
#' @param transfo Transformation to apply to normalized counts
#' before modeling with "Normal" Mixture Models. 
#' It must be : “arcsin”, “logit”, “logMedianRef”, “profile”, “logclr”, 
#'  “clr”, “alr”, “ilr”, or “none”. For "Poisson",
#' no transformation will be used, this argument will be ignored.
#' @param seed seed for random state to ensure reproducible runs
#' @importFrom coseq coseq clusters
#'
#' @return Named list containing the coseq run result as "model", and the cluster membership 
#' for each gene as "membership".
#' @export
#' @examples
#' data("abiotic_stresses")
#' genes <- abiotic_stresses$heat_DEGs
#' clustering <- run_coseq(conds = unique(abiotic_stresses$conditions), 
#' data = abiotic_stresses$normalized_counts, genes = genes, K = 6:9)
run_coseq <- function(conds, genes, data, K = 6:12, transfo = "none",
                      model = "Poisson", seed = NULL) {

  conditions <- colnames(data)[
    stringr::str_split_fixed(colnames(data), '_',2)[,1] %in% conds]
  
  
  if(model == "Poisson")
    transfo = "none"
  
  if(!transfo %in% c("voom", "arcsin", "logit", 
                    "logMedianRef", "profile", 
                    "logclr", "clr", "alr", 
                    "ilr", "none"))
    stop("The required transformation is not known. Please refer
         to function documentation")
  
  if(!model %in% c("Poisson", "Normal"))
    stop("The required model is not known. Please refer
         to function documentation")
  
  groups <- stringr::str_split_fixed(conditions, '_', 2)[, 1]
  dataC <- round(data[genes, conditions], 0)

  clustering_run <-
    suppressMessages(coseq::coseq(
      dataC,
      conds = groups,
      K = K,
      model = model,
      iter = 5,
      transformation = transfo,
      normFactors = "none", # norm already done in DIANE
      parallel = FALSE, #to ensure reproducibility
      GaussianModel = "Gaussian_pk_Lk_Bk", #to avoid singular covariance matrices
      verbose = FALSE,
      seed = seed
    ))
  return(list(
    membership = coseq::clusters(clustering_run),
    model = clustering_run
  ))
}



#' Displays clustering results and quality
#'
#' @param clustering_run result of a coseq run
#' @param plot plot to display, either Integrated Complete Likelihood, or barplots
#' of the posterior probabilities for the clustering. Value must be "ICL" or "barplots".
#' @return plot describing the quality of the clustering process
#'
#' @importFrom coseq plot
#' @export
#' @examples
#' data("abiotic_stresses")
#' genes <- abiotic_stresses$heat_DEGs
#' clustering <- DIANE::run_coseq(conds = c("C", "H", "SH", "MH", "SMH"), 
#' data = abiotic_stresses$normalized_counts, genes = genes, K = 6:9)
#' DIANE::draw_coseq_run(clustering$model, plot = "barplots")
#' DIANE::draw_coseq_run(clustering$model, plot = "ICL")
draw_coseq_run <- function(clustering_run, plot = "ICL") {
  if (plot == "ICL")
    p <- coseq::plot(clustering_run, graphs = c("ICL"))
  if (plot == "barplots") {
    p <- coseq::plot(clustering_run, graphs = c("probapost_barplots"))
  }
  p
}


#' Get genes in a specific cluster
#'
#' @param membership named vector from run_coseq
#' @param cluster desired cluster number
#'
#' @return genes belonging to the desired cluster
#' @export
#' @examples 
#'data("abiotic_stresses")
#'get_genes_in_cluster(abiotic_stresses$heat_DEGs_coseq_membership, 3)
get_genes_in_cluster <- function(membership, cluster) {
  return(names(membership[membership == cluster]))
}


#' Draw expression profiles of a clustering
#'
#' @param data normalized counts with genes as rownames and samples as columns
#' @param membership membership item of the coseq object returned by 
#' \code{run_coseq()} object
#' @param expression if it is set to "profiles" (default), 
#' plots expression/sum(expression). if "counts", plots log(Counts+1)
#' @param k if NULL (default), plots all the clusters. Else, plot the clusters in the vetcor k.
#' @param conds conditions on which to display clustering profiles.
#' Must be a unique vector containing the conditions you want to consider
#' for gene clustering, without the replicate information 
#' (string before the underscore in sample names). Default is all the conditions.
#' @param nrow on how many rows display the cluster profiles if k is NULL
#'
#' @importFrom reshape2 melt
#' @importFrom stringr str_split_fixed
#' @export
#' @examples
#' data("abiotic_stresses")
#' DIANE::draw_profiles(data = abiotic_stresses$normalized_counts, 
#' membership = abiotic_stresses$heat_DEGs_coseq_membership,
#' conds = unique(abiotic_stresses$conditions)) 
#' DIANE::draw_profiles(data = abiotic_stresses$normalized_counts, 
#' membership = abiotic_stresses$heat_DEGs_coseq_membership, 
#' conds = unique(abiotic_stresses$conditions), k = 3) 
draw_profiles <-
  function(data,
           membership,
           conds = unique(stringr::str_split_fixed(colnames(data), '_', 2)[, 1]),
           expression = "profiles",
           k = NULL,
           nrow = 3) {
    clusters <- membership
    
    conditions <- colnames(data)[stringr::str_split_fixed(colnames(data), 
                                                          '_',2)[,1] %in% conds]
    
    
    data <- data[names(clusters), conditions]
    
    if (expression == "profiles") {
      profiles <- data.frame(data / rowSums(data))
      ylab <- "Normalized counts/Mean(Normalized counts)"
    }
    
    if (expression == "counts") {
      profiles <- data.frame(log(as.matrix(data) + 1))
      ylab <- "log(Normalized counts)"
    }
    
    profiles$gene <- rownames(profiles)
    d <- suppressMessages(reshape2::melt(profiles))
    d$group <- stringr::str_split_fixed(d$variable, '_', 2)[, 1]
    d$cluster <- clusters[match(d$gene, names(clusters))]
    d$geneRep <-
      paste0(d$gene, stringr::str_split_fixed(d$variable, '_', 2)[, 2])
    d$replicate <- stringr::str_split_fixed(d$variable, '_', 2)[, 2]
    
    if (is.null(k)) {
      g <-
        ggplot2::ggplot(data = d, ggplot2::aes(x = group, y = value))  +
        ggplot2::facet_wrap(~ cluster, scales = "free")
      leg <- "none"
    }
    else{
      k <- as.vector(k)
      g <-
        ggplot2::ggplot(data = d[d$cluster %in% k,], 
                        ggplot2::aes(x = group, y = value)) +
        ggplot2::facet_wrap(~ cluster, scales = "free")
      leg <- "none"
    }
    if (!is.null(k) && length(k) == 1) {
      g <-
        g + ggplot2::geom_line(
          alpha = 0.12,
          lwd = 0.9,
          ggplot2::aes(group = geneRep, color =replicate)
        )
      leg <- "right"
    }
    else {
      g <-
        g + ggplot2::geom_line(alpha = 0.08,
                               lwd = 0.9,
                               ggplot2::aes(group = geneRep, 
                                            color = as.factor(cluster)))
    }
    g <-
      g + ggplot2::geom_boxplot(
        alpha = 0.4,
        lwd = 1,
        color = "black",
        outlier.color = "black",
        outlier.alpha = 0.1
      ) +
      ggplot2::geom_jitter(width = 0.1, alpha = 0.0015) +
      ggplot2::ggtitle("Normalized expression profiles")
    
    g <-
      g + ggplot2::theme(
        plot.title = ggplot2::element_text(size = 22, face = "bold"),
        strip.text.x = ggplot2::element_text(size = 20),
        legend.position = leg, legend.text = ggplot2::element_text(size = 20),
        legend.title = ggplot2::element_text(size = 20, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 12, angle = 30),
        axis.text.x = ggplot2::element_text(
          size = 12,
          hjust = 0,
          angle = -50,
          colour = "grey50"
        ),
        legend.text.align = 1,
        axis.title = ggplot2::element_text(size = 19)
      ) + ggplot2::xlab("") + ggplot2::ylab(ylab)
    g
  }
