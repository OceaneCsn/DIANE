

#' run_coseq
#'
#' @param conds Condition names to be used for clustering
#' @param genes Genes used as an input for the clustering
#' @param data normalized counts (MUST be normalized!)
#' @param K number of clusters range
#'
#' @importFrom coseq coseq clusters
#'
#' @return named list containing the coseq run result as "model", and the cluster membership for each gene as "membership"
#' @export
#' @examples
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' genes <- topTags$table$genes
#' clustering <- DIANE::run_coseq(conds = unique(demo_data_At$conditions), data = normalized_counts, genes = genes, K = 6:9)
run_coseq <- function(conds, genes, data, K = 6:12) {
  conditions <- unique(grep(paste(conds, collapse = "|"),
                            colnames(data), value = TRUE))
  
  groups <- str_split_fixed(conditions, '_', 2)[, 1]
  dataC <- round(data[genes, conditions], 0)
  run_pois <-
    coseq::coseq(
      dataC,
      conds = groups,
      K = K,
      model = "Poisson",
      iter = 5,
      transformation = "none",
      normFactors = "none",
      parallel = TRUE
    )
  return(list(
    membership = coseq::clusters(run_pois),
    model = run_pois
  ))
}



#' draw_coseq_run : displays the indications of clustering
#'
#' @param run_pois result of a coseq run
#' @param plot plot to display, eather integrated Complete Likelihood, or barplots
#' of the posterior probabilities for the clustering. Value must be "ICL" or "barplots".
#' @return plot describing the quality of the clustering process
#'
#' @importFrom coseq plot
#' @export
#' @examples
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' genes <- topTags$table$genes
#' clustering <- DIANE::run_coseq(conds = unique(demo_data_At$conditions), data = normalized_counts, genes = genes, K = 6:9)
#' DIANE::draw_coseq_run(clustering$model, plot = "barplots")
#' DIANE::draw_coseq_run(clustering$model, plot = "ICL")
draw_coseq_run <- function(run_pois, plot = "ICL") {
  if (plot == "ICL")
    p <- coseq::plot(run_pois, graphs = c("ICL"))
  if (plot == "barplots"){
    p <- coseq::plot(run_pois, graphs = c("probapost_barplots"))
  }
  p
}


#' get_genes_in_cluster
#'
#' @param membership named vector from run_coseq
#' @param cluster desired cluster number
#'
#' @return genes belonging to the desired cluster
#' @export
#'
get_genes_in_cluster <- function(membership, cluster) {
  return(names(membership[membership == cluster]))
}


#' Draw profiles of a clustering
#'
#' @param data normalized counts
#' @param membership membership item of the coseq object returned by run_coseq object
#' @param expression if it is set to "profiles" (default), plots expression/sum(expression). if "counts", plots log(Counts+1)
#' @param k if NULL (default), plot all the clusters. Else, plot the clusters in the vetcor k.
#' @param conds conditions on which to perform clustering, ingoring the others
#' @param nrow on how many rows display the cluster profiles if k is NULL
#'
#' @importFrom reshape2 melt
#' @importFrom stringr str_split_fixed
#' @import ggplot2
#' @return ggplot
#' @export
#' @examples
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' genes <- topTags$table$genes
#' clustering <- DIANE::run_coseq(conds = unique(demo_data_At$conditions), data = normalized_counts, genes = genes, K = 6:9)
#' DIANE::draw_profiles(data = normalized_counts, clustering$membership, conds = unique(demo_data_At$conditions)) 
#' DIANE::draw_profiles(data = normalized_counts, clustering$membership, conds = unique(demo_data_At$conditions), k = 3) 
draw_profiles <-
  function(data,
           membership,
           conds,
           expression = "profiles",
           k = NULL,
           nrow = 3) {
    clusters <- membership
    
    conditions <- unique(grep(paste(conds, collapse = "|"),
                              colnames(data), value = TRUE))
    
    
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
    d <- reshape2::melt(profiles)
    d$group <- stringr::str_split_fixed(d$variable, '_', 2)[, 1]
    d$cluster <- clusters[match(d$gene, names(clusters))]
    d$geneRep <-
      paste0(d$gene, stringr::str_split_fixed(d$variable, '_', 2)[, 2])
    
    
    if (is.null(k)) {
      g <-
        ggplot2::ggplot(data = d, aes(x = group, y = value))  +
        ggplot2::facet_wrap(~ cluster, scales = "free")
    }
    else{
      k <- as.vector(k)
      g <-
        ggplot2::ggplot(data = d[d$cluster %in% k,], aes(x = group, y = value)) +
        ggplot2::facet_wrap(~ cluster, scales = "free")
    }
    if (!is.null(k) && length(k) == 1) {
      g <-
        g + ggplot2::geom_line(
          alpha = 0.12,
          lwd = 0.9,
          color = "darkgreen",
          aes(group = geneRep)
        )
    }
    else {
      g <-
        g + ggplot2::geom_line(alpha = 0.08,
                               lwd = 0.9,
                               aes(group = geneRep, color = as.factor(cluster)))
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
      ggtitle("Expression profiles of the clusters")
    
    g <-
      g + ggplot2::theme(
        plot.title = element_text(size = 22, face = "bold"),
        strip.text.x = element_text(size = 20),
        legend.position = "none",
        axis.text.y = element_text(size = 18, angle = 30),
        axis.text.x = element_text(
          size = 12,
          hjust = 0,
          angle = -50,
          colour = "grey50"
        ),
        legend.text.align = 1,
        axis.title = element_text(size = 24)
      ) + xlab("") + ylab(ylab)
    g
  }