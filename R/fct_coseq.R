



#' run_coseq
#'
#' @param conds Condition names to be used for clustering
#' @param genes Genes used as an input for the clustering
#' @param data normalized counts
#' @param K number of clusters range
#' 
#' @importFrom coseq coseq clusters
#'
#' @return named list containing the coseq run result, and the cluster membership for each gene
#' @export
#'
run_coseq <- function(conds, genes, data, K = 6:12){
  
  conditions <- unique(grep(paste(conds, collapse = "|"), 
              colnames(data), value = TRUE))
  
  groups <- str_split_fixed(conditions, '_', 2)[,1]
  dataC <- round(data[genes,conditions],0)
  run_pois <- coseq::coseq(dataC, conds = groups, K = K, model = "Poisson", iter = 5, transformation = "none", 
                           normFactors = "none", parallel = TRUE)
  return(list(membership = coseq::clusters(run_pois),
              model = run_pois))
}



#' draw_coseq_run
#'
#' @param run_pois result of a coseq run
#'
#' @return plot describing the quality of the clustering process
#' 
#' @importFrom coseq plot
#' @export
#'
draw_coseq_run <- function(run_pois, plot = "ICL"){
  if (plot == "ICL") coseq::plot(run_pois, graphs = c("ICL"))
  else coseq::plot(run_pois, graphs = c("probapost_barplots"))

}


#' Draw profiles of a clustering
#'
#' @param data normalized counts
#' @param membership membership item of the coseq object returned by run_coseq object
#' @param expression if it is set to "profiles" (default), plots expression/sum(expression). if "counts", plots log(Counts+1)
#' @param k if NULL (default), plot all the clusters. Else, plot the clusters in the vetcor k.
#' @param nrow on how many rows display the cluster profiles if k is NULL
#'
#' @importFrom reshape2 melt
#' @importFrom stringr str_split_fixed
#' @import ggplot2
#' @return ggplot
#' @export
#'
draw_profiles <- function(data, membership, expression = "profiles", k = NULL,
                          nrow = 3){
  
  clusters <- membership
  
  
  data <- data[names(clusters),]
  
  
  if(expression == "profiles"){profiles <- data.frame(data/rowSums(data))
  ylab <- "Normalized counts/Mean(Normalized counts)"}
  
  if(expression == "counts") {
    profiles <- data.frame(log(as.matrix(data)+1))
    ylab <- "log(Normalized counts)"
  }
  
  profiles$gene <- rownames(profiles)
  d <- reshape2::melt(profiles)
  d$group <- stringr::str_split_fixed(d$variable, '_', 2)[,1]
  d$cluster <- clusters[match(d$gene, names(clusters))]
  d$geneRep <- paste0(d$gene,stringr::str_split_fixed(d$variable, '_', 2)[,2])
  
  print(head(d))
  
  if (is.null(k)) {
    g <- ggplot2::ggplot(data = d, aes(x = group, y = value))  + ggplot2::facet_wrap(~cluster, scales = "free") 
  }
  else{
    k <- as.vector(k)
    g <- ggplot2::ggplot(data = d[d$cluster %in% k,], aes(x=group, y=value)) + ggplot2::facet_wrap(~cluster, scales = "free") 
  }
  
  g <- g + ggplot2::geom_line(alpha = 0.08,lwd = 0.9, aes(group = geneRep, color = as.factor(cluster)))
  
  g <- g + ggplot2::geom_boxplot(alpha=0.4, lwd=1, color = "black", outlier.color = "black", outlier.alpha = 0.1) + 
    ggplot2::geom_jitter(width = 0.1, alpha = 0.0015) + ggtitle("Expression profiles of the clusters")
  
  g <- g + ggplot2::theme(plot.title = element_text(size=22, face="bold"), strip.text.x = element_text(size = 20), legend.position="none",
                axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 12, hjust = 0, angle = -50, colour = "grey50"),legend.text.align=1,
                axis.title=element_text(size=24)) + xlab("") + ylab(ylab)
  g
}