



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
  if(plot == "ICL") coseq::plot(run_pois, graphs = c("ICL"))
  else coseq::plot(run_pois, graphs = c("probapost_barplots"))
}