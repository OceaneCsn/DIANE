
#' Regulatory importances estimation via Random Forests
#' 
#' @description GENIE3 needs to be given a list of genes, that will be the nodes of the inferred network. 
#' Among those genes, some must be considered as potential regulators. 
#' GENIE3 can determine the influence if every regulators over each input genes, 
#' using their respective expression profiles. You can specify which conditions 
#' you want to be considered for those profiles during the network inference.
#' For each target gene, the methods uses Random Forests to provide a ranking of all 
#' regulators based on their influence on the target expression. This ranking is then merged 
#' across all targets, giving a global regulatory links ranking stored in the result matrix.
#' 
#' @param normalized.count normalized expression matrix containing the regressors and target genes
#' in its rows, and samples a columns
#' @param conds condition names to be used in the inference (not columns names, conditions names
#' before the underscore)
#' @param regressors genes to be taken as regressors during the inference procedures (regulator genes)
#' @param targets genes to be included in the inferred network. Regressors can also be in the targets
#' @param nTrees Number of trees by Random Forest
#' @param nCores Number of CPU cores to use during the procedure. Default is the detected number of cores
#' minus one.
#' @param verbose If set to TRUE, a feedback on the progress of the calculations is given. Default: TRUE
#' @param importance_metric character being either node_purity or MSEincrease_oob. This is the importance type
#' computed for the regulator-gene pairs, as returned by the randomForest package.  Default is node_purity,
#' the metric used in GENIE3. Our improvement of the method uses MSEincrease_oob for consistency reasons regarding
#' to statistical edges testing. The default one is around 4 times fatser, but more sensitive to the number of
#'  regulators and to over-fitting. Too few samples will lead to NA in MSEincrease_oob, so in that case, it is
#'  advised to used GENIE3's default one.
#'
#' @return Matrix filled with regulator-target regulatory weights
#' @export
#' @examples
#' \dontrun{
#' data("abiotic_stresses")
#' data("regulators_per_organism")
#' 
#' aggregated_data <- aggregate_splice_variants(data = abiotic_stresses$normalized_counts)
#' 
#' genes <- get_locus(abiotic_stresses$heat_DEGs)
#' regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])
#' 
#' mat <- network_inference(aggregated_data, conds = abiotic_stresses$conditions, 
#' targets = genes, regressors = regressors, nTrees = 1000, nCores = 4)
#' }

network_inference <- function(normalized.count, conds, regressors, targets, nTrees=1000, 
                              nCores = ifelse(is.na(parallel::detectCores()), 
                                              1, max(parallel::detectCores() - 1, 1)),
                              verbose = TRUE, importance_metric = "node_purity"){

  conditions <- colnames(normalized.count)[stringr::str_split_fixed(
    colnames(normalized.count), '_',2)[,1] %in% conds]
  
  
  if (length(conditions) == 0) {
    stop("The required conditions were not found in the expression data")
  }
  
  if (sum(targets %in% rownames(normalized.count)) == 0) {
    stop("The required target genes were not found in the expression data")
  }
  
  if (!importance_metric %in% c("node_purity", "MSEincrease_oob")) {
    stop("The importance_metric argument must be either node_purity or 
         MSEincrease_oob")
  }
  
  normalized.count <- normalized.count[,conditions]
  
  regressors <- intersect(rownames(normalized.count), regressors)
  
  if (length(regressors) == 0) {
    stop("No regressors were found among the targets")
  }
  
  if (importance_metric == "node_purity")
    mat <- GENIE3::GENIE3(as.matrix(normalized.count), regulators = regressors, targets = targets, 
                treeMethod = "RF", K = "sqrt", nTrees = nTrees, nCores = nCores, verbose = verbose)
  if (importance_metric == "MSEincrease_oob")
    mat <- GENIE3OOB(as.matrix(normalized.count), regulators = regressors, targets = targets,
                     nTrees = nTrees, nCores = nCores, verbose = verbose)
  return(mat)
}


#' Thresholds a non sparse adjascency matrix to return the strongest weights only
#' 
#' @description 
#' Without thresholding, we would obtain a fully connected weighted 
#' graph from GENIE3, with far too many links to be interpretable.
#' In order build a meaningful network, this weighted adjacency matrix between 
#' regulators and targets has to be sparsified, and we have to determine the regulatory 
#' weights that we consider significant.
#' 
#' This method is a nice exploratory way to threshold complete networks,
#' but to get more robust and significant results, consider using 
#' the \code{DIANE::test_edges()} function.
#'
#' @param mat matrix containing the importance values for each target and regulator,
#' as returned by \code{DIANE::network_inference()}
#' @param n_edges number of edges top edges to keep in the final network.
#'
#' @return igraph object representing the Gene Regulatory Network
#' @export
#' @examples
#' \dontrun{
#' data("abiotic_stresses")
#' data("regulators_per_organism")
#' genes <- get_locus(abiotic_stresses$heat_DEGs)
#' # mat was inferred using the function network_inference
#' mat <- abiotic_stresses$heat_DEGs_regulatory_links
#' network <- DIANE::network_thresholding(mat, n_edges = length(genes))
#' }
network_thresholding <- function(mat, n_edges){
  links <- GENIE3::getLinkList(mat, reportMax = n_edges)
  g <- igraph::graph_from_data_frame(links, directed = T)
  return(g)
}


#' Creates data describing a gene regulatory network
#' 
#' 
#' @description Creates dataframe that describe the network
#' compatible with visNetwork, from an igraph object.
#' The degree and community of each node are computed.
#' Communities are defined by the louvain algorithm, that 
#' maximizes local modularity. A type is assigned to
#' each gene, either a regulator or a target gene.
#'
#' @param graph igraph object, as returned by \code{DIANE::network_thresholding()},
#' or \code{DIANE::network_from_tests()}.
#' @param regulators list of regulators, so they can be marked
#' as special nodes in the network
#' @param gene_info dataframe with gene IDs as rownames, containing
#' additional info on genes, in columns named label and description.
#' For known organisms, it can be obtained from \code{DIANE::get_gene_information()}
#'
#' @return Named list of dataframes containing nodes and edges information.
#' @export
#' @examples
#' \dontrun{
#' data("abiotic_stresses")
#' data("regulators_per_organism")
#' genes <- get_locus(abiotic_stresses$heat_DEGs)
#' # mat was inferred using the function network_inference
#' mat <- abiotic_stresses$heat_DEGs_regulatory_links
#' network <- DIANE::network_thresholding(mat, n_edges = length(genes))
#' data <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]])
#' }
network_data <- function(graph, regulators, gene_info = NULL){
  data <- visNetwork::toVisNetworkData(graph)
  
  degree <- igraph::degree(graph)
  
  # degree computaton
  data$nodes$degree <- degree[match(data$nodes$id, names(degree))]
  
  # modules computation
  memberships <- community_structure(graph)
  
  data$nodes$community <- memberships[match(data$nodes$id, names(memberships))]
  
  data$nodes$group <- ifelse(data$nodes$id %in% regulators, "Regulator", 
                             ifelse(grepl("mean_", data$nodes$id), 
                                    "Grouped Regulators", "Target Gene"))
  data$nodes$gene_type <- data$nodes$group
  if("weight" %in% colnames(data$edges))
    data$edges$value <- data$edges$weight
  if("fdr" %in% colnames(data$edges))
    data$edges$value <- 2
  #data$edges$value <- log10(-log10(data$edges$fdr))
  # adding additional infos
  if(!is.null(gene_info)){
    data$nodes[,colnames(gene_info)] <- 
      gene_info[match(data$nodes$id, rownames(gene_info)), ]
    
    # for grouped nodes, concatenates their labels
    if("label" %in% colnames(gene_info)){
      for(id in data$nodes$id){
        if(grepl("mean_", id)){
          ids <- stringr::str_split_fixed(id, '_', 2)[,2]
          ids <- unlist(strsplit(ids, '-'))
          labels <- paste0("mean_", paste(gene_info[match(ids, rownames(gene_info)), "label"], collapse = '-'))
          data$nodes$label[data$nodes$id == id] <- labels
        }
      }
    }
  }
  return(data)
}


#' Interactive network view
#' 
#' Visualize the GRN with regulators as green square nodes,
#' and targets as gray discs.
#'
#' @param nodes dataframe containing nodes, from \code{DIANE::network_data()$nodes}
#' @param edges dataframe containing edges, from \code{DIANE::network_data()$edges}
#' 
#' @import visNetwork
#' @export
#' @examples \dontrun{
#' data("abiotic_stresses")
#' data("regulators_per_organism")
#' genes <- get_locus(abiotic_stresses$heat_DEGs)
#' # mat was inferred using the function network_inference
#' mat <- abiotic_stresses$heat_DEGs_regulatory_links
#' network <- DIANE::network_thresholding(mat, n_edges = length(genes))
#' data <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]])
#' DIANE::draw_network(data$nodes, data$edges)}
draw_network <- function(nodes, edges){
  library("visNetwork")
  visNetwork::visNetwork(nodes = nodes, edges = edges) %>%
    visNetwork::visEdges(smooth = FALSE, arrows = 'to', color = '#333366') %>%
    visPhysics(
      solver = "forceAtlas2Based",
      timestep = 0.6,
      minVelocity = 12,
      maxVelocity = 10,
      stabilization = F
    ) %>%
    visNetwork::visGroups(
      groupname = "Regulator",
      size = 28,
      color = list("background" = "#49A346", "border" = "#FFFFCC"),
      shape = "square"
    ) %>%
    visNetwork::visGroups(
      groupname = "Grouped Regulators",
      size = 45,
      color = list("background" = "#1C5435", "border" = "#FFFFCC"),
      shape = "square"
    ) %>%
    visNetwork::visGroups(groupname = "Target Gene",
              color = list("background" = "#B6B3B3", hover = "grey",
                           "border" = "#96E69A")) %>%
    visNetwork::visNodes(borderWidth = 0.5, font = list("size" = 35))
}
