
#' Genie3 regulatory importances estimation
#' 
#' @description GENIE3 needs to be fed a list of genes, that will be the nodes of the inferred network. 
#' Among those genes, some must be considered as potential regulators. 
#' GENIE3 can determine the influence if every regulators over each input genes, 
#' using their respective expression profiles. You can specify which conditions 
#' you want to be consired for those profiles during the network inference.
#' For each target gene, the methods uses Random Forests to provide a ranking of all 
#' regulators based on their influence on the target expression. This ranking is then merged 
#' across all targets, giving a global regulatory links ranking stored in the result matrix.
#' 
#' @param normalized.count normalized expression matrix containing the regressors and target genes
#' @param conds condition names to be used in the inference
#' @param regressors genes to be taken as regressors during the inference procedures (regulator genes)
#' @param targets genes to be included in the inferred network. Regressors can also be in the targets
#' @param nTrees Number of trees by Random Forest
#' @param nCores Number of CPU cores to use during the procedure
#'
#' @return matrix object
#' @export
#' @examples
#' \dontrun{
#' data("demo_data_At")
#' data("regulators_per_organism")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' genes <- topTags$table$genes
#' 
#' regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])
#' mat <- DIANE::network_inference(normalized_counts, conds = demo_data_At$conditions, 
#' targets = genes, regressors = regressors)
#' }

network_inference <- function(normalized.count, conds, regressors, targets, nTrees=5000, 
                              nCores=ifelse(is.na(parallel::detectCores()), 1, max(parallel::detectCores()-1, 1))){

  conditions <- colnames(normalized.count)[str_split_fixed(colnames(normalized.count), '_',2)[,1] %in% conds]
  
  normalized.count <- normalized.count[,conditions]
  
  regressors <- intersect(rownames(normalized.count), regressors)
  mat <- GENIE3::GENIE3(as.matrix(normalized.count), regulators = regressors, targets = targets, 
                treeMethod = "RF", K = "sqrt", nTrees = nTrees, nCores = 2, verbose = FALSE)
  
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
#' @param mat matrix containing the importance values for each target and regulator
#' @param n_edges number of edges top edges (based on their weight) to keep in the final network
#'
#' @return igraph object
#' @export
#' @examples
#' \dontrun{
#' data("demo_data_At")
#' data("regulators_per_organism")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' genes <- topTags$table$genes
#' 
#' regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])
#' mat <- DIANE::network_inference(normalized_counts, conds = demo_data_At$conditions, 
#' targets = genes, regressors = regressors)
#' network <- DIANE::network_thresholding(mat, n_edges = length(genes))
#' }
network_thresholding <- function(mat, n_edges){
  links <- GENIE3::getLinkList(mat, reportMax = n_edges)
  g <- igraph::graph_from_data_frame(links, directed = T)
  return(g)
}


#' Creates data compatible with visNetwork from an igraph object
#' 
#' @description Creates dataframe that describe the network.
#' The degree and community of each node are computed.
#' Communities are defined by the louvain algorithm, that 
#' maximizes local modularity. The gene type is asssigned to
#' eacg gene, either a regulator or a target gene.
#'
#' @param graph igraph object
#' @param regulators list of regulators
#'
#' @return list of dataframes containing nodes and edges information
#' @export
#' @examples
#' \dontrun{
#' data("demo_data_At")
#' data("regulators_per_organism")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' genes <- topTags$table$genes
#' 
#' regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])
#' mat <- DIANE::network_inference(normalized_counts, conds = demo_data_At$conditions, 
#' targets = genes, regressors = regressors)
#' network <- DIANE::network_thresholding(mat, n_edges = length(genes))
#' data <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]])
#' }
network_data <- function(graph, regulators){
  data <- visNetwork::toVisNetworkData(graph)
  
  degree <- igraph::degree(graph)
  
  # degree computaton
  data$nodes$degree <- degree[match(data$nodes$id, names(degree))]
  
  # modules computation
  memberships <- community_structure(graph)
  
  data$nodes$community <- memberships[match(data$nodes$id, names(memberships))]
  
  data$nodes$group <- ifelse(data$nodes$id %in% regulators, "Regulator", "Target Gene")
  data$nodes$gene_type <- ifelse(data$nodes$id %in% regulators, "Regulator", "Target Gene")
  data$edges$value <- data$edges$weight
  return(data)
}


#' Displays an interactive network view, with regulators as green sqaure nodes.
#'
#' @param nodes dataframe containing nodes
#' @param edges dataframe containing edges
#' 
#' @import visNetwork
#' @export
#' @examples 
#' data("demo_data_At")
#' data("regulators_per_organism")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' genes <- topTags$table$genes
#' 
#' regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])
#' mat <- DIANE::network_inference(normalized_counts, conds = demo_data_At$conditions, 
#' targets = genes, regressors = regressors)
#' network <- DIANE::network_thresholding(mat, n_edges = length(genes))
#' data <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]])
#' # adding common names as label for network visualisation
#' data$nodes$label <- demo_data_At$gene_info[match(data$nodes$id, rownames(demo_data_At$gene_info)), "label"]
#' DIANE::draw_network(data$nodes, data$edges)
draw_network <- function(nodes, edges){
  visNetwork(nodes = nodes, edges = edges) %>%
  visEdges(smooth = FALSE, arrows = 'to', color = '#333366') %>%
    visPhysics(
      solver = "forceAtlas2Based",
      timestep = 0.6,
      minVelocity = 12,
      maxVelocity = 10,
      stabilization = F
    ) %>%
    visGroups(
      groupname = "Regulator",
      size = 28,
      color = list("background" = "#49A346", "border" = "#FFFFCC"),
      shape = "square"
    ) %>%
    visGroups(groupname = "Target Gene",
              color = list("background" = "#B6B3B3", hover = "grey",
                           "border" = "#96E69A")) %>%
    visNodes(borderWidth = 0.5, font = list("size" = 35))
}