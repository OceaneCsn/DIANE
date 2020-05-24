
#' Genie3 regulatory importances estimation
#'
#' @param normalized.count normalized expression matrix containing the regressors and target genes
#' @param conds condition names to be used in the inference
#' @param regressors genes to be taken as regressors during the inference procedures (regulator genes)
#' @param targets genes to be included in the inferred network. Regressors can also be in the targets
#' @param nTrees Number of trees by Random Forest
#' @param nCores Number of CPU cores to use during the procedure
#'
#' @return matrix object

network_inference <- function(normalized.count, conds, regressors, targets, nTrees=5000, nCores=1){
  
  conditions <- unique(grep(paste(conds, collapse = "|"),
                            colnames(normalized.count), value = TRUE))
  
  normalized.count <- normalized.count[,conditions]
  
  regressors <- intersect(rownames(normalized.count), regressors)
  mat <- GENIE3::GENIE3(as.matrix(normalized.count), regulators = regressors, targets = targets, 
                treeMethod = "RF", K = "sqrt", nTrees = nTrees, nCores = nCores, verbose = T)
  
  return(mat)
}


#' Thresholds a non sparse adjascency matrix to return the strongest weights only
#'
#' @param mat matrix containing the importance values for each target and regulator
#' @param n_edges number of edges to keep in the final network
#'
#' @return igraph object

network_thresholding <- function(mat, n_edges){
  links <- GENIE3::getLinkList(mat, reportMax = n_edges)
  g <- igraph::graph_from_data_frame(links, directed = T)
  return(g)
}


#' Get data compatible with visNetwork from igraph object
#'
#' @param graph igraph object
#' @param regulators list of regulators
#'
#' @return list of dataframes containing nodes and edges information

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


#' Displays an interavtive network view
#'
#' @param nodes dataframe containing nodes
#' @param edges dataframe containing edges
#' 
#' @import visNetwork
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
    visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}") %>%
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