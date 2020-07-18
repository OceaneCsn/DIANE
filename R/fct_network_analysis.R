#' Returns modules-communities memberships
#'
#' @param graph igraph object, directed
#' @return (named) vector
community_structure <- function(graph) {
  g <- igraph::as.undirected(graph, mode = "collapse")
  return(igraph::membership(igraph::cluster_louvain(g)))
}

#' Plots the histogram of in and out degrees, betweeness of 
#' regultors and target genes.
#'
#' @param nodes dataframe containing the nodes information
#' @param graph igraph object
#'
#' @export
#' @examples 
#' data("abiotic_stresses")
#' data("regulators_per_organism")
#' 
#' mat <- abiotic_stresses$heat_DEGs_regulatory_links
#' network <- DIANE::network_thresholding(mat, n_edges = length(abiotic_stresses$heat_DEGs))
#' 
#' data <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]])
#' DIANE::draw_network_degrees(data$nodes, network)
draw_network_degrees <- function(nodes, graph) {
  targets <- nodes[nodes$gene_type == "Target Gene", "id"]
  TFs <- nodes[nodes$gene_type == "Regulator" | nodes$gene_type == "Grouped Regulators", "id"]
  
  
  degree_in_targets <-
    igraph::degree(graph, mode = 'in', v = targets)
  degree_in_tfs <- igraph::degree(graph, mode = "in", v = TFs)
  degree_out_tfs <- igraph::degree(graph, mode = "out", v = TFs)
  betweenness <- igraph::betweenness(graph, weights = NA, v = TFs)
  deg_targ = data.frame(degree_in_targets)
  Node_nw_st <-
    data.frame(degree_in_tfs, degree_out_tfs, betweenness)
  deg_in_targ <-
    ggplot2::ggplot(data = deg_targ, ggplot2::aes(x = degree_in_targets)) +
    ggplot2::geom_histogram(fill = "#69b3a2",
                            color = "#e9ecef",
                            alpha = 0.7) +
    ggplot2::ggtitle("In-Degree distribution of target genes") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
  
  deg_in_tfs <-
    ggplot2::ggplot(data = Node_nw_st, ggplot2::aes(x = degree_in_tfs)) +
    ggplot2::geom_histogram(fill = "#69b3a2",
                            color = "#e9ecef",
                            alpha = 0.7) +
    ggplot2::ggtitle("In-Degree distribution of regulators") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
  
  deg_out_tfs <-
    ggplot2::ggplot(data = Node_nw_st, ggplot2::aes(x = degree_out_tfs)) +
    ggplot2::geom_histogram(fill = "#69b322",
                            color = "#e9ecef",
                            alpha = 0.7) + ggplot2::xlim(0, 10) +
    ggplot2::ggtitle("Out-Degree distribution of regulators") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
  
  bet <-
    ggplot2::ggplot(data = Node_nw_st, ggplot2::aes(x = betweenness)) +
    ggplot2::geom_histogram(fill = "#E69F00",
                            color = "#e9ecef",
                            alpha = 0.7) +
    ggplot2::ggtitle("Betweeness distribution of the regulators") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
  return(gridExtra::grid.arrange(deg_in_targ, 
                                 deg_in_tfs, deg_out_tfs, bet, nrow = 2))
  
}

#' Targets and regulators of a gene in a network
#'
#' @param graph igraph object
#' @param node gene to describe
#'
#' @return list with targets and regulators attributes
#' @export
describe_node <- function(graph, node){
  
  targets <- igraph::neighborhood(graph, nodes = node, mode = "out")[[1]]$name
  regulators <- igraph::neighborhood(graph, nodes = node, mode = "in")[[1]]$name
  
  return(list(targets = setdiff(targets, node), 
              regulators = setdiff(regulators, node)))
}
