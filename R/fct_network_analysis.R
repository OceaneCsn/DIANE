

#' Returns modules-communities memberships
#'
#' @param graph igraph object, directed
#'
#' @return (named) vector
#' @export
#'
#' @examples
#' g <- igraph::sample_pa(1000, directed = TRUE)
#' head(community_structure(g))
community_structure <- function(graph){
  
  g <- igraph::as.undirected(graph, mode = "collapse")
  
  return(igraph::membership(igraph::cluster_louvain(g)))
}
