#' Get correlation
#'
#' @param pair character of two gene IDs sepearte by a white space
#' @param normalized.count normalized expression matrix containing 
#' the IDs of the pair
#'
#' @return spearman correlation bewteen the two expression vectors
get_correlation <- function(pair, normalized.count){
  
  tf1 <- str_split_fixed(pair, ' ',2)[1]
  tf2 <- str_split_fixed(pair, ' ',2)[2]
  
  if(tf1 %in% rownames(normalized.count) & tf2 %in% rownames(normalized.count))
  return(cor(as.numeric(normalized.count[tf1,]), as.numeric(normalized.count[tf2,]), method = "spearman"))
}

#' Group correlated regulators
#'
#' During network inference, the importance of regressors is assessed for aech target gene.
#' But in the case of highly correlated regulators, only one of them could be captured by the random 
#' forests, stealing importance from the others, and leave the other ones with no edge, which would be incorrect.
#' To account for this, this function summarizes the expression of the regulators being correlated above
#' a certain threshold in new variables with the mean of their profiles.
#' 
#' To do so, a graph of all regressors correlated (spearman) above the threshold id built, and groups
#' are formed with a community detection algorithm. Each group is then averaged in a single variable.
#' Strongly negatively correlated regressors are added to the group it is correlated to, but their expression
#' is not taken into account in the summary profile.
#' 
#' @param normalized.count normalized expression matrix containing genes and regressors.
#' @param genes target genes that want to be used in the inference process
#' @param regressors regressor genes that want to be used in the inference process.
#' @param corr_thr correlation threshold to be used for regressors grouping.
#'
#' @return a named list.
#' counts : the normalized expression data containing the new summarized variables, in 
#' the format mean_geneID1-geneID2... The individual genes that were grouped are removed.
#' correlated_regressors_graph : visNetwork interactive plot of the correlated regulators
#' grouped_genes : new vector of target genes, with individual correlated genes replaced by groups
#' grouped_regressors : new vector of regressors, with individual correlated genes replaced by groups
#' @export
#'
#' @examples
#' data(abiotic_stresses)
#' aggregated_data <- aggregate_splice_variants(abiotic_stresses$normalized_counts)
#' genes <- get_locus(abiotic_stresses$heat_DEGs)
#' regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])
#' grouping <- DIANE::group_regressors(aggregated_data, genes, regressors)
#' print(names(grouping))
group_regressors <- function(normalized.count, genes, regressors, corr_thr = 0.9){
  #calculating correlations for each TF pairs
  pairs <- data.frame(t(combn(regressors, 2)))
  pairs$cor <- sapply(paste(pairs[,1], pairs[,2]), get_correlation, 
                      normalized.count = normalized.count)
  top <- pairs[pairs$cor > corr_thr,]
  
  
  # graph and communities detection of highly correlated TFs
  net_un <- igraph::graph_from_data_frame(top, directed = FALSE)
  louvain <- igraph::cluster_louvain(net_un)
  groups <- igraph::membership(louvain)
  
  # graph plot creation
  d <- visNetwork::toVisNetworkData(net_un)
  
  d$edges$value <- d$edges$cor
  d$edges$label <- round(d$edges$cor,3)
  d$nodes$group <- groups[match(d$nodes$id, names(groups))]
  
  graph_plot <- d

  other_tfs <- regressors[!regressors %in% names(groups)]
  
  # Builiding the new consensus variables
  new_reg <- c()
  grouped_regs <- data.frame(matrix(nrow = length(unique(groups)), 
                                    ncol = length(colnames(normalized.count))))
  colnames(grouped_regs) <- colnames(normalized.count)
  rownames(grouped_regs) <- unique(groups)
  
  grouped_tfs <- c()
  
  for(group in unique(groups)){
    tfs <- names(groups)[groups==group]
    mean_tf <- colMeans(normalized.count[tfs,])
    grouped_regs[group,] <- mean_tf
    
    # find the negatively correlated regulators to that group, and add them, without
    # using their expression un the mean group expression
    
    for (tf in other_tfs){
      if(cor(as.numeric(normalized.count[tf,]), mean_tf, method = "spearman") < - corr_thr){
        print(paste("adding tf ", tf, " to group ", group, "because correlation of", 
                    cor(as.numeric(normalized.count[tf,]), mean_tf, method = "spearman"), "to mean"))
        tfs <- c(tfs, tf)
        # remove this tf from the list so it is not assigned to another group later
        other_tfs <- other_tfs[other_tfs != tf]
        print(length(other_tfs))
      }
    }
    new_reg <- c(new_reg, paste0("mean_",paste(tfs, collapse = "-")))
    grouped_tfs <- c(grouped_tfs, tfs)
  }
  
  rownames(grouped_regs) <- new_reg
  normalized.count <- rbind.data.frame(normalized.count, grouped_regs)
  
  #remove regressors that are grouped from the data
  normalized.count <- normalized.count[!rownames(normalized.count) %in% grouped_tfs,]
  return(list(counts = normalized.count, correlated_regressors_graph = graph_plot,
              grouped_genes = c(intersect(genes, rownames(normalized.count)), 
                                rownames(normalized.count)[grepl("means_", rownames(normalized.count))]),
              grouped_regressors = c(new_reg, regressors[!regressors %in% grouped_tfs])))
}
