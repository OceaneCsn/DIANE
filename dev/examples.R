
library(DIANE)

data("abiotic_stresses")
data("gene_annotations")

####### go analysis

genes <- convert_from_agi(get_locus(abiotic_stresses$heat_DEGs))
background <- convert_from_agi(get_locus(rownames(abiotic_stresses$normalized_counts)))



go <- enrich_go(genes, background)
#DIANE::draw_enrich_go(go, max_go = 30)

#' Number of shared genes between Two GO terms
#'
#' @param go_pair 
#' @param go_table go as returned by enrich_go
#'
#' @return integer
interGO <- function(go_pair, go_table){
  go1 <- unlist(strsplit(go_pair, ' '))[1]
  go2 <- unlist(strsplit(go_pair, ' '))[2]
  g1 <- unlist(strsplit(go_table[go_table$ID == go1, "geneID"], '/'))
  g1 <- g1[g1 != "NA"]
  g2 <- unlist(strsplit(go_table[go_table$ID == go2, "geneID"], '/'))
  g2 <- g2[g2 != "NA"]
  
  return(length(intersect(g1, g2)))
}

#' Draw GO enrichment map
#' 
#' This plots each enriched GO term, with its pvalue and gene count 
#' as color and size aesthetics.
#' Two GO terms are linked by an edges if they share common genes
#' in the input gene list given for enrichment analysis. The width 
#' of the edge is proportionnal to that number of shared genes.
#'
#' @param go dataframe as returned by \code{DIANE::enrich_go()} function.
#'
#' @export
#'
#' @examples
#' data("abiotic_stresses")
#' data("gene_annotations")
#' genes <- convert_from_agi(get_locus(abiotic_stresses$heat_DEGs))
#' background <- convert_from_agi(get_locus(rownames(abiotic_stresses$normalized_counts)))
#' go <- enrich_go(genes, background)
#' draw_enrich_go_map(go)
draw_enrich_go_map <- function(go){
  gos <- go$ID
  pairs <- data.frame(t(combn(gos, 2)))
  
  pairs$common_genes <- sapply(paste(pairs[,1], pairs[,2]), interGO, go_table = go)
  pairs <- pairs[pairs$common_genes > 0,]
  graph <- igraph::graph_from_data_frame(pairs)
  
  layout <- create_layout(graph, layout = 'kk')
  layout$GeneCount <- go[match(layout$name, go$ID),"Count"]
  layout$adj.p.value <- go[match(layout$name, go$ID),"p.adjust"]
  layout$label <- go[match(layout$name, go$ID),"Description"]
  
  ggraph::ggraph(layout) + 
    ggraph::geom_edge_link(aes(width = common_genes), alpha=0.1) + 
    ggraph::geom_node_point(aes(size = GeneCount, color=adj.p.value)) +
    ggraph::geom_node_label(aes(label = label), repel = TRUE, alpha = 0.75)+
    ggtitle("Gene Ontology enrichment map") + theme(plot.title = element_text(size = 20, face = "bold") ) +
    ggplot2::scale_color_gradient(low="#114223", high="#92D9A2")
}


library(ggraph)

########### tests ACP

draw_pca(normalized_counts)


## demo 

raw_counts <- read.csv("D:/These/DIANE_inputs/demo_counts.csv", sep = ',', h = T, row.names = "Gene")

design <- read.csv("D:/These/DIANE_inputs/abiotic_stresses_design.csv", sep = ',', h = T, row.names = "Condition")

library(DIANE)
data("abiotic_stresses")
data("regulators_per_organism")
data("gene_annotations")


conditions <- stringr::str_split_fixed(colnames(raw_counts), '_', 2)[,1]

gene_annotations[["Arabidopsis thaliana"]] <- gene_annotations[["Arabidopsis thaliana"]][1:dim(gene_annotations[["Arabidopsis thaliana"]])[1]-1,]

tail(gene_annotations[["Arabidopsis thaliana"]])

abiotic_stresses <- list(raw_counts = raw_counts, design = design, conditions = conditions)


############ tests group tfs


genes <- get_locus(abiotic_stresses$heat_DEGs)
regressors <- intersect(genes, regulators_per_organism$`Arabidopsis thaliana`)


data <- aggregate_splice_variants(abiotic_stresses$normalized_counts)
head(data)



r <- DIANE::group_regressors(data, genes, regressors)
tail(r$counts)


r$correlated_regressors_graph


mat <- DIANE::network_inference(r$counts, conds = abiotic_stresses$conditions, targets = r$grouped_genes,
                                 regressors = r$grouped_regressors)

network <- DIANE::network_thresholding(mat, n_edges = 500)

d_GENIE3 <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]])

DIANE::draw_network(d_GENIE3$nodes, d_GENIE3$edges)
DIANE::draw_network_degrees(d_GENIE3$nodes, network)


################# gene expression plot
library(ggplot2)






