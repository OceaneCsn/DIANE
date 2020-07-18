#' For Arabdopsis thaliana, converts TAIR IDs to entrez IDs, 
#' gene symbol or description
#'
#' @param ids vector of AGI gene IDs
#' @param to value in c("entrez", "symbol", "name")
#'
#' @return named vector
#' @export
#'
#' @examples
#' genes <- c("AT2G05940", "AT4G16480", "AT4G04570", "AT2G30130", "AT1G56300")
#' entrez_ids <- convert_from_agi(genes)
#' print(entrez_ids)
#' symbols <- convert_from_agi(genes, to = "symbol")
#' print(symbols)
convert_from_agi <- function(ids, to = "entrez"){
  if (to == "entrez")
    x <- org.At.tair.db::org.At.tairENTREZID
  if (to == "symbol")
    x <- org.At.tair.db::org.At.tairSYMBOL
  if (to == "name")
    x <- org.At.tair.db::org.At.tairGENENAME
  mapped_genes <- AnnotationDbi::mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  return(unlist(xx[as.vector(ids)]))
}


#' For Homo sapiens, converts ensembl IDs to entrez IDs, 
#' symbol or name, 
#'
#' @param ids genes to convert, ensembl
#' @param to value in c("entrez", "symbol", "name")
#'
#' @return named list
#' @export
#' @examples
#' genes <- c("ENSG00000000003", "ENSG00000003989", "ENSG00000005884", "ENSG00000007168")
#' convert_from_ensembl(genes)
#' convert_from_ensembl(genes, to = "symbol")
convert_from_ensembl <- function(ids, to = "entrez"){
  if(to == "entrez"){
    xx <- as.list(org.Hs.eg.db::org.Hs.egENSEMBL2EG)
    return(unlist(xx[ids]))
  }
  else{
    entrez <- convert_from_ensembl(ids, to = "entrez")
    if (to == "symbol")
      x <- org.Hs.eg.db::org.Hs.egSYMBOL
    if (to == "name")
      x <- org.Hs.eg.db::org.Hs.egGENENAME
    mapped_genes <- AnnotationDbi::mappedkeys(x)
    xx <- as.list(x[mapped_genes])
    return(unlist(xx[as.vector(entrez)]))
  }
}


#' Enriched Gene Ontology terms in a set of genes
#' 
#' @description This function returns the enriched biological processes
#' in a set of genes, as compared to a background set of genes.
#' The GO terms are retreived by the clusterProfiler package and the
#' bioconductor organism databases (for example org.At.tair.db).
#' The frequencies comparisons between GO terms in the genes of interest and the 
#' background are performed using Fischer tests, with 0.05 adjusted pvalues
#' threshold. Those significantly enriched GO terms are then simplified 
#' with a semantic similarity threshold that can be chosen. It returns a
#' dataframe containing, for each significant GO, their gene count, gene symbols,
#' adjusted pvalue, description, etc. 
#'
#' @param genes gene set of interest. MUST be entrez IDs.
#' @param background gene set considered as the "universe" against which perfom 
#' the tests. MUST be entrez IDS.
#' @param org organism, in the format of bioconductor organisms databases (e.g 
#' org.xx.xx.db)
#' @param sim_cutoff similarity cutoff to use for pooling similar GO terms
#' @param GO_type character between BP (Biological Process), CC(Cellular Component) 
#' or MF (Molecular Function), depending on the GO subtype
#' to use in the analysis
#' @return data.frame
#' @export
#' @examples 
#' data("abiotic_stresses")
#' genes <- abiotic_stresses$heat_DEGs
#' 
#' genes <- get_locus(genes)
#' background <- get_locus(rownames(abiotic_stresses$normalized_counts))
#' 
#' genes <- convert_from_agi(genes)
#' background <- convert_from_agi(background)
#' 
#' go <- enrich_go(genes, background)
#' head(go)
enrich_go <- function(genes, background,
                      org = org.At.tair.db::org.At.tair.db,
                      sim_cutoff = 0.85, GO_type = "BP"){
  
  if(!GO_type %in% c("BP", "CC", "MF")){
    stop("GO_type must be either BP, CC or MF.")
  }
  
  ego <- clusterProfiler::enrichGO(gene = genes,
                                   OrgDb = org,
                                   ont = GO_type,
                                   universe = background,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable = TRUE)
  if(!is.null(ego)){
    simpOnt <- clusterProfiler::simplify(ego, cutoff=sim_cutoff, 
                                         by="p.adjust", select_fun=min)
    res <- simpOnt@result[order(-simpOnt@result$Count),]
    return(res)
  }
  else{return(NULL)}
  
}

#' Enriched custom Gene Ontology terms in a set of genes
#'
#' @param genes list of gene IDs, as present in the first column od genes_to_GO
#' @param universe list of gene IDs as background, default is set to all the genes 
#' contained in genesto_GO
#' @param genes_to_GO dataframe with gene IDs as first column, and GO IDs as second column
#' @param qvalue qvalue cutoff for enriched GO terms, default to 0.1
#' @param pvalue qvalue cutoff for enriched GO terms, default to 0.05
#'
#' @return dataframe containing enriched go terms and their description
#' @export
enrich_go_custom <- function(genes, universe = genes_to_GO[,1], genes_to_GO, qvalue = 0.1, pvalue = 0.05){
  
  if (length(genes) == 0) {
    stop("Empty list of genes")
  }
  
  if (length(universe) == 0) {
    stop("Empty universe")
  }
   
  if (length(intersect(genes, genes_to_GO[,1])) == 0) {
    stop("The genes are not in the first column of the custom dataframe")
  }
  
  go <- clusterProfiler::enricher(gene = genes, universe = universe, 
                                  TERM2GENE = genes_to_GO[,order(ncol(genes_to_GO):1)])
  go <- go@result
  go <- go[go$qvalue < qvalue & go$p.adjust < pvalue,]
  
  xx <- as.list(GO.db::GOTERM)
  goTerms <- sapply(go$ID, function(id){return(AnnotationDbi::Term(xx[[id]]))})
  
  go$Description <- goTerms[match(go$ID, names(goTerms))]
  return(go)
}

#' Plot the enriched go terms of an enrich_go result
#'
#' @param go_data dataframe returned by the enrich_go function,
#' containing for each significant GO, their gene count, gene symbols,
#' adjusted pvalue, description, etc. 
#' @param max_go maximum number of GO terms to display. Default is the 
#' total significant GO number, but for legibility reasons, you may want 
#' to reduce that number
#' @export
#'
#' @examples \dontrun{
#' data("abiotic_stresses")
#' genes <- abiotic_stresses$heat_DEGs
#' 
#' genes <- get_locus(genes)
#' background <- get_locus(rownames(abiotic_stresses$normalized_counts))
#' 
#' genes <- convert_from_agi(genes)
#' background <- convert_from_agi(background)
#' go <- enrich_go(genes, background)
#' draw_enrich_go(go)
#' draw_enrich_go(go, max_go = 20)
#' }

draw_enrich_go <- function(go_data, max_go = dim(go_data)[1]){
  
  go_data <- go_data[order(go_data$p.adjust),]
  res <- go_data[1:min(dim(go_data)[1],max_go),]
  
  res$Description <- substr(res$Description, 1,60)
  
  plotly::ggplotly(
    ggplot2::ggplot(data = res, ggplot2::aes(x = Count, y = Description, color = p.adjust)) +
      ggplot2::geom_point(ggplot2::aes(size = Count))+
      ggplot2::ylab("") + ggplot2::ggtitle("Enriched ontologies and their gene count") +
      ggplot2::scale_color_gradient(low="#114223", high="#92D9A2"), 
      tooltip = c("x", "y"))
}





#' Gives gene information (common name and description) for a specific organism
#'
#' @param ids vector of genes, AGI for Arabidopsis and ensembl for Human
#' @param organism value in c("Arabidopsis thaliana", "Homo sapiens")
#'
#' @return a dataframe with input genes as rownames, and columns label and desciption
#' @export
#' 
#' @examples 
#' genes <-  c("AT2G05940", "AT4G16480", "AT4G04570", "AT2G30130", "AT1G56300")
#' get_gene_information(genes, organism = "Arabidopsis thaliana")
get_gene_information <- function(ids, organism){
  if(organism == "Arabidopsis thaliana"){
   data("gene_annotations", package = "DIANE")
    d <- gene_annotations$`Arabidopsis thaliana`[
      match(ids, rownames(gene_annotations$`Arabidopsis thaliana`)),]
  }
  if (organism == "Homo sapiens"){
    # handling missing values in entrez ids
    entrez <- convert_from_ensembl(ids)
    entrez <- entrez[match(ids, names(entrez))]
    
    label <- convert_from_ensembl(ids, to = "symbol")
    label <- label[match(entrez, names(label))]
    
    description = convert_from_ensembl(ids, to = "name")
    description = description[match(entrez, names(description))]
    
    d <- data.frame(genes = ids, label = label,
                   description = description)
    
    rownames(d) <- d$genes
  }
  
  return(d[,c("label", "description")])
}


#' Number of shared genes between Two GO terms
#'
#' @param go_pair pair of go terms
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
  pairs <- data.frame(t(utils::combn(gos, 2)))
  
  pairs$common_genes <- sapply(paste(pairs[,1], pairs[,2]), interGO, go_table = go)
  pairs <- pairs[pairs$common_genes > 0,]
  graph <- igraph::graph_from_data_frame(pairs)
  
  layout <- ggraph::create_layout(graph, layout = 'kk')
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
