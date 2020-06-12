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


#' Computes enrich Gene Ontology terms in a set of genes
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
#'
#' @return data.frame
#' @export
#' @examples 
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
#' # interest and background gene sets as entrez ids
#' genes <- convert_from_agi(topTags$table$genes)
#' background <- convert_from_agi(rownames(normalized_counts))
#' go <- enrich_go(genes, background)
#' head(go)
enrich_go <- function(genes, background,
                      org = org.At.tair.db::org.At.tair.db,
                      sim_cutoff = 0.85){
  
  ego <- clusterProfiler::enrichGO(gene = genes,
                                   OrgDb = org,
                                   ont = "BP",
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
#' @examples
#' data("demo_data_At")
#' tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, 
#' iteration = FALSE)
#' threshold = 10*length(demo_data_At$conditions)
#' tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
#' normalized_counts <- TCC::getNormalizedData(tcc_object)
#' fit <- DIANE::estimateDispersion(tcc = tcc_object,
#'  conditions = demo_data_At$conditions)
#' topTags <- DIANE::estimateDEGs(fit, reference = "cNF", 
#' perturbation = "cnF", p.value = 0.01)
#' # interest and background gene sets as entrez ids
#' genes <- convert_from_agi(topTags$table$genes)
#' background <- convert_from_agi(rownames(normalized_counts))
#' 
#' go <- enrich_go(genes, background)
#' draw_enrich_go(go)
#' draw_enrich_go(go, max_go = 20)
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
   data("demo_data_At", package = "DIANE")
    d <- demo_data_At$gene_info[match(ids, rownames(demo_data_At$gene_info)),]
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
