
#data(DIANE::raw_counts_At.RData)

load(system.file("extdata", "raw_counts_At.RData", package = "DIANE"))


raw_counts_At <- raw_data_At
raw_counts_At <- raw_counts_At[,grepl("c", colnames(raw_counts_At))]

library(stringr)
conditions_At <- str_split_fixed(colnames(raw_counts_At), '_', 2)[,1]

load(system.file("extdata", "design_At.RData", package = "DIANE"))
design_At <- design_At[grepl('c', rownames(design_At)), c("Nitrate", "Iron")]

demo_data_At <- list(raw_counts = raw_counts_At, conditions = conditions_At, design = design_At)


library(DIANE)
data("demo_data_At")

# visualize data before normalizing and filtering :
DIANE::draw_distributions(demo_data_At$raw_counts, boxplot = F)

DIANE::draw_heatmap(demo_data_At$raw_counts)
# normalization and low counts removal
tcc_object <- DIANE::normalize(demo_data_At$raw_counts, demo_data_At$conditions, iteration = FALSE)
threshold = 10*length(demo_data_At$conditions)
tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
normalized_counts <- TCC::getNormalizedData(tcc_object)


DIANE::draw_MDS(normalized.count = normalized_counts)
# visualize data after normalizing and filtering :
DIANE::draw_distributions(normalized_counts, boxplot = F)


# differential expression analysis
fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = demo_data_At$conditions)
topTags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 0.01)
topTags$table

# visualize differential expression analysis
tags <- DIANE::estimateDEGs(fit, reference = "cNF", perturbation = "cnF", p.value = 1)
DIANE::plotDEGs(tags, fdr = 0.01, lfc = 1)

# running clustering
genes <- topTags$table$genes
clustering <- DIANE::run_coseq(conds = unique(demo_data_At$conditions), data = normalized_counts, genes = genes, K = 6:9)
DIANE::draw_coseq_run(clustering$model, plot = "barplots")
DIANE::draw_coseq_run(clustering$model, plot = "ICL")


# visualize clustering profiles
DIANE::draw_profiles(data = normalized_counts, clustering$membership, conds = unique(demo_data_At$conditions)) 
DIANE::draw_profiles(data = normalized_counts, clustering$membership, conds = unique(demo_data_At$conditions), k = 3) 

genes_cluster <- DIANE::get_genes_in_cluster(clustering$membership, cluster = 3)
# clustering

DIANE::fit_glm(normalized_counts, genes_cluster, demo_data_At$design)



####### essais viseago

genes <- topTags$table$genes
background <- rownames(normalized_counts)
#Ensembl<-ViSEAGO::Ensembl2GO()

Bioconductor<-ViSEAGO::Bioconductor2GO()
#ViSEAGO::available_organisms(Bioconductor, )

myGENE2GO<-ViSEAGO::annotate(
  "org.At.tair.db",
  Bioconductor
)


convert_from_agi <- function(ids, to = "entrez"){
  if(to =="entrez")
    x <- org.At.tair.db::org.At.tairENTREZID
  if(to =="symbol")
    x <- org.At.tair.db::org.At.tairSYMBOL
  if(to =="name")
    x <- org.At.tair.db::org.At.tairGENENAME
  mapped_genes <- AnnotationDbi::mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  return(unlist(xx[as.vector(ids)]))
}
  
  
entrez <- convert_from_agi(genes)
entre_bg <- convert_from_agi(background)


s <- convert_from_agi(genes, to = "symbol")
head(s)

BP<-ViSEAGO::create_topGOdata(
  geneSel=entrez,
  allGenes=entre_bg,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

classic<-topGO::runTest(
  BP,
  algorithm ="classic",
  statistic = "fisher"
)

classic@score
classic@geneData

debug("GOcount")

body(ViSEAGO::GOcount)


standardGeneric

(BP_sResults)

if(length(xx) > 0) {# Get the Entrez gene IDs for the first five genes
  xx[1:5]# Get the first one
  xx[[1]]
}
