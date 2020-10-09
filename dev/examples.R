
library(DIANE)

data("abiotic_stresses")
data("gene_annotations")
data("regulators_per_organism")

####### go analysis

genes <- convert_from_agi(get_locus(abiotic_stresses$heat_DEGs))
background <- convert_from_agi(get_locus(rownames(abiotic_stresses$normalized_counts)))



go <- enrich_go(genes, background)
#DIANE::draw_enrich_go(go, max_go = 30)


library(plotly)
draw_heatmap(data = abiotic_stresses$normalized_counts)


heatmap(abiotic_stresses$heat_DEGs_regulatory_links, )


library(ggraph)

########### tests ACP

draw_PCA(abiotic_stresses$normalized_counts)


## demo 

raw_counts <- read.csv("D:/These/DIANE_inputs/demo_counts.csv", sep = ',', h = T, row.names = "Gene")

design <- read.csv("D:/These/DIANE_inputs/abiotic_stresses_design.csv", sep = ',', h = T, row.names = "Condition")

library(DIANE)
data("abiotic_stresses")
data("regulators_per_organism")
data("gene_annotations")


genes <- sample(abiotic_stresses$heat_DEGs, 4)
ggplotly(draw_expression_levels(abiotic_stresses$normalized_counts, genes = genes))%>%
  layout(legend = list(font = list(size = 15)))

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

network <- DIANE::network_thresholding(mat, n_edges = 150)

d_GENIE3 <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]])

library(visNetwork)
DIANE::draw_network(d_GENIE3$nodes, d_GENIE3$edges) %>% visNodes(font = list("size" = 0) )
DIANE::draw_network_degrees(d_GENIE3$nodes, network)


################# custom GO

library(clusterProfiler)

data <- read.csv("D:/These/DIANE_inputs/lupin_gene_count_matrix.csv", header = TRUE, row.names = "Gene")

tcc <- DIANE::normalize(data, conditions = str_split_fixed(colnames(data), '_', 2)[,1])
tcc <- DIANE::filter_low_counts(tcc, 10*ncol(data))

normalized_counts <- TCC::getNormalizedData(tcc)
DIANE::draw_distributions(normalized_counts)


fit <- DIANE::estimateDispersion(tcc )

degs <- DIANE::estimateDEGs(fit, reference = "S0", perturbation = "S6", p.value = 0.01, lfc = 3)
genes <- degs$table$genes

GOs <- read.table("D:/These/DIANE_inputs/lupin_golist.txt", header = TRUE, sep = ',')

universe <- intersect(rownames(normalized_counts), GOs[,1])

go <- enrich_go_custom(genes, universe, GOs)


DIANE::draw_enrich_go_map(go)


########## R to cytoscape

load("~/Documents/network_CO2genes12Conds.RData")
ig <- igraph::graph_from_data_frame(net_data$edges, directed=TRUE, vertices = net_data$nodes)
igraph::plot.igraph(ig, vertex.size = 2, vertex.cex = 2)
RCy3::createNetworkFromIgraph(ig,"myIgraph")

########## sig genie3 tests and abiotic stresses values for mat and tests and grouped
library(DIANE)

data("abiotic_stresses")
data("gene_annotations")
data("regulators_per_organism")


genes <- get_locus(abiotic_stresses$heat_DEGs)
regressors <- intersect(genes, 
                        regulators_per_organism$`Arabidopsis thaliana`)

data <- aggregate_splice_variants(abiotic_stresses$normalized_counts)

r <- DIANE::group_regressors(data, genes, regressors)



mat <- DIANE::network_inference(r$counts, 
                                conds = abiotic_stresses$conditions, 
                                targets = r$grouped_genes,
                                regressors = r$grouped_regressors, 
                                importance_metric = "MSEincrease_oob", 
                                verbose = TRUE)

library(tictoc)

tic("test edges")
res <- DIANE::test_edges(mat, normalized_counts = r$counts, density = 0.03,
                         nGenes = length(r$grouped_genes), 
                         nRegulators = length(r$grouped_regressors), 
                         nTrees = 1000, verbose = TRUE)
toc()

# abiotic_stresses$heat_edge_tests = res
# abiotic_stresses$heat_DEGs_regulatory_links <- mat

########## verif de la correlation des fdrs d'un run à l'autre des tests

links0 <- abiotic_stresses$heat_edge_tests$links
links <- res$links

links$pair <- paste(links$regulatoryGene, links$targetGene)
links0$pair <- paste(links0$regulatoryGene, links0$targetGene)

common <- intersect(links$pair, links0$pair)
commonlinks <- links[links$pair %in% common,]

plot(commonlinks$fdr, links0[match(commonlinks$pair, links0$pair),]$fdr)
cor(commonlinks$fdr, links0[match(commonlinks$pair, links0$pair),]$fdr)


res$fdr_nEdges_curve

res$pvalues_distributions + xlim(0,0.1)

net <- DIANE::network_from_tests(res$links, fdr = 0.01)


net_data <- network_data(net, regulators_per_organism[["Arabidopsis thaliana"]], gene_annotations$`Arabidopsis thaliana`)

DIANE::draw_network(net_data$nodes, net_data$edges)

abiotic_stresses$heat_grouped_data <- r



# mat2 <- DIANE::network_inference(r$counts, conds = abiotic_stresses$conditions, targets = r$grouped_genes,
#                                 regressors = r$grouped_regressors, 
#                                 verbose = TRUE)

# res <- data.frame(density = seq(0.001, 0.1, length.out = 20),
#                   nEdges = sapply(seq(0.001, 0.1, length.out = 20), 
#                                   get_nEdges, nGenes, nRegulators))
# ggplot(res, aes(x = density, y = nEdges)) + geom_line(size = 1) + ggtitle()


######## Mm TFs 

tfs <- read.csv('D:/These/DIANE_inputs/Mus_musculus_TF.txt', sep = '\t')
tfs <- tfs$Ensembl
write.table(tfs, row.names = F, quote = F, file = "D:/These/DIANE_inputs/Mus_musculus_TF.txt")
regulators_per_organism[["Mus musculus"]] <- tfs





################ lupine data



annot <- read.csv("D:/These/Thesis/DIANE_inputs/genesLupinAnnot.csv", sep = '\t', row.names = 1, header = FALSE)
colnames(annot) <- c("description")
annot$label <- rownames(annot)

library(DIANE)
data("regulators_per_organism")

regulators_per_organism[["Lupinus albus"]] <- annot[stringr::str_detect(annot$description, "transcription"), 'label']
gos <- read.csv("D:/These/Thesis/DIANE_inputs/lupin_golist.txt", sep = '\t')
lupine <- list(annotation = annot, go_list = gos)



######### labels TAIT/ensembl


tair <- read.csv("~/Documents/Seafile/Thesis/DIANE_inputs/col_brm_lbd_David_TAIR10.csv", sep = '\t', row.names = "Gene", h = TRUE)
ens <- read.csv("~/Documents/Seafile/Thesis/DIANE_inputs/col_brm_lbd_David_Ensembl.csv", sep = '\t', row.names = "Gene", h = TRUE)


DIANE::check_IDs(rownames(tair), organism = "Arabidopsis thaliana")
DIANE::check_IDs(rownames(ens), organism = "Arabidopsis thaliana")
length(intersect(rownames(tair), rownames(ens)))

library(DIANE)
data("abiotic_stresses")
data("gene_annotations")
data("regulators_per_organism")


conditions <- stringr::str_split_fixed(colnames(ens), '_', 2)[,1]

annot <- gene_annotations$`Arabidopsis thaliana`


which(!stringr::str_detect( rownames(ens), pattern = "\\.[[:digit:]]+$"))

loci.ens <- get_locus(rownames(ens))

DIANE::

gene_annotations$`Arabidopsis thaliana`[sample(rownames(ens), 20),]
gene_annotations$`Arabidopsis thaliana`[sample(rownames(tair), 20),]


tcc_object <- DIANE::normalize(abiotic_stresses$raw_counts, abiotic_stresses$conditions, iteration = FALSE)
threshold = 10*length(conditions)
tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
normalized_counts <- TCC::getNormalizedData(tcc_object)

############### new coseq tests

library(DIANE)
data("abiotic_stresses")
data("gene_annotations")
data("regulators_per_organism")


genes <- abiotic_stresses$heat_DEGs


clustering <- DIANE::run_coseq(
  conds = unique(abiotic_stresses$conditions), 
  data = abiotic_stresses$normalized_counts, genes = genes, K = 6:9)

normSin <- run_coseq(
  conds = unique(abiotic_stresses$conditions), transfo = "arcsin", model = "Normal",
  data = abiotic_stresses$normalized_counts, genes = genes, K = 2:9)

normLogit <- run_coseq(
  conds = unique(abiotic_stresses$conditions), transfo = "logit", model = "Poisson",
  data = abiotic_stresses$normalized_counts, genes = genes, K = 2:9)


DIANE::draw_coseq_run(normSin$model, plot = "barplots")
DIANE::draw_coseq_run(normSin$model, plot = "ICL")
DIANE::draw_coseq_run(clustering$model, plot = "ICL")


draw_profiles(abiotic_stresses$normalized_counts, normLogit$membership, abiotic_stresses$conditions)
coseq::compareICL(list(normLogit$model, normSin$model))

summary(normLogit)




########## sig genie time estimation
library(DIANE)

data("abiotic_stresses")
data("gene_annotations")
data("regulators_per_organism")


genes <- get_locus(abiotic_stresses$heat_DEGs)
regressors <- intersect(genes, 
                        regulators_per_organism$`Arabidopsis thaliana`)

data <- aggregate_splice_variants(abiotic_stresses$normalized_counts)

r <- DIANE::group_regressors(data, genes, regressors)



mat <- DIANE::network_inference(r$counts, 
                                conds = abiotic_stresses$conditions, 
                                targets = r$grouped_genes,
                                regressors = r$grouped_regressors, 
                                importance_metric = "MSEincrease_oob", 
                                verbose = TRUE)

library(tictoc)

tic("test edges")
res <- DIANE::test_edges(mat, normalized_counts = r$counts, density = 0.02,
                         nGenes = length(r$grouped_genes), 
                         nRegulators = length(r$grouped_regressors), 
                         nTrees = 1000, verbose = TRUE)

toc()

tic("test")
a <- toc()
elapsed <- a$toc - a$tic




t <- estimate_test_edges_time(mat, normalized_counts = r$counts, density = 0.02,
                         nGenes = length(r$grouped_genes), 
                         nRegulators = length(r$grouped_regressors), 
                         nTrees = 1000, verbose = TRUE)
t/60


g2 <- sample(rownames(abiotic_stresses$normalized_counts), size = 4000)
g1 <- sample(rownames(abiotic_stresses$normalized_counts), size = 4040)

venn <- ggVennDiagram::ggVennDiagram(list(g1 = g1, g2 = g2), color = "#888888")
gene_lists <- list(g1 = g1, g2 = g2)

library(rlist)

length(setdiff(g1, g2))


list.count(gene_list, )


get_specific <- function(gene_lists, l ){
  common <- Reduce(intersect, gene_lists)
  return(setdiff(gene_lists[[l]], common))
}



length(get_specific(gene_lists, "g1"))



### tfs droso

tfs <- read.table("../Drosophilia_TFs.txt", h = F)
regulators_per_organism[["Drosophilia melanogaster"]] <- tfs$V1
