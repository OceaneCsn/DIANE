
library(DIANE)

data("abiotic_stresses")
data("gene_annotations")

####### go analysis

genes <- convert_from_agi(get_locus(abiotic_stresses$heat_DEGs))
background <- convert_from_agi(get_locus(rownames(abiotic_stresses$normalized_counts)))



go <- enrich_go(genes, background)
#DIANE::draw_enrich_go(go, max_go = 30)




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






