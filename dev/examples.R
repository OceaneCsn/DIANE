
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
data("regulators_per_organism")

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
DIANE::draw_DEGs(tags, fdr = 0.005, lfc = 1)



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




############### network inference



genes <- topTags$table$genes

regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])

mat <- DIANE::network_inference(normalized_counts, conds = demo_data_At$conditions, targets = genes,
                         regressors = regressors)
network <- DIANE::network_thresholding(mat, n_edges = length(genes))

data <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]])

# adding common names as label for network visualisation
data$nodes$label <- demo_data_At$gene_info[match(data$nodes$id, rownames(demo_data_At$gene_info)), "label"]

DIANE::draw_network(data$nodes, data$edges)


DIANE::draw_network_degrees(data$nodes, network)

head(mat)

####### go analysis

genes <- convert_from_agi(topTags$table$genes)
background <- convert_from_agi(rownames(normalized_counts))



go <- enrich_go(genes, background)
DIANE::draw_enrich_go(go, max_go = 30)




