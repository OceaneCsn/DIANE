################# seed tests for DIANE

# Every step of the pipeline is supposed to return 
# the same results from two same runs


# clustering (initialization of mixtures) OK


#set.seed(12)
data("abiotic_stresses")
genes <- abiotic_stresses$heat_DEGs
clustering <- run_coseq(conds = unique(abiotic_stresses$conditions), 
                        data = abiotic_stresses$normalized_counts, genes = genes, K = 6:8,
                        seed = 123)


conds <- abiotic_stresses$conditions
DIANE::draw_profiles(data = abiotic_stresses$normalized_counts, 
                     membership = clustering$membership,
                     conds = unique(abiotic_stresses$conditions)) 

data("abiotic_stresses")
genes <- abiotic_stresses$heat_DEGs
clustering <- run_coseq(conds = unique(abiotic_stresses$conditions), 
                        data = abiotic_stresses$normalized_counts, genes = genes, K = 6:8, seed = 123)


conds <- abiotic_stresses$conditions
DIANE::draw_profiles(data = abiotic_stresses$normalized_counts, 
                     membership = clustering$membership,
                     conds = unique(abiotic_stresses$conditions)) 




## network inference

# grouping des TFs correlés avec Louvain OK
# apparently it is deterministic, or the seed is enforced in 
# the package, as a have reproducible results directly

aggregated_data <- aggregate_splice_variants(abiotic_stresses$normalized_counts)
genes <- get_locus(abiotic_stresses$heat_DEGs)
regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])
grouping1 <- group_regressors(aggregated_data, genes, regressors)
grouping2 <- group_regressors(aggregated_data, genes, regressors)

length(intersect(grouping1$grouped_genes, grouping2$grouped_genes))






# standard importance metric (calls c++ code) OK

aggregated_data <- aggregate_splice_variants(data = abiotic_stresses$normalized_counts)

genes <- get_locus(abiotic_stresses$heat_DEGs)
regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])

set.seed(123)
mat1 <- network_inference(aggregated_data, conds = abiotic_stresses$conditions,
targets = genes, regressors = regressors, nTrees = 1000, nCores = 4)
set.seed(123)
mat2 <- network_inference(aggregated_data, conds = abiotic_stresses$conditions,
                          targets = genes, regressors = regressors, nTrees = 1000, nCores = 4)

cor(c(mat1[regressors, genes]), c(mat2[regressors, genes]))


# oob prediction loss metric (calls R code) OK


aggregated_data <- aggregate_splice_variants(data = abiotic_stresses$normalized_counts)

genes <- get_locus(abiotic_stresses$heat_DEGs)
regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])

set.seed(123)
mat1 <- network_inference(aggregated_data, conds = abiotic_stresses$conditions,
                          targets = genes, regressors = regressors, nTrees = 1000, nCores = 4,
                          importance_metric = "MSEincrease_oob")
set.seed(123)
mat2 <- network_inference(aggregated_data, conds = abiotic_stresses$conditions,
                          targets = genes, regressors = regressors, nTrees = 1000, nCores = 4,
                          importance_metric = "MSEincrease_oob")

cor(c(mat1[regressors, genes]), c(mat2[regressors, genes]))






# network thresholding
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
get_nEdges(0.001, length(r$grouped_genes), length(r$grouped_regressors))

set.seed(123)
res <- DIANE::test_edges(mat, normalized_counts = r$counts, density = 0.001,
                        nGenes = length(r$grouped_genes),
                        nRegulators = length(r$grouped_regressors),
                        nTrees = 1000, verbose = TRUE)
set.seed(123)
res2 <- DIANE::test_edges(mat, normalized_counts = r$counts, density = 0.001,
                          nGenes = length(r$grouped_genes),
                          nRegulators = length(r$grouped_regressors),
                          nTrees = 1000, verbose = TRUE)

l1 <- res$links
l2 <- res2$links
# louvain module detection OK

