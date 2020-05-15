# Clustering quality descriptors



The [coseq package](https://www.bioconductor.org/packages/release/bioc/vignettes/coseq/inst/doc/coseq.html) tests a range of different clusters in order to give the best fit to the data.

For each number of cluster, the ICL (Integrated Completed Likelihood) is computed. It combines two elements : 

+ The global **likelihood** of the clustering. It quantifies how accurate the clustering seems, regarding the posterior probability of each element to belong the its predicted cluster. It can be computed using the Poisson probability densities resulting from the proposed clustering, for all the genes.

+ The **number of clusters**. As the likelihood tends to grow monotonously with the number of clusters, resulting in a very big number of groups, that would not be very informative for the user. Thus, the ICL penalizes the clustering quality criteria with the number of clusters.

This is why the maximal value of ICL can be interpreted as an approximation of the ideal number of clusters. 