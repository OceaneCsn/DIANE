
# Poisson mixture models for clustering


The [coseq package](https://www.bioconductor.org/packages/release/bioc/vignettes/coseq/inst/doc/coseq.html) performs clustering on genes depending on their expression profiles. 

It relies on the statistical framework of **mixture models**. Each cluster of genes is represented by a Poisson distribution, which parameter is estimated using Expectation-Maximisation algorithms. This algorithm ends up with the **Poisson** distributions parameters that are the best fit to the data, as they maximises the complete likelihood of the clustering.

One EM algorithm is launched for each number of clusters specified in the input range. In the end of the procedure, the number of clusters that maximizes the clustering quality is chosen. The tab "clustering quality" will give more details about this criterion.

The input genes of a clustering can't be all the genes of the data, we use the output of differential expression analysis instead. Indeed, clustering genes that have similar expression accros conditions would not be very informative, and would increase computation time for no reason.

You can also specify a subset of the conditions to be used during the clustering if you're not interested in all the conditions, but a minimum of three different conditions recommended to expect an interesting result.

