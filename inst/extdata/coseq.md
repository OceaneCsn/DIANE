
# Poisson mixture models for clustering

---

The [coseq package](https://www.bioconductor.org/packages/release/bioc/vignettes/coseq/inst/doc/coseq.html) performs clustering on genes depending on their expression profiles. 


Coseq package relies on the statistical framework of mixture models. Each cluster of genes is represented by a Poisson or Gaussian distribution, which parameters is estimated using Expectation-Maximisation algorithms. For Gaussian mixtures, a prior transformation of the count data is required to bring it closer to a normally distributed measure.

One EM algorithm is launched for each number of clusters specified in the input range. 
In the end of the procedure, the number of clusters that maximizes the clustering quality is chosen.  The tab "clustering quality" will give more details about this criterion.

The input genes of a clustering can't be all the genes of the data, we use the output of differential expression analysis instead.

You can also specify a subset of the conditions to be used during the clustering if we're not interested in all the conditions.