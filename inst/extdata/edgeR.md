# Differential expression methods

---

To detect the genes that have significant changes in their expression caused by experimental perturbations, we use the [edgeR package](http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf#section.2.10), based on negative binomals models.

The idea is to model gene expression using Poisson ditributions, well suited for count data. Because biological replicates induce an overdispersion of gene counts, the preferred model is the Negative binomial. This distribution has 2 parameters, its expected value, and an overdispersion paramater.

In order to detect differential expression, these parameters have to be estimated first.

The first step is thus to estimate the gene dispersions, which is acheived by pooling genes with similar expression level, and using empirical Bayes stragtegies.

Then, the expected value of each gene is estimate using genewise generalized linear models, as a combination of the coefficients associeted to each factor of the experiment (or perturbation).

Likelihood ratio tests can then be conducted on the fitted model coefficients to determine weather or not they are satistically different from each other, and thus conclude on differential expression.

In DIANE, the dispersion estimation and model fitting is done once, and statistical tests for differentially expressed genes can be done for different contrasts.
To do so, select the conditions that you want to compare for differential expression, and they will be the one tested againt one another via likelihood ratio tests.

The results are presented in a dataframe, ordered by adjusted pvalues (FDR). The dataframe contains the log fold changes (logFC), the average expression (logCPM) for each genes which FDR is lower than the specified adjusted p-value threshold.
You can also choose to to select one genes having an absolute log fold change over a certain constant.