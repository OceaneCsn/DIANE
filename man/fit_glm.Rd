% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_glm.R
\name{fit_glm}
\alias{fit_glm}
\title{Fits a Poisson generalized linear model to a set of genes}
\usage{
fit_glm(normalized_counts, genes, design, factors = colnames(design))
}
\arguments{
\item{normalized_counts}{normalized counts}

\item{genes}{genes belonging to a specific expression-based clusters}

\item{design}{experimental design as a dataframe}

\item{factors}{factors to use for the fit (defalut is
all the factors of the design)}
}
\value{
glm object
}
\description{
The idea is to extract the importance and effect of each factor.
To do so, the expression of each gene is modeled as a Poisson distribution.
The log of its parameter (the expected value) is approximated by a linear
combination of the factors in the experiment. The coefficients associated to 
each factors are estimated to fit gene expression, and can be insightful to
characterize genes behavior in a particular cluster.
The model with interactions is considered. It your design in not a
complete crossed design, the interaction term will be null.
}
\section{Note }{

Note that we can only apply a glm fit to a set of genes that have very close expression 
profiles accros conditions, else we would have to introduce a new variable related to the genes
themselves.
}

\examples{
data("abiotic_stresses")
genes_cluster <- DIANE::get_genes_in_cluster(
abiotic_stresses$heat_DEGs_coseq_membership, cluster = 3)
glm <- DIANE::fit_glm(abiotic_stresses$normalized_counts, genes_cluster, 
abiotic_stresses$design)
summary(glm)
}
