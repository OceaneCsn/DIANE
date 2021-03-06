% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_geneIDs.R
\name{check_IDs}
\alias{check_IDs}
\title{Check compatibility between gene IDs and an organism}
\usage{
check_IDs(ids, organism)
}
\arguments{
\item{ids}{character vector of gene identifiers to be tested}

\item{organism}{organism, should be betwwen "Arabidopsis thaliana", "Homo sapiens", "Mus musculus", 
"Caenorhabditis elegans", "Escherichia coli", "Drosophilia melanogaster",
"Lupinus albus"}
}
\value{
boolean, TRUE if all of the gene IDs match the expected regex for
the specified organism.
}
\description{
Check compatibility between gene IDs and an organism
}
\examples{
data("abiotic_stresses")
check_IDs(rownames(abiotic_stresses$raw_counts),
organism = "Arabidopsis thaliana")
check_IDs(rownames(abiotic_stresses$raw_counts),
organism = "Homo sapiens")
}
