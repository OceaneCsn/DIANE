% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regulators_per_organism.R
\docType{data}
\name{regulators_per_organism}
\alias{regulators_per_organism}
\title{List of regulator genes for a number of model organisms}
\format{
A named list with the following elements:
\describe{
 \item{Arabidopsis thaliana}{Character vector of the known transcriptional regulators in Arabidopsis thaliana}
 \item{Homo sapiens}{Character vector of the known transcription factors in Human}
 \item{Mus musculus}{Character vector of the known transcription factors in Mouse}
 \item{Drosophilia melanogaster}{Character vector of the known transcription factors in Drosophilia melanogaster}
 \item{Caenorhabditis elegans}{Character vector of the known transcription factors in C elegans}
 \item{Escherichia coli}{Character vector of the known transcription factors and sigma elements in E coli K12}
}
}
\source{
\describe{
 \item{For Arabidopsis thaliana}{\url{http://plntfdb.bio.uni-potsdam.de/v3.0/index.php?sp_id=ATH}}
 \item{For Human}{\url{http://humantfs.ccbr.utoronto.ca/download.php}}
 \item{Mus musculus}{\url{http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download}}
 \item{Drosophilia melanogaster}{\url{https://flybase.org/reports/FBgg0000745.html}}
 \item{C elegans}{\url{http://www.wormbook.org/chapters/www_transcriptionalregulation.2/transregulate.html#sec5-1}}
 \item{E coli}{\url{http://regulondb.ccg.unam.mx/menu/download/datasets/index.jsp}}
 }
}
\usage{
regulators_per_organism
}
\description{
Those particular genes are known to be transcriptionnal actors,
and are needed to infer gene regulatory networks in DIANE.
}
\examples{
{
 print(head(regulators_per_organism[["Arabidopsis thaliana"]]))
 print(head(regulators_per_organism[["Homo sapiens"]]))
 print(head(regulators_per_organism[["Mus musculus"]]))
 print(head(regulators_per_organism[["Drosophilia melanogaster"]]))
 print(head(regulators_per_organism[["Caenorhabditis elegans"]]))
 print(head(regulators_per_organism[["Escherichia coli"]]))
}
}
\keyword{datasets}
