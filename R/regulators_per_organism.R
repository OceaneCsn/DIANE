#' List of regulator genes for a number of model organisms
#'
#'Those particular genes are known to be transcriptionnal actors,
#'and are needed to infer gene regulatory networks in DIANE.
#'
#' @source 
#' \describe{
#'  \item{For Arabidopsis thaliana}{\url{http://plntfdb.bio.uni-potsdam.de/v3.0/index.php?sp_id=ATH}}
#'  \item{For Human}{\url{http://humantfs.ccbr.utoronto.ca/download.php}}
#'  \item{Mus musculus}{\url{http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download}}
#'  }
#' @format A named list with the following elements:
#' \describe{
#'  \item{Arabidopsis thaliana}{Character vector of the known regulators in Arabidopsis thaliana}
#'  \item{Homo sapiens}{Character vector of the known transcription factors in Human}
#'  \item{Mus musculus}{Character vector of the known transcription factors in Mouse}
#' }
#' @examples
#' {
#'  print(head(regulators_per_organism[["Arabidopsis thaliana"]]))
#'  print(head(regulators_per_organism[["Homo sapiens"]]))
#'  print(head(regulators_per_organism[["Mus musculus"]]))
#' }
"regulators_per_organism"