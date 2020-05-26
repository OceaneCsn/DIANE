#' List of regulator genes for a number of model organisms
#'
#' @source For Arabidopsis thaliana : http://plntfdb.bio.uni-potsdam.de/v3.0/index.php?sp_id=ATH
#' For Human : http://humantfs.ccbr.utoronto.ca/download.php
#' @format A named list with the following elements:
#' \describe{
#'  \item{Arabidopsis thaliana}{character vector of the known regulators in Arabidopsis thaliana}
#'  \item{Homo sapiens}{character vector of the known transcription factors in Human}
#' }
#' @examples
#' {
#'  head(regulators_per_organism[["Arabidopsis thaliana"]])
#'  head(regulators_per_organism[["Homo sapiens"]])
#' }
"regulators_per_organism"