#' Expression matrix from Arabidopsis thaliana.
#'
#' @source BPMP Institute, Supagro Montpellier
#' @format A list with the following elements:
#' \describe{
#'  \item{raw_counts}{raw transcript aboundances as obtained after mapping and quantifying RNASeq reads (dataframe of integer positive values).}
#'  \item{conditions}{Experimental conditions of each sample, regardless of biological replicates. (character vector)}
#'  \item{design}{For each condition, the level of the corresponding factors.}
#' }
#' @examples
#' {
#'  head(demo_data_At$raw_counts)
#'  demo_data_At$conditions
#'  demo_data_At$design
#' }
"demo_data_At"