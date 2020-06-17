#' Combined abiotic stresses
#' 
#' 
#' This dataset contians the transcriptome of Arabidopsis thaliana plants exposed to global warming 
#' induced conditions. It was generated in the article "Molecular plant responses to combined abiotic 
#' stresses put a spotlight on unknown and abundant genes", by Sewelam et Al. in Journal of experimental
#' Botany, 2020.
#' The experimental perturbations studied are high tempreature, hight salinity and osmotic changes in 
#' the soil. Each factors has two levels, one of them considered as the reference, and the other one as 
#' the stress level.
#' 
#'
#' @source 
#' 
#' \describe{
#'  \item{Data download : }
#'  {\href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146206}{GEO accession}}
#'  \item{Publication :  }
#'  {\href{https://academic.oup.com/jxb/advance-article/doi/10.1093/jxb/eraa250/5842162#204457392}{Molecular 
#'  plant responses to combined abiotic stresses put a spotlight on unknown and abundant genes}}
#' }
#' 
#' @format A named list with the following elements:
#' \describe{
#'  \item{raw_counts}{Dataframe of positive values. Raw transcript aboundances 
#'  as obtained after mapping and quantifying RNASeq reads. Rows are transcripts, and
#'  columns are experimental triplicate conditions}
#'  \item{design}{Dataframe. Describes for each condition, the level of each factor.}
#'  \item{conditions}{Character vector. Gives the condition corresponding to each 
#'  column of the raw_counts element}
#' }
#' @examples
#' {
#'  head(abiotic_stresses$raw_counts)
#'  abiotic_stresses$design
#'  abiotic_stresses$conditions
#' }
"abiotic_stresses"