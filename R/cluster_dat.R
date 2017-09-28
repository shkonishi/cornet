#' The example of cluster_mat output
#'
#' @description The list of data frame which result of 'cluster_mat.R'.
#'     These of dataframe objects cutting of hieralchicalclustering.
#'
#' @docType data
#' @usage data(cluster_dat)
#' @format An object of list with 4 length.
#' @keywords datasets
#' @examples
#' data(cluster_dat)
#' lapply(cluster_dat, dim)
#' lapply(cluster_dat, function(x){x[1:6,1:6]})
#'
#'
"cluster_dat"
