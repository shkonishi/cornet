#' The example of dycutdf output
#'
#' @description The list of data frame which result of 'dycutdf.R'.
#'     These of dataframe objects cutting of hieralchicalclustering.
#'
#' @docType data
#' @usage data(cl_dat)
#' @format An object of list with 10 length.
#' @keywords datasets
#' @examples
#' data(cl_dat)
#' lapply(cl_dat, dim)
#' lapply(cl_dat, function(x){x[1:3,1:3]})
#'
#' # search target genes conataining cluster ----
#' sct <- c("gene1", "gene2")
#' sapply(cl_dat, function(x)match(sct, names(x)))
"cl_dat"
