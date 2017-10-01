#' Construction of a network graph from correlation matxix
#' @description  Construction of a correlation network graph from correlation matrix through Random Matrix Theory-based methods.
#' @usage corgraph(mat, it = seq(0.30,0.99,by=0.01))
#' @return A list consists from one igraph object and two dataframes. The one is a weighted edge list, the another one is
#'     a result from ks-test of difference eigen sequence based on thresh correlation matrix.
#' @param mat correlation matrix
#' @param it a vector of threshold sequence as iteration. The initial value as default, 'it =seq(0.30,0.99,by=0.01)'
#' @examples
#' # sample data
#' data(cluster_dat)
#' cormat <- cor(cluster_dat[[2]])
#' res <- corgraph(mat = cormat)
#' res[[1]]
#' head(res[[2]])
#' head(res[[3]])
#' @importFrom igraph graph.edgelist
#' @export
corgraph <- function(mat, it = seq(0.30,0.99,by=0.01)){
  diag(mat) <- 0
  mtx <- abs(mat)
  res.ks <- lapply(it,
                   function(i){
                     mtx_tmp <- ifelse(mtx > i, mtx,0)
                     ei <- eigen(mtx_tmp) # eigenvalue decomposition
                     diff_eigenvalue <- diff(sort(ei$values))
                     mean_diff_eigenvalue <- mean(diff_eigenvalue)
                     if(mean_diff_eigenvalue > 0){
                       diff_eigenvalue <- diff_eigenvalue / mean_diff_eigenvalue
                       # calc. Kolmogorov-Smirnov distance
                       res.ks <- suppressWarnings(stats::ks.test(diff_eigenvalue,"pexp"))
                       unlist(res.ks)[1:2]
                     }else{ NA }
                   }
  )
  res.ks.dat <- data.frame(thresh = it,
                           ks_d = as.numeric(sapply(res.ks, "[",1)),
                           ks_p = as.numeric(sapply(res.ks, "[",2)),
                           stringsAsFactors = F)
  th <- res.ks.dat$thresh[which.min(res.ks.dat$ks_d)]

  mtx_opt <- ifelse(mat > th | mat < -th, mat, 0)
  edgel <- cornet::matoedge(mtx_opt)
  g <- igraph::graph.edgelist(as.matrix(edgel[edgel$value != 0, 1:2]), directed = F )
  return(list(undir.graph = g, edge.list=edgel, res.ks.text=res.ks.dat))
}
