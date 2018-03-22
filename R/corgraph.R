#' Construction of a network graph from correlation matxix
#' @description  Construction of a correlation network graph from correlation matrix through Random Matrix Theory-based methods.
#'    if negetive correlation exists, these information were added to result of igraph object.
#' @usage corgraph(mat, thresh = seq(0.30,0.99,by=0.01))
#' @return A list consists from one igraph object and two dataframes. The one is a weighted edge list, the another one is
#'     a result from ks-test of difference eigen sequence based on thresh correlation matrix.
#' @param mat correlation matrix
#' @param thresh a vector of threshold sequence as iteration. The initial value as default, 'thresh =seq(0.30,0.99,by=0.01)'
#'     if using an arbitral threshold, set a threshold value like as "thresh = 0.85".
#' @examples
#' # # arguments
#' # nfpkm <- as.data.frame(rskodat::nfpkm)
#' # res_cld <- cornet::dycutdf(dat = nfpkm, column = -1:-4)$cluster_dat
#' # cormat <- cor(res_cld[["3"]])
#' # # create cor graph with a selected cluster
#' # res <- corgraph(mat = cormat)
#' # res[[1]]
#' # head(res[[2]])
#' # head(res[[3]])
#' @importFrom igraph get.edge.ids graph.adjacency
#' @export
corgraph <- function(mat, thresh = seq(0.30,0.99,by=0.01)){
  diag(mat) <- 0
  mtx <- abs(mat) # absolute mat

  # threshold using random matrix -----
  if(length(thresh) == 1){
    th <- thresh
    res.ks.dat <- NULL

  } else {
    res.ks <- lapply(thresh,
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
    res.ks.dat <- data.frame(thresh = thresh,
                             ks_d = as.numeric(sapply(res.ks, "[",1)),
                             ks_p = as.numeric(sapply(res.ks, "[",2)),
                             stringsAsFactors = F)
    th <- res.ks.dat$thresh[which.min(res.ks.dat$ks_d)]
    pv <- res.ks.dat$ks_p[which.min(res.ks.dat$ks_d)]
    print(paste0("thresh:", th,  "  p:", pv))
  }

  # adjmat for igraph
  mtx_opt <- ifelse(mtx > th, 1, 0)
  # create igraph object
  g <- igraph::graph.adjacency(adjmatrix = mtx_opt, mode = "undirected", diag = F)

  # negative correlation edge to matrix for searching edge from igraph object
  nedgel <- matoedge(ifelse(mat > th | mat < -th, mat, 0), format = "df")
  #nedge.mt <-  # negative correlation edge mat
  nedge.pos <- igraph::get.edge.ids(g, as.vector(t(nedgel[nedgel$value < 0, 1:2])))

  igraph::E(g)$e.c <-
    ifelse(igraph::E(g) %in% igraph::E(g)[nedge.pos], "deepskyblue", "grey80")
  return(list(undir.graph = g, edge.list = nedgel, res.ks.text=res.ks.dat))
}
