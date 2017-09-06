#' Construction of MIC graph
#' @usage micgraph(dat, thresh)
#' @param dat data frame or matrix
#' @param thresh threshold of correlation or MIC
#' @return  data frame of edge list of
#' @examples

micgraph <- function(dat, threshold){
  library(minerva)
  micg <- minerva::mine(t(d))
  micg <- micg^2
  diag(micg) <- 0
  micg[which(micg >= threshold )  ] <- 1
  micg[which(micg <  threshold )  ] <- 0
  g <- igraph::graph.adjacency( micg, mode="undirected", weighted=TRUE)
  g <- igraph::simplify(g, remove.multiple=T, remove.loops=T)
  cat(paste("node :", vcount(g), "\n",
            "edge :", ecount(g), "\n",
            "transitivity :", transitivity(g), "\n",
            "degree :", sum(igraph::degree(g))/vcount(g), "\n",
            "density :", graph.density(g), "\n",sep=""))
  return(g)
}
