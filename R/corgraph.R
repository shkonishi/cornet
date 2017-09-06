#' Construction of correlation network from multivariate data
#' @usage micgraph(dat, thresh)
#' @param dat data frame or matrix
#' @param thresh threshold of correlation or MIC
#' @return  data frame of edge list of
#' @examples
#'
corgraph <- function(d, cor.method, threshold){
  # test
  d = t(mtcars) # dynamic cut data set
  cor.method = "spearman"

  cor_mat <- cor(d, method=cor.method)
  diag(cor_mat) <- 0

  hist(cor_mat, breaks = seq(0,1, 0.01))

  g <- igraph::graph.adjacency(cor_mat, weighted=TRUE, mode="lower")

  cord[which(cord >= threshold )  ] <- 1
  cord[which(cord <  threshold )  ] <- 0
  g <- graph.adjacency( cord, mode="undirected", weighted=TRUE)
  g <- igraph::simplify(g, remove.multiple=T, remove.loops=T)
  cat(paste("node :", vcount(g), "\n",
            "edge :", ecount(g), "\n",
            "transitivity :", transitivity(g), "\n",
            "degree :", sum(igraph::degree(g))/vcount(g), "\n",
            "density :", graph.density(g), "\n",sep=""))
  return(g)
}
