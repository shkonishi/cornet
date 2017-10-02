#' simple plot of igraph object
#' @description several options has initial value, and
#' @usage igplot(ig, lay, v.l, v.l.c, v.f.c, v.s, v.c, e.c, e.w, ...)
#' @param ig igraph object
#' @param lay graph layout
#' @param v.c,v.f.c,v.l,v.l.c,v.s vertex.color, vertex.frame.color, vertex.label, vertex.label.color, vertex.size
#' @param e.c,e.w edge.color, edge.width
#' @param ... other arguments of plot.igraph
#' @examples
#' # sample data
#' dat <- data.frame(
#'  S1 =c(43.26, 166.6, 12.53, 28.77, 114.7, 119.1, 118.9, 3.76, 32.73, 17.46),
#'  S2=c(40.89, 41.87, 39.55, 191.92, 79.7, 80.57, 156.69, 2.48, 11.99, 56.11),
#'  S3=c(5.05, 136.65, 42.09, 236.56, 99.76, 114.59, 186.95, 136.78, 118.8, 21.41)
#'  )
#' rownames(dat) <- paste0("G", 1:10)
#'
#' # correlation matrix
#' cormat <- round(cor(t(dat)),2)
#'
#' # threshold graph
#' res <- cornet::corgraph(mat=cormat)
#' g1 <- res[[1]]
#' igplot(ig=g1, v.s=15)
#'
#' # complete graph
#' g2 <- cornet::matoedge(cormat)
#'  ewid <- abs(igraph::E(g2)$weight)
#'  ecol <-  ifelse(igraph::E(g2)$weight < 0 , "red", "grey80")
#' igplot(ig = g2, v.s = 15, e.c = ecol, e.w = ewid*4)
#'
#' @importFrom igraph layout_nicely V E
#' @export
igplot <- function(ig, lay=igraph::layout_nicely, v.l=igraph::V(ig)$name,
                        v.l.c = "white", v.f.c="white", v.s = 5, v.c = "turquoise3",
                        e.c = "grey80", e.w = 1, ...){
  par(mar=c(0,0,0,0))
  plot(ig,
       layout=lay,
       vertex.label = v.l,
       vertex.label.color = v.l.c,
       vertex.frame.color=v.f.c,
       vertex.size = v.s,
       vertex.color = v.c,
       edge.color = e.c,
       edge.width = e.w,
       ...
  )
}
