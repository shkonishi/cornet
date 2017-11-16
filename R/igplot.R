#' simple plot of igraph object
#' @description Shortening several 'plot.igraph' options, and these has initial value.
#' @usage igplot(ig, lay, connected, v.l, v.l.c, v.f.c, v.s, v.c, e.c, e.w, e.lty, ...)
#' @param ig igraph object
#' @param lay layout function of igraph, if 'lay="all"', all layout functions were performed.
#' @param connected logical: default is true, connected vertex was shown.
#' @param v.c,v.f.c,v.l,v.l.c,v.s vertex parameters, v.c(vertex.color), v.f(vertex.frame.color), v.l(vertex.label), v.l.c(vertex.label.color), v.s(vertex.size)
#' @param e.c,e.w,e.lty edge parameters, e.c(edge.color), e.w(edge.width), e.lty(edge.lty)
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
#'  ecol <-  ifelse(igraph::E(g2)$weight < 0 , "steelblue3", "grey80")
#' igplot(ig = g2, lay=igraph::layout.circle, v.s = 15, e.c = ecol, e.w = ewid*4)
#' @importFrom igraph graph.union decompose.graph components layout_nicely V E layout.auto layout.bipartite layout.circle layout.davidson.harel layout.drl layout.fruchterman.reingold layout.fruchterman.reingold.grid layout.gem layout.graphopt layout.grid layout.grid.3d layout.kamada.kawai layout.lgl layout.mds layout.merge layout.norm layout.random layout.reingold.tilford layout.sphere layout.spring layout.star layout.sugiyama layout.svd
#' @importFrom utils lsf.str
#' @importFrom graphics par plot
#' @export
igplot <- function(ig, lay = igraph::layout_nicely, connected = TRUE,
                   v.l = igraph::V(ig)$name, v.l.c = "white", v.f.c = "white",
                   v.s = 5, v.c = "#8B887880",
                   e.c = "grey80", e.w = 1, e.lty = 1, ...){

  # if igraph object has attibutes of vertices and edge, intial value was replaced.
  v.c <- if(!is.null(V(ig)$color)){V(ig)$color}else{v.c}
  e.c <- if(!is.null(E(ig)$color)){E(ig)$color}else{e.c}

  # Only connected vertex was shown
  if (connected == TRUE){
    ig <- delete.vertices(ig, V(ig)$name[degree(ig)==0])
  }

  if (!is.function(lay)){
    if (lay == "all"){
      lays <- grep("^layout\\.", unclass(utils::lsf.str(envir = asNamespace("igraph"), all = T)), value = T)
      suppressWarnings(
        for (x in lays) {
          #ERROR HANDLING
          possibleError <- tryCatch({
            coords <- do.call(x, list(ig))
            par(mar=c(0,0,1,0))
            igraph::plot.igraph(ig, layout=coords,
                                vertex.label = v.l, vertex.label.color = v.l.c, vertex.frame.color=v.f.c,
                                vertex.size = v.s, vertex.color = v.c, edge.color = e.c,
                                edge.width = e.w, edge.lty=e.lty, main = x, edge.lty=e.lty, ...)
          },
          error=function(e) {
            e
            print(paste("Error in layout ",x,sep = ""))
          }
          )
          if (inherits(possibleError, "error")) next
        }
      )
    }
  } else {
    par(mar=c(0,0,1,0))
    plot(ig, layout=lay, vertex.label = v.l, vertex.label.color = v.l.c,
         vertex.frame.color=v.f.c, vertex.size = v.s, vertex.color = v.c,
         edge.color = e.c, edge.width = e.w, edge.lty=e.lty, ...
    )
  }
}
