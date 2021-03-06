#' simple plot of igraph object
#' @description Shortening several 'plot.igraph' options, and these has initial value.
#' @usage igplot(g, lay, connected, v.l, v.l.c, v.l.cx, v.f.c, v.s, v.c, v.l.f, v.shp,
#'     e.c, e.w, e.lty, ...)
#' @param g igraph object
#' @param lay layout function of igraph, if 'lay="all"', all layout functions were performed.
#' @param connected logical: default is true, connected vertex was shown.
#' @param v.c,v.f.c,v.l,v.l.c,v.l.cx,v.s,v.l.f,v.shp, vertex parameters, v.c(vertex.color), v.f.c(vertex.frame.color), v.l(vertex.label), v.l.c(vertex.label.color), v.l.cx(vertex.label.cex), v.s(vertex.size), v.l.f(vertex.label.family), v.shp(vertex.shape)
#' @param e.c,e.w,e.lty edge parameters, e.c(edge.color), e.w(edge.width), e.lty(edge.lty)
#' @param ... other arguments of plot.igraph. E.g. margin, frame, and main
#' @examples
#' # sample data
#' dat <- data.frame(
#'  S1 = c(43.26, 166.6, 12.53, 28.77, 114.7, 119.1, 118.9, 3.76, 32.73, 17.46),
#'  S2 = c(40.89, 41.87, 39.55, 191.92, 79.7, 80.57, 156.69, 2.48, 11.99, 56.11),
#'  S3 = c(5.05, 136.65, 42.09, 236.56, 99.76, 114.59, 186.95, 136.78, 118.8, 21.41)
#'  )
#' rownames(dat) <- paste0("G", 1:10)
#'
#' # correlation matrix
#' cormat <- round(cor(t(dat)),2)
#'
#' # threshold graph
#' res <- cornet::corgraph(mat=cormat)
#' g1 <- res[[1]]
#' cornet::igplot(g=g1, v.s=15)
#'
#' # complete graph
#' g2 <- cornet::matoedge(cormat)
#'  ewid <- abs(igraph::E(g2)$weight)
#'  ecol <-  ifelse(igraph::E(g2)$weight < 0 , "steelblue3", "grey80")
#' cornet::igplot(g = g2, lay=igraph::layout.circle, v.s = 15, e.c = ecol, e.w = ewid*4)
#' @importFrom igraph graph.union decompose.graph components layout_nicely V E layout.auto layout.bipartite layout.circle layout.davidson.harel layout.drl layout.fruchterman.reingold layout.fruchterman.reingold.grid layout.gem layout.graphopt layout.grid layout.grid.3d layout.kamada.kawai layout.lgl layout.mds layout.merge layout.norm layout.random layout.reingold.tilford layout.sphere layout.spring layout.star layout.sugiyama layout.svd
#' @importFrom utils lsf.str
#' @importFrom graphics par plot
#' @export
igplot <- function(g, lay = igraph::layout_nicely, connected = TRUE,
                   v.l = NULL, v.l.c = NULL, v.l.cx = 0.8, v.f.c = NA,
                   v.s = 5, v.c = "#8B887880", v.l.f ="Helvetica", v.shp = NULL,
                   e.c = "grey30", e.w = 0.5, e.lty = 1, ...){

  # if igraph object has attibutes of vertices and edge, initial value was replaced.
  v.c <- if(all(v.c =="#8B887880") & !is.null(igraph::V(g)$v.c)){igraph::V(g)$v.c}else{v.c}
  v.l <- if(is.null(v.l) & !is.null(igraph::V(g)$v.l)){igraph::V(g)$v.l}else{v.l}
  v.l.c <- if(is.null(v.l.c) & !is.null(igraph::V(g)$v.l.c)){igraph::V(g)$v.l.c}else{v.l.c}
  v.l.cx <- if(all(v.l.cx == 0.8) & !is.null(igraph::V(g)$v.l.cx)){igraph::V(g)$v.l.cx}else{v.l.cx}
  v.s <- if(all(v.s == 5) & !is.null(igraph::V(g)$v.s)){igraph::V(g)$v.s}else{v.s}
  v.f.c <- if(all(is.na(v.f.c)) & !is.null(igraph::V(g)$v.f.c)){igraph::V(g)$v.f.c}else{v.f.c}
  v.shp <- if(is.null(v.shp) & !is.null(igraph::V(g)$v.shp)){igraph::V(g)$v.shp}else{v.shp}
  e.c <- if(all(e.c == "grey80") & !is.null(igraph::E(g)$e.c)){igraph::E(g)$e.c}else{e.c}
  e.w <- if(all(e.w == 0.5) & !is.null(igraph::E(g)$e.w)){igraph::E(g)$e.w}else{e.w}
  e.lty <- if(all(e.lty == 1) & !is.null(igraph::E(g)$e.lty)){igraph::E(g)$e.lty}else{e.lty}

  # Only connected vertex was shown
  if (connected == TRUE){
    g <- igraph::delete.vertices(g, igraph::V(g)$name[igraph::degree(g)==0])
  }

  # all layout function
  lays <- grep("^layout_.", unclass(utils::lsf.str(envir = asNamespace("igraph"), all = T)), value = T)

  # draw graph
  if (!is.function(lay) & !is.matrix(lay)){
    if (lay == "all"){
      suppressWarnings(
        for (x in lays) {
          #ERROR HANDLING
          possibleError <- tryCatch({
            coords <- eval(parse(text= paste0("igraph::", x, "(g)")))
            par(mar=c(0,0,1,0))
            igraph::plot.igraph(g, layout=coords,
                                vertex.label = v.l, vertex.label.color = v.l.c, vertex.frame.color=v.f.c,
                                vertex.label.cex = v.l.cx, vertex.size = v.s, vertex.color = v.c,
                                vertex.label.family = v.l.f, vertex.shape = v.shp,
                                edge.color = e.c, edge.width = e.w, edge.lty = e.lty, edge.lty=e.lty,
                                main = x, ...)
          },
          error = function(e) {
            print(paste("Error in layout ", x, sep = ""))
          })
          if (inherits(possibleError, "error")) next
        }
      )
    }

  } else if (is.function(lay) | is.matrix(lay)){
    par(mar=c(0,0,1,0))
    plot(g, layout=lay, vertex.label = v.l, vertex.label.color = v.l.c, vertex.label.cex = v.l.cx,
         vertex.frame.color=v.f.c, vertex.size = v.s, vertex.color = v.c,
         vertex.label.family = v.l.f, vertex.shape = v.shp,
         edge.color = e.c, edge.width = e.w, edge.lty=e.lty, ...
    )
  } else {
    stop("'lay' must to be a function of layout, or the return value of layout function.")
  }
}
