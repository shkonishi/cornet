#' The reshape of mine output
#' @description The result of mine is multi matrix of all node pairs, and that is too many memory usage using big data.
#'    this program retuns data frame as edge list and atrributes of these edge.
#' @usage cluster_mine(cl_dat)
#' @return A list of data frames, result of minerva::mine.
#' @param cl_dat A list of data frames, result of 'cornet::cluster_mat'
#' @examples
#' # sample data, result of 'cluster_mat'
#' @importFrom minerva mine
#' @importFrom foreach foreach %dopar%
#' @export
cluster_mine <- function(cl_dat){
  # argument check: cl_dat
  if(class(cl_dat)=="list"){
    loops <- seq_along(cl_dat)
  } else if (class(cl_dat)=="data.frame"){
    loops <- 1
  } else {
    stop("'cl_dat' is a list of multi data frame, or a data frame")
  }

  # parse result of minerva::mine ----
  f <- function(x){
    res.mine <-  minerva::mine(x, alpha = 0.6, n.cores = 2) # mine
    # create edge list
    edge <- cornet::matoedge(res.mine$MIC)[1:2]
    mic <- cornet::matoedge(res.mine$MIC)[,3] # mic
    mas <- cornet::matoedge(res.mine$MAS)[,3] # mas
    mev <- cornet::matoedge(res.mine$MEV)[,3] # mev
    mcn <- cornet::matoedge(res.mine$MCN)[,3] # mcn
    micr2 <- cornet::matoedge(res.mine$MICR2)[,3] # micr2
    gmic <- cornet::matoedge(res.mine$GMIC)[,3] # gmic
    tic <- cornet::matoedge(res.mine$TIC)[,3] # tic
    dat <- data.frame(edge, mic, mas, mev, mcn, micr2, gmic, tic)
    dat <- dat[order(dat$tic, decreasing = T),]
  }
  x <- NULL
  sum_mine <- foreach::foreach(x = cl_dat) %dopar% {f(x)}

  names(sum_mine) <- names(cl_dat)
  return(sum_mine)
}
