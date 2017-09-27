#' The utility of minerva package
#' @description The result of mine is multi matrix of all node pairs, and that is too many memory usage using big data.
#'    this program retuns data frame as edge list and atrributes of these edge.
#' @usage cluster_mine(cl_dat)
#' @return A list of data frames, result of minerva::mine.
#' @param cl_dat A list of data frames, result of 'cornet::cluster_mat'
#' @examples
#' # sample data, result of 'cluster_mat'
#' @importFrom minerva mine
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

  sum_mine <-
    lapply(loops,
           function(i){
             x <- cl_dat[[i]] # data frame of respectiv clusters
             res.mine <-  minerva::mine(x, alpha = 0.6, n.cores = 2) # mine
             mic <- res.mine$MIC; mic <- mic[lower.tri(mic)]
             mas <- res.mine$MAS; mas <- mas[lower.tri(mas)]
             mev <- res.mine$MEV; mev <- mev[lower.tri(mev)]
             mcn <- res.mine$MCN; mcn <- mcn[lower.tri(mcn)]
             micr2 <- res.mine$MICR2; micr2 <- micr2[lower.tri(micr2)]
             gmic <- res.mine$GMIC; gmic <- gmic[lower.tri(gmic)]
             tic <- res.mine$TIC; tic <- tic[lower.tri(tic)]
             edge <- mat_to_edge(res.mine$MIC)
             dat <- data.frame(edge, mic, mas, mev, mcn, micr2, gmic, tic)
           }
    )
  names(sum_mine) <- names(cl_dat)
  return(sum_mine)
}
