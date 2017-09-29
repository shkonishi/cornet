#' The reshape of mine output
#' @description The result of mine is multi matrix of all node pairs, and that is too many memory usage using big data.
#'    this program retuns dataframe as edge list and atrributes of these edge.
#' @usage cluster_mine(cl_dat)
#' @return A list of dataframes, result of minerva::mine soreted by 'TIC'.
#' @param cl_dat A list of dataframes, result of 'cornet::cluster_mat'
#' @examples
#' # sample data, result of 'cluster_mat'
#' data("cluster_dat")
#' res <- cluster_mine(cl_dat=cluster_dat[[3]])
#'
#' @importFrom minerva mine
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @export
cluster_mine <- function(cl_dat){
  # argument check ----
  if(class(cl_dat) != "list" & class(cl_dat) != "data.frame"){
    stop("'cl_dat' is a dataframe or a list of multiple datafrme.")
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
    r <- cornet::matoedge(cor(x, method="pearson"))[,3]
    rho <- cornet::matoedge(cor(x, method="spearman"))[,3]
    dat <- data.frame(edge, mic, mas, mev, mcn, micr2, gmic, tic, pearson=r, spearman=rho)
    dat <- dat[order(dat$tic, decreasing = T),]
  }

  if (class(cl_dat) == "list" & all(sapply(cl_dat, class)== "data.frame")){
    cl <- parallel::makeCluster(parallel::detectCores())
    doParallel::registerDoParallel(cl)
    x <- NULL
    sum_mine <- foreach::foreach(x = cl_dat) %dopar% {f(x)}
    parallel::stopCluster(cl)
    names(sum_mine) <- names(cl_dat)
  } else {
    sum_mine <- f(cl_dat)
  }

  return(sum_mine)
}
