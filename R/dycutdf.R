#' Hierarchical clustering using amp, and split dataframe by cluster, with dynamic cut.
#' @description dycutdf
#' @return list of results of clustering, objects was hclust, dendrogram, and list of splitted data.frame
#' @usage dycutdf(dat, distm, clm, column, method_dycut, ...)
#' @param dat data frame or matrix
#' @param distm distance measure from amap::Dist options. default is "pearson"
#' @param clm clustering method select from hclust options default is "average"
#' @param column column which containing expression data default 'column=1:ncol(dat)'
#' @param method_dycut method of dynamic cut, otherwise k number of cuttree function.
#' @param ... additional dynamicTreeCut::cutreeDynamic options
#' @examples
#' d <- data.frame(t(iris[-5]))
#' res <- dycutdf(dat = d, distm = "euclidean", clm = "average",
#'                column = 1:150, method_dycut = "hybrid")
#' res[[1]] # hclust object
#' res[[2]] # result of dynamic cut object
#' res[[3]] # list of splitted data.frame
#'
#' ##  The following examples takes too much times
#' # fp <- system.file("extdata/nfpkm_200.txt", package = "cornet")
#' # nfpkm <- read.table(fp, header=TRUE, stringsAsFactors = FALSE)
#' # res1 <- cornet::dycutdf(dat = nfpkm, distm = "abscorrelation", clm = "average",
#' #     column = 5:ncol(nfpkm), method_dycut = "tree")
#' @import dplyr
#' @import ggplot2
#' @importFrom amap Dist
#' @importFrom dendextend hang.dendrogram
#' @importFrom stats hclust as.dendrogram median cutree
#' @export
#'
dycutdf <- function(dat, distm="pearson", clm="average", column=1:ncol(dat), method_dycut="tree", ...){
  # argument check: dat ----
  if (class(dat)=="data.frame"){
    mat <- dat[column]
  } else if (class(dat) == "matrix"){
    mat <- data.frame(dat)
    if(identical(dimnames(dat)[[2]], NULL)){
      names(mat) <- 1:ncol(mat)
    }
  }

  # hierarchical clustering ----
  distms <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman", "kendall")
  if (any(distms %in% distm)){
    dis.mat <- amap::Dist(t(mat), method = distm)
    r_hcl <- stats::hclust(dis.mat, method = clm)
  } else if (class(distm)=="dist"){
    dis.mat <- distm
    r_hcl <- stats::hclust(dis.mat, method = clm)
  } else {
    stop('"distm" is select from "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman", "kendall", \n
         or your dist object.')
  }
  # dendrogram object ----
  r_den <- dendextend::hang.dendrogram(as.dendrogram(r_hcl))
  #### plot(r_den, leaflab="none")

  # cutreeDynamic  ----
  if (method_dycut=="tree"){
    resdyct <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, method = "tree", ... )
  } else if (method_dycut =="hybrid") {
    resdyct <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, distM = as.matrix(dis.mat), method = "hybrid", ... )
  } else if (is.integer(method_dycut)){
    resdyct <- cutree(r_hcl, k=method_dycut)
  }

  # leaf label (leaf order) ----
  names(resdyct) <- r_hcl$labels
  resdyct.lo <- resdyct[r_hcl$order]

  # split cluster and data frame corresponding to cluster, and   ----
  cl_dat <- lapply(split(names(mat), factor(resdyct)), function(x)mat[x])

  # result data shaping ----
  return(list(hclust_obj = r_hcl, dend_obj = r_den, dynamic_cut=resdyct, cluster_dat = cl_dat))

}
