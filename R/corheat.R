#' Gene clustering and cutting cluster with dynamic cut
#' @description Gene clustering for RNA-seq data and cutting cluster with dynamic cut
#' @usage corheat(dat, distm, clm, method_dycut, draw)
#' @return draw heatmap and results of dynamic cut
#' @param dat RNA-seq count data(row:sample, column:gene)
#' @param distm select from c("pearson", "abspearson", "squarepearson", "spearman")
#' @param clm clustering method of hclust
#' @param method_dycut method of 'cutreeDynamic' select from "tree", or "hybrid
#' @param draw logical draw heatmap or nod.
#' @importFrom amap Dist
#' @importFrom grDevices hcl heat.colors
#' @importFrom gplots heatmap.2
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom dendextend hang.dendrogram
#' @importFrom stats as.dist cor
#' @examples
#' # data: normalized fpkm
#' fp <- system.file("extdata/nfpkm_rnsq.txt", package = "rsko")
#' nfpkm <- read.table(fp, header=TRUE, stringsAsFactors = FALSE)
#'
#' # remove row count
#' dat <- nfpkm[-1:-3][colSums(nfpkm[-1:-3]) != 0]
#'
#' # corheat with dynamic tree cut
#' res <- corheat(dat=dat, distm="spearman", clm="average", method_dycut="tree", draw=TRUE)
#' @export
corheat <- function(dat, distm, clm, method_dycut, draw){
# minerva::mine(dat, alpha=1.0)
  # argument check: dat ----
  if(class(dat)=="data.frame"|class(dat)=="matrix"){
    if(ncol(dat) > 1000){
      cat(" In case of 'ncol(dat)' over 1000, no drawing heat map of all genes, and just return hclust object.")
      smp_genes <- sample(1:ncol(dat), 1000)
    }
  } else {
    stop("'dat' class is data frame or matrix")
  }

  # argument check: dism ----
  distms <- c("pearson","abspearson", "squarepearson", "spearman")
  if(!any(distms %in% distm)){
    stop(paste('Select a distance measure form', paste(distms, collapse = ", ")))
  }

  # argument check: clm ----
  clms <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")
  if(!any(clms %in% clm)){
    stop(paste('Select a clustering meathod form', paste(clms, collapse = ", ")))
  }

  # clustering ----
  dat <- as.matrix(dat)
  if (distm=="squarepearson"){
    cor.dat <- stats::cor(dat)^2
    dis.mat <- as.dist(1-cor.dat)
    r_hcl <- hclust(dis.mat, method = clm)
    r_den <- dendextend::hang.dendrogram(as.dendrogram(r_hcl))

  } else if (distm=="pearson" ){
    Distm <- "correlation"
    cor.dat <- stats::cor(dat, method = distm)
    dis.mat <- stats::as.dist(1-cor.dat)
    r_hcl <- hclust(amap::Dist(t(as.matrix(dat)), method = Distm), method = clm)
    r_den <- dendextend::hang.dendrogram(as.dendrogram(r_hcl))

  } else if (distm=="abspearson"){
    cor.dat <- abs(stats::cor(dat))
    dis.mat <- stats::as.dist(1-cor.dat)
    r_hcl <- hclust(amap::Dist(t(dat), method = distm), method = clm)
    r_den <- dendextend::hang.dendrogram(as.dendrogram(r_hcl))

  } else if (distm=="spearman"){
    cor.dat <- stats::cor(dat, method = distm)
    dis.mat <- stats::as.dist(1-cor.dat)
    r_hcl <- hclust(amap::Dist(t(dat), method = distm), method = clm)
    r_den <- dendextend::hang.dendrogram(as.dendrogram(r_hcl))

  } else {
    stop('Select a distance measure form c(\"squarepearson\", \"pearson\", \"abspearson\") ')
  }


  # dynamic cut  method "tree", "hybrid" ----
  dycmethods <- c("tree", "hybrid", "height")
  if(!any(dycmethods %in% method_dycut)){
    stop("'method_dycut' select from 'tree','hybrid', and 'height'.")
  }
  # method tree  method hybrid
  if (method_dycut=="tree"){
    dyct <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, method = method_dycut )
  } else if (method_dycut=="hybrid"){
    dyct <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, distM = as.matrix(dis.mat), method = method_dycut )
  } else {
    stop("cut height still be deveropping")
  }

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  ncl <- length(unique(dyct))
  vcol <- gg_color_hue(ncl)
  cols.ccl <- vcol[as.factor(dyct)]

  # leaf label ----
  if(any(length(labels(r_hcl)))){
    lab <- r_hcl$labels
  }else{
    lab <- as.character(1:ncol(dat))
  }

  names(dyct) <- lab

  # create leaf colour with dynamiccut cluster
  leaf_col <- labels(r_den)
  llab <- split(lab, factor(dyct))
  invisible(lapply(seq_along(llab),
                   function(i){leaf_col[leaf_col %in% llab[[i]]] <<- vcol[i]}))

  # add edge colour ----
  f <- function(x, y){
    r_den <<- dendextend::set(r_den, "by_labels_branches_col", value = x, TF_values = y)
  }
  invisible(mapply(f, llab, vcol))


  # set branch width, leaf labels color and labels cex
  r_den <- r_den %>%
    dendextend::set("branches_lwd", value = 0.5) %>%
    dendextend::set("labels_col", value = leaf_col) %>%
    dendextend::set("labels_cex", 0.5)

  # method hybrid
  # dyct <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, distM = as.matrix(dis.mat), method = method_dycut )

  # draw heat map -----
  if(draw==TRUE){

    gplots::heatmap.2(
      as.matrix(cor.dat),
      symm = T,
      #revC = T,
      col=rev(heat.colors(256)), # bluered(256), rev(heat.colors(256))
      scale="none", #"row", "Column", "both", "none"
      dendrogram = "both", #"none",#"col", #"both",
      Colv=(r_den),
      Rowv=r_den,
      ColSideColors = cols.ccl,
      RowSideColors = cols.ccl,
      key=TRUE,
      keysize=1,
      symkey=FALSE,
      density.info="none",
      trace="none",
      margin=c(6,5),
      cexCol=0.8,
      labRow = FALSE,
      labCol = FALSE
    )
    return(list(cormat=cor.dat, r_hcl=r_hcl, cl_with_dynamiccut = dyct))
  } else {
    return(list(cormat=cor.dat, r_hcl=r_hcl, cl_with_dynamiccut = dyct))
  }
}
