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
#' # sampling 200 genes
#' dat <- nfpkm[,sample(5:ncol(nfpkm), 200)]
#'
#' # corheat with dynamic tree cut
#' res <- corheat(dat=dat, distm="spearman", clm="average", method_dycut="tree", draw=TRUE)
#'
#' @export
corheat <- function(dat, distm, clm, method_dycut, draw){
  # distm="spearman"; clm="average"; method_dycut="tree"; draw=TRUE

  # argument check: dat ----
  if(class(dat)=="data.frame"|class(dat)=="matrix"){
    if(ncol(dat) > 1000){
      cat(" In case of 'ncol(dat)' over 1000, no drawing heat map of all genes, and just return hclust object.")
      dat <- dat[,sample(1:ncol(dat), 1000)]
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

  # dynamicTree cutmethod "tree" or "hybrid" ----
  ## leaf label(tree order) ----
  if(any(length(labels(r_hcl)))){
    lab <- r_hcl$labels[r_hcl$order]
  }else{
    lab <- as.character(1:ncol(dat))
  }
  ## cutreeDynamic
  if (method_dycut=="tree"){
    dyct <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, method = method_dycut )
    dyct <- dyct[r_hcl$order]
    names(dyct) <- lab

  } else if (method_dycut=="hybrid"){
    dyct <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, distM = as.matrix(dis.mat), method = method_dycut )
    dyct <- dyct[r_hcl$order]
    names(dyct) <- lab
  } else {
    stop("cut height still be deveropping")
  }

  ## create colour(tree order) -----
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  vcol <- gg_color_hue(length(unique(dyct))) # colours corresponding to clusters
  cols.ccl <- vcol[as.factor(dyct)] # vector of leaf colours
  llab <- split(lab, factor(dyct))


  # modify dendrogram object with 'dendextend' ----
  ## add edge colour concatenated to leaf belongs to same cluster ----
  f <- function(x, y){
    r_den <<- dendextend::set(r_den, "by_labels_branches_col", value = x, TF_values = y)
  }
  invisible(mapply(f, llab, vcol))


  ## set parameter of these branch width, leaf labels color and labels cex to dendrogram object. ----
  r_den <- r_den %>%
    dendextend::set("branches_lwd", value = 0.5) %>%
    dendextend::set("labels_col", value = cols.ccl) %>%
    dendextend::set("labels_cex", 0.5)


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
