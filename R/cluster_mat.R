#' Cluster_mat
#' @description Cluster_mat
#' @usage cluster_mat(dat, distm, clm, column, method_dycut, x_fctr, y_fctr, rep_fctr, ...)
#' @param dat data frame or matrix
#' @param distm distance measure from amap::Dist options. default is "pearson"
#' @param clm clustering method select from hclust options default is "average"
#' @param column column which containing expression data default 'column=1:ncol(dat)'
#' @param method_dycut method of dynamic cut
#' @param x_fctr x-axis factor for ggplot object like matplot
#' @param y_fctr optional default is NA
#' @param rep_fctr optional default is NA
#' @param ... additional dynamicTreeCut::cutreeDynamic options
#' @examples
#' d <- data.frame(t(iris[-5]))
#' xfctr <- factor(names(iris[-5]), levels = names(iris[-5]))
#' res <- cluster_mat(dat = d, distm = "euclidean", clm = "average",
#'        column = 1:150, method_dycut = "hybrid", x_fctr = xfctr)
#' res[[1]] # hclust object
#' res[[2]] # result of dynamic cut object
#' res[[3]] # ggplot object of matplot of all clusters respectively.
#'
#' ##  The following examples takes too much times
#' # sample data
#' # fp <- system.file("extdata/nfpkm_rnsq.txt", package = "cornet")
#' # nfpkm <- read.table(fp, header=TRUE, stringsAsFactors = FALSE)
#'
#' ## dat1:
#' # dat1 <- nfpkm
#' # res1 <- cluster_mat(dat = dat1, distm = "spearman", clm = "average",
#' #     column = 5:ncol(dat1), method_dycut = "tree",
#' #     y_fctr = dat1$runs, x_fctr = dat1$days, rep_fctr = dat1$reps)
#'
#' # dat2: no replicate and y_fctr
#' # dat2 <- nfpkm[which(nfpkm$runs==1)[seq(1,36, 3)] ,]
#' # dat2 <- data.frame(dat2[2], dat2[-1:-3][colSums(dat2[-1:-3]) != 0])
#'
#' # res2 <- cluster_mat(dat = dat2, distm = "spearman", clm = "average",
#' #     column = 2:ncol(dat2), method_dycut = "tree", x_fctr=dat2$days)
#'
#' ## dat3: no y_fctr
#' # dat3 <- nfpkm[nfpkm$runs == 1, -1]
#' # dat3 <- data.frame(dat3[1:2], dat3[-1:-2][colSums(dat3[-1:-2]) != 0])
#' # res3 <- cluster_mat(dat = dat3, distm = "spearman", clm = "average",
#' #     column = 3:ncol(dat3), method_dycut = "tree", x_fctr=dat3$days, rep_fctr=dat3$reps)
#' @import dplyr
#' @import ggplot2
#' @importFrom amap Dist
#' @importFrom tidyr gather
#' @importFrom graphics barplot layout layout.show par plot text
#' @importFrom stats sd hclust as.dendrogram median
#' @export
#'
cluster_mat <- function(dat, distm="pearson", clm="average", column=1:ncol(dat), method_dycut="tree", x_fctr, y_fctr=NA, rep_fctr=NA, ...){
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
  r_den <- dendextend::hang.dendrogram(as.dendrogram(r_hcl))

  # cutreeDynamic  ----
  dyct.tr <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, method = "tree", ... )
  dyct.hb <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, distM = as.matrix(dis.mat), method = "hybrid", ... )

  # leaf label (leaf order) ----
  lab <- r_hcl$labels
  names(dyct.hb) <- lab; names(dyct.tr) <- lab
  dyct.hb.lo <- dyct.hb[r_hcl$order];  dyct.tr.lo <- dyct.tr[r_hcl$order]

  # draw dendrogram with cluster
  ## cluster colour, cluster number ----
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  if(any(unique(dyct.hb) == 0)){
    hb.vcol =c("grey50",gg_color_hue(length(unique(dyct.hb)) -1))
    hb.cl_levels <- c(0, unique(dyct.hb.lo)[!unique(dyct.hb.lo) %in% 0])
  }else{
    hb.vcol <- gg_color_hue(length(unique(dyct.hb)))
    hb.cl_levels <- unique(dyct.hb.lo)
  }
  if(any(unique(dyct.tr) == 0)){
    tr.vcol =c("grey50",gg_color_hue(length(unique(dyct.tr)) -1))
    tr.cl_levels <- c(0, unique(dyct.tr.lo)[!unique(dyct.tr.lo) %in% 0])
  }else{
    tr.vcol <- gg_color_hue(length(unique(dyct.tr)))
    tr.cl_levels <- unique(dyct.tr.lo)
  }

  ## leaf colour(sample order for side bar colour) ----
  tr.leaf_col_so <- tr.vcol[factor(dyct.tr, levels=tr.cl_levels)] # leaf col (sample order)
  tr.leaf_col_lo <- tr.leaf_col_so[r_hcl$order] # leaf col(leaf order)
  hb.leaf_col_so <- hb.vcol[factor(dyct.hb, levels=hb.cl_levels)] # leaf col (sample order)
  hb.leaf_col_lo <- hb.leaf_col_so[r_hcl$order] # leaf col(leaf order)

  ## split cluster and data frame corresponding to cluster, and  add edge colour to dendrogram ----
  f <- function(x, y){
    r_den <<- dendextend::set(r_den, "by_labels_branches_col", value = x, TF_values = y)
  }
  if (method_dycut == "tree"){
    llab <- split(lab[r_hcl$order], factor(dyct.tr.lo, levels=tr.cl_levels))
    cl_levels <- tr.cl_levels
    invisible(mapply(f, llab, tr.vcol))
    dyct <- dyct.tr
  } else if (method_dycut == "hybrid"){
    llab <- split(lab[r_hcl$order], factor(dyct.hb.lo, levels=hb.cl_levels))
    cl_levels <- hb.cl_levels
    invisible(mapply(f, llab, hb.vcol))
    dyct <- dyct.hb
  }

  cl_dat <- lapply(llab, function(x)mat[x])


  # barplot as cluster side bar ----
  ## graph layout ----
  def.par <- par(no.readonly = TRUE) # default
  gmat <- matrix(c(1,2,3), nrow=3, byrow = TRUE)
  lay <- graphics::layout(gmat, widths=c(10,10,10), heights=c(4,1,1), respect = T)
  graphics::par(mar=c(0,2,0,0))
  graphics::plot(r_den, leaflab="none")

  ## cluster side bar ----
  par(mar=c(0,2,2,0))
  ## number of cluster elements
  ### hybrid bar
  hb.cl_txt <- if(any(hb.cl_levels==0)){
    as.character(hb.cl_levels[hb.cl_levels!=0])
  }else{
    as.character(hb.cl_levels)
  }
  ### tree bar
  tr.cl_txt <- if(any(tr.cl_levels==0)){
    as.character(tr.cl_levels[tr.cl_levels!=0])
  }else{
    as.character(tr.cl_levels)
  }

  hb.cl_num <- table(dyct.hb.lo)[hb.cl_txt] # cluster number table
  hb.cl_txt_al <- rep(NA, length(dyct.hb.lo)) # text on color side bar
  hb.cl_txt_al[match(hb.cl_txt, dyct.hb.lo) + ceiling(hb.cl_num/2)] <- hb.cl_txt

  tr.cl_num <- table(dyct.tr.lo)[tr.cl_txt] # cluster number table
  tr.cl_txt_al <- rep(NA, length(dyct.tr.lo)) # text on color side bar
  tr.cl_txt_al[match(tr.cl_txt, dyct.tr.lo) + ceiling(tr.cl_num/2)] <- tr.cl_txt

  if (method_dycut=="tree"){
    bp <- graphics::barplot(rep(1, length(tr.leaf_col_lo)), yaxt="n", ylab="tree", border = tr.leaf_col_lo, col=tr.leaf_col_lo)
    #graphics::text(bp, y = 0.5, labels = tr.cl_txt_al, col="white")
    graphics::mtext(side = 2, outer = 0, text = "tree")

    bp <- graphics::barplot(rep(1, length(hb.leaf_col_lo)), yaxt="n", ylab="hybrid", border = hb.leaf_col_lo, col=hb.leaf_col_lo)
    #graphics::text(bp, y = 0.5, labels = hb.cl_txt_al, col="white")
    graphics::mtext(side = 2, outer = 0, text = "hybrid")

  }else if (method_dycut=="hybrid"){
    bp <- graphics::barplot(rep(1, length(hb.leaf_col_lo)), yaxt="n", border = hb.leaf_col_lo, col=hb.leaf_col_lo)
    #graphics::text(bp, y = 0.5, labels = hb.cl_txt_al, col="white")
    graphics::mtext(side = 2, outer = 0, text = "hybrid")
    bp <- graphics::barplot(rep(1, length(tr.leaf_col_lo)), yaxt="n", ylab="tree", border = tr.leaf_col_lo, col=tr.leaf_col_lo)
    #graphics::text(bp, y = 0.5, labels = tr.cl_txt_al, col="white")
    graphics::mtext(side = 2, outer = 0, text = "tree")

  }
  par(def.par) # reset to default






  # matplot ----
  ## z-conversion ----
  zconverter <- function(dat, colums, mrgn){
    dat[,colums] <- sweep(sweep(dat[,colums], mrgn, apply(dat[,colums], mrgn, mean), "-"), mrgn,
                          apply(dat[,colums], mrgn, stats::sd), "/")
    return(dat)
  }
  zdat <- zconverter(dat = dat, colums = column, mrgn = 2)

  ## insert cluster number to gg_dat, and tidy data of ggplots. ----
  genes <- NULL; value <- NULL; cl <- NULL;

  if (missing(y_fctr) & !missing(rep_fctr)){ # x- and replicate factor, and no y_fctr
    if (!is.factor(x_fctr)){
      x_fctr <- factor(as.character(x_fctr), unique(as.character(x_fctr)))
    }
    if (!is.factor(rep_fctr)){
      rep_fctr <- factor(as.character(rep_fctr), unique(as.character(rep_fctr)))
    }

    ## add cluster number and mean of replicate. ----
    gg_dat <- data.frame(x_fctr, rep_fctr, zdat[column]) %>%
      tidyr::gather(key="genes", value="value", -1:-2, factor_key = TRUE) %>%
      mutate(cl=factor(dyct[factor(genes, levels=names(zdat[column]))], levels=cl_levels)) %>%
      group_by(rep_fctr, genes) %>%
      mutate(value = mean(value)) %>%
      dplyr::slice(1) %>%
      arrange(genes)

    ## median of all genes corresponding to same cluster. ----
    gg_dat_med <- gg_dat %>%
      ungroup() %>%
      group_by(x_fctr, cl) %>%
      mutate(value = median(value)) %>%
      dplyr::slice(1) %>%
      arrange(cl, x_fctr)

    ## plot all genes belogs to a cluster. ----
    ggmat <- ggplot2::ggplot(gg_dat, ggplot2::aes(x=x_fctr, y=value, group=genes)) +
      ggplot2::geom_line(alpha=0.3, linetype=3) +
      ggplot2::facet_wrap(~cl, ncol=3)

    ## plot median of all genes belogs to a cluster. ----
    ggmat2 <- ggplot2::ggplot(gg_dat_med, ggplot2::aes(x=x_fctr, y=value, group=genes)) +
      ggplot2::geom_line(size = 2, alpha = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~cl, ncol=3)

  } else if (missing(y_fctr) & missing(rep_fctr)){ # no rep_fctr and y_fctr
    if (!is.factor(x_fctr)){
      x_fctr <- factor(as.character(x_fctr), unique(as.character(x_fctr)))
    }

    ## add cluster number ----
    gg_dat <- data.frame(x_fctr, zdat[column]) %>%
      tidyr::gather(key="genes", value="value", -1, factor_key = TRUE) %>%
      mutate(cl=factor(dyct[factor(genes, levels=names(zdat[column]))], levels=cl_levels))

    ## median of all genes corresponding to same cluster. ----
    gg_dat_med <- gg_dat %>%
      ungroup() %>%
      group_by(x_fctr, cl) %>%
      mutate(value = median(value)) %>%
      dplyr::slice(1) %>%
      arrange(genes)

    ## plot all genes belogs to a cluster ----
    ggmat <- ggplot2::ggplot(gg_dat, ggplot2::aes(x=x_fctr, y=value, group=genes)) +
      ggplot2::geom_line(alpha=0.5, linetype=3) +
      ggplot2::facet_wrap(~cl, ncol=3)

    ## plot median of all genes belogs to a cluster. ----
    ggmat2 <- ggplot2::ggplot(gg_dat_med, ggplot2::aes(x=x_fctr, y=value,  group=genes)) +
      ggplot2::geom_line(size = 2, alpha = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~cl, ncol=3)


  }else if (!missing(y_fctr) & missing(rep_fctr)){ # no rep_fctr
    if (!is.factor(x_fctr)){
      x_fctr <- factor(as.character(x_fctr), unique(as.character(x_fctr)))
    }
    if (!is.factor(y_fctr)){
      y_fctr <- factor(as.character(y_fctr), unique(as.character(y_fctr)))
    }

    ## add cluster number ----
    gg_dat <- data.frame(y_fctr, x_fctr, zdat[column]) %>%
      tidyr::gather(key="genes", value="value", -1:-2, factor_key = TRUE) %>%
      mutate(cl=factor(dyct[factor(genes, levels=names(zdat[column]))], levels=cl_levels))

    ## median of all genes corresponding to same cluster. ----
    gg_dat_med <- gg_dat %>%
      ungroup() %>%
      group_by(x_fctr, cl) %>%
      mutate(value = median(value)) %>%
      dplyr::slice(1) %>%
      arrange(genes)

    ## plot all genes belogs to a cluster ----
    ggmat <- ggplot2::ggplot(gg_dat, ggplot2::aes(x=x_fctr, y=value, group=genes)) +
      ggplot2::geom_line(alpha=0.5, linetype=3) +
      ggplot2::facet_wrap(~cl, ncol=3)

    ## plot median of all genes belogs to a cluster. ----
    ggmat2 <- ggplot2::ggplot(gg_dat_med, ggplot2::aes(x=x_fctr, y=value,  group=genes)) +
      ggplot2::geom_line(size = 2, alpha = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~cl, ncol=3)


  } else if (!missing(y_fctr) & !missing(rep_fctr)){ # y_fctr and rep_fctr
    if (!is.factor(x_fctr)){
      x_fctr <- factor(as.character(x_fctr), unique(as.character(x_fctr)))
    }
    if (!is.factor(rep_fctr)){
      rep_fctr <- factor(as.character(rep_fctr), unique(as.character(rep_fctr)))
    }
    if (!is.factor(y_fctr)){
      y_fctr <- factor(as.character(y_fctr), unique(as.character(y_fctr)))
    }

    ## add cluster number and mean of replicate. ----
    gg_dat <- data.frame(y_fctr, x_fctr, rep_fctr, zdat[column]) %>%
      tidyr::gather(key="genes", value="value", -1:-3, factor_key = TRUE) %>%
      mutate(cl=factor(dyct[factor(genes, levels=names(zdat[column]))], levels=cl_levels)) %>%
      group_by(y_fctr, rep_fctr, genes) %>%
      mutate(value=mean(value)) %>%
      dplyr::slice(1) %>%
      arrange(genes)

    ## median of all genes corresponding to same cluster. ----
    gg_dat_med <- gg_dat %>%
      ungroup() %>%
      group_by(y_fctr, x_fctr, cl) %>%
      mutate(value = median(value)) %>%
      dplyr::slice(1) %>%
      arrange(genes)

    ## plot all genes belogs to a cluster.----
    if(length(cl_levels)%%2 == 0){
      fct_col <- 2
    }else{
      fct_col <- 3
    }
    ggmat <- ggplot2::ggplot(gg_dat,
                             ggplot2::aes(x=x_fctr, y=value, colour=y_fctr, group=interaction(genes, y_fctr))) +
      ggplot2::geom_line(alpha=0.5, linetype=3) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~cl, ncol=fct_col)

    ## plot median of all genes belogs to a cluster. ----
    ggmat2 <- ggplot2::ggplot(gg_dat_med,
                              ggplot2::aes(x=x_fctr, y=value, colour=y_fctr, group=interaction(genes, y_fctr))) +
      ggplot2::geom_line(size = 2, alpha = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~cl, ncol=fct_col)

  }else{
    stop("at least 'x_fctr' is needed, 'y_fctr' and 'rep_fctr' are optional")
  }

  # data shaping
  return(list(hclust_obj = r_hcl, dend_obj = r_den, dynamic_cut=dyct, cluster_dat = cl_dat,
              gg_mat_all = ggmat, gg_mat_med = ggmat2))

}
