#' Cluster_mat
#' @description Cluster_mat
#' @usage cluster_mat(dat, distm, clm, column, method_dycut, y_fctr, x_fctr, rep_fctr, ...)
#' @param dat data frame or matrix
#' @param distm distance measure from amap::Dist
#' @param clm hclust methods
#' @param column column which containing expression data
#' @param method_dycut method of dynamic cut
#' @param x_fctr x-axis factor for ggplot object like matplot
#' @param y_fctr optional
#' @param rep_fctr optional
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
cluster_mat <- function(dat, distm, clm, column, method_dycut, y_fctr, x_fctr, rep_fctr, ...){

  if (class(dat)=="data.frame"){
    mat <- dat[column]
  } else if (class(dat) == "matrix"){
    mat <- data.frame(dat)
    if(identical(dimnames(dat)[[2]], NULL)){
      names(mat) <- 1:ncol(mat)
    }
  }
  dis.mat <- amap::Dist(t(mat), method = distm)
  r_hcl <- stats::hclust(dis.mat, method = clm)
  r_den <- dendextend::hang.dendrogram(as.dendrogram(r_hcl))

  # method tree  method hybrid ----
  if (method_dycut=="tree"){
    dyct <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, method = method_dycut, ... )
  } else if (method_dycut=="hybrid"){
    dyct <- dynamicTreeCut::cutreeDynamic(dendro = r_hcl, distM = as.matrix(dis.mat), method = method_dycut, ... )
  } else {
    stop("cut height still be deveropping")
  }

  # leaf label (leaf order) ----
  lab <- r_hcl$labels
  names(dyct) <- lab
  dyct.lo <- dyct[r_hcl$order]

  # draw dendrogram with cluster -----
  ## cluster colour, cluster number, leaf colour
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  vcol <- gg_color_hue(length(unique(dyct)))

  ## cluster level
  cl_levels <- unique(dyct.lo)

  ## leaf colour(sample order for side bar colour)  &
  leaf_col_so <- vcol[factor(dyct, levels=cl_levels)] # leaf col (sample order)
  leaf_col_lo <- leaf_col_so[r_hcl$order] # leaf col(leaf order)

  ## cluster id
  llab <- split(lab[r_hcl$order], factor(dyct.lo, levels=cl_levels))

  ## split data frame corresponding to cluster.
  cl_dat <- lapply(llab, function(x)mat[x])

  ## add edge colour ----
  f <- function(x, y){
    r_den <<- dendextend::set(r_den, "by_labels_branches_col", value = x, TF_values = y)
  }
  invisible(mapply(f, llab, vcol))

  # barplot cluster side bar ----
  ## graph layout ----
  gmat <- matrix(c(1,2), nrow=2, byrow = TRUE)
  lay <- graphics::layout(gmat, widths=c(10,10), heights=c(4,1), respect = T)
  graphics::par(mar=c(0,2,0,0))
  graphics::plot(r_den, leaflab="none")

  ## cluster side bar ----
  par(mar=c(0,2,2,0))
  cl_txt <- as.character(cl_levels) # number of cluster elements
  cl_num <- table(dyct.lo)[cl_txt] # cluster number table
  cl_txt_al <- rep(NA, length(dyct)) # text on color side bar
  cl_txt_al[ceiling(cumsum(cl_num)-cl_num/2)] <- cl_txt
  bp <- graphics::barplot(rep(1, length(leaf_col_lo)), yaxt="n", border = leaf_col_lo, col=leaf_col_lo)
  graphics::text(bp, y = 0.5, labels = cl_txt_al, col="white")

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
      ggplot2::geom_line(alpha=0.3) +
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
      ggplot2::geom_line(alpha=0.2) +
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
      ggplot2::geom_line(alpha=0.2) +
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
    if(length(cl_num)%%2 == 0){
      fct_col <- 2
    }else{
        fct_col <- 3
      }
    ggmat <- ggplot2::ggplot(gg_dat,
                ggplot2::aes(x=x_fctr, y=value, colour=y_fctr, group=interaction(genes, y_fctr))) +
      ggplot2::geom_line(alpha=0.2) +
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
