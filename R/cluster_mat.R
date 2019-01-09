#' matplot per cluster
#' @description split count data of transcriptome with result of gene clustering, such as "cutree" or "cutreeDynamic".
#' @usage cluster_mat(dat, res_dycut, fcdat, fctrl_col, facet_ncol, facet_nrow)
#' @param dat a data frame of count data for transcriptome, which columns are genes and row are samples.
#' @param res_dycut result of "dynamicTreeCut::cutreeDynamic", or "cutree"
#' @param fcdat data.frame of one or two factor for samples, which first column is a factor for Y-axis and a second column is a factor for X-axis.
#' If samples factor is not unique and replicate, median of replicate was expression value.
#' @param fctrl_col colour code for colour factor. The default value is NULL, then default ggplot colours.
#' @param facet_ncol ncol argument of facet_wrap. The default value is 3.
#' @param facet_nrow nrow argument of facet_wrap. The default value is 4
#' @examples#
#' # # argument1: data containing replicate or not.
#' # nfpkm <- rskodat::nfpkm
#' # nfpkm_norep <- nfpkm %>% group_by_at(1:4) %>% summarise_at(vars(-1:-4), funs(median))
#' #
#' # # argument2: result of dycutdf
#' # res_dyct1 <- cornet::dycutdf(nfpkm, distm = "spearman", column = -1:-4)
#' # res_dyct2 <- cornet::dycutdf(nfpkm_norep, distm = "spearman", column = -1:-4)
#' #
#' # cl1 <- res_dyct1$dynamic_cut
#' # cl2 <- res_dyct2$dynamic_cut
#' #
#' # # argument3: factorial data.frame one(x) or two(y,x)
#' # fcdat1 <- data.frame(runs=nfpkm$runs, days=nfpkm$days)
#' # fcdat2 <- data.frame(runs=nfpkm_norep$runs, days=nfpkm_norep$days)
#' #
#' # # execution
#' # res1 <- cornet::cluster_mat(dat = nfpkm[-1:-4], res_dycut = cl1, fcdat = fcdat1)
#' # res2 <- cornet::cluster_mat(dat = nfpkm_norep[-1:-4], res_dycut = cl2, fcdat = fcdat2,
#' #                            fctrl_col = c(1,2,4))
#' #
#' # # drawing selected panels
#' # do.call(gridExtra::grid.arrange, c(res2[c(2,5)], list(ncol=2)))
#'
#' @importFrom dplyr %>% group_by_at mutate arrange slice ungroup
#' @importFrom ggplot2 geom_line aes facet_wrap theme_bw labs aes_string
#' @importFrom plyr .
#' @importFrom tidyr gather
#' @export
cluster_mat <- function(dat, res_dycut, fcdat, fctrl_col=NULL, facet_ncol=3, facet_nrow=4){
  # argument check: dat and fcdat ----
  if (nrow(dat) != nrow(fcdat)){
    stop("nrow dat and nrow fcdat must to be same.")
  }
  # argument check: length of res_dycut and ncol of dat ----
  if (!identical(names(dat), names(res_dycut))){
    stop("res_dycut is result of gene clustering of dat")
  }

  # data.frame of samples factor ----
  labs <- paste0("cl.",names(table(res_dycut)), "(",table(res_dycut), ")")
  cllabs <- data.frame(cluster = as.numeric(names(table(res_dycut))),
                       lab = factor(labs, levels = labs))

  # if samples factor is not unique and replicate, median of replicate was expression value ----
  value = NULL
  tmp <- dat %>%
      apply(., 2, scale) %>%
      data.frame(fcdat, ., stringsAsFactors=F, check.names=F) %>%
      tidyr::gather(., key="gene", value="value", -1:-ncol(fcdat), factor_key=T) %>%
      dplyr::group_by_at(1:(ncol(fcdat)+1)) %>%
      dplyr::mutate(value=median(value)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      merge(., data.frame(gene=names(res_dycut), res_dycut), by="gene", all=T) %>%
      merge(., cllabs, by.x="res_dycut", by.y="cluster") %>%
      dplyr::arrange(.[[2]], .[[3]], .[[4]])


  # multiple pages arguments ----
  panels <- levels(tmp$lab)
  n_panels <- length(panels)
  n_p <- facet_ncol*facet_nrow
  n_pages <- ceiling(n_panels/n_p)
  layout_pages <- factor(findInterval(seq_along(levels(tmp$lab)), n_p*seq(n_pages)[1:(n_pages-1)]+1)+1)
  panels_per_page <- split(panels, layout_pages)

  # multiple pages ----
  for(i in seq_along(panels_per_page)){
    ggdi <- tmp %>% filter(tmp$lab %in% panels_per_page[[i]])
    if (ncol(fcdat) == 2){
      if (is.null(fctrl_col)){
        ggi <- ggplot2::ggplot(ggdi,ggplot2::aes_string(x=names(ggdi)[4], y=names(ggdi)[5], colour=names(ggdi)[3],
                                                        group=interaction(ggdi[,2], ggdi[,3]))) +
          ggplot2::geom_line(size=0.1, alpha=0.5) +
          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size = 5))) +
          ggplot2::theme_bw() +
          ggplot2::facet_wrap(~lab, ncol=facet_ncol, nrow=facet_nrow)
        plot(ggi)
      } else {
        ggi <- ggplot2::ggplot(ggdi,ggplot2::aes_string(x=names(ggdi)[4], y=names(ggdi)[5], colour=names(ggdi)[3],
                                                        group=interaction(ggdi[,2], ggdi[,3]))) +
          ggplot2::geom_line(size=0.1, alpha=0.5) +
          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size = 5))) +
          ggplot2::theme_bw() +
          ggplot2::scale_colour_manual(values = fctrl_col) +
          ggplot2::facet_wrap(~lab, ncol=facet_ncol, nrow=facet_nrow)
        plot(ggi)

      }

    }else{
      ggi <- ggplot2::ggplot(tmp,ggplot2::aes_string(x=names(tmp)[4], y=names(tmp)[5],colour=names(tmp)[3],
                                                     group=names(tmp)[2])) +
        ggplot2::geom_line(size=0.1, alpha=0.5) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size = 5))) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~lab, ncol=4)
      plot(ggi)
    }
  }

  # all ggplot objects into 'gg_list'
  gg_list <- vector("list", length = n_panels)
  gg_dats <- split(tmp, f = tmp$lab)
  if(is.null(fctrl_col)){
    gg_list <- lapply(gg_dats, function(x){
      ggplot2::ggplot(x,ggplot2::aes_string(x=names(x)[4], y=names(x)[5],
                                            colour=names(x)[3], group=interaction(x[,2], x[,3]))) +
        ggplot2::geom_line(size=0.1, alpha=0.5) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size = 5))) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~lab, ncol=4)
    })

  }else{
    gg_list <- lapply(gg_dats, function(x){
      ggplot2::ggplot(x,ggplot2::aes_string(x=names(x)[4], y=names(x)[5],
                                            colour=names(x)[3], group=interaction(x[,2], x[,3]))) +
        ggplot2::geom_line(size=0.1, alpha=0.5) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size = 5))) +
        ggplot2::scale_colour_manual(values = fctrl_col) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~lab, ncol=4)
    })

  }

  return(gg_list)
}

