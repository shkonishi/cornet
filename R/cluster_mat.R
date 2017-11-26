#' matplot per cluster
#' @description split count data of transcriptome with result of gene clustering, such as "cutree" or "cutreeDynamic".
#' @usage cluster_mat(dat, res_dycut, fcdat)
#' @param dat a data frame of count data for transcriptome, which columns are genes and row are samples.
#' @param res_dycut result of "dynamicTreeCut::cutreeDynamic", or "cutree"
#' @param fcdat data.frame of one or two factor for samples, which first column is a factor for Y-axis and a second column is a factor for X-axis
#' @examples
#' # # argument1: data
#' # nfpkm <- rskodat::nfpkm
#'
#' # # argument2: result of dycutdf
#' # res_dycut <- dycutdf(nfpkm, distm = "abscorrelation", column = -1:-4)
#'
#' # # argument3: factorial data.frame one(x) or two(y,x)
#' # fcdat <- data.frame(runs=nfpkm$runs, days=nfpkm$days)
#' # cluster_mat(nfpkm[-1:-4], res_dycut$dynamic_cut, fcdat)
#'
#' @importFrom dplyr %>% group_by_at mutate arrange slice ungroup
#' @importFrom ggplot2 geom_line aes facet_wrap theme_bw labs aes_string
#' @importFrom plyr .
#' @importFrom tidyr gather
#' @export

cluster_mat <- function(dat, res_dycut, fcdat){
  # argument check: dat and fcdat ----
  if (nrow(dat) != nrow(fcdat)){
    stop("nrow dat and nrow fcdat must to be same.")
  }
  # argument check: length of res_dycut and ncol of dat ----
  if (!identical(names(dat), names(res_dycut))){
    stop("res_dycut is result of gene clustering of dat")
  }

  # data.frame of samples factor ----
  cllabs <- data.frame(cluster = as.numeric(names(table(res_dycut))),
                       lab = factor(paste0(names(table(res_dycut)), "(",table(res_dycut), ")"),
                                    levels=unique(paste0(names(table(res_dycut)), "(",table(res_dycut), ")"))),
                       stringsAsFactors=F)

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

  if (ncol(fcdat) == 2){
    ggplot2::ggplot(tmp,ggplot2::aes_string(x=names(tmp)[4], y=names(tmp)[5], colour=names(tmp)[3],
                                     group=interaction(tmp[,2], tmp[,3]))) +
      ggplot2::geom_line(size=0.1, alpha=0.5) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size = 5))) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~lab, ncol=4)
  } else {
    ggplot2::ggplot(tmp,ggplot2::aes_string(x=names(tmp)[4], y=names(tmp)[5],
                                     colour=names(tmp)[3], group=names(tmp)[2])) +
      ggplot2::geom_line(size=0.1, alpha=0.5) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size = 5))) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~lab, ncol=4)
  }
}

