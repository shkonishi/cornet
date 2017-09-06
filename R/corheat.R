# dat  cor matrix
# pearson, abspearson, squarepearson
# dat <- temp[-1]; d.method <- "pearson"; c.method <- "average"; corr.col="blured"
# n_sample <- 500
corheat <- function(dat, d.method, c.method, corr.col, n_sample){
  dat <- as.matrix(dat)
  if (d.method=="squarepearson"){
    cor.dat <- cor(t(dat))^2
    cl <- hclust(as.dist(1-(cor(t(dat))^2)), method = c.method)
    rden <- as.dendrogram(cl)

  } else if (d.method=="pearson" ){
    cor.dat <- cor(t(dat))
    cl <- hclust(amap::Dist(as.matrix(dat), method = d.method), method = c.method)
    rden <- as.dendrogram(cl)

  } else if (d.method=="abspearson"){
    cor.dat <- abs(cor(t(dat)))
    cl <- hclust(amap::Dist(dat, method = d.method), method = c.method)
    rden <- as.dendrogram(cl)
  } else if (d.method=="spearman"){
    cor.dat <- cor(t(dat), method = d.method)
    cl <- hclust(amap::Dist(dat, method = d.method), method = c.method)
    rden <- as.dendrogram(cl)

  } else {
    stop('clustering method is different. select one form c(\"squarepearson\", \"pearson\", \"abspearson\") ')
  }

  if (corr.col=="bluered"){
    corr_col <- bluered(256)
  } else if(corr.col=="heat.colors"){
    corr_col <- rev(heat.colors(256))
  } else {
    stop('Select one from c(\"bluered\", \"head.colors\")')
  }

  #pdf(pdfn)
  gplots::heatmap.2(
    as.matrix(cor.dat),
    symm = T,
    revC = T,
    col=corr_col, # bluered(256), rev(heat.colors(256))
    scale="none", #"row", "Column"
    dendrogram = "both", #"none",#"col", #"both",
    Colv=rden,
    Rowv=rden,
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
  #dev.off()
  return(list(cormat =cor.dat, res_hcl=cl))
}
