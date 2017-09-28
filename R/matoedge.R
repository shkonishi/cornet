#' Convert a correlation matrix to a weighted edge list
#' @description Convert a correlation matrix to a weighted edge list.
#'     if this matrix is named matrix, these names as node id. Only lower triangled data is used.
#' @usage matoedge(mat)
#' @param mat A squared matrix, and lower and upper triangled data to be a same.
#' @examples
#' mat <- cor(iris[-5])
#' matoedge(mat)
#' @export
matoedge <- function(mat){
  # argument check: mat is a squared matrix
  if(nrow(mat) != ncol(mat)){
    stop("This matrix is not squared matrix")
  }

  bit <- cbind(unlist(lapply(seq(nrow(mat)-1), function(i)rep(i,(nrow(mat)-i)))),
               unlist(lapply(seq(nrow(mat))[-1], function(i)i:(nrow(mat))))
  )
  if(identical(rownames(mat), NULL)){
    x_id=bit[,1]; y_id=bit[,2]
  } else {
    x_id=rownames(mat)[bit[,1]]; y_id=rownames(mat)[bit[,2]]
  }

  data.frame(x_id, y_id, value = mat[lower.tri(mat)])
}
