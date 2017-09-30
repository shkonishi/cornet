#' Convert a correlation matrix to a weighted edge list
#' @description Convert a correlation matrix to a weighted edge list from lower triangled data.
#'     if this matrix is named matrix, these names as node id.
#' @usage matoedge(mat, diag)
#' @param mat A squared matrix, and lower and upper triangled data to be a same.
#' @param diag default, FALSE
#' @examples
#' mat <- cor(iris[-5])
#' matoedge(mat = mat)
#' matoedge(mat = mat, diag = TRUE)
#' @export
matoedge <- function(mat, diag = FALSE){
  # argument check: mat is a squared matrix
  if(nrow(mat) != ncol(mat)){
    stop("This matrix is not squared matrix")
  }

  if(diag == FALSE){
    bit <- cbind(rep(1:(nrow(mat)-1), (nrow(mat)-1):1),
                 unlist(lapply(seq(nrow(mat))[-1], function(i)i:(nrow(mat)))) )

    if(identical(rownames(mat), NULL)){
      x_id=bit[,1]; y_id=bit[,2]
    } else {
      x_id=rownames(mat)[bit[,1]]; y_id=rownames(mat)[bit[,2]]
    }
    d <- data.frame(x_id, y_id, value = mat[lower.tri(mat)])
    d[d$value != 0,]

  }else{
    bit <- cbind(rep(1:nrow(mat), nrow(mat):1),
                 unlist(lapply(1:ncol(mat), function(i)seq(i, ncol(mat)))))

    if(identical(rownames(mat), NULL)){
      x_id=bit[,1]; y_id=bit[,2]
    } else {
      x_id=rownames(mat)[bit[,1]]; y_id=rownames(mat)[bit[,2]]
    }
    d <- data.frame(x_id, y_id, value = mat[lower.tri(mat, diag = T)])
    d[d$value != 0,]

  }


}
