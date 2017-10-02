#' Convert a correlation matrix to a weighted edge list or igraph object
#' @description Convert a correlation matrix to a weighted edge list from lower triangled data.
#'     if this matrix is named matrix, these names as node id.
#' @usage matoedge(mat, diag, zero.weight, format)
#' @param mat A squared and synmetric matrix, like a correlation matirx.
#' @param format select "igraph"(default) or "df" .
#' @param diag default, FALSE
#' @param zero.weight default, FALSE
#' @examples
#' mat <- cor(iris[-5])
#' matoedge(mat = mat)
#' matoedge(mat = mat, diag = TRUE, format = "df")
#' @export
matoedge <- function(mat, diag = FALSE, zero.weight = FALSE, format="igraph"){
  # argument check: mat is a squared matrix
  if(nrow(mat) != ncol(mat)){
    stop("This matrix is not squared matrix")
  }

  if (format == "igraph"){
    if (diag == FALSE){
      igraph::graph.adjacency(mat, mode="undirected", weighted=TRUE, diag = FALSE)
    } else {
      igraph::graph.adjacency(mat, mode="undirected", weighted=TRUE)
    }

  } else if (format =="df"){
    if (diag == FALSE){
      bit <- cbind(rep(1:(nrow(mat)-1), (nrow(mat)-1):1),
                   unlist(lapply(seq(nrow(mat))[-1], function(i)i:(nrow(mat)))) )

      if(identical(rownames(mat), NULL)){
        x_id=bit[,1]; y_id=bit[,2]
      } else {
        x_id=rownames(mat)[bit[,1]]; y_id=rownames(mat)[bit[,2]]
      }
      d <- data.frame(x_id, y_id, value = mat[lower.tri(mat)])
      if(zero.weight==FALSE){ d[d$value != 0,] } else if (zero.weight==TRUE){d}

    }else{
      bit <- cbind(rep(1:nrow(mat), nrow(mat):1),
                   unlist(lapply(1:ncol(mat), function(i)seq(i, ncol(mat)))))

      if(identical(rownames(mat), NULL)){
        x_id=bit[,1]; y_id=bit[,2]
      } else {
        x_id=rownames(mat)[bit[,1]]; y_id=rownames(mat)[bit[,2]]
      }
      d <- data.frame(x_id, y_id, value = mat[lower.tri(mat, diag = T)])
      if(zero.weight==FALSE){ d[d$value != 0,] } else if (zero.weight==TRUE){d}

    }
  }
}
