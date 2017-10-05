#' node centrality
#' @description hub do.
#' @usage gethub(g, com_fun)
#' @param g igraph object
#' @param com_fun character: function name of igraph for community detection
#' @examples
#' ##
#' @importFrom igraph delete.vertices degree fastgreedy.community cluster_louvain delete.vertices V E vcount degree betweenness
#' @importFrom grDevices adjustcolor
#' @importFrom RColorBrewer brewer.pal
#' @export
gethub <- function(g, com_fun){
  com.algs <- c("fastgreedy.community", "cluster_louvain")
  if(any(com.algs %in% com_fun)){
    commu <- eval(parse(text=paste0(com_fun, "(g)")))
    mem <- commu$membership # module membership
    num_mod <- max(mem)  # number of max module
    num_nodes <- igraph::vcount(g) # numbero of vertex
    deg <- igraph::degree(g) # degree
  }else{
    stop('select from "fastgreedy.community" or "cluster_louvain"')
  }

  # Calculation of the within-module degree
  z_score <- numeric(num_nodes)
  for(s in 1:num_mod){
    v_seq <- subset(igraph::V(g),mem == s)
    g_sub <- igraph::delete.vertices(g,subset(igraph::V(g), mem != s))
    deg_sub <- igraph::degree(g_sub)
    z_score[v_seq] <- (deg_sub - mean(deg_sub)) / sd(deg_sub)
  }

  # Calculation of the participation coefficient
  participation_coeff <- numeric(num_nodes) + 1.0
  for(i in 1:num_nodes){
    mem_nei <- mem[igraph::neighbors(g,i)]
    for(s in 1:num_mod){
      deg_to_s <- length(subset(mem_nei,mem_nei == s))
      participation_coeff[[i]] <- participation_coeff[[i]] - (deg_to_s / deg[[i]]) ** 2.0
    }
  }

  # Classification
  role <- vector(mode = "numeric", length = num_nodes)

  pc <- as.data.frame(cbind(z_score,participation_coeff))
  names(pc) <- c("within_module_degree","participation_coefficient")

  # R1: Ultra-peripheral nodes
  v_seq <- which(pc$within_module_degree<2.5 & pc$participation_coeff<0.05)
  role[v_seq] <- "R1: Ultra-peripheral nodes"

  # R2: Peripheral nodes
  v_seq <- which(pc$within_module_degree<2.5 & pc$participation_coeff>=0.05 & pc$participation_coeff<0.625)
  role[v_seq] <- "R2: Peripheral nodes"

  # R3: Non-hub connectors
  v_seq <- which(pc$within_module_degree<2.5 & pc$participation_coeff>=0.625 & pc$participation_coeff<0.8)
  role[v_seq] <- "R3: Non-hub connectors"

  # R4: Non-hub kinless nodes
  v_seq <- which(pc$within_module_degree<2.5 & pc$participation_coeff>=0.8)
  role[v_seq] <- "R4: Non-hub kinless nodes"

  # R5: Provincial hubs
  v_seq <- which(pc$within_module_degree>=2.5 & pc$participation_coeff<0.3)
  role[v_seq] <- "R5: Provincial hubs"

  # R6: Connector hubs
  v_seq <- which(pc$within_module_degree>=2.5 & pc$participation_coeff>=0.3 & pc$participation_coeff<0.75)
  role[v_seq] <- "R6: Connector hubs"

  # R7: Kinless hubs
  v_seq <- which(pc$within_module_degree>=2.5 & pc$participation_coeff>=0.75)
  role[v_seq] <- "R7: Kinless hubs"

  # result dataframe
  #pc <- cbind(pc,role=role)
  result <- data.frame(node=igraph::V(g)$name, between=igraph::betweenness(g),  degree=deg, ndegree=deg/vcount(g), pc, role)
  sortresult <- result[order(result$between, decreasing = T),]

  # igplot with roll
  roles <- sapply(strsplit(as.character(result$role), ":"), "[", 1)
  cong <- delete.vertices(g, V(g)$name[degree(g)==0])
  roles <- sapply(strsplit(as.character(result$role), ":"), "[", 1)
  vrole <- c("R7","R6","R5","R4","R3")
  vcol <- RColorBrewer::brewer.pal(5, "Set1")
  vcols <- adjustcolor(rep("cornsilk4", vcount(cong)), alpha.f = 0.5)

  invisible(sapply(seq_along(vrole), function(i){
    n <- as.character(result$node[roles %in% vrole[i]])
    if(!identical(n, character(0))){
      vcols[which(V(cong)$name %in% n)] <<- vcol[i]
    }
  }
  ))
  cornet::igplot(cong, v.c = vcols, v.l = NA, v.s = degree(cong)/vcount(cong)*100)
  return(sortresult)
}

