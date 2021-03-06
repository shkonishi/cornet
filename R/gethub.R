#' Node centrality
#' @description get hub node.
#' @usage gethub(g, com_fun)
#' @param g igraph object
#' @param com_fun character: function name of igraph for community detection, "fastgreedy.community" or "cluster_louvain"
#' @examples
#' # sample data
#' data(cl_dat)
#' cormat <- cor(cl_dat[["3"]])
#' res <- cornet::corgraph(mat = cormat)
#' g <- res$undir.graph
#' reshub <- cornet::gethub(g=g, com_fun="fastgreedy.community")
#' #
#' @importFrom igraph delete.vertices degree fastgreedy.community cluster_louvain delete.vertices V E vcount degree betweenness
#' @importFrom grDevices adjustcolor
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats sd setNames
#' @importFrom graphics legend
#' @export
gethub <- function(g, com_fun="cluster_louvain"){

  # argument check: 'g' is an igraph object ----
  if (!class(g) == "igraph"){
    stop("'g' must be an igraph object.")
  }
  # argument check: com_fun ----
  com.algs <- c("fastgreedy.community", "cluster_louvain")
  if(any(com.algs %in% com_fun)){
    commu <- eval(parse(text=paste0("igraph::", com_fun, "(g)")))
    mem <- commu$membership # module membership
    num_mod <- max(mem)  # number of max module
    num_nodes <- igraph::vcount(g) # number of vertex
    deg <- igraph::degree(g) # degree

    # print modularity
    print(paste0("modularity:", round(commu$modularity, 3)))

  }else{
    stop('select from "fastgreedy.community" or "cluster_louvain"')
  }


  # Calculation of the within-module degree ----
  z_score <- numeric(num_nodes)
  for(s in 1:num_mod){
    v_seq <- subset(igraph::V(g),mem == s)
    g_sub <- igraph::delete.vertices(g,subset(igraph::V(g), mem != s))
    deg_sub <- igraph::degree(g_sub)
    z_score[v_seq] <- (deg_sub - mean(deg_sub)) / sd(deg_sub)
  }

  # Calculation of the participation coefficient ----
  participation_coeff <- numeric(num_nodes) + 1.0
  for(i in 1:num_nodes){
    mem_nei <- mem[igraph::neighbors(g,i)]
    for(s in 1:num_mod){
      deg_to_s <- length(subset(mem_nei,mem_nei == s))
      participation_coeff[[i]] <- participation_coeff[[i]] - (deg_to_s / deg[[i]]) ** 2.0
    }
  }

  # Classification ----
  role <- vector(mode = "numeric", length = num_nodes)
  pc <- setNames(as.data.frame(cbind(z_score,participation_coeff)),
                 c("within_module_degree","participation_coefficient"))


  # R1: Ultra-peripheral nodes
  v_seq <- which(pc$within_module_degree < 2.5 & pc$participation_coeff < 0.05)
  role[v_seq] <- "R1: Ultra-peripheral nodes"

  # R2: Peripheral nodes
  v_seq <- which(pc$within_module_degree < 2.5 & pc$participation_coeff >= 0.05 & pc$participation_coeff < 0.625)
  role[v_seq] <- "R2: Peripheral nodes"

  # R3: Non-hub connectors
  v_seq <- which(pc$within_module_degree < 2.5 & pc$participation_coeff >= 0.625 & pc$participation_coeff < 0.8)
  role[v_seq] <- "R3: Non-hub connectors"

  # R4: Non-hub kinless nodes
  v_seq <- which(pc$within_module_degree < 2.5 & pc$participation_coeff >= 0.8)
  role[v_seq] <- "R4: Non-hub kinless nodes"

  # R5: Provincial hubs
  v_seq <- which(pc$within_module_degree >= 2.5 & pc$participation_coeff < 0.3)
  role[v_seq] <- "R5: Provincial hubs"

  # R6: Connector hubs
  v_seq <- which(pc$within_module_degree >= 2.5 & pc$participation_coeff >= 0.3 & pc$participation_coeff < 0.75)
  role[v_seq] <- "R6: Connector hubs"

  # R7: Kinless hubs
  v_seq <- which(pc$within_module_degree >= 2.5 & pc$participation_coeff >= 0.75)
  role[v_seq] <- "R7: Kinless hubs"

  # result dataframe----
  #pc <- cbind(pc,role=role)
  result <- data.frame(node=igraph::V(g)$name,
                       module_member = mem,
                       num_of_members = sapply(seq_along(mem), function(i)sum(mem %in% mem[i])),
                       between=igraph::betweenness(g),
                       degree=deg,
                       ndegree=deg/igraph::vcount(g),
                       pc,
                       role,
                       stringsAsFactors = F)
  sortresult <- result[order(result$between, decreasing = T),]

  # igplot with role
  ## delete non connective node and connected graph ----
  cong <- igraph::delete.vertices(g, igraph::V(g)$name[igraph::degree(g)==0])

  ## facter of roles ----
  roles <- sapply(strsplit(as.character(sortresult$role[match(igraph::V(cong)$name, sortresult$node)]), "\\:"), "[", 1)
  role_fctr <- factor(roles, levels=c("R7","R6","R5","R4","R3","R2","R1","0"))

  ## node colour ----
  vcol <- adjustcolor(c(RColorBrewer::brewer.pal(5, "Set1"), rep("cornsilk4", 3)), alpha.f = 0.5)
  vcols <- vcol[role_fctr]

  ## node size ----
  vsiz <- c(7,7,5,5,5,2,2,2)
  vsizs <- vsiz[role_fctr]

  ## node label ----
  vlabs <- ifelse(roles %in% c("R7","R6","R5","R4","R3"), igraph::V(cong)$name, NA)

  ## modified igraph object with attribuete ----
  igraph::V(cong)$v.c <- vcols
  igraph::V(cong)$v.s <- vsizs
  igraph::V(cong)$v.l <- vlabs

  ## igplot ----
  # cornet::igplot(cong, v.c = vcols, v.s = vsizs, e.w = 0.5, v.l=vlabs, v.l.c = "black")
  # graphics::legend("topleft",pt.cex = 3, border = F,
  #                  legend = c("R7:Kinless hubs", "R6:Connector hubs", "R5:Provincial hubs",
  #                             "R4:Non-hub kinless nodes", "R3:Non-hub connectors"),
  #                  bg = "transparent",pch = 20, col = vcol, cex = 0.8)

  # return modified igraph object and data.frame of hub extranction. -----
  return(setNames(list(cong, sortresult), c("hubbed_graph", "res_hub")))
}
