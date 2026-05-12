# Function to check basic descriptives of simulated network lists

summariseNetworks <- function(net_list) {
  
  data.frame(
    #target_density = 0.0066,
    mean_density = mean(sapply(net_list, network::network.density)),
    
    #target_triangles = 7817,
    mean_triangles = mean(sapply(net_list, function(g)
      sna::triad.census(g, mode = "graph")[4])),
    
    #target_degree = 6.79,
    mean_degree = mean(sapply(net_list, function(g)
      mean(sna::degree(g, gmode = "graph")))),
    
    #target_max = 316,
    mean_max_degree = mean(sapply(net_list, function(g)
      max(sna::degree(g, gmode = "graph")))),
    
    #target_comp = 1,
    mean_components = mean(sapply(net_list, function(g)
      length(sna::component.dist(g)$csize))),
    
    mean_component_coverage = mean(sapply(net_list, function(g)
      (max(sna::component.dist(g)$csize)/network.size(g))*100)),
    
    mean_centralisation = mean(sapply(net_list, function(g)
      igraph::centr_degree(intergraph::asIgraph(g), normalized = TRUE)$centralization)
    ))
}

# Function to turn generated degree sequences into network objects

netFromDegSeq <- function(degree_sequences) {
  
  g_list <- lapply(degree_sequences, function(x) {
    igraph::realize_degseq(
      x,
      allowed.edge.types = "simple",
      method = "smallest"
    )
  })
  
  net_list <- lapply(g_list, intergraph::asNetwork)
  
  return(net_list)
}

# Function to simply plot a list of simulated networks

plotSimNetworks <- function(net_list) {
  
  obj_name <- deparse(substitute(net_list))
  
  # flatten one level if needed
  if (is.list(net_list[[1]]) && !inherits(net_list[[1]], "network")) {
    net_list <- unlist(net_list, recursive = FALSE)
  }
  
  for (i in seq_along(net_list)) {
    plot(
      net_list[[i]],
      main = paste0(obj_name, " ", i)
    )
  }
}