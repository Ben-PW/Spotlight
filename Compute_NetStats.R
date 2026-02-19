#####################################################################################################

# Helper functions for simulation script

####################################################################################################

#library(igraph)

# Coerce directed networks to symmetrical

undirect <- function(graph_list) {
  lapply(graph_list, function(g) {
    if (igraph::is.directed(g)) {
      g <- igraph::as.undirected(g, mode = "collapse")
    }
    g
  })
}

# Assign IDs for testing between node removal or addition situations

IDNodes <- function(graph_list){
  lapply(graph_list, function(g){
    V(g)$NodeID <- seq_len(igraph::vcount(g))
    g
  })
}

##### Compute network level metrics #####
# Very easy to add more metric here if needed
# Note to self
# I think I need to address the rounding here?

computeMetrics <- function(graph_list, name) {
  data.frame(
    id = paste0(name,"_", seq_along(graph_list)),
    
    density = sapply(graph_list, function(g) {
      igraph::edge_density(g, loops = FALSE)
    }),
    
    dcent = sapply(graph_list, function(g) {
      cent <- igraph::centr_degree(g, mode = "all", normalized = TRUE)
      cent$centralization
    }),
    
    clustering = sapply(graph_list, function(g) {
      igraph::transitivity(g, type = "global")
    }),
    
    size = sapply(graph_list, igraph::vcount),
    
    APL = sapply(graph_list, function(g) {
      # all‐pairs shortest paths; exclude Infs
      dist_mat <- igraph::distances(g, mode = "all")
      mean(dist_mat[is.finite(dist_mat)], na.rm = TRUE)
    })
  )
}


computeCentralityDf <- function(graph_list,
                                network_label,
                                alpha,
                                b,
                                spotlight_pct,
                                miss_level,
                                normalized = FALSE) {
  
  purrr::imap_dfr(graph_list, function(g, net_id) {
    
    Degree <- igraph::degree(g, mode = "all", normalized = normalized)
    Betweenness <- igraph::betweenness(
      g,
      directed = igraph::is_directed(g), # this is handy if interested in directed later
      normalized = normalized
    )
    Closeness <- igraph::closeness(
      g,
      mode = "all",
      normalized = normalized # normalised for now but that's cosmetic
    )
    Eigenvector <- igraph::eigen_centrality(
      g,
      directed = igraph::is_directed(g)
    )$vector
    
    # Write results to tibble
    tibble::tibble(
      net_id = net_id,
      network_label = network_label,
      alpha = alpha,
      b = b,
      spotlight_pct = spotlight_pct,
      miss_level = miss_level,
      NodeID = as.integer(igraph::V(g)$NodeID),
      Spotlight = as.integer(igraph::V(g)$Spotlight),
      Degree = Degree,
      Betweenness = Betweenness,
      Closeness = Closeness,
      Eigenvector = Eigenvector
    )
  })
}

#detach(package:igraph)