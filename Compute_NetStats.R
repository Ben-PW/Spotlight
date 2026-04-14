#####################################################################################################

# Helper functions for simulation script

####################################################################################################

#library(igraph)

# Coerce directed networks to symmetrical

undirect <- function(graph_list) {
  lapply(graph_list, function(g) {
    if (igraph::is_directed(g)) {
      g <- igraph::as.undirected(g, mode = "collapse")
    }
    g
  })
}

# Assign IDs for testing between node removal or addition situations

IDNodes <- function(graph_list){
  lapply(graph_list, function(g){
    igraph::V(g)$NodeID <- seq_len(igraph::vcount(g))
    g
  })
}

computeMetrics <- function(graph_list) {
  purrr::map_dfr(graph_list, function(g) {
    tibble::tibble(
      dataset = igraph::graph_attr(g, "dataset"),
      replicate_id = igraph::graph_attr(g, "replicate_id"),
      source = igraph::graph_attr(g, "source"),
      alpha = igraph::graph_attr(g, "alpha"),
      spotlight_pct = igraph::graph_attr(g, "spotlight_pct"),
      b = igraph::graph_attr(g, "b"),
      miss_level = igraph::graph_attr(g, "miss_level"),
      #stage = igraph::graph_attr(g, "stage"),
      graph_id = igraph::graph_attr(g, "graph_id"),
      
      density = igraph::edge_density(g, loops = FALSE),
      dcent = igraph::centr_degree(g, mode = "all", normalized = TRUE)$centralization,
      clustering = igraph::transitivity(g, type = "global"),
      size = igraph::vcount(g),
      APL = {
        dist_mat <- igraph::distances(g, mode = "all")
        mean(dist_mat[is.finite(dist_mat)], na.rm = TRUE)
      }
    )
  })
}

computeCentralityDf <- function(graph_list, normalized = FALSE) {
  purrr::map_dfr(graph_list, function(g) {
    tibble::tibble(
      dataset = igraph::graph_attr(g, "dataset"),
      replicate_id = igraph::graph_attr(g, "replicate_id"),
      source = igraph::graph_attr(g, "source"),
      alpha = igraph::graph_attr(g, "alpha"),
      spotlight_pct = igraph::graph_attr(g, "spotlight_pct"),
      b = igraph::graph_attr(g, "b"),
      miss_level = igraph::graph_attr(g, "miss_level"),
      graph_id = igraph::graph_attr(g, "graph_id"),
      
      NodeID = as.integer(igraph::V(g)$NodeID),
      Spotlight = as.integer(igraph::V(g)$Spotlight),
      # ADD THE ATTRIBUTE LOGIC HERE
      
      Degree = igraph::degree(g, mode = "all", normalized = normalized),
      Betweenness = igraph::betweenness(
        g,
        directed = igraph::is_directed(g),
        normalized = normalized
      ),
      Closeness = igraph::closeness(
        g,
        mode = "all",
        normalized = normalized
      ),
      Eigenvector = igraph::eigen_centrality(
        g,
        directed = igraph::is_directed(g)
      )$vector
    )
  })
}

#detach(package:igraph)