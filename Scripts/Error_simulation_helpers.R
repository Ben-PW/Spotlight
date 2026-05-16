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

# Below function no longer required, kept for now in case of any 
# unforseen downstream effects fo removing

#tagGraphs_init <- function(graph_list,
#                           dataset,
#                           source,
#                           alpha = NA_real_,
#                           spotlight_pct = NA_real_,
#                           b = NA_real_,
#                           miss_level = NA_real_) {
#  purrr::imap(graph_list, function(g, replicate_id) {
#    igraph::graph_attr(g, "dataset") <- dataset
#    igraph::graph_attr(g, "replicate_id") <- as.integer(replicate_id)
#    igraph::graph_attr(g, "source") <- source
#    igraph::graph_attr(g, "alpha") <- alpha
#    igraph::graph_attr(g, "spotlight_pct") <- spotlight_pct
#    igraph::graph_attr(g, "b") <- b
#    igraph::graph_attr(g, "miss_level") <- miss_level
#    g
#  })
#}


# Function to create a unique graph id within the spotlight loop
makeGraphID <- function(dataset, replicate_id, source,
                        alpha = NA_real_,
                        spotlight_pct = NA_real_,
                        b = NA_real_,
                        miss_level = NA_real_) {
  paste(
    paste0("ds_", dataset),
    paste0("rep_", replicate_id),
    paste0("src_", source),
    paste0("a_", alpha),
    paste0("sp_", spotlight_pct),
    paste0("b_", b),
    paste0("ml_", miss_level),
    sep = "__"
  )
}

# Function to create and assign graph metadata, including graph id, to the 
# graph as a graph attribute

tagGraphs <- function(graph_list,
                      dataset,
                      source,
                      alpha = NA_real_,
                      spotlight_pct = NA_real_,
                      b = NA_real_,
                      miss_level = NA_real_) {
  purrr::imap(graph_list, function(g, replicate_id) {
    replicate_id <- as.integer(replicate_id)
    
    igraph::graph_attr(g, "dataset") <- dataset
    igraph::graph_attr(g, "replicate_id") <- replicate_id # this is just graph position in list
    igraph::graph_attr(g, "source") <- source
    igraph::graph_attr(g, "alpha") <- alpha
    igraph::graph_attr(g, "spotlight_pct") <- spotlight_pct
    igraph::graph_attr(g, "b") <- b
    igraph::graph_attr(g, "miss_level") <- miss_level
    igraph::graph_attr(g, "graph_id") <- makeGraphID(
      dataset = dataset,
      replicate_id = replicate_id,
      source = source,
      alpha = alpha,
      spotlight_pct = spotlight_pct,
      b = b,
      miss_level = miss_level
    )
    
    g
  })
}


# Compute network level metrics over a list of networks and store as tibble
# with graph metadata extracted from graph attributes

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
      graph_id = igraph::graph_attr(g, "graph_id"),
      
      density = igraph::edge_density(g, loops = FALSE),
      dcent = igraph::centr_degree(g, mode = "all", normalized = TRUE)$centralization,
      clustering = igraph::transitivity(g, type = "global"),
      size = igraph::vcount(g),
      APL = {
        dist_mat <- igraph::distances(g, mode = "all")
        mean(dist_mat[is.finite(dist_mat)], na.rm = TRUE)
      },
      components = igraph::count_components(g)
    )
  })
}

# Compute node level metrics over a list of networks and store as tibble, with
# graph metadata extrcated from graph attributes

computeCentralityDf <- function(graph_list) {
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
      
      Degree_norm = igraph::degree(g, mode = "all", normalized = TRUE),
      Degree_raw = igraph::degree(g, mode = "all", normalized = FALSE),
      
      Betweenness_norm = igraph::betweenness(
        g,
        directed = igraph::is_directed(g),
        normalized = TRUE
      ),
      Betweenness_raw = igraph::betweenness(
        g,
        directed = igraph::is_directed(g),
        normalized = FALSE
      ),
      
      Closeness_norm = igraph::closeness(
        g,
        mode = "all",
        normalized = TRUE
      ),
      Closeness_raw = igraph::closeness(
        g,
        mode = "all",
        normalized = FALSE
      ),
      
      Eigenvector = igraph::eigen_centrality(
        g,
        directed = igraph::is_directed(g)
      )$vector
    )
  })
}

#detach(package:igraph)