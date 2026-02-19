#################################################################################################

# Spotlight function now uses exponential weighting

####################################################################################

# library(netUtils) - probs don't need anymore tbf

# Function to assign spotlight conditional on spotlight degree
# It's a bit hacky but for attribute related spotlight, I'm just going to manually
# assign 'Degree' as a proxy for attribute

assignSpotlight <- function(graph_list, spotlight_pct, alpha = 0) {
  lapply(graph_list, function(g) {
    
    n <- igraph::vcount(g)
    k <- max(1L, round(n * spotlight_pct))
    
    Spotlight <- integer(n)
    
    deg <- igraph::degree(g)
    w <- (deg + 1) ^ alpha
    
    idx <- sample.int(n, size = k, replace = FALSE, prob = w)
    Spotlight[idx] <- 1L
    
    # Assign to nodes
    igraph::V(g)$Spotlight <- Spotlight
    
    # Assign to edges
      ends_mat <- igraph::ends(g, igraph::E(g), names = FALSE)
      igraph::E(g)$Spotlight <-
        as.integer(Spotlight[ends_mat[,1]] | Spotlight[ends_mat[,2]])
  
    
    g
  })
}


# Function to sample nodes conditionally on spotlight (b = 1 means equal)
# b > 1 biases towards non-spotlit
# might want to change the weighting system later but not difficult

sampleSpotlight <- function(graph_list, miss_level, b = 1) {
  lapply(graph_list, function(g) {
    
    m <- igraph::ecount(g)
    k <- round(m * miss_level)
    
    sp <- igraph::E(g)$Spotlight
    
    w <- ifelse(sp == 0, b, 1) # if no spotlight, assign weight b
    
    idx_drop <- sample.int(m, size = k, replace = FALSE, prob = w)
    
    igraph::delete_edges(g, igraph::E(g)[idx_drop])
  })
}

#iFlo <- assignSpotlight(datasets$Flo, 0.1, 1)
#iFlo <- sampleSpotlight(datasets$iFlo, 0.1, 1)
