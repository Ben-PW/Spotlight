###################################################################################

# This script is going to be for the various functions to transform the simulation
# data into something useable. Empty for now, as I'm focussing on testing the 
# simulation output to make sure it's valid, but eventually these will be 
# specified here and integrated into Spotlight_main.R, along with dedicated plotting
# functions

##################################################################################

tagGraphs_init <- function(graph_list,
                       dataset,
                       source,
                       alpha = NA_real_,
                       spotlight_pct = NA_real_,
                       b = NA_real_,
                       miss_level = NA_real_) {
  purrr::imap(graph_list, function(g, replicate_id) {
    igraph::graph_attr(g, "dataset") <- dataset
    igraph::graph_attr(g, "replicate_id") <- as.integer(replicate_id)
    igraph::graph_attr(g, "source") <- source
    igraph::graph_attr(g, "alpha") <- alpha
    igraph::graph_attr(g, "spotlight_pct") <- spotlight_pct
    igraph::graph_attr(g, "b") <- b
    igraph::graph_attr(g, "miss_level") <- miss_level
    g
  })
}

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

tagGraphs <- function(graph_list,
                       dataset,
                       source,
                       alpha = NA_real_,
                       spotlight_pct = NA_real_,
                       b = NA_real_,
                       miss_level = NA_real_) {
  purrr::imap(graph_list, function(g, replicate_id) {
    igraph::graph_attr(g, "dataset") <- dataset
    igraph::graph_attr(g, "replicate_id") <- as.integer(replicate_id)
    igraph::graph_attr(g, "source") <- source
    igraph::graph_attr(g, "alpha") <- alpha
    igraph::graph_attr(g, "spotlight_pct") <- spotlight_pct
    igraph::graph_attr(g, "b") <- b
    igraph::graph_attr(g, "miss_level") <- miss_level
    g
  })
}

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
    igraph::graph_attr(g, "replicate_id") <- replicate_id
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

