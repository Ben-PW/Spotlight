##################################################################################################

# Main for spotlight simulation. Huge refactor from Combine5.R
# Old code generated and stored separate networks for every condition, this code
# generates then dynamically and deletes when no longer required. Much better
# on memory, not sure about how to parallelise yet.

########################################### Begin ##############################################

#### Importing data and functions ####
# This needs to be done before loading packages as some scripts load and detach packages

source('Data_process.R') # ERGMs fit to data
source('Compute_NetStats.R') # Assorted helper functions
source('Spotlight_function.R') # Spotlight relevant functions
source('Data_handlers.R')

#library(dplyr)

####################################### Simulate networks #######################################

# Toy networks for now

iFlo <- stats::simulate(flo1, nsim = 100) # 100 networks per condition
iAc1 <- stats::simulate(acct1, nsim = 100)
iAc2 <- stats::simulate(acct2, nsim = 100)

rm(flo1, acct1, acct2)

##### Define datasets #####

datasets <- list(
 Flo = iFlo,
#  Ac1 = iAc1,
 Ac2 = iAc2
)

##### Convert all to igraph objects #####

datasets <- lapply(datasets, function(netlist) {
  lapply(netlist, intergraph::asIgraph)
})

##### Ensure all are undirected #####

datasets <- lapply(datasets, undirect)

##### Assign ID to nodes #####

datasets <- lapply(datasets, IDNodes)

##### Compute Ground Truth metrics #####

GTAll <- Map(computeMetrics, datasets, names(datasets)) |>
  dplyr::bind_rows()

#### Spotlight params ####

# Basic params for testing
spotlight_pcts <- c(0.01, 0.10) # % nodes spotlit
miss_levels    <- c(0.10, 0.20) # missingness levels
alphas         <- c(0, 0.5, 1, 2) # exponential degree bias
bs             <- c(1, 2, 4) # weights for non-spotlit ties

#### Begin sim ####

global_rows <- list()
node_rows   <- list()
kg <- 1L # metrics list counter
kn <- 1L # nodes counter

set.seed(123)

for (ds in names(datasets)) {
  
  base_list <- datasets[[ds]]
  
  for (a in alphas) {
    for (sp in spotlight_pcts) {
      
      # Do spotlight here so it's consistent across ml and b levels
      sp_list <- assignSpotlight(base_list, spotlight_pct = sp, alpha = a)
      
      for (bv in bs) {
        for (ml in miss_levels) {
          
          # apply spotlight
          obs_list <- sampleSpotlight(sp_list, miss_level = ml, b = bv)
          
          # Global metrics
          global_rows[[kg]] <- computeMetrics(obs_list, name = ds) |>
            dplyr::mutate(
              network_label = ds,
              alpha = a,
              spotlight_pct = sp,
              b = bv,
              miss_level = ml
            )
          kg <- kg + 1L
          
          # Node centralities
          node_rows[[kn]] <- computeCentralityDf(
            graph_list = obs_list,
            network_label  = ds,
            alpha = a,
            b = bv,
            spotlight_pct = sp,
            miss_level = ml
          )
          kn <- kn + 1L
          
          rm(obs_list)
          gc(FALSE) # immediately bin from memory
        }
      }
      
      rm(sp_list)
      gc(FALSE)
    }
  }
}

testResults <- dplyr::bind_rows(global_rows)
nodeResults <- dplyr::bind_rows(node_rows)

#### Get baseline node centralities
# Baseline don't have alpha/b/miss_level so set to NA

nodeGT <- purrr::imap_dfr(datasets, function(base_list, ds) {
  
  # Temp spotlight variable as computeCentralityDf expects it
  base_tmp <- lapply(base_list, function(g) { igraph::V(g)$Spotlight <- 0L; g })
  
  computeCentralityDf(
    graph_list     = base_tmp,
    network_label  = ds,
    alpha          = NA_real_,
    b              = NA_real_,
    spotlight_pct  = NA_real_,
    miss_level     = NA_real_
  ) |> dplyr::mutate(source = "true")
})

# Tag observed nodes
nodeResults <- nodeResults |> dplyr::mutate(source = "observed")




