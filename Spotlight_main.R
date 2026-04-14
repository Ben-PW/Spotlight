##################################################################################################

# Main for spotlight simulation. Huge refactor from Combine5.R
# Old code generated and stored separate networks for every condition, this code
# generates then dynamically and deletes when no longer required. Much better
# on memory, not sure about how to parallelise yet.

########################################### Begin ##############################################

source('Data_process.R') # ERGMs fit to data
source('Compute_NetStats.R') # Assorted helper functions
source('Spotlight_function.R') # Spotlight relevant functions
source('Data_handlers.R') # Network ID system and associated handlers

################################### Check file paths ###########################################

node_gt_dir <- here::here("Results", "node_gt")
net_gt_dir <- here::here("Results", "net_gt")
node_out_dir <- here::here("Results", "node_results")
net_out_dir <- here::here("Results", "network_results")

# create if missing

if (!dir.exists(node_gt_dir)) {
  dir.create(node_gt_dir, recursive = TRUE)
}

if (!dir.exists(net_gt_dir)) {
  dir.create(net_gt_dir, recursive = TRUE)
}

if (!dir.exists(node_out_dir)) {
  dir.create(node_out_dir, recursive = TRUE)
}

if (!dir.exists(net_out_dir)) {
  dir.create(net_out_dir, recursive = TRUE)
}

# wipe folder before sim if necessary to avoid duplicates

nodefiles <- length(list.files(node_out_dir, full.names = TRUE))
nodefilesgt <- length(list.files(node_gt_dir, full.names = TRUE))
netfiles <- length(list.files(net_out_dir, full.names = TRUE))
netfilesgt <- length(list.files(net_gt_dir, full.names = TRUE))

if (nodefiles + nodefilesgt + netfiles + netfilesgt > 0) {
  ans <- readline(prompt = "Warning: Files detected in results folders. Y = delete and continue:")
  if (tolower(ans) %in% "y") {
    file.remove(list.files(node_out_dir, full.names = TRUE))
    file.remove(list.files(net_out_dir, full.names = TRUE))
    file.remove(list.files(node_gt_dir, full.names = TRUE))
    file.remove(list.files(net_gt_dir, full.names = TRUE))
  } else {
    stop("Exiting")
  }
}

################################## Generate networks ######################################

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

##### Assign ID to networks #####

datasets <- purrr::imap(datasets, function(graph_list, ds) {
  tagGraphs_init(
    graph_list = graph_list,
    dataset = ds,
    source = "true"
  )
})

##### Assign ID to nodes #####

datasets <- lapply(datasets, IDNodes)

##### Compute ground truth metrics #####

graphGT <- purrr::map_dfr(datasets, computeMetrics)

nodeGT <- purrr::map_dfr(datasets, function(graph_list) {
  graph_list <- lapply(graph_list, function(g) {
    if (is.null(igraph::V(g)$Spotlight)) { #compute centrality expects spotlight
      igraph::V(g)$Spotlight <- 0L
    }
    g
  })
  computeCentralityDf(graph_list)
})

# save gt results to disk 

saveRDS(graphGT,
        file = here::here(
          "Results", "net_gt", "network_results_gt.rds"
        ))

saveRDS(nodeGT,
        file = here::here(
          "Results", "node_gt", "node_results_gt.rds"
        ))


#################################### Begin sim ##################################

#### Spotlight params ####

# Basic params for testing
spotlight_pcts <- c(0.01, 0.025, 0.05, 0.075, 0.10) # % nodes spotlit
miss_levels <- c(0.10, 0.20, 0.30, 0.40, 0.50) # missingness levels
alphas <- c(0, 0.5, 1, 2, 3) # exponential degree bias
bs <- c(1, 2, 3, 4) # weights for non-spotlit ties

#### Loop setup ####

global_rows <- list()
node_rows <- list()
kg <- 1L # metrics list counter
kn <- 1L # nodes counter
flush <- 50L # save increments
node_batch_id <- 1L # node batch ids
network_batch_id <- 1L # network batch ids

set.seed(123)

##### Begin loop #####

################ IMPORTANT
################ IMPORTANT
################ LOGIC FOR EXTRACTING "att" VARIABLE IN SIMULATED NETWORKS NEEDS ADDING
################ RELEVANT SECTION IN Compute_NetStats.R HAS BEEN COMMENTED

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
          
          # store simulation parameters
          obs_list <- tagGraphs(
            graph_list = obs_list,
            dataset = ds,
            source = "observed",
            alpha = a,
            spotlight_pct = sp,
            b = bv,
            miss_level = ml
          )
          
          # Calculate network level metrics
          global_rows[[kg]] <- computeMetrics(obs_list) 
          
          kg <- kg + 1L # increment counter
          
          # Batch save network calculations to disk every 50 loops and 
          # re-initialise required variables
          
          if ((kg - 1L) >= flush) {
            network_batch <- dplyr::bind_rows(global_rows)
            
            saveRDS(
              network_batch,
              file = here::here(
                "Results", "network_results",
                paste0("network_results_batch_", network_batch_id, ".rds")
              )
            )
            
            global_rows <- list() # refresh list for future loops
            kg <- 1L # reset node metrics loop counter
            network_batch_id <- network_batch_id + 1L # increment batch counter
            rm(network_batch)
            gc(FALSE)
          }
          
          # Calculate node level centralities
          node_rows[[kn]] <- computeCentralityDf(obs_list)
          
          kn <- kn + 1L # increment counter
          
          # Batch save node calculations to disk every 50 loops and 
          # re-initialise required variables
          
          if ((kn - 1L) >= flush) {
            node_batch <- dplyr::bind_rows(node_rows)
            
            saveRDS(
              node_batch,
              file = here::here(
                "Results", "node_results",
                paste0("node_results_batch_", node_batch_id, ".rds")
              )
            )
            
            node_rows <- list() # refresh list for future loops
            kn <- 1L # reset node metrics loop counter
            node_batch_id <- node_batch_id + 1L # increment batch counter
            rm(node_batch)
            gc(FALSE)
          }
          
          rm(obs_list)
          #gc(FALSE) # immediately bin from memory
        }
      }
      
      rm(sp_list)
      gc(FALSE)
    }
  }
}

##### End of loop cleanup #####

# End of loop flush if results remain but kg < flush
if (length(global_rows) > 0) {
  network_batch <- dplyr::bind_rows(global_rows)
  saveRDS(
    network_batch,
    file = here::here(
      "Results", "network_results",
      paste0("network_results_batch_", network_batch_id, ".rds")
    )
  )
  rm(network_batch)
}

# End of loop flush if results remain but kn < flush
if (length(node_rows) > 0) {
  node_batch <- dplyr::bind_rows(node_rows)
  saveRDS(
    node_batch,
    file = here::here(
      "Results", "node_results",
      paste0("node_results_batch_", node_batch_id, ".rds")
    )
  )
  rm(node_batch)
}

# Cleanup

rm(a, 
   alphas, 
   ans,
   bs, 
   bv, 
   ds, 
   flush, 
   kg, 
   kn, 
   miss_levels, 
   ml, 
   node_batch_id, 
   network_batch_id,
   node_out_dir,
   net_out_dir,
   sp,
   spotlight_pcts)


print("Loop completed and environment cleaned")

################################# Simulation Complete ###############################





