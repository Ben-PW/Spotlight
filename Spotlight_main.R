##################################################################################################

# Main for spotlight simulation. Huge refactor from Combine5.R
# Old code generated and stored separate networks for every condition, this code
# generates then dynamically and deletes when no longer required. Much better
# on memory, not sure about how to parallelise yet.

########################################### Begin ##############################################

source('Data_process.R') # ERGMs fit to data
source('Compute_NetStats.R') # Assorted helper functions
source('Spotlight_function.R') # Spotlight relevant functions
source('Data_handlers.R') # Not written yet

#### Simulate networks ####

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

#################################### Begin sim ##################################

#### Check file paths ####

out_dir <- here::here("Results", "node_results")

# create if missing
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# wipe folder before sim
if (length(list.files(out_dir, full.names = TRUE)) > 0) {
  ans <- readline(prompt = "Warning: Files detected in node results folder. Y = delete and continue")
  if (tolower(ans) %in% "y") {
    file.remove(list.files(out_dir, full.names = TRUE))
  } else {
    stop("Exiting")
  }
}

#### Spotlight params ####

# Basic params for testing
spotlight_pcts <- c(0.01, 0.10) # % nodes spotlit
miss_levels <- c(0.10, 0.20) # missingness levels
alphas <- c(0, 0.5, 1, 2) # exponential degree bias
bs <- c(1, 2, 4) # weights for non-spotlit ties

#### Loop setup ####

global_rows <- list()
node_rows <- list()
kg <- 1L # metrics list counter
kn <- 1L # nodes counter
flush <- 50L # save increments
node_batch_id <- 1L # batch ids

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
          
          # Calculate network level metrics
          global_rows[[kg]] <- computeMetrics(obs_list, name = ds) |>
            dplyr::mutate(
              network_label = ds,
              alpha = a,
              spotlight_pct = sp,
              b = bv,
              miss_level = ml
            ) |>
            dplyr::mutate(source = "observed") # tag as error networks
          
          kg <- kg + 1L # increment counter
          
          # Calculate node level centralities
          node_rows[[kn]] <- computeCentralityDf(
            graph_list = obs_list,
            network_label  = ds,
            alpha = a,
            b = bv,
            spotlight_pct = sp,
            miss_level = ml
            ) |>
            dplyr::mutate(source = "observed") # tag as error networks
          
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
            gc(FALSE)
          }
          
          rm(obs_list)
          gc(FALSE) # immediately bin from memory
        }
      }
      
      rm(sp_list)
      gc(FALSE)
    }
  }
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
}

graphResults <- dplyr::bind_rows(global_rows)
#nodeResults <- dplyr::bind_rows(node_rows)

#### Get baseline node centralities ####
# Baseline don't have alpha/b/miss_level so set to NA

nodeGT <- purrr::imap_dfr(datasets, function(base_list, ds) {
  
  # Temp spotlight variable as computeCentralityDf expects it
  base_tmp <- lapply(base_list, function(g) { igraph::V(g)$Spotlight <- 0L; g })
  
  computeCentralityDf(
    graph_list = base_tmp,
    network_label = ds,
    alpha = NA_real_,
    b = NA_real_,
    spotlight_pct = NA_real_,
    miss_level = NA_real_
  ) |> 
    dplyr::mutate(source = "true")
})






