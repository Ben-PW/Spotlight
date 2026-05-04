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

################################### Database setup ###########################################

# All file management is now done with respect to current wd

here::here()

# Check for local results folder and create if missing

if(!dir.exists(here::here("Results"))){
  dir.create(here::here("Results"), recursive = TRUE)
}

# .duckdb file is a database file which stores results as tables, not files

db_path <- here::here("Results", "spotlight_results.duckdb")

# check for exisiting db so don't append simulation to any test results

if (file.exists(db_path)) {
  ans <- readline(prompt = "Existing DuckDB detected. Y = delete and continue:")
  if (tolower(ans) == "y") {
    file.remove(db_path)
    file.remove(paste0(db_path, ".wal")) # remove any log files as well
  } else {
    stop("Exiting")
  }
}

# connect to db

con <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = db_path
)

on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

################################## Generate networks ######################################

# Toy networks for now

### THIS SHOULD BE ergm::simulate_formula()
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
      igraph::V(g)$Spotlight <- NA_integer_ # tag as NA for gt graphs
    }
    g
  })
  computeCentralityDf(graph_list)
})

# save gt results to db
DBI::dbWriteTable(con, "network_results_gt", graphGT, overwrite = TRUE)
DBI::dbWriteTable(con, "node_results_gt", nodeGT, overwrite = TRUE)


#################################### Begin sim ##################################

#### Spotlight params ####

# Basic params for testing
spotlight_pcts <- c(0.01, 0.05, 0.10) # % nodes spotlit
miss_levels <- c(0.10, 0.30, 0.50) # missingness levels
alphas <- c(0, 10) # exponential degree bias ##### TEST PARAMETERS HERE ####
bs <- c(1, 2, 4) # weights for non-spotlit ties

#### Loop setup ####

global_rows <- list()
node_rows <- list()
kg <- 1L # metrics list counter
kn <- 1L # nodes counter
flush <- 50L # save increments
node_batch_id <- 1L # node batch ids
network_batch_id <- 1L # network batch ids

#### Progress setup ####

start_time <- Sys.time()

ds_names <- names(datasets)
ds_total <- length(ds_names)
ds_counter <- 0L

set.seed(123)

##### Begin loop #####

################ IMPORTANT
################ IMPORTANT
################ LOGIC FOR EXTRACTING "att" VARIABLE IN SIMULATED NETWORKS NEEDS ADDING
################ RELEVANT SECTION IN Compute_NetStats.R HAS BEEN COMMENTED

tryCatch(
  
  # Main loop for most of error simulation
  
  expr = {
    
    for(ds in names(datasets)){
      
      # Count datasets for basic progress bar
      ds_counter <- ds_counter + 1L
    
      base_list <- datasets[[ds]]
    
      # Iterate over degree bias
      for (a in alphas) {
      
        # Iterate over spotlight percentages
        for (sp in spotlight_pcts) {
        
          # Do spotlight here so it's consistent across ml and b levels
          sp_list <- assignSpotlight(base_list, spotlight_pct = sp, alpha = a)
        
          # Iterate over non-spotlit tie weight
          for (bv in bs) {
          
            # Iterate over missingness levels
            for (ml in miss_levels) {
            
              #### Main spotlight and data storage section below ####
            
              tryCatch({
              
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
              
                kg <- kg + 1L
              
                # Append network calculations to db every 50 loops
                if ((kg - 1L) >= flush) {
                  network_batch <- dplyr::bind_rows(global_rows)
                
                 DBI::dbWriteTable(
                   con,
                   "network_results",
                   network_batch,
                   append = DBI::dbExistsTable(con, "network_results")
                 )
                
                  global_rows <- list()
                  kg <- 1L
                  network_batch_id <- network_batch_id + 1L
                  rm(network_batch)
                  gc(FALSE)
                }
              
                # Calculate node level centralities
                node_rows[[kn]] <- computeCentralityDf(obs_list)
              
                kn <- kn + 1L
              
                # Append node calculations to db every 50 loops
                if ((kn - 1L) >= flush) {
                  node_batch <- dplyr::bind_rows(node_rows)
                
                  DBI::dbWriteTable(
                    con,
                    "node_results",
                    node_batch,
                    append = DBI::dbExistsTable(con, "node_results")
                  )
                
                  node_rows <- list()
                  kn <- 1L
                  node_batch_id <- node_batch_id + 1L
                  rm(node_batch)
                  gc(FALSE)
                }
              
                rm(obs_list)
              
              }, 
            
              # Error handler for any issues writing the results to the database
              # very basic, just lists the dataset and error conditions
            
              error = function(e) {
              
                msg <- paste0(
                  "ERROR | dataset=", ds,
                  " | alpha=", a,
                  " | spotlight_pct=", sp,
                  " | b=", bv,
                  " | miss_level=", ml,
                  " | message=", e$message
                )
              
                message(msg)
              
                write(
                  paste(Sys.time(), msg, sep = " | "),
                  file = here::here("Results", "error_log.txt"),
                  append = TRUE
                )
              
                rm(list = "obs_list")
                gc(FALSE)
              
              })
            }
          }
        
        rm(sp_list)
        gc(FALSE)
        
        }
      }
      
      # End of for(ds in datasets) block, this is a basic progress tracker to 
      # print how far along the datasets we are
      
      message(
        paste0(
          "Dataset ", ds_counter, "/", ds_total,
          " completed: ", ds,
          " | elapsed: ",
          round(difftime(Sys.time(), start_time, units = "mins"), 2),
          " mins"
        )
      )
    }
  },
  
  # Final flush of remaining values and disconnect from db
  
  finally = {
    
    if (length(global_rows) > 0) {
      network_batch <- dplyr::bind_rows(global_rows)
      
      DBI::dbWriteTable(
        con,
        "network_results",
        network_batch,
        append = DBI::dbExistsTable(con, "network_results")
      )
      
      rm(network_batch)
    }
    
    if (length(node_rows) > 0) {
      node_batch <- dplyr::bind_rows(node_rows)
      
      DBI::dbWriteTable(
        con,
        "node_results",
        node_batch,
        append = DBI::dbExistsTable(con, "node_results")
      )
      
      rm(node_batch)
    }
    
    if (exists("con") && DBI::dbIsValid(con)) {
      DBI::dbDisconnect(con, shutdown = TRUE)
    }
  }
)

##### End of loop cleanup #####

rm(list = intersect(
  c("a", "alphas", "ans", "bs", "bv", "ds", "flush", "kg", "kn",
    "miss_levels", "ml", "node_batch_id", "network_batch_id",
    "sp", "spotlight_pcts"),
  ls()
))

# Disconnect db, not necessarily required but there for completeness and I am new
# to this method

if (DBI::dbIsValid(con)) {
  DBI::dbDisconnect(con, shutdown = TRUE)
}

print("Loop completed and environment cleaned")

################################# Simulation Complete ###############################





