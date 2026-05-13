#####################################################################################################

#Data preprocessing script. Candidate networks for error simulation created here

#####################################################################################################

# requires
source(here::here("Scripts", "Data_simulation_helpers.R"))
source(here::here("Scripts", "Degree_sampler.R"))
source(here::here("Scripts", "ERGM_simulator.R"))

############################## Generate degree sequences #############################

# Create parameter grid 

basis_grid <- tidyr::expand_grid(
  size = c(30, 60, 120),
  average_degree = c(3, 6),
  freeman_centralisation = c(0.1, 0.3, 0.5)
) |>
  dplyr::mutate(
    average_degree_tolerance = calc_ad_tol(size), # allows 0.01 density variation
    name = purrr::pmap_chr(
      list(size, average_degree, freeman_centralisation),
      make_basis_name
    ),
    seed = 123 + dplyr::row_number()
  )

# If below is taking ages, change the max steps argument in the degree sampler
# file, currently set to 500,000

basis_list <- basis_grid |>
  dplyr::mutate(
    result = purrr::pmap(
      list(
        name,
        size,
        average_degree,
        average_degree_tolerance,
        freeman_centralisation,
        seed
      ),
      function(name,
               size, 
               average_degree, 
               average_degree_tolerance, 
               freeman_centralisation,
               seed) {
        
        message("Sampling degseq for ", name, "...")
        
        out <- makeNetworkBasis(
          nsim = 500,
          size = size,
          average_degree = average_degree,
          average_degree_tolerance = average_degree_tolerance,
          freeman_centralisation = freeman_centralisation,
          tolerance = 0.05,
          min_degree = 1,
          seed = seed,
          verbose = FALSE,
          cut_breaks = 4,
          slice_n = 3
        )
        
        message("Finished ", name)
        
        out
      }
    )
  ) |>
  dplyr::select(name, result) |>
  tibble::deframe()

################### n30_ad3_c01 #############

n30_ad3_c01 <- makeNetworkBasis(nsim = 500,
                                size = 30,
                                average_degree = 3,
                                average_degree_tolerance = 0.3,
                                freeman_centralisation = 0.1,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 3)

################### n30_ad6_c01 ##################

n30_ad6_c01 <- makeNetworkBasis(nsim = 500,
                                size = 30,
                                average_degree = 6,
                                average_degree_tolerance = 0.3,
                                freeman_centralisation = 0.1,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 3)

################### n30_ad3_c05 ################################


n30_ad3_c05 <- makeNetworkBasis(nsim = 500,
                                size = 30,
                                average_degree = 3,
                                average_degree_tolerance = 0.3,
                                freeman_centralisation = 0.5,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 3)

################### n30_ad6_c05 #########

n30_ad3_c05 <- makeNetworkBasis(nsim = 500,
                                size = 30,
                                average_degree = 6,
                                average_degree_tolerance = 0.3,
                                freeman_centralisation = 0.5,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 3)


################### n60_ad3_c01 #####################

n60_ad3_c01 <- makeNetworkBasis(nsim = 500,
                                size = 60,
                                average_degree = 3,
                                average_degree_tolerance = 0.59,
                                freeman_centralisation = 0.1,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 3)



################### n60_ad6_c01 ################

n60_ad6_c01 <- makeNetworkBasis(nsim = 500,
                                size = 60,
                                average_degree = 6,
                                average_degree_tolerance = 0.59,
                                freeman_centralisation = 0.1,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 3)

################### n60_ad3_c05 ##########

n60_ad3_c05 <- makeNetworkBasis(nsim = 500,
                                size = 60,
                                average_degree = 3,
                                average_degree_tolerance = 0.59,
                                freeman_centralisation = 0.5,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 3)

################### n60_ad6_c05 ############

n60_ad6_c05 <- makeNetworkBasis(nsim = 500,
                                size = 60,
                                average_degree = 6,
                                average_degree_tolerance = 0.59,
                                freeman_centralisation = 0.5,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 3)

################### n120_ad3_c01 ############

n120_ad3_c01 <- makeNetworkBasis(nsim = 500,
                                size = 120,
                                average_degree = 3,
                                average_degree_tolerance = 1.19,
                                freeman_centralisation = 0.1,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 3)

################### n120_ad6_c01 ###########

n120_ad6_c01 <- makeNetworkBasis(nsim = 500,
                                 size = 120,
                                 average_degree = 6,
                                 average_degree_tolerance = 1.19,
                                 freeman_centralisation = 0.1,
                                 tolerance = 0.05,
                                 min_degree = 1,
                                 seed = 123,
                                 verbose = FALSE,
                                 cut_breaks = 4,
                                 slice_n = 3)

################### n120_ad3_c05 ############

n120_ad3_c05 <- makeNetworkBasis(nsim = 500,
                                 size = 120,
                                 average_degree = 3,
                                 average_degree_tolerance = 1.19,
                                 freeman_centralisation = 0.5,
                                 tolerance = 0.05,
                                 min_degree = 1,
                                 seed = 123,
                                 verbose = FALSE,
                                 cut_breaks = 4,
                                 slice_n = 3)

################### n120_ad6_c05 ########## 

n120_ad6_c01 <- makeNetworkBasis(nsim = 500,
                                 size = 120,
                                 average_degree = 6,
                                 average_degree_tolerance = 1.19,
                                 freeman_centralisation = 0.5,
                                 tolerance = 0.05,
                                 min_degree = 1,
                                 seed = 123,
                                 verbose = FALSE,
                                 cut_breaks = 4,
                                 slice_n = 3)

################################# Testing n120 c04

n120_ad3_c04 <- makeNetworkBasis(nsim = 500,
                                size = 120,
                                average_degree = 3,
                                average_degree_tolerance = 1.19,
                                freeman_centralisation = 0.4,
                                tolerance = 0.05,
                                min_degree = 1,
                                seed = 123,
                                verbose = FALSE,
                                cut_breaks = 4,
                                slice_n = 2)



tgt <- ceiling(200/length(n120_ad3_c04$selected_degree_sequences))

n120_ad3_c04_test <- simulateNetworks(n120_ad3_c04$networks,
                                     nmAtt = 1,
                                     gwdeg = 1,
                                     gwesp = 0.4,
                                     gwdsp = -0.025,
                                     target_connected = tgt,
                                     max_attempts = 300)

summariseNetworks(n120_ad3_c04_test$networks)

par(mfrow = c(2, 2), mar = c(0.1, 0.1, 0.6, 0.1))
plotSimNetworks(n120_ad3_c05_test$networks[1:4])
plotSimNetworks(n120_ad3_c04_test$networks[1:16])

############################# Testing n120 c03

n120_ad3_c03 <- makeNetworkBasis(nsim = 500,
                                 size = 120,
                                 average_degree = 3,
                                 average_degree_tolerance = 1.19,
                                 freeman_centralisation = 0.3,
                                 tolerance = 0.05,
                                 min_degree = 1,
                                 seed = 123,
                                 verbose = FALSE,
                                 cut_breaks = 4,
                                 slice_n = 1)

tgt <- ceiling(200/length(n120_ad3_c03$selected_degree_sequences))

n120_ad3_c03_test <- simulateNetworks(n120_ad3_c03$networks,
                                      nmAtt = 1,
                                      gwdeg = 1,
                                      gwesp = 0.4,
                                      gwdsp = -0.025,
                                      target_connected = tgt,
                                      max_attempts = 300)


###################################### Sim params ################################

# Empirical data finds param values ranging between below
# gwdegree: -1.24 : 3.93
# nodematch: -2.65 : 3.66
# nodefactor: -2.7 : 1.9
# gwesp: 0.27 : 2.47
# edges: -12.45 : 3.21
# gwdsp: -0.69 : 0.66
# size: 6 : 1036
# density: 0.0066 : 0.47

################################################################################

############################ ERGM Simulation Stage #########################

################################################################################

###### Test: n30_ad3_c01 #####

# Good performance across all degree sequences
tgt <- ceiling(100/length(n30_ad3_c01_degs))

n30_ad3_c01_test <- simulateNetworks(n30_ad3_c01_degs,
                                     nmAtt = 1,
                                     gwdeg = 1,
                                     gwesp = 0.4,
                                     gwdsp = -0.025,
                                     target_connected = tgt)

n30_ad3_c01_test$diagnostics
par(mfrow = c(5,5))
plotSimNetworks(n30_ad3_c01_test$networks)

###### Test: n30_ad6_c01 #####

# Good performance across all degree sequences
tgt <- ceiling(100/length(n30_ad6_c01_degs))

n30_ad6_c01_test <- simulateNetworks(n30_ad6_c01_degs,
                                     nmAtt = 1,
                                     gwdeg = 1,
                                     gwesp = 0.4,
                                     gwdsp = -0.025,
                                     target_connected = tgt)

n30_ad6_c01_test$diagnostics
par(mfrow = c(5,5))
plotSimNetworks(n30_ad6_c01_test$networks)

###### Test: n30_ad3_c05 #####

# Good performance across all degree sequences
tgt <- ceiling(100/length(n30_ad3_c05_degs))

n30_ad3_c05_test <- simulateNetworks(n30_ad3_c05_degs,
                                     nmAtt = 1,
                                     gwdeg = 1,
                                     gwesp = 0.4,
                                     gwdsp = -0.025,
                                     target_connected = tgt)

n30_ad3_c05_test$diagnostics
par(mfrow = c(5,5))
plotSimNetworks(n30_ad3_c05_test$networks)

###### Test: n30_ad6_c05 #####

# Good performance across all degree sequences
tgt <- ceiling(100/length(n30_ad6_c05_degs))

n30_ad6_c05_test <- simulateNetworks(n30_ad6_c05_degs,
                                     nmAtt = 1,
                                     gwdeg = 1,
                                     gwesp = 0.4,
                                     gwdsp = -0.025,
                                     target_connected = tgt)

n30_ad6_c05_test$diagnostics
par(mfrow = c(5,5))
plotSimNetworks(n30_ad6_c05_test$networks)

summariseNetworks(n30_ad6_c05_test$networks)

###### Test: n60_ad3_c01 #####

# Poorer performance, increase attempts
tgt <- ceiling(100/length(n60_ad3_c01_degs))

n60_ad3_c01_test <- simulateNetworks(n60_ad3_c01_degs,
                                     nmAtt = 1,
                                     gwdeg = 1,
                                     gwesp = 0.4,
                                     gwdsp = -0.025,
                                     target_connected = tgt,
                                     max_attempts = 300)

n60_ad3_c01_test$diagnostics
par(mfrow = c(5,5))
plotSimNetworks(n60_ad3_c01_test$networks)

###### Test: n60_ad3_c05 #####

# Poorer performance, increase attempts
tgt <- ceiling(100/length(n60_ad3_c05_degs))

n60_ad3_c05_test <- simulateNetworks(n60_ad3_c05_degs,
                                     nmAtt = 1,
                                     gwdeg = 1,
                                     gwesp = 0.4,
                                     gwdsp = -0.025,
                                     target_connected = tgt,
                                     max_attempts = 300)

n60_ad3_c05_test$diagnostics
par(mfrow = c(5,5))
plotSimNetworks(n60_ad3_c05_test$networks)

###### Test: n60_ad6_c01 #####

# Good performance
tgt <- ceiling(100/length(n60_ad6_c01_degs))

n60_ad6_c01_test <- simulateNetworks(n60_ad6_c01_degs,
                                     nmAtt = 1,
                                     gwdeg = 1,
                                     gwesp = 0.4,
                                     gwdsp = -0.025,
                                     target_connected = tgt,
                                     max_attempts = 100)

n60_ad6_c01_test$diagnostics
par(mfrow = c(5,5))
plotSimNetworks(n60_ad6_c01_test$networks)

###### Test: n60_ad6_c05 #####

# Good performance
tgt <- ceiling(100/length(n60_ad6_c05_degs))

n60_ad6_c05_test <- simulateNetworks(n60_ad6_c05_degs,
                                     nmAtt = 1,
                                     gwdeg = 1,
                                     gwesp = 0.4,
                                     gwdsp = -0.025,
                                     target_connected = tgt,
                                     max_attempts = 100)

n60_ad6_c05_test$diagnostics
par(mfrow = c(5,5))
plotSimNetworks(n60_ad6_c05_test$networks)

###### Test: n120_ad3_c01 #####

# This one is problematic
# Increase number of attempts
tgt <- ceiling(100/length(n120_ad3_c01_degs))

n120_ad3_c01_test <- simulateNetworks(n120_ad3_c01_degs,
                                      nmAtt = 1,
                                      gwdeg = 1,
                                      gwesp = 0.4,
                                      gwdsp = -0.025,
                                      target_connected = tgt,
                                      max_attempts = 500)

n120_ad3_c01_test$diagnostics
sum(n120_ad3_c01_test$diagnostics$accepted)
par(mfrow = c(5,5))
plotSimNetworks(n120_ad3_c01_test$networks)

###### Test: n120_ad3_c05 #####

# Good lord this one is rough, will need to up target accepted, oversample
# degree sequences and increase max attempts
tgt <- ceiling(200/length(n120_ad3_c05_degs))

n120_ad3_c05_test <- simulateNetworks(n120_ad3_c05_degs,
                                      nmAtt = 1,
                                      gwdeg = 1,
                                      gwesp = 0.4,
                                      gwdsp = -0.025,
                                      target_connected = tgt,
                                      max_attempts = 500)

n120_ad3_c05_test$diagnostics
sum(n120_ad3_c05_test$diagnostics$accepted)
par(mfrow = c(2,2))
par(mfrow = c(3, 3), mar = c(0.1, 0.1, 0.6, 0.1))
plotSimNetworks(n120_ad3_c05_test$networks)
summariseNetworks(n120_ad3_c05_test$networks)

###### Test: n120_ad6_c01 #####

# Good performance
tgt <- ceiling(100/length(n120_ad6_c01_degs))

n120_ad6_c01_test <- simulateNetworks(n120_ad6_c01_degs,
                                      nmAtt = 1,
                                      gwdeg = 1,
                                      gwesp = 0.4,
                                      gwdsp = -0.025,
                                      target_connected = tgt)

n120_ad6_c01_test$diagnostics
par(mfrow = c(4,4))
plotSimNetworks(n120_ad6_c01_test$networks)

###### Test: n120_ad6_c05 #####

tgt <- ceiling(100/length(n120_ad6_c05_degs))

n120_ad6_c05_test <- simulateNetworks(n120_ad6_c05_degs,
                                      nmAtt = 1,
                                      gwdeg = 1,
                                      gwesp = 0.4,
                                      gwdsp = -0.025,
                                      target_connected = tgt)

n120_ad6_c05_test$diagnostics
par(mfrow = c(4,4))
plotSimNetworks(n120_ad6_c05_test$networks)

################################################################################

################################ Collect datasets ###############################

test_batch_1 <- list(
  n30_ad3_c01 = n30_ad3_c01_test$networks,
  n30_ad6_c01 = n30_ad6_c01_test$networks,
  n30_ad3_c05 = n30_ad3_c05_test$networks,
  n30_ad6_c05 = n30_ad6_c05_test$networks,
  
  n60_ad3_c01 = n60_ad3_c01_test$networks,
  n60_ad6_c01 = n60_ad6_c01_test$networks,
  n60_ad3_c05 = n60_ad3_c05_test$networks,
  n60_ad6_c05 = n60_ad6_c05_test$networks,
  
  n120_ad3_c01 = n120_ad3_c01_test$networks,
  n120_ad6_c01 = n120_ad6_c01_test$networks,
  n120_ad3_c05 = n120_ad3_c05_test$networks,
  n120_ad6_c05 = n120_ad6_c05_test$networks
)

saveRDS(
  test_batch_1,
  file = here::here("Data", "Test1", "test_batch_1.rds")
)
