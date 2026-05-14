#####################################################################################################

#Data preprocessing script. Candidate networks for error simulation created here

#####################################################################################################

# requires
source(here::here("Scripts", "Data_simulation_helpers.R"))
source(here::here("Scripts", "Degree_sampler.R"))
source(here::here("Scripts", "ERGM_simulator.R"))

############################## Generate degree sequences #############################

############# Create parameter grid ############

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


############### Sample degree sequences ##########

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


################## Simulate from degree sequences ##########

datasets <- purrr::imap(
  basis_list,
  function(basis, name) {
    message("Simulating networks from basis ", name, "...")
    
    out <- simulateFromBasis(basis, verbose = TRUE)
    
    message("Completed simulation for basis ", name)
    
    out
  }
)


#saveRDS(
#  test_batch_1,
#  file = here::here("Data", "Test1", "test_batch_1.rds")
#)
