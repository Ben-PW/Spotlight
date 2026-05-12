#####################################################################################################

#Data preprocessing script. Candidate networks for error simulation created here

#####################################################################################################

# requires
source(here::here("Scripts", "Data_simulation_helpers.R"))
source(here::here("Scripts", "Degree_sampler.R"))
source(here::here("Scripts", "ERGM_simulator.R"))

############################## Generate degree sequences #############################

################################## Low size #############################

#################### Low avdeg ###################

########## Low c ###################

n30_ad3_c01 <- sampleDegSeq(
  nsim = 500, # this argument is now defunct due to sampler changes
  size = 30,
  average_degree = 3,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  verbose = FALSE
)

n30_ad3_c01_tab <- summariseDegSeq(n30_ad3_c01)

n30_ad3_c01_tab %>%
  summarise(
    n = n(),
    min_ad = min(average_degree),
    max_ad = max(average_degree),
    mean_ad = mean(average_degree), # 3.03
    min_c = min(centralisation),
    mean_c = mean(centralisation), # 0.112
    max_c = max(centralisation),
    min_dmax = min(dmax),
    max_dmax = max(dmax)
  )

selected_ids <- n30_ad3_c01_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4),    
    c_bin  = cut(centralisation, breaks = 4),    
    iqr_bin = cut(degree_iqr, breaks = 4),       
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 2) %>% # sample 2 as there are fewer sequences in this cell
  ungroup() %>%
  pull(seq_id)

n30_ad3_c01_degs <- n30_ad3_c01$degree_sequences[selected_ids]

n30_ad3_c01_degs <- netFromDegSeq(n30_ad3_c01_degs)

network_summary(n30_ad3_c01_degs) # av deg = 3.032
# av cent = 0.102

################################# Low size ###############################

#################### High avdeg ##############

########## Low C ################

n30_ad6_c01 <- sampleDegSeq(
  nsim = 50,
  size = 30,
  average_degree = 6,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  verbose = FALSE
)

n30_ad6_c01_tab <- summariseDegSeq(n30_ad6_c01)

n30_ad6_c01_tab %>%
  summarise(
    n = n(),
    min_ad = min(average_degree),
    max_ad = max(average_degree),
    mean_ad = mean(average_degree), # 6.26
    min_c = min(centralisation),
    mean_c = mean(centralisation), # 0.126
    max_c = max(centralisation),
    min_dmax = min(dmax),
    max_dmax = max(dmax)
  )

selected_ids <- n30_ad6_c01_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4),   
    c_bin  = cut(centralisation, breaks = 4),    
    iqr_bin = cut(degree_iqr, breaks = 4),       
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 2) %>%
  ungroup() %>%
  pull(seq_id)

n30_ad6_c01_degs <- n30_ad6_c01$degree_sequences[selected_ids]

n30_ad6_c01_degs <- netFromDegSeq(n30_ad6_c01_degs)

network_summary(n30_ad6_c01_degs) # avdeg = 6.019
# avcent = 0.115

################################# Low size ################################

#################### Low avdeg ######################

########## High C #######################

n30_ad3_c05 <- sampleDegSeq(
  nsim = 50,
  size = 30,
  average_degree = 3,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n30_ad3_c05_tab <- summariseDegSeq(n30_ad3_c05)

n30_ad3_c05_tab %>%
  summarise(
    n = n(),
    min_ad = min(average_degree),
    max_ad = max(average_degree),
    mean_ad = mean(average_degree), # 3.03
    min_c = min(centralisation),
    mean_c = mean(centralisation), # 0.5
    max_c = max(centralisation),
    min_dmax = min(dmax),
    max_dmax = max(dmax)
  )

selected_ids <- n30_ad3_c05_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4), 
    c_bin  = cut(centralisation, breaks = 4),   
    iqr_bin = cut(degree_iqr, breaks = 4),  
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, 
           c_bin, 
           iqr_bin, 
           dmax_bin
  ) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)

n30_ad3_c05_degs <- n30_ad3_c05$degree_sequences[selected_ids]

n30_ad3_c05_degs <- netFromDegSeq(n30_ad3_c05_degs)

network_summary(n30_ad3_c05_degs) # av deg = 3.018
# av cent = 0.468

################################# Low size #########
#################### High avdeg ############
########## High C #########

n30_ad6_c05 <- sampleDegSeq(
  nsim = 50,
  size = 30,
  average_degree = 6,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n30_ad6_c05_tab <- summariseDegSeq(n30_ad6_c05)

n30_ad6_c05_tab %>%
  summarise(
    n = n(),
    min_ad = min(average_degree),
    max_ad = max(average_degree),
    mean_ad = mean(average_degree), # 3.03
    min_c = min(centralisation),
    mean_c = mean(centralisation), # 0.5
    max_c = max(centralisation),
    min_dmax = min(dmax),
    max_dmax = max(dmax)
  )

selected_ids <- n30_ad6_c05_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4), 
    c_bin  = cut(centralisation, breaks = 4),   
    iqr_bin = cut(degree_iqr, breaks = 4),  
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, 
           c_bin, 
           iqr_bin, 
           dmax_bin
  ) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)

n30_ad6_c05_degs <- n30_ad6_c05$degree_sequences[selected_ids]

n30_ad6_c05_degs <- netFromDegSeq(n30_ad6_c05_degs)

network_summary(n30_ad6_c05_degs)

################################# Mid size #####################
#################### Low avdeg #########
########## Low C #########

n60_ad3_c01 <- sampleDegSeq(
  nsim = 1000,
  size = 60,
  average_degree = 3,
  average_degree_tolerance = 0.59,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n60_ad3_c01_tab <- summariseDegSeq(n60_ad3_c01)

selected_ids <- n60_ad3_c01_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4),    
    c_bin  = cut(centralisation, breaks = 4),   
    iqr_bin = cut(degree_iqr, breaks = 4),      
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)

n60_ad3_c01_degs <- n60_ad3_c01$degree_sequences[selected_ids]

n60_ad3_c01_degs <- netFromDegSeq(n60_ad3_c01_degs)

network_summary(n60_ad3_c01_degs) # av deg = 2.971
# av cent = 0.103

################################# Mid size ################
#################### High avdeg ##########
########## Low C ##########

n60_ad6_c01 <- sampleDegSeq(
  nsim = 50,
  size = 60,
  average_degree = 6,
  average_degree_tolerance = 0.59,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n60_ad6_c01_tab <- summariseDegSeq(n60_ad6_c01)

selected_ids <- n60_ad6_c01_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4),   
    c_bin  = cut(centralisation, breaks = 4), 
    iqr_bin = cut(degree_iqr, breaks = 4),      
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)

n60_ad6_c01_degs <- n60_ad6_c01$degree_sequences[selected_ids]

n60_ad6_c01_degs <- netFromDegSeq(n60_ad6_c01_degs)

network_summary(n60_ad6_c01_degs) # av deg = 5.989
# av cent = 0.120

################################# Mid size ##########
#################### Low avdeg ########### 
########## High C #########

n60_ad3_c05 <- sampleDegSeq(
  nsim = 50,
  size = 60,
  average_degree = 3,
  average_degree_tolerance = 0.59,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n60_ad3_c05_tab <- summariseDegSeq(n60_ad3_c05)

selected_ids <- n60_ad3_c05_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4), 
    c_bin  = cut(centralisation, breaks = 4),   
    iqr_bin = cut(degree_iqr, breaks = 4),  
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)


n60_ad3_c05_degs <- n60_ad3_c05$degree_sequences[selected_ids]

n60_ad3_c05_degs <- netFromDegSeq(n60_ad3_c05_degs)

network_summary(n60_ad3_c05_degs) # av deg = 3.118
# av cent = 0.483

################################# Mid size ############
#################### High avdeg ############
########## High C ############

n60_ad6_c05 <- sampleDegSeq(
  nsim = 50,
  size = 60,
  average_degree = 6,
  average_degree_tolerance = 0.59,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n60_ad6_c05_tab <- summariseDegSeq(n60_ad6_c05)

selected_ids <- n60_ad6_c05_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4),    
    c_bin  = cut(centralisation, breaks = 4),    
    iqr_bin = cut(degree_iqr, breaks = 4),       
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)

n60_ad6_c05_degs <- n60_ad6_c05$degree_sequences[selected_ids]

n60_ad6_c05_degs <- netFromDegSeq(n60_ad6_c05_degs)

network_summary(n60_ad6_c05_degs) # av deg = 6.051
# av cent = 0.483

################################# High size ############
#################### Low avdeg ##########
########## Low C ############

n120_ad3_c01 <- sampleDegSeq(
  nsim = 500,
  size = 120,
  average_degree = 3,
  average_degree_tolerance = 1.19,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n120_ad3_c01_tab <- summariseDegSeq(n120_ad3_c01)

n120_ad3_c01_tab %>%
  summarise(
    n = n(),
    min_ad = min(average_degree),
    max_ad = max(average_degree),
    mean_ad = mean(average_degree),
    min_c = min(centralisation),
    max_c = max(centralisation),
    min_dmax = min(dmax),
    max_dmax = max(dmax)
  )

selected_ids <- n120_ad3_c01_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4),   
    c_bin  = cut(centralisation, breaks = 4),    
    iqr_bin = cut(degree_iqr, breaks = 4),       
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)

n120_ad3_c01_degs <- n120_ad3_c01$degree_sequences[selected_ids]

n120_ad3_c01_degs <- netFromDegSeq(n120_ad3_c01_degs)

network_summary(n120_ad3_c01_degs) # av deg = 3.317
# av cent = 0.094

################################# High size ###########
#################### High avdeg ########## 
########## Low C ###########

n120_ad6_c01 <- sampleDegSeq(
  nsim = 500,
  size = 120,
  average_degree = 6,
  average_degree_tolerance = 1.19,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n120_ad6_c01_tab <- summariseDegSeq(n120_ad6_c01)

n120_ad6_c01_tab %>%
  summarise(
    n = n(),
    min_ad = min(average_degree),
    max_ad = max(average_degree),
    mean_ad = mean(average_degree),
    min_c = min(centralisation),
    max_c = max(centralisation),
    min_dmax = min(dmax),
    max_dmax = max(dmax)
  )

selected_ids <- n120_ad6_c01_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4), 
    c_bin  = cut(centralisation, breaks = 4),    
    iqr_bin = cut(degree_iqr, breaks = 4),       
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)

n120_ad6_c01_degs <- n120_ad6_c01$degree_sequences[selected_ids]

n120_ad6_c01_degs <- netFromDegSeq(n120_ad6_c01_degs)

network_summary(n120_ad6_c01_degs) # av deg = 5.84
# av cent = 0.112

################################# High size ############
#################### Low avdeg #############
########## High C ##########

n120_ad3_c05 <- sampleDegSeq(
  nsim = 500,
  size = 120,
  average_degree = 3,
  average_degree_tolerance = 1.19,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n120_ad3_c05_tab <- summariseDegSeq(n120_ad3_c05)

n120_ad3_c05_tab %>%
  summarise(
    n = n(),
    min_ad = min(average_degree),
    max_ad = max(average_degree),
    mean_ad = mean(average_degree),
    min_c = min(centralisation),
    max_c = max(centralisation),
    min_dmax = min(dmax),
    max_dmax = max(dmax)
  )

selected_ids <- n120_ad3_c05_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4), 
    c_bin  = cut(centralisation, breaks = 4),    
    iqr_bin = cut(degree_iqr, breaks = 4),       
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)

n120_ad3_c05_degs <- n120_ad3_c05$degree_sequences[selected_ids]

n120_ad3_c05_degs <- netFromDegSeq(n120_ad3_c05_degs)

network_summary(n120_ad3_c05_degs)

################################# High size ########## 
#################### High avdeg ########### 
########## High C #######

n120_ad6_c05 <- sampleDegSeq(
  nsim = 500,
  size = 120,
  average_degree = 6,
  average_degree_tolerance = 1.19,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

n120_ad6_c05_tab <- summariseDegSeq(n120_ad6_c05)

n120_ad6_c05_tab %>%
  summarise(
    n = n(),
    min_ad = min(average_degree),
    max_ad = max(average_degree),
    mean_ad = mean(average_degree),
    min_c = min(centralisation),
    max_c = max(centralisation),
    min_dmax = min(dmax),
    max_dmax = max(dmax)
  )

selected_ids <- n120_ad6_c05_tab %>%
  mutate(
    ad_bin = cut(average_degree, breaks = 4), 
    c_bin  = cut(centralisation, breaks = 4),   
    iqr_bin = cut(degree_iqr, breaks = 4),      
    dmax_bin = cut(dmax, breaks = 4)
  ) %>%
  group_by(ad_bin, c_bin, iqr_bin, dmax_bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(seq_id)


n120_ad6_c05_degs <- n120_ad6_c05$degree_sequences[selected_ids]

n120_ad6_c05_degs <- netFromDegSeq(n120_ad6_c05_degs)

network_summary(n120_ad6_c05_degs) # av deg = 6.155
# av cent = 0.491



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

tgt <- ceiling(100/length(n120_ad3_c05_degs))

n120_ad3_c05_test <- simulateNetworks(n120_ad3_c05_degs,
                                      nmAtt = 1,
                                      gwdeg = 1,
                                      gwesp = 0.4,
                                      gwdsp = -0.025,
                                      target_connected = tgt,
                                      max_attempts = 500)

n120_ad3_c05_test$diagnostics
sum(n120_ad3_c05_test$diagnostics$accepted)
par(mfrow = c(5,5))
plotSimNetworks(n120_ad3_c05_test$networks)

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

