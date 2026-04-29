#####################################################################################################

#Data preprocessing script. Candidate networks for error simulation created here

#####################################################################################################

source("Degree_sequence_sampler2.R")

here::here()

# Function to check descriptives of simulated networks
network_summary <- function(net_list) {
  
  data.frame(
    #target_density = 0.0066,
    mean_density = mean(sapply(net_list, network::network.density)),
    
    #target_triangles = 7817,
    mean_triangles = mean(sapply(net_list, function(g)
      sna::triad.census(g, mode = "graph")[4])),
    
    #target_degree = 6.79,
    mean_degree = mean(sapply(net_list, function(g)
      mean(sna::degree(g, gmode = "graph")))),
    
    #target_max = 316,
    mean_max_degree = mean(sapply(net_list, function(g)
      max(sna::degree(g, gmode = "graph")))),
    
    #target_comp = 1,
    mean_components = mean(sapply(net_list, function(g)
      length(sna::component.dist(g)$csize))),
    
    mean_component_coverage = mean(sapply(net_list, function(g)
      (max(sna::component.dist(g)$csize)/network.size(g))*100)),
    
    mean_centralisation = mean(sapply(net_list, function(g)
      igraph::centr_degree(intergraph::asIgraph(g), normalized = TRUE)$centralization)
    ))
}

simulateNetworks <- function(net_list, 
                             nsim = 1,
                             nfAtt = 0,
                             nmAtt = 0,
                             gwdeg = 0.5,
                             gwesp = 0.5,
                             gwdsp = -0.025) {
  
  all_sims <- list()
  counter <- 1
  
  for (i in seq_along(net_list)) {
    
    net <- net_list[[i]]
    
    n <- network::network.size(net)
    
    # Assign attributes ONCE per basis network
    network::set.vertex.attribute(
      net,
      attrname = "att",
      value = sample(c("A", "B"), n, replace = TRUE, prob = c(3, 1))
    )
    
    form <- net ~
      nodefactor("att") +
      nodematch("att") +
      gwdegree(0.3, fixed = TRUE) +
      gwesp(0.3, fixed = TRUE) +
      gwdsp(0.3, fixed = TRUE)
    
    coefs <- c(
      nodefactor.att.B = nfAtt,
      nodematch.att = nmAtt,
      gwdeg.fixed = gwdeg,
      gwesp.fixed = gwesp,
      gwdsp.fixed = gwdsp
    )
    
    sim <- ergm::simulate_formula(
      form,
      constraints = ~degreedist,
      coef = coefs,
      nsim = nsim,
      output = "network"
    )
    
    # Flatten into one list
    for (j in seq_len(nsim)) {
      all_sims[[counter]] <- sim[[j]]
      counter <- counter + 1
    }
  }
  
  return(all_sims)
}

netFromDegSeq <- function(degree_sequences) {
  
  g_list <- lapply(degree_sequences, function(x) {
    igraph::realize_degseq(
      x,
      allowed.edge.types = "simple",
      method = "smallest"
    )
  })
  
  net_list <- lapply(g_list, intergraph::asNetwork)
  
  return(net_list)
}

plotSimNetworks <- function(net_list) {
  
  obj_name <- deparse(substitute(net_list))
  
  # flatten one level if needed
  if (is.list(net_list[[1]]) && !inherits(net_list[[1]], "network")) {
    net_list <- unlist(net_list, recursive = FALSE)
  }
  
  for (i in seq_along(net_list)) {
    plot(
      net_list[[i]],
      main = paste0(obj_name, " ", i)
    )
  }
}
############################## Generate degree sequences #############################

######################## n50 ####

############## ad3 ####

##### c01 ####

n50ad3c01 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 50,
  average_degree = 3,
  freeman_centralisation = 0.1,
  tolerance = 0.01,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = TRUE,
  verbose = TRUE,
  store_trace = TRUE
)

##### c05 ####

n50ad3c05 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 50,
  average_degree = 3,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

############## ad6 ####

##### c01 ####

n50ad6c01 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 50,
  average_degree = 6,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

##### c05 ####

n50ad6c05 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 50,
  average_degree = 6,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

######################## n100 ####

############## ad3 ####

##### c01 ####

n100ad3c01 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

##### c05 ####

n100ad3c05 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

############## ad6 ####

##### c01 ####

n100ad6c01 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 100,
  average_degree = 6,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)
##### c05 ####

n100ad6c05 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 100,
  average_degree = 6,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

######################## n150 ####

############## ad3 ####

##### c01 ####

n150ad3c01 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 150,
  average_degree = 3,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

##### c05 ####

n150ad3c05 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 150,
  average_degree = 3,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

############## ad6 ####

##### c01 ####

n150ad6c01 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 150,
  average_degree = 6,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

##### c05 ####

n150ad6c05 <- degree_sequence_sample_mcmc(
  nsim = 10,
  size = 150,
  average_degree = 6,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)



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

#### Supervision feedback
# Keep attributes out of primary design
# Focus:
# Density, centralisation
# Density: High - Low
# Centralisation: High - Low
# Attribute: No effect - High effect
# Core periph: dense, centralised
# Cell struct: sparser, decentralised
# Don't model archetypes directly, relate variable levels to whichever is
# closest
# Model whether attribute is related to spotlight or not later
# Adding an option to save networks would be useful for replication

############################## IMPORTANT
# gwdegree coefficient is interpreted REVERSE
# https://environmentalpolicy.ucdavis.edu/blog/shiny-app-help-interpret-gw-degree-estimates-ergms
# https://joss.theoj.org/papers/10.21105/joss.00036
# https://www.sciencedirect.com/science/article/pii/S0378873306000396
# It is interpreted as an anti-perferential attachment term apparently

# Set up environment for some simulation param experiments

####################### Trial network parameters
# Here i need to trial parameters that are representative of real world networks
# but also produce useable data for the simulation

sapply(net_list, function(g) {
  igraph::centr_degree(
    intergraph::asIgraph(g),
    normalized = TRUE
  )$centralization
})

par(mfrow = c(4, 5), mar = c(0.2, 0.2, 1, 0.2))

n100ad3c01_list <- netFromDegSeq(n100ad3c01$degree_sequences)

n100ad3c01sim <- simulateNetworks(n100ad3c01_list,
                                  gwdeg = 1,
                                  gwesp = 1,
                                  gwdsp = 1,
                                  nsim = 2)

network_summary(n100ad3c01sim)
plotSimNetworks(n100ad3c01sim)

sapply(n100ad3c01sim, function(g) {
  igraph::centr_degree(
    intergraph::asIgraph(g),
    normalized = TRUE
  )$centralization
})

n100ad3c05_list <- netFromDegSeq(n100ad3c05$degree_sequences)

n100ad3c05sim <- simulateNetworks(n100ad3c05_list,
                                  gwdeg = 1,
                                  gwesp = 0.1,
                                  gwdsp = 0.1,
                                  nsim = 2)

network_summary(n100ad3c05sim)
plotSimNetworks(n100ad3c05sim)

for (i in 1:20) {
  plot(trial[[i]], main = paste0("sim ", i))
}












































































simulateNetwork <- function(size,
                            avdeg,
                            prob = c(3,1),
                            nfAtt = 0,
                            nmAtt = 0,
                            gwdeg = 1,
                            gwesp = 1,
                            gwdsp = -0.025,
                            nsim = 20){
  
  net <- network::network.initialize(size, directed = FALSE)
  
  network::set.vertex.attribute(
    net,
    attrname = "att",
    value = sample(c("A", "B"), size, replace = TRUE, prob = prob)
  )
  
  net <- setDensity(net, av_deg = avdeg)
  
  form <- net ~
    nodefactor("att") +
    nodematch("att") +
    gwdegree(0.3, fixed = TRUE) +
    gwesp(0.3, fixed = TRUE) +
    gwdsp(0.3, fixed = TRUE)
  
  coefs <- c(
    nodefactor.att.B = nfAtt,
    nodematch.att = nmAtt,
    gwdeg.fixed = gwdeg,
    gwesp.fixed = gwesp,
    gwdsp.fixed = gwdsp
  )
  
  sim <- stats::simulate(
    form,
    constraints = ~edges,
    coef = coefs,
    nsim = nsim,
    output = "network"
  )
  
  sim
  
}

############################### LESSONS FROM TEST 1
# Don't include nodefactor, it results in networks which are far too centralised and
# induces a huge number of isolates. It likely interacts significantly with negative
# gwdeg parameter values, which just gets messy.
# For next tests, likely switch to nodematch, just as a structural control

test1 <- simulateNetwork(size = 500,
                        avdeg = 3,
                        prob = c(2,1),
                        nfAtt = 0,
                        gwdeg = 1,
                        nsim = 20)

test2 <- simulateNetwork(size = 500,
                         avdeg = 3,
                         prob = c(2,1),
                         nfAtt = 0,
                         gwdeg = -1,
                         nsim = 20)

test3 <- simulateNetwork(size = 500,
                         avdeg = 3,
                         prob = c(2,1),
                         nfAtt = 1,
                         gwdeg = -1,
                         nsim = 20)

test4 <- simulateNetwork(size = 500,
                         avdeg = 3,
                         prob = c(2,1),
                         nfAtt = 1,
                         gwdeg = 1,
                         nsim = 20)

test5 <- simulateNetwork(size = 100,
                         avdeg = 3,
                         prob = c(2,1),
                         nfAtt = 0,
                         gwdeg = 1,
                         nsim = 20)

test6 <- simulateNetwork(size = 100,
                         avdeg = 3,
                         prob = c(2,1),
                         nfAtt = 0,
                         gwdeg = -1,
                         nsim = 20)

test7 <- simulateNetwork(size = 100,
                         avdeg = 3,
                         prob = c(2,1),
                         nfAtt = 1,
                         gwdeg = -1,
                         nsim = 20)

test8 <- simulateNetwork(size = 100,
                         avdeg = 3,
                         prob = c(2,1),
                         nfAtt = 1,
                         gwdeg = 1,
                         nsim = 20)

testgwdeg0 <- simulateNetwork(size = 100,
                              avdeg = 4,
                              prob = c(1,1),
                              nfAtt = 0,
                              nmAtt = 1,
                              gwdeg = 1,
                              nsim = 20)

par(mfrow = c(4, 5), mar = c(0.2, 0.2, 1, 0.2))


network::set.vertex.attribute(
    net,
    attrname = "att",
    value = sample(c("A", "B"), size, replace = TRUE, prob = prob)
  )
  
  net <- setDensity(net, av_deg = avdeg)
  
  form <- net ~
    nodefactor("att") +
    nodematch("att") +
    gwdegree(0.3, fixed = TRUE) +
    gwesp(0.3, fixed = TRUE) +
    gwdsp(0.3, fixed = TRUE)
  
  coefs <- c(
    nodefactor.att.B = nfAtt,
    nodematch.att = nmAtt,
    gwdeg.fixed = gwdeg,
    gwesp.fixed = gwesp,
    gwdsp.fixed = gwdsp
  )
  
  sim <- stats::simulate(
    form,
    constraints = ~edges,
    coef = coefs,
    nsim = nsim,
    output = "network"
  )
  
  sim
  
for (i in 1:20) {
  plot(test2[[i]], main = paste0("test_2 ", i))
}
for (i in 1:20) {
  plot(test3[[i]], main = paste0("test_3 ", i))
}
for (i in 1:20) {
  plot(test4[[i]], main = paste0("test_4 ", i))
}
for (i in 1:20) {
  plot(test5[[i]], main = paste0("test_5 ", i))
}
for (i in 1:20) {
  plot(test6[[i]], main = paste0("test_6 ", i))
}
for (i in 1:20) {
  plot(test7[[i]], main = paste0("test_7 ", i))
}
for (i in 1:20) {
  plot(test8[[i]], main = paste0("test_8 ", i))
}
for (i in 1:20) {
  plot(testgwdeg0[[i]], main = paste0("test_nm1 ", i))
}
summary(test1 ~ edges)

network_summary(test1)
network_summary(test2)
network_summary(test3)
network_summary(test4)
network_summary(test5)
network_summary(test6)
network_summary(test7)
network_summary(test8)
network_summary(testgwdeg0)

################################# LESSONS FROM TEST2
# These parameters are producing very decentralised networks
# gwesp may be producing eqalitarian structures?
# likely need to use nodefactor to encourage hubs

par(mfrow = c(4, 5), mar = c(0.2, 0.2, 1, 0.2))

simulateNetwork <- function(size,
                            avdeg,
                            prob = c(3,1),
                            nfAtt = 0,
                            nmAtt = 0,
                            gwdeg = 1,
                            #gwesp = 1,
                            #gwdsp = -0.025,
                            nsim = 20){
  
  net <- network::network.initialize(size, directed = FALSE)
  
  network::set.vertex.attribute(
    net,
    attrname = "att",
    value = sample(c("A", "B"), size, replace = TRUE, prob = prob)
  )
  
  net <- setDensity(net, av_deg = avdeg)
  
  form <- net ~
    nodefactor("att") +
    nodematch("att") +
    gwdegree(0.3, fixed = TRUE) #+
    #gwesp(0.3, fixed = TRUE) #+
    #gwdsp(0.3, fixed = TRUE)
  
  coefs <- c(
    nodefactor.att.B = nfAtt,
    nodematch.att = nmAtt,
    gwdeg.fixed = gwdeg#,
    #gwesp.fixed = gwesp,
    #gwdsp.fixed = gwdsp
  )
  
  sim <- stats::simulate(
    form,
    constraints = ~edges,
    coef = coefs,
    nsim = nsim,
    output = "network"
  )
  
  sim
  
}

test2 <- simulateNetwork(size = 100,
                         avdeg = 3,
                         nfAtt = 1,
                         gwdeg = -0.5,
                         nsim = 20)

for (i in 1:20) {
  plot(test2[[i]], main = paste0("test_2 ", i))
}

########################## PLACEHOLDER NETWORKS TO TEST PIPELINE #######################
########################################## Florentine Families ######################################
# 16x16 undirected, low density, isolates

data("flo", package = "network")

flo <- network::network(flo, directed = F)

network::set.vertex.attribute(flo,
                              "wealth",
                              c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3))
#flo %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)

flo1 <- ergm::ergm(flo ~ edges + absdiff("wealth"))
#rm(flo)

####################################### Books network ##################################

# For some reason I can't access these datasets any more
# data(books)

###############################  ############################
# 4 x 24 x 24 weighted adjacency matrices

data("covert_17", package = "networkdata")

ac1 <- intergraph::asNetwork(covert_17)
ac2 <- intergraph::asNetwork(covert_17) # just using the same data for now


# plot(ac1)
# plot(ac2)

acct1 <- ergm::ergm(ac1 ~ edges) 
summary(acct1)

acct2 <- ergm::ergm(ac2 ~ edges)
summary(acct2)

#rm(ac1, ac2)


