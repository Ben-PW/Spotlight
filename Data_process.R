#####################################################################################################

#Data preprocessing script. Import, pre-process, fit ERGMs and check GoFs here. 

#####################################################################################################

#### IMPORTANT !!! ####
# If you are trying to run the demo section, you must uncomment below

# install.packages("remotes")
# remotes::install_github("schochastics/networkdata")

here::here()

######################################### Simulation #################################

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
      length(sna::component.dist(g)$csize)))
  )
}

# Uses specified average degree to determine number of edges added to network
getEdgeNo <- function(nodes, av_deg) {
  round((av_deg * nodes) / 2)
}

# Adds random ties to network with number specified by average degree
setDensity <- function(net, av_deg) {
  nodes <- network.size(net)
  possible <- utils::combn(nodes, 2)
  m <- getEdgeNo(nodes = nodes, av_deg = av_deg)
  idx <- sample(seq_len(ncol(possible)), m, replace = FALSE)
  chosen <- possible[, idx, drop = FALSE]
  
  for (j in seq_len(ncol(chosen))) {
    network::add.edge(net, chosen[1, j], chosen[2, j])
  }
  net
}


library(ergm)
library(network)



###################################### Sim params ################################

# gwdegree: -1.24 : 3.93
# nodematch: -2.65 : 3.66
# nodefactor: -2.7 : 1.9
# gwesp: 0.27 : 2.47
# edges: -12.45 : 3.21
# gwdsp: -0.69 : 0.66
# size: 6 : 1036
# density: 0.0066 : 0.47

# fixing density and abandoning the edges parameter will be better

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

################################## 'Low trust' archetype ############################

# In theory, a low trust network could take two forms, a cell structure with sparse
# inter-group connections, or perhaps a more centralised network with low triadic 
# closure
# The expanded 9/11 terrorist network might be a good example here, a dense cluster
# in the hijacking cell, with very sparse connections outside it, with some highly
# central figures
# In a larger network the PIRA stage 3 example holds, significant within-brigade
# clustering, with sparse connections outside that and some very highly centralised 
# figures

# This proves quite difficult to specify manually, however...

# Update: Below specification is useful for producing cell like structures, slightly
# negative gwdsp means inter-group connections are sparse, but triadic closure
# encourages clustering within cells

#form <- n3 ~
#  edges +
#  nodematch("role") +
#  nodefactor("stat") +
#  gwdegree(0.3, fixed = TRUE) +
#  gwesp(0.3, fixed = TRUE) +
#  gwdsp(0.3, fixed = TRUE) # would expect gwdsp to be significant in this case

#coefs <- c(
#  edges = -5.7,
#  nodematch.role = 1.5,
#  nodefactor.stat.B = 0,
#  gwdeg.fixed = 2,
#  gwesp.fixed = 2,
#  gwdsp.fixed = -0.03 # very light touch with this otherwise no LCC
#)

# The below specification resulted in quite nice looking cell structure
# Groups were constrained to be equal sizes, produced more pronounced
# /visual clustering

#n3 %v% "role" <- sample(c("A", "B", "C", "D"),
#                        n, TRUE, prob = c(8, 8, 8, 8))
#n3 %v% "stat" <- sample(c("A", "B"),
#                        n, TRUE, prob = c(6, 3))

#form <- n3 ~
#  edges +
#  nodematch("role") +
#  nodefactor("stat") +
#  gwdegree(0.3, fixed = TRUE) +
#  gwesp(0.3, fixed = TRUE) +
#  gwdsp(0.3, fixed = TRUE) # would expect gwdsp to be significant in this case

#coefs <- c(
#  edges = -5.7,
#  nodematch.role = 2,
#  nodefactor.stat.B = 0,
#  gwdeg.fixed = 2,
#  gwesp.fixed = 2,
#  gwdsp.fixed = -0.02
#)

# Below specification is nice for slightly sparser networks
#coefs <- c(
#  edges = -5.3,
#  nodematch.role = 1.7,
#  nodefactor.stat.B = 0,
#  gwdeg.fixed = 1.7,
#  gwesp.fixed = 1.5,
#  gwdsp.fixed = -0.02
#)

############################## IMPORTANT
# gwdegree coefficient is interpreted REVERSE
# https://environmentalpolicy.ucdavis.edu/blog/shiny-app-help-interpret-gw-degree-estimates-ergms
# https://joss.theoj.org/papers/10.21105/joss.00036
# https://www.sciencedirect.com/science/article/pii/S0378873306000396
# It is interpreted as an anti-perferential attachment term apparently

n <- 60

n3 <- network.initialize(n,
                         directed = FALSE)

n3 %v% "role" <- sample(c("A", "B", "C"),
                        n, TRUE, prob = c(8, 8, 8))
n3 %v% "stat" <- sample(c("A", "B"),
                        n, TRUE, prob = c(6, 3))


n3 <- setDensity(n3, av_deg = 3)

form <- n3 ~
  nodematch("role") +
  nodefactor("stat") +
  gwdegree(0.3, fixed = TRUE) +
  gwesp(0.3, fixed = TRUE) +
  gwdsp(0.3, fixed = TRUE) # w

coefs <- c(
  nodematch.role = 1.7,
  nodefactor.stat.B = 0,
  gwdeg.fixed = 1.4,
  gwesp.fixed = 1.5,
  gwdsp.fixed = -0.025
)

sim33 <- simulate(
  form,
  constraints = ~edges,
  coef = coefs,
  nsim = 20,
  output = "network"
  )


coefs <- c(
  edges = -4.9,
  nodematch.role = 1.4,
  nodefactor.stat.B = 0.5,
  gwdeg.fixed = 1,
  gwesp.fixed = 1,
  gwdsp.fixed = -0.025
)

sim4 <- simulate(
  form,
  coef = coefs,
  nsim = 20,
  output = "network"
)


cols <- as.factor()
par(mfrow = c(4, 5), mar = c(0.2, 0.2, 1, 0.2))

for (i in 1:20) {
  plot(sim33[[i]], main = paste0("sim ", i))
}

summary(sim33 ~ edges)

network_summary(sim3)

############################# Fully artificial specifications ########################

n <- 100

n1 <- network.initialize(n,
                         directed = FALSE)

n1 %v% "A" <- sample(c("A", "B", "C", "D"),
                        n, TRUE, prob = c(8, 32, 15, 8))
n1 %v% "B" <- sample(c("A", "B"),
                        n, TRUE, prob = c(6, 3))

form <- n1 ~
  edges +
  nodematch("A") +
  nodefactor("B") +
  gwdegree(0.1, fixed = TRUE) +
  gwesp(0.1, fixed = TRUE) +
  gwdsp(0, fixed = TRUE) 

coefs <- c(
  edges = -6,
  nodematch.A = 3,
  nodefactor.B.B = 0,
  gwdeg.fixed = -2,
  gwesp.fixed = 2,
  gwdsp.fixed = 0.3
)


sim <- simulate(
  form,
  coef = coefs,
  nsim = 20,
  output = "network",
  control = control.simulate.formula(
    MCMC.maxedges = 1000
  )
)

plot(sim[[4]])

network_summary(sim)


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


