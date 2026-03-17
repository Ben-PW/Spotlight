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

##### Smith & Papachristos 2016 based network #####

# Approximate terms of criminal network from the paper
# Alt triangles

library(ergm)
library(network)

n <- 1030
n1 <- network.initialize(n,
                         directed = FALSE)

n1 %v% "kingpin" <- as.integer(seq_len(n) == sample(n, 1))

n1 %v% "eth" <- sample(c("English","German","Irish","Italian","Jewish","Other"),
                        n, TRUE, prob=c(.23,.10,.22,.27,.06,.13))

form <- n1 ~ isolates + 
  edges + 
  nodematch("eth") + 
  nodefactor("kingpin") +
  gwesp(0.5, fixed = TRUE) + 
  gwdegree(0.5, fixed = TRUE) +
  gwdsp(0.5, fixed=TRUE)

coefs <- c(
  isolates = -10,
  edges = -5.5,                 
  nodematch.eth = 0.1,
  nodefactor.kingpin = 6,
  gwesp.fixed = 2,
  gwdeg.fixed = 0.5,
  gwdsp.fixed = -0.1
)

sim1 <- simulate(
  form,
  coef = coefs,
  nsim = 20,
  output = "network"
)

plot(sim1[[1]])

network_summary(sim1)

form2 <- n1 ~ #components + 
  edges + 
  nodematch("eth") + 
  nodefactor("kingpin") +
  gwesp(0.7, fixed = TRUE) + 
  gwdegree(0.2, fixed = TRUE) #+
  #gwdsp(0.5, fixed=TRUE)

coefs2 <- c(
  #components = -1,
  edges = -8,                 
  nodematch.eth = 0.2,
  nodefactor.kingpin = 8,
  gwesp.fixed = 4,
  gwdeg.fixed = 3#,
  #gwdsp.fixed = -0.05
)

sim2 <- simulate(
  form2,
  coef = coefs2,
  nsim = 20,
  output = "network"
)

plot(sim2[[1]])

network_summary(sim2)


# ERGM formula
form <- n1 ~ edges +
  nodematch("group") +
  gwesp(0.5, fixed = TRUE) +
  gwdegree(0.5, fixed = TRUE) +
  nodefactor("group")

coefs <- c(
  edges = -7,         
  nodematch.group = 1,    
  gwesp.fixed = 2.6,
  gwdeg.fixed = 2.6,
  nodefactor.group.B = -0.5
)

sim1 <- simulate(
  form,
  coef = coefs,
  nsim = 100,
  output = "network"
)

plot(sim1[[1]])


############################# Grund and Densley 2014 ###############################

n <- 48
n1 <- network.initialize(n,
                         directed = FALSE)

n1 %v% "eth" <- sample(c("D","A","B","C"),
                       n, TRUE, prob=c(24, 11, 10, 6))

form <- n1 ~ edges +
  nodematch("eth") +
  nodefactor("eth", levels = c("A", "B", "C")) +
  gwesp(0, fixed = TRUE)

coefs <- c(
  edges = -4.36,                 
  nodematch.eth = 1.12,
  #nodefactor.eth.D = 0,
  nodefactor.eth.A = 0.5,
  nodefactor.eth.B = 0.49,
  nodefactor.eth.C = 0.61,
  gwesp.fixed = 1.16
)

sim1 <- simulate(
  form,
  coef = coefs,
  nsim = 20,
  output = "network"
)
par(mfrow = c(4, 5), mar = c(0.2, 0.2, 1, 0.2))

for (i in 1:20) {
  plot(sim1[[i]], main = paste0("sim ", i))
}

network_summary(sim1)

##################################### Berlusconi 2021 ################################

# Phase 1 base

# This one does not seem to produce valid networks, likely due to missing multiplex tie 
# information and other attribute related info?

# Yes, ERGM terms compete with each other. By removing some terms and keeping the rest
# constant, it can result in very different looking networks

# This is going to have to resemble Robins et al. 2004 in approach more than I thought

n <- 63

n2 <- network.initialize(n,
                         directed = FALSE)

n2 %v% "role" <- sample(c("A", "B", "C", "D"),
                        n, TRUE, prob = c(8, 32, 15, 8))
n2 %v% "stat" <- sample(c("A", "B"),
                        n, TRUE, prob = c(6, 3))

form <- n2 ~
  edges +
  nodematch("role") +
  nodefactor("stat") +
  gwdegree(0.3, fixed = TRUE) +
  gwesp(0.3, fixed = TRUE) #+
  #gwdsp(0.3, fixed = TRUE) # including gwdsp results in highly centralised
  
coefs <- c(
  edges = -5.5,
  nodematch.role = 0.6,
  nodefactor.stat.B = 0.5,
  gwdeg.fixed = 2.87,
  gwesp.fixed = 1.81 #,
  #gwdsp.fixed = 0.2
)

sim2 <- simulate(
  form,
  coef = coefs,
  nsim = 20,
  output = "network"
)

plot(sim2[[1]])

network_summary(sim2)

################################# Diviak et al., 2019 #############################

n <- 32

n4 <- network.initialize(n,
                         directed = FALSE)

n4 %v% "stat" <- sample(c("A", "B"),
                   n, TRUE, prob = c(1, 1))

form <- n4 ~
  edges +
  gwdegree(0, fixed = TRUE) +
  gwesp(0, fixed = TRUE) +
  gwdsp(0, fixed = TRUE) +
  nodefactor("stat") +
  nodematch("stat")

coefs <- c(
  edges = -3,
  gwdeg.fixed = -1,
  gwesp.fixed = 0.7,
  gwdsp.fixed = 0.1,
  nodefactor.stat.B = 0.8,
  nodematch.stat = -0.5
)

sim4 <- simulate(
  form,
  coef = coefs,
  nsim = 20,
  output = "network"
)

plot(sim4[[1]])

network_summary(sim4)

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


n <- 60

n3 <- network.initialize(n,
                         directed = FALSE)

n3 %v% "role" <- sample(c("A", "B", "C"),
                        n, TRUE, prob = c(8, 8, 8))
n3 %v% "stat" <- sample(c("A", "B"),
                        n, TRUE, prob = c(6, 3))

form <- n3 ~
  edges +
  nodematch("role") +
  nodefactor("stat") +
  gwdegree(0.3, fixed = TRUE) +
  gwesp(0.3, fixed = TRUE) +
  gwdsp(0.3, fixed = TRUE) # w

coefs <- c(
  edges = -5.1,
  nodematch.role = 1.7,
  nodefactor.stat.B = 0,
  gwdeg.fixed = 1.4,
  gwesp.fixed = 1.5,
  gwdsp.fixed = -0.025
)

coefs <- c(
  edges = -4.9,
  nodematch.role = 1.4,
  nodefactor.stat.B = 0.5,
  gwdeg.fixed = 1,
  gwesp.fixed = 1,
  gwdsp.fixed = -0.025
)

sim3 <- simulate(
  form,
  coef = coefs,
  nsim = 20,
  output = "network"
)

simulate.ergm()

cols <- as.factor()
par(mfrow = c(4, 5), mar = c(0.2, 0.2, 1, 0.2))

for (i in 1:20) {
  plot(sim3[[i]], main = paste0("sim ", i))
}

for (i in 1:20) {
  sna::gplot.layout.kamadakawai(sim3[[i]])
}



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
  nodematch.role = 3,
  nodefactor.stat.B = 0,
  gwdeg.fixed = -2,
  gwesp.fixed = 2,
  gwdsp.fixed = 0.3
)

sim <- simulate(
  form,
  coef = coefs,
  nsim = 20,
  output = "network"
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
rm(flo)

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

rm(ac1, ac2)


