#####################################################################################################

#Data preprocessing script. Import, pre-process, fit ERGMs and check GoFs here. 

#####################################################################################################

#### IMPORTANT !!! ####
# If you are trying to run the demo section, you must run the below commented
# out code

# install.packages("remotes")
# remotes::install_github("schochastics/networkdata")

library(statnet)
library(networkdata)
library(intergraph)

########################################## Florentine Families ######################################
# 16x16 undirected, low density, isolates

data(flo)

flo <- network(flo, directed = F)

flo %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)

flo1 <- ergm(flo ~ edges + absdiff("wealth"))
rm(flo)

####################################### Books network ##################################

# For some reason I can't access these datasets any more
# data(books)

###############################  ############################
# 4 x 24 x 24 weighted adjacency matrices

data(covert_17)

ac1 <- intergraph::asNetwork(covert_17)
ac2 <- intergraph::asNetwork(covert_17) # just using the same data for now


plot(ac1)
plot(ac2)



acct1 <- ergm(ac1 ~ edges)
summary(acct1)

acct2 <- ergm(ac2 ~ edges)
summary(acct2)

rm(ac1, ac2)

#################################### Detach data package ###################################

detach(package:networkdata)