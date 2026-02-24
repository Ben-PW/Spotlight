#####################################################################################################

#Data preprocessing script. Import, pre-process, fit ERGMs and check GoFs here. 

#####################################################################################################

#### IMPORTANT !!! ####
# If you are trying to run the demo section, you must uncomment below

# install.packages("remotes")
# remotes::install_github("schochastics/networkdata")

here::here()

############################# Provisional Irish Republican Army 
######################### Load network #########################

pira_df <- as.matrix(
  read.csv(
    here::here("Data", "PIRA", "60_PERIOD3_NET.csv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

PIRA3 <- network::network(
  pira_df,
  directed = FALSE,
  loops = FALSE,
  matrix.type = "adjacency"
)

######################### Load attributes #########################

pira_att <- read.csv(
  here::here("Data", "PIRA", "60_PERIOD3_ATT.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 99999 assumed to be placeholder for missing data
att_mat <- as.matrix(pira_att[, -1, drop = FALSE])
att_mat[att_mat == 99999] <- "Unknown"
pira_att[, -1] <- att_mat

# Remove Age due to high missingness 
if ("Age at Recruitment" %in% names(pira_att)) {
  pira_att[["Age at Recruitment"]] <- NULL
}

# Attach all attributes to network
for (col in names(pira_att)[-1]) {
  network::set.vertex.attribute(PIRA3, col, pira_att[[col]])
}


######################### Brigade collapse #########################

# Brigade membership stored as multiple binary dummies so collapse to one categorical
brig_vars <- c(
  "Antrim Brigade", "Armagh Brigade", "Derry Brigade", "Down Brigade",
  "Fermanagh Brigade", "Tyrone Brigade"
)

# Pull dummies from attached attributes 
M <- sapply(brig_vars, function(v) {
  as.integer(network::get.vertex.attribute(PIRA3, v) == 1)
})

table(rowSums(M)) # No multibrigade membership

# Unknown brigade if not exactly one brigade flagged
rs <- rowSums(M)
BrigadeUnknown <- as.integer(rs != 1)
Brigade <- ifelse(rs == 1, brig_vars[max.col(M, ties.method = "first")], "Unknown")

# Set attributes
network::set.vertex.attribute(PIRA3, "BrigadeUnknown", BrigadeUnknown)
network::set.vertex.attribute(PIRA3, "Brigade", Brigade)

# Store known brigades
known_brig <- setdiff(unique(Brigade), "Unknown")

table(rowSums(M))

######################### Cleanup #########################

rm(pira_df, pira_att, att_mat, node_ids, brig_vars, M, rs, BrigadeUnknown, Brigade)



#### Specify ERGM ####

#ctrl <- ergm::control.ergm(
#  MCMC.burnin = 2e5,
#  MCMC.interval = 4096,
#  MCMC.samplesize = 12000,
#  MCMLE.maxit = 60
#)

m_test <- ergm::ergm(PIRA3 ~ edges + 
                       gwdegree(0.5, fixed = TRUE) +
                       #gwesp(0.5, fixed = TRUE) +
                       nodefactor("Marital Status") +
                       nodematch("Marital Status", 
                                         levels = c(0,1)) +
                       nodematch("Gender") +
                       nodematch("University") +
                       nodefactor("Brigade") + 
                       nodematch("Brigade", levels = known_brig),
                     #control = ctrl
                     )

summary(m_test)

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


