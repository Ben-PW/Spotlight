#####################################################################################################

#Data preprocessing script. Import, pre-process, fit ERGMs and check GoFs here. 

#####################################################################################################

#### IMPORTANT !!! ####
# If you are trying to run the demo section, you must uncomment below

# install.packages("remotes")
# remotes::install_github("schochastics/networkdata")

here::here()

############################# Provisional Irish Republican Army 

###################################### PIRA 1 #####################################

#### Trying with alternative attribute sets ####

pira1_df <- as.matrix(
  read.csv(
    here::here("Data", "PIRA", "60_PERIOD1_NET.csv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

PIRA1 <- network::network(
  pira1_df,
  directed = FALSE,
  loops = FALSE,
  matrix.type = "adjacency"
)

#### Load attributes ####

pira_att <- read.csv(
  here::here("Data", "PIRA", "60_PERIOD1_ATTT.csv"),
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
  network::set.vertex.attribute(PIRA1, col, pira_att[[col]])
}


#### Brigade collapse ####

# Brigade membership stored as multiple binary dummies so collapse to one categorical
brig_vars <- c(
  "Antrim Brigade", "Armagh Brigade", "Derry Brigade", "Down Brigade",
  "Fermanagh Brigade", "Tyrone Brigade"
)

# Pull dummies from attached attributes 
M <- sapply(brig_vars, function(v) {
  as.integer(network::get.vertex.attribute(PIRA1, v) == 1)
})

table(rowSums(M)) # No multibrigade membership

# Unknown brigade if not exactly one brigade flagged
rs <- rowSums(M)
BrigadeUnknown <- as.integer(rs != 1)
Brigade <- ifelse(rs == 1, brig_vars[max.col(M, ties.method = "first")], "Unknown")

# Set attributes
network::set.vertex.attribute(PIRA1, "BrigadeUnknown", BrigadeUnknown)
network::set.vertex.attribute(PIRA1, "Brigade", Brigade)

# Store known brigades
known_brig <- setdiff(unique(Brigade), "Unknown")

table(rowSums(M))

# Cleanup 

rm(pira_df, 
   pira_att, 
   att_mat, 
   node_ids, 
   brig_vars, 
   M, 
   rs, 
   BrigadeUnknown, 
   Brigade)



#### Specify ERGM ####
par(mfrow = c(1,1))
plot(PIRA1) # highly disconnected, coerce to LCC

######################################### PIRA 2 ###################################

pira2_df <- as.matrix(
  read.csv(
    here::here("Data", "PIRA", "60_PERIOD2_NET.csv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

PIRA2 <- network::network(
  pira2_df,
  directed = FALSE,
  loops = FALSE,
  matrix.type = "adjacency"
)

#### Load attributes ####

pira_att <- read.csv(
  here::here("Data", "PIRA", "60_PERIOD2_ATT.csv"),
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
  network::set.vertex.attribute(PIRA2, col, pira_att[[col]])
}


#### Brigade collapse ####

# Brigade membership stored as multiple binary dummies so collapse to one categorical
brig_vars <- c(
  "Antrim Brigade", "Armagh Brigade", "Derry Brigade", "Down Brigade",
  "Fermanagh Brigade", "Tyrone Brigade"
)

# Pull dummies from attached attributes 
M <- sapply(brig_vars, function(v) {
  as.integer(network::get.vertex.attribute(PIRA2, v) == 1)
})

table(rowSums(M)) # No multibrigade membership

# Unknown brigade if not exactly one brigade flagged
rs <- rowSums(M)
BrigadeUnknown <- as.integer(rs != 1)
Brigade <- ifelse(rs == 1, brig_vars[max.col(M, ties.method = "first")], "Unknown")

# Set attributes
network::set.vertex.attribute(PIRA2, "BrigadeUnknown", BrigadeUnknown)
network::set.vertex.attribute(PIRA2, "Brigade", Brigade)

# Store known brigades
known_brig <- setdiff(unique(Brigade), "Unknown")

table(rowSums(M))

# Cleanup 

rm(pira_df, 
   pira_att, 
   att_mat, 
   node_ids, 
   brig_vars, 
   M, 
   rs, 
   BrigadeUnknown, 
   Brigade)



#### Specify ERGM ####
par(mfrow = c(1,1))
plot(PIRA2) # highly disconnected, coerce to LCC

####################################### PIRA 5 ####################################

pira5_df <- as.matrix(
  read.csv(
    here::here("Data", "PIRA", "60_PERIOD4_5_NET.csv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

PIRA5 <- network::network(
  pira5_df,
  directed = FALSE,
  loops = FALSE,
  matrix.type = "adjacency"
)

#### Load attributes ####

pira_att <- read.csv(
  here::here("Data", "PIRA", "60_PERIOD4_5_ATT.csv"),
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
  network::set.vertex.attribute(PIRA5, col, pira_att[[col]])
}


#### Brigade collapse ####

# Brigade membership stored as multiple binary dummies so collapse to one categorical
brig_vars <- c(
  "Antrim Brigade", "Armagh Brigade", "Derry Brigade", "Down Brigade",
  "Fermanagh Brigade", "Tyrone Brigade"
)

# Pull dummies from attached attributes 
M <- sapply(brig_vars, function(v) {
  as.integer(network::get.vertex.attribute(PIRA5, v) == 1)
})

table(rowSums(M)) # No multibrigade membership

# Unknown brigade if not exactly one brigade flagged
rs <- rowSums(M)
BrigadeUnknown <- as.integer(rs != 1)
Brigade <- ifelse(rs == 1, brig_vars[max.col(M, ties.method = "first")], "Unknown")

# Set attributes
network::set.vertex.attribute(PIRA5, "BrigadeUnknown", BrigadeUnknown)
network::set.vertex.attribute(PIRA5, "Brigade", Brigade)

# Store known brigades
known_brig <- setdiff(unique(Brigade), "Unknown")

table(rowSums(M))

# Cleanup 

rm(pira_df, 
   pira_att, 
   att_mat, 
   node_ids, 
   brig_vars, 
   M, 
   rs, 
   BrigadeUnknown, 
   Brigade)



#### Specify ERGM ####
par(mfrow = c(1,1))
plot(PIRA5) # highly disconnected, coerce to LCC



################################# PIRA3 ####################################

source("PIRA3_process.R")

#### Specify ERGM ####
par(mfrow = c(1,1))
plot(PIRA3) # highly disconnected, coerce to LCC

# Coerce to LCC
comp <- sna::component.dist(PIRA3)
largest_comp_id <- which.max(table(comp$membership))
keep_nodes <- which(comp$membership == largest_comp_id)
PIRA3 <- network::get.inducedSubgraph(PIRA3, keep_nodes)

plot(PIRA3)

# ERGM with missing data included

ctrl <- ergm::control.ergm(
  seed = 1234
)

# IMPORTANT FOR BELOW !!
# The following linear dependence has been detected among the model statistics:
#edges = nodematch.Gender + `nodefactor.Brigade.Down Brigade`
#edges = nodefactor.Period3NonVio.1 + nodematch.Period3NonVio
#2 * edges + `nodematch.Marital Status` + 2 * nodefactor.Period3NonVio.1 = nodematch.Gender + nodematch.University + nodefactor.Period3VForOp.1
#nodematch.Gender + nodematch.University = edges + `nodematch.Marital Status` + nodematch.Period3VForOp
#`nodematch.Marital Status` + nodematch.Gender + nodematch.University + 1/2 * nodefactor.Period3Vio.1 + 1/2 * nodematch.Period3Gun = 2 * edges + 1/2 * nodematch.Period3Vio + nodefactor.Period3NonVio.1 + 1/2 * nodefactor.Period3Gun.1 + nodefactor.Period3IED_C.1
#3 * edges + 1/2 * nodematch.Period3Vio + nodefactor.Period3NonVio.1 + 1/2 * nodefactor.Period3Gun.1 = `nodematch.Marital Status` + nodematch.Gender + nodematch.University + 1/2 * nodefactor.Period3Vio.1 + 1/2 * nodematch.Period3Gun + nodematch.Period3IED_C
#2 * edges + 1/2  [... truncated]

pira3 <- ergm::ergm(PIRA3 ~ edges + 
                       #gwdegree(0.5, fixed = TRUE) +
                       gwesp(0.5, fixed = TRUE) +
                       nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                  levels = 1) +
                       nodematch("Marital Status", 
                                         levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                       nodematch("Gender") +
                       nodematch("University") +
                       nodefactor("Brigade", 
                                  levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                       nodematch("Brigade", 
                                 levels = known_brig) + # match only between known brigades
                       nodefactor("Period3Vio") +
                       nodematch("Period3Vio") +
                       nodefactor("Period3NonVio") +
                       nodematch("Period3NonVio") +
                       nodefactor("Period3VForOp") +
                       nodematch("Period3VForOp") +
                       nodefactor("Period3Senior") +
                       #nodematch("Period3Senior") + excluded as linear combination with nodefactor. Nodefactor had higher coefficient in base study so was retained
                       nodefactor("Period3Gun") +
                       nodematch("Period3Gun") +
                       nodefactor("Period3IED_C") +
                       nodematch("Period3IED_C") +
                       nodefactor("Period3IED_P") +
                       nodematch("Period3IED_P") +
                       nodefactor("Period3ForOp") +
                       nodematch("Period3ForOp") +
                       nodefactor("Period3Rob") +
                       nodematch("Period3Rob"),
                     control = ctrl
                     )

summary(pira3)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3)# looks ok, decent chain mixing, some distributions have heavier left tail
pira3gof <- ergm::gof(pira3, GOF = ~ distance + espartners + triadcensus + degree)
pira3gof
par(mfrow = c(2,2))
plot(pira3gof, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)

rm(ctrl)

pira3.1 <- ergm::ergm(PIRA3 ~ edges + 
                      gwdegree(0.3, fixed = TRUE) +
                      gwesp(0.3, fixed = TRUE) +
                      #nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                      #           levels = 1) +
                      nodematch("Marital Status", 
                                levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                      nodematch("Gender") +
                      nodematch("University") +
                      #nodefactor("Brigade", 
                      #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                      nodematch("Brigade") + # match only between known brigades
                      #nodefactor("Period3Vio") +
                      nodematch("Period3Vio") +
                      #nodefactor("Period3NonVio") +
                      nodematch("Period3NonVio") +
                      #nodefactor("Period3VForOp") +
                      nodematch("Period3VForOp") +
                      #nodefactor("Period3Senior") +
                      #nodematch("Period3Senior") + excluded as linear combination with nodefactor. Nodefactor had higher coefficient in base study so was retained
                      #nodefactor("Period3Gun") +
                      nodematch("Period3Gun") +
                      #nodefactor("Period3IED_C") +
                      nodematch("Period3IED_C") +
                      #nodefactor("Period3IED_P") +
                      nodematch("Period3IED_P") +
                      #nodefactor("Period3ForOp") +
                      nodematch("Period3ForOp") +
                      #nodefactor("Period3Rob") +
                      nodematch("Period3Rob"),
                    control = ctrl
)

summary(pira3.1)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.1) 
pira3.1gof <- ergm::gof(pira3.1, GOF = ~ distance + espartners + triadcensus + degree)
pira3.1gof
par(mfrow = c(2,2))
plot(pira3.1gof, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)


pira3.1 <- ergm::ergm(PIRA3 ~ edges + 
                        gwdegree(0.3, fixed = TRUE) +
                        gwesp(0.3, fixed = TRUE) +
                        nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                   levels = 1) +
                        nodematch("Marital Status", 
                                  levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                        nodematch("Gender") +
                        nodefactor("Gender") +
                        nodematch("University") +
                        nodefactor("Brigade", 
                                  levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                        nodematch("Brigade") + 
                        nodefactor("ViolenceStatus") +
                        nodematch("ViolenceStatus"),
                      control = ctrl
)

pira3.1_LCC <- ergm::ergm(PIRA3_LCC ~ edges + 
                        gwdegree(0.3, fixed = TRUE) +
                        gwesp(0.3, fixed = TRUE) +
                        nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                   levels = 1) +
                        nodematch("Marital Status", 
                                  levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                        nodematch("Gender") +
                        nodefactor("Gender") +
                        nodematch("University") +
                        nodefactor("Brigade", 
                                   levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                        nodematch("Brigade") + 
                        nodefactor("ViolenceStatus") +
                        nodematch("ViolenceStatus"),
                      control = ctrl
)

summary(pira3.1_LCC)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.1_LCC) 
pira3.1gof_LCC <- ergm::gof(pira3.1_LCC, GOF = ~ distance + espartners + triadcensus + degree)
pira3.1gof_LCC
par(mfrow = c(2,2))
plot(pira3.1gof_LCC, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)


pira3.2_LCC <- ergm::ergm(PIRA3_LCC ~ edges + 
                            gwdegree(0.4, fixed = TRUE) +
                            gwesp(0.3, fixed = TRUE) +
                            nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                       levels = 1) +
                            nodematch("Marital Status", 
                                      levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                            nodematch("Gender") +
                            nodefactor("Gender") +
                            #nodematch("University") +
                            #nodefactor("Brigade", 
                            #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                            nodematch("Brigade",
                                      levels = known_brig) + 
                            nodefactor("ViolenceStatus"),
                            #nodematch("ViolenceStatus"),
                          control = ctrl
)

summary(pira3.2_LCC)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.2_LCC) 
pira3.2gof_LCC <- ergm::gof(pira3.2_LCC, GOF = ~ distance + espartners + triadcensus + degree)
pira3.2gof_LCC
par(mfrow = c(2,2))
plot(pira3.2gof_LCC, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)

pira3.2_LCC <- ergm::ergm(PIRA3_LCC ~ edges + 
                            gwdegree(0.3, fixed = TRUE) +
                            gwesp(0.35, fixed = TRUE) +
                            nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                       levels = 1) +
                            nodematch("Marital Status", 
                                      levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                            nodematch("Gender") +
                            nodefactor("Gender") +
                            #nodematch("University") +
                            #nodefactor("Brigade", 
                            #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                            nodematch("Brigade",
                                      levels = known_brig) + 
                            nodefactor("ViolenceStatus"),
                          #nodematch("ViolenceStatus"),
                          control = ctrl
)

summary(pira3.2_LCC)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.2_LCC) 
pira3.2gof_LCC <- ergm::gof(pira3.2_LCC, GOF = ~ distance + espartners + triadcensus + degree)
pira3.2gof_LCC
par(mfrow = c(2,2))
plot(pira3.2gof_LCC, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)

# Trying with core periphery coding
pira3.2_LCC <- ergm::ergm(PIRA3_LCC ~ edges + 
                            gwdegree(0.3, fixed = TRUE) +
                            gwesp(0.3, fixed = TRUE) +
                            nodefactor("Core") +
                            nodematch("Core") +
                            nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                       levels = 1) +
                            nodematch("Marital Status", 
                                      levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                            nodematch("Gender") +
                            nodefactor("Gender") +
                            #nodematch("University") +
                            #nodefactor("Brigade", 
                            #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                            nodematch("Brigade",
                                      levels = known_brig) + 
                            nodefactor("ViolenceStatus"),
                          #nodematch("ViolenceStatus"),
                          control = ctrl
)

summary(pira3.2_LCC)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.2_LCC) 
pira3.2gof_LCC <- ergm::gof(pira3.2_LCC, GOF = ~ distance + espartners + triadcensus + degree)
pira3.2gof_LCC
par(mfrow = c(2,2))
plot(pira3.2gof_LCC, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)

# Trying with core periphery coding
pira3.2_LCC <- ergm::ergm(PIRA3_LCC ~ edges + 
                            gwdegree(0.3, fixed = TRUE) +
                            gwesp(0.4, fixed = TRUE) +
                            nodefactor("Core") +
                            nodematch("Core") +
                            nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                       levels = 1) +
                            nodematch("Marital Status", 
                                      levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                            nodematch("Gender") +
                            nodefactor("Gender") +
                            #nodematch("University") +
                            #nodefactor("Brigade", 
                            #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                            nodematch("Brigade",
                                      levels = known_brig) + 
                            nodefactor("ViolenceStatus"),
                          #nodematch("ViolenceStatus"),
                          control = ctrl
)

summary(pira3.2_LCC)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.2_LCC) 
pira3.2gof_LCC <- ergm::gof(pira3.2_LCC, GOF = ~ distance + espartners + triadcensus + degree)
pira3.2gof_LCC
par(mfrow = c(2,2))
plot(pira3.2gof_LCC, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)

# Trying with core periphery coding
pira3.3_LCC <- ergm::ergm(PIRA3_LCC ~ edges + 
                            gwdegree(0.4, fixed = TRUE) +
                            gwesp(0.45, fixed = TRUE) +
                            nodefactor("Core") +
                            nodematch("Core") +
                            nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                       levels = 1) +
                            nodematch("Marital Status", 
                                      levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                            nodematch("Gender") +
                            nodefactor("Gender") +
                            #nodematch("University") +
                            #nodefactor("Brigade", 
                            #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                            nodematch("Brigade",
                                      levels = known_brig),
                            #nodefactor("ViolenceStatus"),
                          #nodematch("ViolenceStatus"),
                          control = ctrl
)

summary(pira3.3_LCC)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.3_LCC) 
pira3.3gof_LCC <- ergm::gof(pira3.3_LCC, GOF = ~ distance + espartners + triadcensus + degree )
pira3.3gof_LCC
par(mfrow = c(2,2))
plot(pira3.3gof_LCC, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)

# This one might actually be good, chain mixing might not be ideal though
pira3.4_LCC <- ergm::ergm(PIRA3_LCC ~ edges + 
                            #dsp(0) +
                            gwdegree(0.5, fixed = TRUE) +
                            gwesp(0.5, fixed = TRUE) +
                            nodefactor("Core") +
                            nodematch("Core") +
                            #nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                            #           levels = 1) +
                            nodematch("Marital Status", 
                                      levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                            nodematch("Gender") +
                            nodefactor("Gender") +
                            #nodematch("University") +
                            #nodefactor("Brigade", 
                            #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                            nodematch("Brigade",
                                      levels = known_brig),
                          #nodefactor("ViolenceStatus"),
                          #nodematch("ViolenceStatus"),
                          control = ctrl
)

summary(pira3.4_LCC)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.4_LCC) 
pira3.4gof_LCC <- ergm::gof(pira3.4_LCC, GOF = ~ distance + espartners + triadcensus + degree )
pira3.4gof_LCC
par(mfrow = c(2,2))
plot(pira3.4gof_LCC, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)

ctrl <- ergm::control.ergm(
  seed = 1234,
  MCMLE.maxit = 200,
  MCMC.burnin = 1000000
)

pira3.5_LCC <- ergm::ergm(PIRA3_LCC ~ edges + 
                            #dsp(0) +
                            gwdegree(0.5, fixed = TRUE) +
                            gwesp(0.4, fixed = TRUE ) +
                            nodefactor("Core") +
                            nodematch("Core") +
                            nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                      levels = 1) +
                            nodematch("Marital Status", 
                                      levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                            nodematch("Gender") +
                            nodefactor("Gender") +
                            #nodematch("University") +
                            nodefactor("Brigade", 
                                       levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                            nodematch("Brigade",
                                      levels = known_brig) +
                          #nodefactor("ViolenceStatus"),
                          nodematch("ViolenceStatus"),
                          control = ctrl
)

summary(pira3.5_LCC)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.5_LCC) 
pira3.5gof_LCC <- ergm::gof(pira3.5_LCC, GOF = ~ distance + espartners + triadcensus + degree )
pira3.5gof_LCC
par(mfrow = c(2,2))
plot(pira3.5gof_LCC, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)

ctrl <- ergm::control.ergm(
  seed = 1234,
  MCMLE.maxit = 200,
  MCMC.burnin = 1000000
)

pira3.6_LCC <- ergm::ergm(PIRA3_LCC ~ edges + 
                            #dsp(0) +
                            gwdegree(0.3, fixed = TRUE) +
                            gwesp(0.3, fixed = TRUE ) +
                            gwdsp(0.3, fixed = TRUE) +
                            nodefactor("Core") +
                            nodematch("Core") +
                            nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                                       levels = 1) +
                            nodematch("Marital Status", 
                                      levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                            nodematch("Gender") +
                            nodefactor("Gender") +
                            #nodematch("University") +
                            #nodefactor("Brigade", 
                            #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                            nodematch("Brigade",
                                      levels = known_brig) +
                            #nodemix("Brigade") + 
                            #nodefactor("ViolenceStatus"),
                            nodematch("ViolenceStatus"),
                          control = ctrl
)

summary(pira3.5_LCC)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira3.5_LCC) 
pira3.5gof_LCC <- ergm::gof(pira3.5_LCC, GOF = ~ distance + espartners + triadcensus + degree )
pira3.5gof_LCC
par(mfrow = c(2,2))
plot(pira3.5gof_LCC, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)


ctrl <- ergm::control.ergm(
  seed = 1234
)

pira4 <- ergm::ergm(PIRA3 ~ edges + 
                      gwdegree(0.2, fixed = TRUE) +
                      gwesp(0.2, fixed = TRUE) +
                      #nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                      #           levels = 1) +
                      nodematch("Marital Status", 
                                levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                      nodematch("Gender") +
                      nodematch("University") +
                      #nodefactor("Brigade", 
                      #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                      nodematch("Brigade", 
                                levels = known_brig) + # match only between known brigades
                      nodefactor("Period3Vio") +
                      nodematch("Period3Vio") +
                      nodefactor("Period3NonVio") +
                      nodematch("Period3NonVio") +
                      nodefactor("Period3VForOp") +
                      nodematch("Period3VForOp") +
                      nodefactor("Period3Senior") +
                      #nodematch("Period3Senior") + excluded as linear combination with nodefactor. Nodefactor had higher coefficient in base study so was retained
                      nodefactor("Period3Gun") +
                      nodematch("Period3Gun") +
                      nodefactor("Period3IED_C") +
                      nodematch("Period3IED_C") +
                      nodefactor("Period3IED_P") +
                      nodematch("Period3IED_P") +
                      nodefactor("Period3ForOp") +
                      nodematch("Period3ForOp") +
                      nodefactor("Period3Rob") +
                      nodematch("Period3Rob"),
                    control = ctrl
)

summary(pira4)
par(mar = c(2, 2, 2, 2))
ergm::mcmc.diagnostics(pira4) # decent chain mixing, some distributions have heavier right tail

rm(ctrl)

ergm::gof(pira3, GOF = ~ model)
pira3gof <- ergm::gof(pira3, GOF = ~ distance + espartners + triadcensus + degree)
pira3gof
par(mfrow = c(2,2))
plot(pira3gof, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)
plot(pira3gof) # bad

ergm::gof(pira4, GOF = ~ model)
pira4gof <- ergm::gof(pira4, GOF = ~ distance + espartners + triadcensus + degree)
pira4gof
par(mfrow = c(2,2))
plot(pira4gof, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)
plot(pira4gof) # bad

################################ PIRA 4 model rebuild ##############################

pira4 <- ergm::ergm(PIRA3 ~ edges + 
                      altkstar(1, fixed = TRUE) +
                      gwesp(0.2, fixed = TRUE),
                      #nodefactor("Marital Status", # estimating coefficient for married vs not married + unknown
                      #           levels = 1) +
                      #nodematch("Marital Status", 
                      #          levels = c(0,1)) + # estimating only coefficients for married and unmarried, exclude unknown
                      #nodematch("Gender") +
                      #nodematch("University") +
                      #nodefactor("Brigade", 
                      #           levels = known_brig) + #estimate coefficients for real brigades and leave unknowns as reference
                      #nodematch("Brigade", 
                      #          levels = known_brig) + # match only between known brigades
                      #nodefactor("Period3Vio") +
                      #nodematch("Period3Vio") +
                      #nodefactor("Period3NonVio") +
                      #nodematch("Period3NonVio") +
                      #nodefactor("Period3VForOp") +
                      #nodematch("Period3VForOp") +
                      #nodefactor("Period3Senior") +
                      #nodematch("Period3Senior") + excluded as linear combination with nodefactor. Nodefactor had higher coefficient in base study so was retained
                      #nodefactor("Period3Gun") +
                      #nodematch("Period3Gun") +
                      #nodefactor("Period3IED_C") +
                      #nodematch("Period3IED_C") +
                      #nodefactor("Period3IED_P") +
                      #nodematch("Period3IED_P") +
                      #nodefactor("Period3ForOp") +
                      #nodematch("Period3ForOp") +
                      #nodefactor("Period3Rob") +
                      #nodematch("Period3Rob"),
                    control = ctrl
)

ergm::gof(pira4, GOF = ~ model)
pira4gof <- ergm::gof(pira4, GOF = ~ distance + espartners + triadcensus + degree)
pira4gof
par(mfrow = c(2,2))
plot(pira4gof, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)
plot(pira4gof) # bad

pira4.1 <- ergm::ergm(PIRA3 ~ edges + # pass
                      gwdegree(0.1, fixed = TRUE) +
                      gwesp(0.1, fixed = TRUE),
                    control = ctrl
)

pira4.2 <- ergm::ergm(PIRA3 ~ edges + # pass
                        gwdegree(0.2, fixed = TRUE) +
                        gwesp(0.2, fixed = TRUE),
                      control = ctrl
)

pira4.3 <- ergm::ergm(PIRA3 ~ edges + # fail
                        gwdegree(0.3, fixed = TRUE) +
                        gwesp(0.3, fixed = TRUE),
                      control = ctrl
)

pira4.4 <- ergm::ergm(PIRA3 ~ edges + # fail
                        gwdegree(0.4, fixed = TRUE) +
                        gwesp(0.2, fixed = TRUE),
                      control = ctrl
)

pira4.4 <- ergm::ergm(PIRA3 ~ edges + # fail
                        gwdegree(0.2, fixed = TRUE) +
                        gwesp(0.4, fixed = TRUE),
                      control = ctrl
)

pira4.5 <- ergm::ergm(PIRA3 ~ edges + # fail
                        gwdegree(0.25, fixed = TRUE) +
                        gwesp(0.25, fixed = TRUE),
                      control = ctrl
)

pira4.6 <- ergm::ergm(PIRA3 ~ edges + #fail
                        degree(1) +
                        degree(2) +
                        gwdegree(1.2, fixed = TRUE) +
                        gwesp(2, fixed = TRUE),
                      control = ctrl
)

pira4.7 <- ergm::ergm(PIRA3 ~ edges + # fail
                        gwdegree(fixed = FALSE, cutoff = 200) +
                        gwesp(0.2, fixed = TRUE),
                      control = ctrl
)

ergm::gof(pira4, GOF = ~ model)
pira4gof <- ergm::gof(pira4, GOF = ~ distance + espartners + triadcensus + degree)
pira4gof
par(mfrow = c(2,2))
plot(pira4gof, cex.lab=1.6, cex.axis=1.6, 
     plotlogodds = TRUE)
plot(pira4gof) # bad


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


