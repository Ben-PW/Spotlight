here::here()

pira_df <- as.matrix(
  read.csv(
    here::here("Data", "PIRA", "60_PERIOD3_NET.csv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

range(pira_df) # checking for non binary vals

PIRA3 <- network::network(
  pira_df,
  directed = FALSE,
  loops = FALSE,
  matrix.type = "adjacency"
)

#### Load attributes ####

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

ids_net <- network::get.vertex.attribute(PIRA3, "vertex.names")
ids_att <- as.character(pira_att[[1]])

pira_att2 <- pira_att[match(ids_net, ids_att), , drop=FALSE]
stopifnot(all(as.character(pira_att2[[1]]) == ids_net))

# Collapse Brigade
brig_vars <- c(
  "Antrim Brigade", "Derry Brigade", "Armagh Brigade",
  "Down Brigade", "Tyrone Brigade", "Fermanagh Brigade"
)


is_unknown <- apply(pira_att2[brig_vars] == "Unknown", 1, all)

member_matrix <- pira_att2[brig_vars] == 1
brig_name <- brig_vars[max.col(member_matrix, ties.method = "first")]

pira_att2$Brigade <- ifelse(is_unknown, "Unknown", brig_name)

# Use this to specify parameter levels later
known_brig <- setdiff(unique(pira_att2$Brigade), "Unknown")

# Collapse violence status (multicolinearity issues)

pira_att2$ViolenceStatus <- ifelse(
  pira_att2$Period3Vio == 1, "violent",
  ifelse(pira_att2$Period3NonVio == 1, "nonviolent", "other")
)

# Attach all attributes to network
for (col in names(pira_att2)[-1]) {
  network::set.vertex.attribute(PIRA3, col, pira_att2[[col]])
}

par(mfrow = c(1,1))
plot(PIRA3) # highly disconnected, coerce to LCC

# Coerce to LCC
comp <- sna::component.dist(PIRA3)
largest_comp_id <- which.max(table(comp$membership))
keep_nodes <- which(comp$membership == largest_comp_id)
PIRA3_LCC <- network::get.inducedSubgraph(PIRA3, keep_nodes)

plot(PIRA3_LCC)

# Significant core periphery structure
deg <- sna::degree(PIRA3_LCC, gmode = "graph")
hist(deg, breaks = 49, freq = TRUE)
?hist
