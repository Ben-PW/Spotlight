##################################################################################

# Degree sequence v 3
# Changed the sampler so it starts by sampling from admissable dmax values, then
# constructs degree sequences from there. One issue here is this is not proportional
# sampling, so extreme dmax values (both high and low) are going to be overrepresented
# This could be considered a feature though
# Another problem is the way the sampler works inherently smoothes the curve of
# the degree distribution. This isn't necessarily an issue as we can argue that
# ERGM specifications would also have resulted in smoothed curved
# Higher centralisation specifications will necessarily have to avoid smoothed
# distributions as well

freeman_from_degree <- function(deg) {
  n <- length(deg)
  dmax <- max(deg)
  sum(dmax - deg) / ((n - 1) * (n - 2))
}

is_graphical_safe <- function(deg) {
  if ("is_graphical" %in% getNamespaceExports("igraph")) {
    return(igraph::is_graphical(deg, allowed.edge.types = "simple"))
  }
  
  if ("is.graphical.degree.sequence" %in% getNamespaceExports("igraph")) {
    return(igraph::is.graphical.degree.sequence(deg))
  }
  
  stop("No graphical degree-sequence checker found in igraph.")
}

get_total_degree <- function(size, average_degree) {
  total_degree <- round(size * average_degree)
  
  # Undirected simple graph requires even total degree
  if ((total_degree %% 2) != 0) {
    total_degree <- total_degree + 1L
  }
  
  total_degree
}

admissible_dmax_values <- function(size,
                                   average_degree,
                                   freeman_centralisation,
                                   tolerance = 0.01,
                                   min_degree = 1L) {
  
  n <- as.integer(size)
  total_degree <- get_total_degree(n, average_degree)
  
  if (n < 3) stop("size must be at least 3.")
  if (min_degree < 0 || min_degree > (n - 1)) {
    stop("min_degree must be between 0 and n - 1.")
  }
  
  dmax_vals <- seq.int(min_degree, n - 1)
  
  cent_vals <- (n * dmax_vals - total_degree) / ((n - 1) * (n - 2))
  
  in_band <- abs(cent_vals - freeman_centralisation) <= tolerance
  
  # Basic feasibility filters beyond the centralisation band:
  # 1. one node fixed at dmax, others at least min_degree
  #    so total degree must be at least dmax + (n - 1) * min_degree
  lower_feasible <- total_degree >= (dmax_vals + (n - 1) * min_degree)
  
  # 2. if dmax is the maximum degree, all n nodes can contribute at most dmax each
  upper_feasible <- total_degree <= (n * dmax_vals)
  
  keep <- in_band & lower_feasible & upper_feasible
  
  data.frame(
    dmax = dmax_vals[keep],
    implied_centralisation = cent_vals[keep]
  )
}

construct_degseq_given_dmax <- function(size,
                                        total_degree,
                                        dmax,
                                        min_degree = 1L,
                                        max_tries = 5000) {
  
  n <- as.integer(size)
  
  if (dmax < min_degree || dmax > (n - 1)) {
    stop("dmax must be between min_degree and n - 1.")
  }
  
  # Start with node 1 fixed at dmax, all others at min_degree
  base <- c(dmax, rep.int(min_degree, n - 1))
  remaining <- total_degree - sum(base)
  
  if (remaining < 0) {
    stop("Chosen dmax is incompatible with total_degree and min_degree.")
  }
  
  # To preserve dmax as the maximum, other nodes cannot exceed dmax
  other_cap <- dmax - min_degree
  
  if (remaining > (n - 1) * other_cap) {
    stop("Not enough capacity in remaining nodes for this dmax.")
  }
  
  for (attempt in seq_len(max_tries)) {
    deg <- base
    rem <- remaining
    
    # Randomly distribute the remaining degree mass over nodes 2:n
    while (rem > 0) {
      eligible <- which(deg[-1] < dmax) + 1L
      if (length(eligible) == 0) break
      
      i <- sample(eligible, 1)
      deg[i] <- deg[i] + 1L
      rem <- rem - 1L
    }
    
    if (rem != 0) next
    
    deg <- sort(deg, decreasing = TRUE)
    
    # Sanity checks
    if (max(deg) != dmax) next
    if (min(deg) < min_degree) next
    if (sum(deg) != total_degree) next
    
    if (is_graphical_safe(deg)) {
      return(deg)
    }
  }
  
  NULL
}

degree_sequence_sample_dmax <- function(nsim = 20,
                                        size,
                                        average_degree,
                                        freeman_centralisation,
                                        tolerance = 0.01,
                                        min_degree = 1L,
                                        max_attempts_per_sim = 500,
                                        seed = NULL,
                                        unique_sequences = FALSE,
                                        verbose = FALSE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- as.integer(size)
  total_degree <- get_total_degree(n, average_degree)
  
  admissible <- admissible_dmax_values(
    size = n,
    average_degree = average_degree,
    freeman_centralisation = freeman_centralisation,
    tolerance = tolerance,
    min_degree = min_degree
  )
  
  if (nrow(admissible) == 0) {
    stop("No admissible dmax values found for these parameters.")
  }
  
  results <- vector("list", nsim)
  seen <- character(0)
  n_found <- 0L
  total_attempts <- 0L
  
  while (n_found < nsim) {
    sim_attempts <- 0L
    found_this_sim <- FALSE
    
    while (sim_attempts < max_attempts_per_sim && !found_this_sim) {
      sim_attempts <- sim_attempts + 1L
      total_attempts <- total_attempts + 1L
      
      # Sample an admissible dmax
      chosen_row <- sample.int(nrow(admissible), 1)
      chosen_dmax <- admissible$dmax[chosen_row]
      
      deg <- construct_degseq_given_dmax(
        size = n,
        total_degree = total_degree,
        dmax = chosen_dmax,
        min_degree = min_degree
      )
      
      if (is.null(deg)) next
      
      key <- paste(deg, collapse = "-")
      
      if (unique_sequences && key %in% seen) next
      
      n_found <- n_found + 1L
      results[[n_found]] <- list(
        degree_sequence = deg,
        dmax = max(deg),
        realised_centralisation = freeman_from_degree(deg),
        realised_average_degree = mean(deg)
      )
      
      if (unique_sequences) {
        seen <- c(seen, key)
      }
      
      found_this_sim <- TRUE
      
      if (verbose) {
        message(
          "Found ", n_found, "/", nsim,
          " | dmax = ", max(deg),
          " | C = ", round(freeman_from_degree(deg), 4)
        )
      }
    }
    
    if (!found_this_sim) {
      warning(
        "Stopped after failing to find a sequence for one requested draw. ",
        "Try increasing tolerance or max_attempts_per_sim."
      )
      break
    }
  }
  
  results <- results[seq_len(n_found)]
  
  list(
    degree_sequences = lapply(results, `[[`, "degree_sequence"),
    dmax = sapply(results, `[[`, "dmax"),
    realised_centralisation = sapply(results, `[[`, "realised_centralisation"),
    realised_average_degree = sapply(results, `[[`, "realised_average_degree"),
    admissible_dmax = admissible,
    n_found = n_found,
    nsim_requested = nsim,
    total_attempts = total_attempts
  )
}

out <- degree_sequence_sample_dmax(
  nsim = 10,
  size = 500,
  average_degree = 6,
  freeman_centralisation = 0.7,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  unique_sequences = FALSE,
  verbose = TRUE
)

out$dmax
out$realised_centralisation
head(out$admissible_dmax)

out <- degree_sequence_sample_dmax(
  nsim = 10,
  size = 500,
  average_degree = 6,
  freeman_centralisation = 0.2,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  unique_sequences = FALSE,
  verbose = TRUE
)

# New degree sequence/base network generation pipeline #

out <- degree_sequence_sample_dmax(
  nsim = 10,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.2,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  unique_sequences = FALSE,
  verbose = TRUE
)

deg <- lapply(out$degree_sequences,
              function(x){
              igraph::realize_degseq(x, allowed.edge.types = "simple",
                                     method = "smallest")
                }
              )

net <- lapply(deg,
              function(x){
                intergraph::asNetwork(x)
              })

# End of new network generation pipeline #

# Assume:
# net = list of network objects (your basis networks)
# nsim = number of simulations per basis network

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
    
    sim <- stats::simulate(
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

trial <- simulateNetworks(net_list = net,
                          nfAtt = 0.5,
                          nmAtt = 0,
                          gwdeg = -0.5,
                          gwesp = 0.5,
                          gwdsp = 0.5,
                          nsim = 2)

summary(trial[[2]] ~ edges)

for (i in 1:20) {
  plot(trial[[i]], main = paste0("sim ", i))
}

############################## Examples for supervision ############################

# PLOT 1 - rationale for passing degree sequence into ERGM, realistic structures
s100ad3fc0.2 <- degree_sequence_sample_dmax(
  nsim = 10,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.2,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  unique_sequences = FALSE,
  verbose = TRUE
)

s100ad3fc0.2 <- lapply(s100ad3fc0.2$degree_sequences,
              function(x){
                igraph::realize_degseq(x, allowed.edge.types = "simple",
                                       method = "smallest")
              }
)

s100ad3fc0.2 <- lapply(s100ad3fc0.2,
              function(x){
                intergraph::asNetwork(x)
              })

s100ad3fc0.2 <- simulateNetworks(net_list = s100ad3fc0.2,
                          nfAtt = 0.5,
                          nmAtt = 0,
                          gwdeg = -0.5,
                          gwesp = 0.5,
                          gwdsp = 0.5,
                          nsim = 2)

for (i in 1:20) {
  plot(s100ad3fc0.2[[i]], main = paste0("LCLD ", i))
}

# PLOT 2 - density and centralisation parameter control

LCLD100 <- degree_sequence_sample_dmax(
  nsim = 5,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.2,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  unique_sequences = FALSE,
  verbose = TRUE
)

HCLD100 <- degree_sequence_sample_dmax(
  nsim = 5,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.7,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  unique_sequences = FALSE,
  verbose = TRUE
)

LCHD100 <- degree_sequence_sample_dmax(
  nsim = 5,
  size = 100,
  average_degree = 6,
  freeman_centralisation = 0.2,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  unique_sequences = FALSE,
  verbose = TRUE
)

HCHD100 <- degree_sequence_sample_dmax(
  nsim = 5,
  size = 100,
  average_degree = 6,
  freeman_centralisation = 0.7,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123,
  unique_sequences = FALSE,
  verbose = TRUE
)



HCLD100 <- lapply(HCLD100$degree_sequences,
                       function(x){
                         igraph::realize_degseq(x, allowed.edge.types = "simple",
                                                method = "smallest")
                       }
)

LCLD100 <- lapply(LCLD100$degree_sequences,
                  function(x){
                    igraph::realize_degseq(x, allowed.edge.types = "simple",
                                           method = "smallest")
                  }
)

HCHD100 <- lapply(HCHD100$degree_sequences,
                  function(x){
                    igraph::realize_degseq(x, allowed.edge.types = "simple",
                                           method = "smallest")
                  }
)

LCHD100 <- lapply(LCHD100$degree_sequences,
                  function(x){
                    igraph::realize_degseq(x, allowed.edge.types = "simple",
                                           method = "smallest")
                  }
)

HCLD100 <- lapply(HCLD100,
                  function(x){
                    intergraph::asNetwork(x)
                  })
  
LCLD100 <- lapply(LCLD100,
                  function(x){
                    intergraph::asNetwork(x)
                  })

HCHD100 <- lapply(HCHD100,
                  function(x){
                    intergraph::asNetwork(x)
                  })

LCHD100 <- lapply(LCHD100,
                  function(x){
                    intergraph::asNetwork(x)
                  })

HCLD100sim <- simulateNetworks(net_list = HCLD100,
                                 nfAtt = 0.5,
                                 nmAtt = 0,
                                 gwdeg = -0.5,
                                 gwesp = 0.5,
                                 gwdsp = 0.5,
                                 nsim = 2)

LCLD100sim <- simulateNetworks(net_list = LCLD100,
                              nfAtt = 0.5,
                              nmAtt = 0,
                              gwdeg = -0.5,
                              gwesp = 0.5,
                              gwdsp = 0.5,
                              nsim = 2)

HCHD100sim <- simulateNetworks(net_list = HCHD100,
                               nfAtt = 0.5,
                               nmAtt = 0,
                               gwdeg = -0.5,
                               gwesp = 0.5,
                               gwdsp = 0.5,
                               nsim = 2)

LCHD100sim <- simulateNetworks(net_list = LCHD100,
                               nfAtt = 0.5,
                               nmAtt = 0,
                               gwdeg = -0.5,
                               gwesp = 0.5,
                               gwdsp = 0.5,
                               nsim = 2)

for (i in 1:5) plot(HCLD100sim[[i]], main = paste0("HCLD ", i))
for (i in 1:5) plot(LCLD100sim[[i]], main = paste0("LCLD ", i))
for (i in 1:5) plot(HCHD100sim[[i]], main = paste0("HCHD ", i))
for (i in 1:5) plot(LCHD100sim[[i]], main = paste0("LCHD ", i))
