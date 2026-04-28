################################################################################
# Degree-sequence MCMC sampler with mixed move sizes
# - fixed size
# - fixed average degree (hence fixed total degree)
# - target Freeman centralisation within tolerance
# - graphicality enforced
# - mixed proposal moves for better exploration
################################################################################

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
  
  # Undirected graph requires even total degree
  if ((total_degree %% 2) != 0) {
    total_degree <- total_degree + 1L
  }
  
  total_degree
}

construct_initial_degseq <- function(size,
                                     total_degree,
                                     min_degree = 1L,
                                     max_tries = 10000) {
  
  n <- as.integer(size)
  
  min_sum <- n * min_degree
  max_sum <- n * (n - 1L)
  
  if (total_degree < min_sum || total_degree > max_sum) {
    stop("total_degree incompatible with size and min_degree.")
  }
  
  for (attempt in seq_len(max_tries)) {
    deg <- rep.int(min_degree, n)
    rem <- total_degree - sum(deg)
    
    while (rem > 0) {
      eligible <- which(deg < (n - 1L))
      if (length(eligible) == 0) break
      
      i <- sample(eligible, 1)
      deg[i] <- deg[i] + 1L
      rem <- rem - 1L
    }
    
    if (rem != 0) next
    
    deg <- sample(deg)
    deg <- sort(deg, decreasing = TRUE)
    
    if (sum(deg) != total_degree) next
    if (min(deg) < min_degree) next
    if (!is_graphical_safe(deg)) next
    
    return(deg)
  }
  
  NULL
}


# Proposal moves


propose_move_1 <- function(deg, min_degree = 1L) {
  # Move 1 degree from one node to another
  n <- length(deg)
  
  donors <- which(deg > min_degree)
  recipients <- which(deg < (n - 1L))
  
  if (length(donors) == 0 || length(recipients) == 0) return(NULL)
  
  i <- sample(donors, 1)
  j <- sample(recipients, 1)
  
  if (i == j) return(NULL)
  
  prop <- deg
  prop[i] <- prop[i] - 1L
  prop[j] <- prop[j] + 1L
  
  sort(prop, decreasing = TRUE)
}

propose_move_2split <- function(deg, min_degree = 1L) {
  # Move 2 degrees from one node to two different nodes
  n <- length(deg)
  
  donors <- which(deg >= (min_degree + 2L))
  recipients <- which(deg < (n - 1L))
  
  if (length(donors) == 0 || length(recipients) < 2) return(NULL)
  
  i <- sample(donors, 1)
  recips <- setdiff(recipients, i)
  
  if (length(recips) < 2) return(NULL)
  
  js <- sample(recips, 2, replace = FALSE)
  
  prop <- deg
  prop[i] <- prop[i] - 2L
  prop[js[1]] <- prop[js[1]] + 1L
  prop[js[2]] <- prop[js[2]] + 1L
  
  sort(prop, decreasing = TRUE)
}

propose_move_2merge <- function(deg, min_degree = 1L) {
  # Move 1 degree each from two nodes to one node
  n <- length(deg)
  
  donors <- which(deg > min_degree)
  recipients <- which(deg <= (n - 3L))  # needs room for +2
  
  if (length(donors) < 2 || length(recipients) == 0) return(NULL)
  
  j <- sample(recipients, 1)
  dons <- setdiff(donors, j)
  
  if (length(dons) < 2) return(NULL)
  
  is <- sample(dons, 2, replace = FALSE)
  
  prop <- deg
  prop[is[1]] <- prop[is[1]] - 1L
  prop[is[2]] <- prop[is[2]] - 1L
  prop[j] <- prop[j] + 2L
  
  sort(prop, decreasing = TRUE)
}

propose_degseq_move_mixed <- function(deg,
                                      min_degree = 1L,
                                      move_probs = c(0.5, 0.25, 0.25)) {
  move_type <- sample(c("move1", "split2", "merge2"), size = 1, prob = move_probs)
  
  if (move_type == "move1") {
    prop <- propose_move_1(deg, min_degree = min_degree)
  } else if (move_type == "split2") {
    prop <- propose_move_2split(deg, min_degree = min_degree)
  } else {
    prop <- propose_move_2merge(deg, min_degree = min_degree)
  }
  
  list(prop = prop, move_type = move_type)
}


# Main sampler


degree_sequence_sample_mcmc <- function(nsim = 20,
                                        size,
                                        average_degree,
                                        freeman_centralisation,
                                        tolerance = 0.01,
                                        min_degree = 1L,
                                        burnin = 5000,
                                        thin = 500,
                                        max_start_tries = 10000,
                                        max_initial_search_steps = 50000,
                                        max_total_steps = 500000,
                                        move_probs = c(0.5, 0.25, 0.25),
                                        unique_sequences = FALSE,
                                        seed = NULL,
                                        verbose = FALSE,
                                        store_trace = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- as.integer(size)
  total_degree <- get_total_degree(n, average_degree)
  
  # Build an initial graphical sequence
  current <- construct_initial_degseq(
    size = n,
    total_degree = total_degree,
    min_degree = min_degree,
    max_tries = max_start_tries
  )
  
  if (is.null(current)) {
    stop("Could not construct an initial graphical degree sequence.")
  }
  
  current_c <- freeman_from_degree(current)
  
  # Search for an initial state inside the target centralisation band
  init_steps <- 0L
  init_accept <- 0L
  
  while (abs(current_c - freeman_centralisation) > tolerance) {
    init_steps <- init_steps + 1L
    
    if (init_steps > max_initial_search_steps) {
      stop("Could not find an initial sequence inside the target centralisation band.")
    }
    
    move <- propose_degseq_move_mixed(
      current,
      min_degree = min_degree,
      move_probs = move_probs
    )
    
    prop <- move$prop
    if (is.null(prop)) next
    
    prop_c <- freeman_from_degree(prop)
    
    # Greedy phase: only accept if closer to target and graphical
    if (abs(prop_c - freeman_centralisation) < abs(current_c - freeman_centralisation) &&
        is_graphical_safe(prop)) {
      current <- prop
      current_c <- prop_c
      init_accept <- init_accept + 1L
    }
  }
  
  if (verbose) {
    message("Initial state found | C = ", round(current_c, 4),
            " | dmax = ", max(current),
            " | init_steps = ", init_steps)
  }
  
  out <- vector("list", nsim)
  seen <- character(0)
  saved <- 0L
  total_steps <- 0L
  accepted_moves <- 0L
  attempted_moves <- 0L
  
  move_attempts <- c(move1 = 0L, split2 = 0L, merge2 = 0L)
  move_accepts  <- c(move1 = 0L, split2 = 0L, merge2 = 0L)
  
  if (store_trace) {
    trace_df <- data.frame(
      step = integer(0),
      centralisation = numeric(0),
      dmax = integer(0),
      accepted = logical(0),
      move_type = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    trace_df <- NULL
  }
  
  while (saved < nsim && total_steps < max_total_steps) {
    total_steps <- total_steps + 1L
    
    move <- propose_degseq_move_mixed(
      current,
      min_degree = min_degree,
      move_probs = move_probs
    )
    
    prop <- move$prop
    move_type <- move$move_type
    move_attempts[move_type] <- move_attempts[move_type] + 1L
    
    accepted <- FALSE
    
    if (!is.null(prop)) {
      attempted_moves <- attempted_moves + 1L
      prop_c <- freeman_from_degree(prop)
      
      # Hard constraints
      if (abs(prop_c - freeman_centralisation) <= tolerance &&
          is_graphical_safe(prop)) {
        current <- prop
        current_c <- prop_c
        accepted_moves <- accepted_moves + 1L
        move_accepts[move_type] <- move_accepts[move_type] + 1L
        accepted <- TRUE
      }
    }
    
    if (store_trace) {
      trace_df[nrow(trace_df) + 1L, ] <- list(
        step = total_steps,
        centralisation = current_c,
        dmax = max(current),
        accepted = accepted,
        move_type = move_type
      )
    }
    
    # Save after burn-in, every thin steps
    if (total_steps > burnin && ((total_steps - burnin) %% thin == 0)) {
      key <- paste(current, collapse = "-")
      
      if (unique_sequences && key %in% seen) next
      
      saved <- saved + 1L
      
      out[[saved]] <- list(
        degree_sequence = current,
        dmax = max(current),
        realised_centralisation = current_c,
        realised_average_degree = mean(current)
      )
      
      if (unique_sequences) {
        seen <- c(seen, key)
      }
      
      if (verbose) {
        message(
          "Saved ", saved, "/", nsim,
          " | step = ", total_steps,
          " | dmax = ", max(current),
          " | C = ", round(current_c, 4)
        )
      }
    }
  }
  
  out <- out[seq_len(saved)]
  
  list(
    degree_sequences = lapply(out, `[[`, "degree_sequence"),
    dmax = sapply(out, `[[`, "dmax"),
    realised_centralisation = sapply(out, `[[`, "realised_centralisation"),
    realised_average_degree = sapply(out, `[[`, "realised_average_degree"),
    n_found = saved,
    nsim_requested = nsim,
    total_degree = total_degree,
    size = n,
    average_degree_target = average_degree,
    freeman_target = freeman_centralisation,
    tolerance = tolerance,
    burnin = burnin,
    thin = thin,
    total_steps = total_steps,
    attempted_moves = attempted_moves,
    accepted_moves = accepted_moves,
    acceptance_rate = if (attempted_moves > 0) accepted_moves / attempted_moves else NA_real_,
    init_steps = init_steps,
    init_accepts = init_accept,
    move_attempts = move_attempts,
    move_accepts = move_accepts,
    trace = trace_df
  )
}

# ------------------------------------------------------------------------------
# Convenience plotting / summary helpers
# ------------------------------------------------------------------------------

summarise_degseq_sample <- function(out) {
  data.frame(
    seq_id = seq_along(out$degree_sequences),
    dmax = out$dmax,
    centralisation = out$realised_centralisation,
    average_degree = out$realised_average_degree
  )
}

plot_degseq_trace <- function(out) {
  if (is.null(out$trace)) stop("No trace stored in output.")
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(1, 2))
  
  plot(
    out$trace$step,
    out$trace$centralisation,
    type = "l",
    xlab = "Step",
    ylab = "Freeman centralisation",
    main = "Trace: centralisation"
  )
  abline(h = out$freeman_target, lty = 2)
  abline(h = out$freeman_target + out$tolerance, lty = 3)
  abline(h = out$freeman_target - out$tolerance, lty = 3)
  
  plot(
    out$trace$step,
    out$trace$dmax,
    type = "l",
    xlab = "Step",
    ylab = "Maximum degree",
    main = "Trace: dmax"
  )
}

# ------------------------------------------------------------------------------
# Example usage
# ------------------------------------------------------------------------------

# Example 1: lower centralisation
out_lowC <- degree_sequence_sample_mcmc(
  nsim = 20,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.15,
  tolerance = 0.01,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

out_lowC$degree_sequences

summarise_degseq_sample(out_lowC)
plot_degseq_trace(out_lowC)

# Turn degree sequences into igraph objects
g_lowC <- lapply(out_lowC$degree_sequences, function(x) {
  igraph::realize_degseq(
    x,
    allowed.edge.types = "simple",
    method = "smallest"
  )
})

# Convert to network objects if needed
net_lowC <- lapply(g_lowC, intergraph::asNetwork)

# Plot them
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 3))
for (i in seq_along(net_lowC)) {
  plot(net_lowC[[i]], main = paste0("Low C ", i))
}
par(oldpar)

# Example 2: higher centralisation
out_highC <- degree_sequence_sample_mcmc(
  nsim = 5,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 50000,
  thin = 5000,
  seed = 123,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

summarise_degseq_sample(out_highC)
plot_degseq_trace(out_highC)

g_highC <- lapply(out_highC$degree_sequences, function(x) {
  igraph::realize_degseq(
    x,
    allowed.edge.types = "simple",
    method = "smallest"
  )
})

net_highC <- lapply(g_highC, intergraph::asNetwork)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 3))
for (i in seq_along(net_highC)) {
  plot(net_highC[[i]], main = paste0("High C ", i))
}
par(oldpar)

# Test high density

# Example 2: higher centralisation
out_highCD <- degree_sequence_sample_mcmc(
  nsim = 5,
  size = 100,
  average_degree = 5,
  freeman_centralisation = 0.7,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 125,
  unique_sequences = FALSE,
  verbose = TRUE,
  store_trace = TRUE
)

summarise_degseq_sample(out_highCD)
plot_degseq_trace(out_highCD)

g_highC <- lapply(out_highCD$degree_sequences, function(x) {
  igraph::realize_degseq(
    x,
    allowed.edge.types = "simple",
    method = "smallest"
  )
})

net_highC <- lapply(g_highC, intergraph::asNetwork)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 3))
for (i in seq_along(net_highC)) {
  plot(net_highC[[i]], main = paste0("High C ", i))
}
par(oldpar)

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

library(ergm)
trial <- simulateNetworks(net_list = net_highC,
                          nfAtt = 0,
                          nmAtt = 0.5,
                          gwdeg = 0.5,
                          gwesp = 0.5,
                          gwdsp = -0.02,
                          nsim = 2)

summary(trial[[1]] ~ edges)

net_highC

par(mfrow = c(4, 5), mar = c(0.2, 0.2, 1, 0.2))

for (i in 1:20) {
  plot(trial[[i]], main = paste0("sim ", i))
}
