################################################################################
# Degree-sequence MCMC sampler with mixed move sizes
# - fixed size
# - target average degree WITH optional tolerance
# - target Freeman centralisation within tolerance
# - minimum degree enforced
# - graphicality enforced
# - mixed proposal moves, including total-degree-changing moves
################################################################################

############################## Utility functions for testing ###########################

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

#################################### MCMC Chain ###############################################

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
  
  if ((total_degree %% 2) != 0) {
    total_degree <- total_degree + 1L
  }
  
  total_degree
}

within_avg_degree_band <- function(deg,
                                   average_degree,
                                   average_degree_tolerance = 0) {
  abs(mean(deg) - average_degree) <= average_degree_tolerance
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

################################################################################
# Proposal moves
################################################################################

propose_move_1 <- function(deg, min_degree = 1L) {
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
  n <- length(deg)
  
  donors <- which(deg > min_degree)
  recipients <- which(deg <= (n - 3L))
  
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

propose_add_2 <- function(deg) {
  n <- length(deg)
  
  recipients <- which(deg <= (n - 3L))
  if (length(recipients) == 0) return(NULL)
  
  j <- sample(recipients, 1)
  
  prop <- deg
  prop[j] <- prop[j] + 2L
  
  if (prop[j] > (n - 1L)) return(NULL)
  
  sort(prop, decreasing = TRUE)
}

propose_drop_2 <- function(deg, min_degree = 1L) {
  donors <- which(deg >= (min_degree + 2L))
  if (length(donors) == 0) return(NULL)
  
  i <- sample(donors, 1)
  
  prop <- deg
  prop[i] <- prop[i] - 2L
  
  if (min(prop) < min_degree) return(NULL)
  
  sort(prop, decreasing = TRUE)
}

propose_degseq_move_mixed <- function(deg,
                                      min_degree = 1L,
                                      move_probs = c(
                                        move1 = 0.40,
                                        split2 = 0.20,
                                        merge2 = 0.20,
                                        add2 = 0.10,
                                        drop2 = 0.10
                                      )) {
  
  move_type <- sample(names(move_probs), size = 1, prob = move_probs)
  
  if (move_type == "move1") {
    prop <- propose_move_1(deg, min_degree = min_degree)
  } else if (move_type == "split2") {
    prop <- propose_move_2split(deg, min_degree = min_degree)
  } else if (move_type == "merge2") {
    prop <- propose_move_2merge(deg, min_degree = min_degree)
  } else if (move_type == "add2") {
    prop <- propose_add_2(deg)
  } else if (move_type == "drop2") {
    prop <- propose_drop_2(deg, min_degree = min_degree)
  } else {
    stop("Unknown move type.")
  }
  
  list(prop = prop, move_type = move_type)
}

################################################################################
# Main sampler
################################################################################

degree_sequence_sample_mcmc <- function(nsim = 20,
                                        size,
                                        average_degree,
                                        average_degree_tolerance = 0,
                                        freeman_centralisation,
                                        tolerance = 0.01,
                                        min_degree = 1L,
                                        burnin = 5000,
                                        thin = 500,
                                        max_start_tries = 10000,
                                        max_initial_search_steps = 50000,
                                        max_total_steps = 500000,
                                        move_probs = c(
                                          move1 = 0.40,
                                          split2 = 0.20,
                                          merge2 = 0.20,
                                          add2 = 0.10,
                                          drop2 = 0.10
                                        ),
                                        unique_sequences = FALSE,
                                        seed = NULL,
                                        verbose = FALSE,
                                        store_trace = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- as.integer(size)
  total_degree <- get_total_degree(n, average_degree)
  
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
  
  init_steps <- 0L
  init_accept <- 0L
  
  while (
    abs(current_c - freeman_centralisation) > tolerance ||
    !within_avg_degree_band(current, average_degree, average_degree_tolerance)
  ) {
    init_steps <- init_steps + 1L
    
    if (init_steps > max_initial_search_steps) {
      stop("Could not find an initial sequence inside the target bands.")
    }
    
    move <- propose_degseq_move_mixed(
      current,
      min_degree = min_degree,
      move_probs = move_probs
    )
    
    prop <- move$prop
    
    if (is.null(prop)) next
    
    prop_c <- freeman_from_degree(prop)
    
    current_dist <- abs(current_c - freeman_centralisation) +
      abs(mean(current) - average_degree)
    
    prop_dist <- abs(prop_c - freeman_centralisation) +
      abs(mean(prop) - average_degree)
    
    if (
      prop_dist < current_dist &&
      within_avg_degree_band(prop, average_degree, average_degree_tolerance) &&
      is_graphical_safe(prop)
    ) {
      current <- prop
      current_c <- prop_c
      init_accept <- init_accept + 1L
    }
  }
  
  if (verbose) {
    message(
      "Initial state found | C = ", round(current_c, 4),
      " | avg degree = ", round(mean(current), 4),
      " | dmax = ", max(current),
      " | init_steps = ", init_steps
    )
  }
  
  out <- vector("list", nsim)
  seen <- character(0)
  saved <- 0L
  total_steps <- 0L
  accepted_moves <- 0L
  attempted_moves <- 0L
  seen_dmax <- integer(0)
  
  move_attempts <- setNames(rep(0L, length(move_probs)), names(move_probs))
  move_accepts <- setNames(rep(0L, length(move_probs)), names(move_probs))
  
  if (store_trace) {
    trace_df <- data.frame(
      step = integer(0),
      centralisation = numeric(0),
      average_degree = numeric(0),
      total_degree = integer(0),
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
      prop_avg <- mean(prop)
      
      ############################ CHANGED THIS SECTION
      ################ Immediately returns any new value of dmax
      
      if (
        abs(prop_c - freeman_centralisation) <= tolerance &&
        abs(prop_avg - average_degree) <= average_degree_tolerance &&
        min(prop) >= min_degree &&
        is_graphical_safe(prop)
      ) {
        
        current <- prop
        current_c <- prop_c
        accepted_moves <- accepted_moves + 1L
        move_accepts[move_type] <- move_accepts[move_type] + 1L
        accepted <- TRUE
        
        # Save immediately if this is a new dmax value
        this_dmax <- max(current)
        
        if (!(this_dmax %in% seen_dmax) && saved < nsim) {
          
          seen_dmax <- c(seen_dmax, this_dmax)
          
          key <- paste(current, collapse = "-")
          
          if (!unique_sequences || !(key %in% seen)) {
            saved <- saved + 1L
            
            out[[saved]] <- list(
              degree_sequence = current,
              dmax = this_dmax,
              realised_centralisation = current_c,
              realised_average_degree = mean(current),
              realised_total_degree = sum(current)
            )
            
            
            
            if (unique_sequences) {
              seen <- c(seen, key)
            }
            
            if (verbose) {
              message(
                "Saved ", saved, "/", nsim,
                " | step = ", total_steps,
                " | dmax = ", max(current),
                " | C = ", round(current_c, 4),
                " | avg degree = ", round(mean(current), 4),
                " | total degree = ", sum(current)
              )
            }
          }
        }
      }
    }
    
    if (store_trace) {
      trace_df[nrow(trace_df) + 1L, ] <- list(
        step = total_steps,
        centralisation = current_c,
        average_degree = mean(current),
        total_degree = sum(current),
        dmax = max(current),
        accepted = accepted,
        move_type = move_type
      )
    }
    
    if (total_steps > burnin && ((total_steps - burnin) %% thin == 0)) {
      
      current_c <- freeman_from_degree(current)
      
      if (
        abs(current_c - freeman_centralisation) > tolerance ||
        abs(mean(current) - average_degree) > average_degree_tolerance ||
        min(current) < min_degree ||
        !is_graphical_safe(current)
      ) {
        next
      }
      
      key <- paste(current, collapse = "-")
      
      if (unique_sequences && key %in% seen) next
      
      saved <- saved + 1L
      
      out[[saved]] <- list(
        degree_sequence = current,
        dmax = max(current),
        realised_centralisation = current_c,
        realised_average_degree = mean(current),
        realised_total_degree = sum(current)
      )
      
      if (unique_sequences) {
        seen <- c(seen, key)
      }
      
      if (verbose) {
        message(
          "Saved ", saved, "/", nsim,
          " | step = ", total_steps,
          " | dmax = ", max(current),
          " | C = ", round(current_c, 4),
          " | avg degree = ", round(mean(current), 4),
          " | total degree = ", sum(current)
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
    realised_total_degree = sapply(out, `[[`, "realised_total_degree"),
    n_found = saved,
    nsim_requested = nsim,
    initial_total_degree = total_degree,
    size = n,
    average_degree_target = average_degree,
    average_degree_tolerance = average_degree_tolerance,
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

################################################################################
# Convenience helpers
################################################################################

summarise_degseq_sample <- function(out) {
  data.frame(
    seq_id = seq_along(out$degree_sequences),
    dmax = out$dmax,
    centralisation = out$realised_centralisation,
    average_degree = out$realised_average_degree,
    total_degree = out$realised_total_degree
  )
}




################################################################################
# Test example
################################################################################

# testsim <- simulateNetworks(test_list,          # These parameters resulted in an 
#                            nmAtt = 0.5,         # average component number of 1.56
#                            gwdeg = 1,           # reducing gwesp and adding in nmAtt
#                            gwesp = 0.3,         # seemed to help things
#                            gwdsp = -0.025,
#                            nsim = 2)


out <- degree_sequence_sample_mcmc(
  nsim = 50,
  size = 30,
  average_degree = 3,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 123,
  unique_sequences = TRUE,
  verbose = TRUE,
  store_trace = TRUE
)


unique(summarise_degseq_sample(out)[2:5])[3]

unique(summarise_degseq_sample(out)[2:5])
plot_degseq_trace(out)

seq_tab <- data.frame(
  seq_id = seq_along(out$degree_sequences),
  dmax = out$dmax,
  centralisation = out$realised_centralisation,
  average_degree = out$realised_average_degree,
  total_degree = out$realised_total_degree,
  sd_degree = sapply(out$degree_sequences, sd),
  degree_iqr = sapply(out$degree_sequences, IQR),
  n_degree_values = sapply(out$degree_sequences, function(x) length(unique(x)))
)

unique(seq_tab[, c(
  "dmax", "centralisation", "average_degree",
  "total_degree", "sd_degree", "degree_iqr", "n_degree_values"
)])

out2 <- degree_sequence_sample_mcmc(
  nsim = 50,
  size = 30,
  average_degree = 5,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.15,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 123,
  unique_sequences = TRUE,
  verbose = TRUE,
  store_trace = TRUE
)

unique(summarise_degseq_sample(out2)[2:5])

tab <- summarise_degseq_sample(out2)

tab_round <- tab %>%
  mutate(
    centralisation = round(centralisation, 4),
    average_degree = round(average_degree, 4)
  )

uniq_tab <- unique(tab_round[, c("dmax", "centralisation", "average_degree", "total_degree")])

uniq_tab
keep_ids <- tab_round %>%
  mutate(seq_id = row_number()) %>%
  inner_join(uniq_tab,
             by = c("dmax", "centralisation", "average_degree", "total_degree")) %>%
  pull(seq_id)

selected_degseqs <- out$degree_sequences[keep_ids]

test_list <- netFromDegSeq(selected_degseqs)

testsim <- simulateNetworks(test_list,
                                  nmAtt = 0.5,
                                  gwdeg = 1,
                                  gwesp = 0.3,
                                  gwdsp = -0.025,
                                  nsim = 2)

network_summary(testsim)
plotSimNetworks(testsim)
plot_degseq_trace(out2)

seq_tab <- data.frame(
  seq_id = seq_along(out2$degree_sequences),
  dmax = out2$dmax,
  centralisation = out2$realised_centralisation,
  average_degree = out2$realised_average_degree,
  total_degree = out2$realised_total_degree,
  sd_degree = sapply(out2$degree_sequences, sd),
  degree_iqr = sapply(out2$degree_sequences, IQR),
  n_degree_values = sapply(out2$degree_sequences, function(x) length(unique(x)))
)

unique(seq_tab[, c(
  "dmax", "centralisation", "average_degree",
  "total_degree", "sd_degree", "degree_iqr", "n_degree_values"
)])


##############################                               ###########################
############################## Simulating connected networks ###########################
##############################                               ###########################

is_connected_network <- function(net) {
  comp <- sna::component.dist(net)
  length(comp$csize) == 1L
}

simulateNetworks <- function(net_list, 
                                       target_connected = 1,
                                       max_attempts = 100,
                                       require_connected = TRUE,
                                       nfAtt = 0,
                                       nmAtt = 0,
                                       gwdeg = 0.5,
                                       gwesp = 0.5,
                                       gwdsp = -0.025,
                                       att_prob = c(3, 1),
                                       verbose = TRUE) {
  
  all_sims <- list()
  diagnostics <- list()
  counter <- 1L
  
  for (i in seq_along(net_list)) {
    
    net <- net_list[[i]]
    n <- network::network.size(net)
    
    deg <- sna::degree(net, gmode = "graph")
    
    accepted <- 0L
    attempted <- 0L
    rejected <- 0L
    
    network::set.vertex.attribute(
      net,
      attrname = "att",
      value = sample(c("A", "B"), n, replace = TRUE, prob = att_prob)
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
    
    while (accepted < target_connected && attempted < max_attempts) {
      
      attempted <- attempted + 1L
      
      sim <- suppressMessages(
        ergm::simulate_formula(
        form,
        constraints = ~degreedist,
        coef = coefs,
        nsim = 1,
        output = "network")
      )
      
      if (inherits(sim, "network")) {
        sim_net <- sim
      } else {
        sim_net <- sim[[1]]
      }
      
      connected <- is_connected_network(sim_net)
      
      if (!require_connected || connected) {
        
        accepted <- accepted + 1L
        
        network::set.network.attribute(sim_net, "basis_id", i)
        network::set.network.attribute(sim_net, "attempt_id", attempted)
        network::set.network.attribute(sim_net, "accepted_id", accepted)
        network::set.network.attribute(sim_net, "connected", connected)
        
        network::set.network.attribute(sim_net, "basis_dmax", max(deg))
        network::set.network.attribute(sim_net, "basis_average_degree", mean(deg))
        network::set.network.attribute(sim_net, "basis_total_degree", sum(deg))
        network::set.network.attribute(sim_net, "basis_sd_degree", stats::sd(deg))
        network::set.network.attribute(sim_net, "basis_degree_iqr", stats::IQR(deg))
        network::set.network.attribute(sim_net, "basis_n_degree_values", length(unique(deg)))
        
        all_sims[[counter]] <- sim_net
        counter <- counter + 1L
        
      } else {
        rejected <- rejected + 1L
      }
    }
    
    diagnostics[[i]] <- data.frame(
      basis_id = i,
      attempted = attempted,
      accepted = accepted,
      rejected = attempted - accepted,
      success_rate = accepted / attempted,
      target_connected = target_connected,
      target_met = accepted >= target_connected,
      dmax = max(deg),
      average_degree = mean(deg),
      total_degree = sum(deg),
      sd_degree = stats::sd(deg),
      degree_iqr = stats::IQR(deg),
      n_degree_values = length(unique(deg))
    )
    
    if (verbose) {
      message(
        "Basis ", i, "/", length(net_list),
        " | accepted ", accepted, "/", target_connected,
        " | attempts = ", attempted,
        " | success rate = ", round(accepted / attempted, 3)
      )
    }
  }
  
  list(
    networks = all_sims,
    diagnostics = dplyr::bind_rows(diagnostics)
  )
}

newTest <- simulateNetworks(test_list,
                           target_connected = 5,
                           nmAtt = 0.5,
                           gwdeg = 1,
                           gwesp = 0.3,
                           gwdsp = -0.025)

newTest$diagnostics

diag <- newTest$diagnostics

diag[diag$target_met == FALSE, ]

summary(diag$success_rate)

diag[order(diag$success_rate), ]

plotSimNetworks(newTest$networks)

################################################################################

# Testing new parameter tolerances

# For AD
# N = 30, AD tol = 0.3
# N = 60, AD tol = 0.59
# N = 120, AD tol = 1.19

# I am slightly concerned with AD3 and N120 as this is quite sparse,
# there are cases like this in the literature, however

# I think shifting to C0.4 is a good move

summarise_degseq_features <- function(out) {
  data.frame(
    seq_id = seq_along(out$degree_sequences),
    dmax = out$dmax,
    centralisation = round(out$realised_centralisation, 3),
    average_degree = round(out$realised_average_degree, 3),
    total_degree = sapply(out$degree_sequences, sum),
    sd_degree = round(sapply(out$degree_sequences, sd), 3),
    degree_iqr = sapply(out$degree_sequences, IQR),
    n_degree_values = sapply(out$degree_sequences, function(x) length(unique(x)))
  )
}

##################### Stress test: Low N, High C

n30ad3c0.4 <- degree_sequence_sample_mcmc(
  nsim = 50,
  size = 30,
  average_degree = 3,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.4,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 123,
  unique_sequences = TRUE,
  verbose = TRUE,
  store_trace = TRUE
)

summarise_degseq_features(n30ad3c0.4)

n30ad3c04_ids <- summarise_degseq_features(n30ad3c0.4) %>%
  distinct(across(c("dmax", 
                    #"centralisation", 
                    "average_degree", 
                    "degree_iqr"#, "n_degree_values"
                    )),
           .keep_all = TRUE) %>%
  pull(seq_id)

n30ad3c0.4$degree_sequences[n30ad3c04_ids]

############################# Stress test, Low N, Low C

n30ad3c0.4 <- degree_sequence_sample_mcmc(
  nsim = 50,
  size = 30,
  average_degree = 3,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 123,
  unique_sequences = TRUE,
  verbose = TRUE,
  store_trace = TRUE
)

summarise_degseq_features(n30ad3c0.4)

summarise_degseq_features(n30ad3c0.4) %>%
  distinct(across(c("dmax", 
                    "centralisation", 
                    "average_degree" #, 
                    #"degree_iqr"#, "n_degree_values"
  )),
  .keep_all = TRUE)


summarise_degseq_features(n30ad3c0.4) %>%
  distinct(across(c("dmax", 
                    "centralisation", 
                    #"average_degree" #, 
                    "degree_iqr"#, "n_degree_values"
  )),
  .keep_all = TRUE)

############################# Stress test, Low N, Low C, High AD

n30ad3c01 <- degree_sequence_sample_mcmc(
  nsim = 50,
  size = 30,
  average_degree = 6,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  burnin = 20000,
  thin = 2000,
  seed = 123,
  unique_sequences = TRUE,
  verbose = TRUE,
  store_trace = TRUE
)

summarise_degseq_features(n30ad3c0.4) %>%
  distinct(across(c("dmax", 
                    "centralisation", 
                    "average_degree" #, 
                    #"degree_iqr"#, "n_degree_values"
  )),
  .keep_all = TRUE)


summarise_degseq_features(n30ad3c0.4) %>%
  distinct(across(c("dmax", 
                    "centralisation", 
                    #"average_degree" #, 
                    "degree_iqr"#, "n_degree_values"
