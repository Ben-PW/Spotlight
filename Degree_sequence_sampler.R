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

make_random_initial_degseq <- function(n, total_degree, min_degree = 1L, max_tries = 5000) {
  if (total_degree < n * min_degree) {
    stop("Total degree too small for requested minimum degree.")
  }
  
  for (try in seq_len(max_tries)) {
    deg <- rep.int(min_degree, n)
    remaining <- total_degree - sum(deg)
    
    while (remaining > 0) {
      eligible <- which(deg < (n - 1))
      if (length(eligible) == 0) break
      
      i <- sample(eligible, 1)
      deg[i] <- deg[i] + 1L
      remaining <- remaining - 1L
    }
    
    deg <- sort(deg, decreasing = TRUE)
    
    if (is_graphical_safe(deg)) {
      return(deg)
    }
  }
  
  stop("Could not generate a graphical random initial sequence.")
}

degree_sequence_target_centralisation <- function(
    size,
    average_degree,
    freeman_centralisation,
    tolerance = 0.01,
    min_degree = 1L,
    max_iter = 50000,
    n_restarts = 20,
    seed = NULL,
    verbose = FALSE,
    allow_sideways = TRUE,
    p_worse = 0.02
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- as.integer(size)
  if (n < 3) stop("size must be at least 3.")
  if (average_degree < 0 || average_degree > (n - 1)) {
    stop("average_degree must be between 0 and n - 1.")
  }
  if (freeman_centralisation < 0 || freeman_centralisation > 1) {
    stop("freeman_centralisation must be between 0 and 1.")
  }
  
  target_total_degree <- round(n * average_degree)
  if ((target_total_degree %% 2) != 0) {
    target_total_degree <- target_total_degree + 1L
  }
  
  if (target_total_degree < n * min_degree) {
    stop("Requested average degree is too small for min_degree.")
  }
  
  best_deg <- NULL
  best_diff <- Inf
  best_cent <- NA_real_
  
  for (restart in seq_len(n_restarts)) {
    
    deg <- make_random_initial_degseq(
      n = n,
      total_degree = target_total_degree,
      min_degree = min_degree
    )
    
    current_cent <- freeman_from_degree(deg)
    current_diff <- abs(current_cent - freeman_centralisation)
    
    if (current_diff < best_diff) {
      best_deg <- deg
      best_diff <- current_diff
      best_cent <- current_cent
    }
    
    if (current_diff <= tolerance) {
      return(list(
        degree_sequence = deg,
        realised_centralisation = current_cent,
        converged = TRUE
      ))
    }
    
    for (iter in seq_len(max_iter)) {
      proposal <- deg
      
      donors <- which(proposal > min_degree)
      receivers <- which(proposal < (n - 1))
      
      if (length(donors) == 0 || length(receivers) == 0) break
      
      from <- sample(donors, 1)
      to   <- sample(receivers, 1)
      
      if (from == to) next
      
      proposal[from] <- proposal[from] - 1L
      proposal[to]   <- proposal[to] + 1L
      proposal <- sort(proposal, decreasing = TRUE)
      
      if (any(proposal < min_degree)) next
      if (!is_graphical_safe(proposal)) next
      
      prop_cent <- freeman_from_degree(proposal)
      prop_diff <- abs(prop_cent - freeman_centralisation)
      
      accept <- FALSE
      
      if (prop_diff < current_diff) {
        accept <- TRUE
      } else if (allow_sideways && prop_diff == current_diff) {
        accept <- TRUE
      } else if (runif(1) < p_worse) {
        accept <- TRUE
      }
      
      if (accept) {
        deg <- proposal
        current_cent <- prop_cent
        current_diff <- prop_diff
        
        if (current_diff < best_diff) {
          best_deg <- deg
          best_diff <- current_diff
          best_cent <- current_cent
        }
        
        if (current_diff <= tolerance) {
          return(list(
            degree_sequence = deg,
            realised_centralisation = current_cent,
            converged = TRUE
          ))
        }
      }
    }
  }
  
  list(
    degree_sequence = best_deg,
    realised_centralisation = best_cent,
    best_abs_diff = best_diff,
    converged = FALSE
  )
}

degree_sequence_sample <- function(
    nsim = 20,
    size,
    average_degree,
    freeman_centralisation,
    tolerance = 0.01,
    min_degree = 1L,
    max_iter = 50000,
    n_restarts = 20,
    max_attempts = 500,
    seed = NULL,
    verbose = FALSE
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  results <- list()
  seen <- character(0)
  n_found <- 0
  attempts <- 0
  
  while (n_found < nsim && attempts < max_attempts) {
    attempts <- attempts + 1
    
    res <- degree_sequence_target_centralisation(
      size = size,
      average_degree = average_degree,
      freeman_centralisation = freeman_centralisation,
      tolerance = tolerance,
      min_degree = min_degree,
      max_iter = max_iter,
      n_restarts = n_restarts,
      seed = sample.int(.Machine$integer.max, 1),
      verbose = FALSE
    )
    
    if (!isTRUE(res$converged)) next
    
    key <- paste(res$degree_sequence, collapse = "-")
    
    if (!(key %in% seen)) {
      n_found <- n_found + 1
      results[[n_found]] <- res$degree_sequence
      seen <- c(seen, key)
      
      if (verbose) {
        message("Found unique sequence ", n_found, "/", nsim,
                " at attempt ", attempts)
      }
    }
  }
  
  if (n_found < nsim) {
    warning("Only found ", n_found, " unique valid sequences out of requested ", nsim)
  }
  
  list(
    degree_sequences = results,
    n_found = n_found,
    attempts = attempts,
    success_rate = n_found / attempts
  )
}


out <- degree_sequence_sample(
  nsim = 20,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.5,
  tolerance = 0.05,
  seed = 123
)

length(out$degree_sequences)   # how many you actually got
out$success_rate

out

deg <- out$degree_sequences[[1]]

g <- igraph::realize_degseq(
  deg,
  allowed.edge.types = "simple",
  method = "smallest"
)

plot(g)

net <- intergraph::asNetwork(g)

network::set.vertex.attribute(
  net,
  attrname = "att",
  value = sample(c("A", "B"), 100, replace = TRUE, prob = c(3,1))
)

form <- net ~
  nodefactor("att") +
  nodematch("att") +
  gwdegree(0.3, fixed = TRUE) +
  gwesp(0.3, fixed = TRUE) +
  gwdsp(0.3, fixed = TRUE)

coefs <- c(
  nodefactor.att.B = 1,
  nodematch.att = 0,
  gwdeg.fixed = 1,
  gwesp.fixed = 2,
  gwdsp.fixed = -0.025
)

sim <- stats::simulate(
  form,
  constraints = ~degreedist,
  coef = coefs,
  nsim = 20,
  output = "network"
)

sim

summary(sim ~ edges)

for (i in 1:20) {
  plot(sim[[i]], main = paste0("sim ", i))
}


out2 <- degree_sequence_sample(
  nsim = 3,
  size = 100,
  average_degree = 6,
  freeman_centralisation = 0.05,
  tolerance = 0.01,
  seed = 123
)

out2$degree_sequences[[1]]

################################################################################
# Degree sequence generator v2

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

make_random_initial_degseq <- function(n, total_degree, min_degree = 1L, max_tries = 5000) {
  if (total_degree < n * min_degree) {
    stop("Total degree too small for requested minimum degree.")
  }
  
  for (try in seq_len(max_tries)) {
    deg <- rep.int(min_degree, n)
    remaining <- total_degree - sum(deg)
    
    while (remaining > 0) {
      eligible <- which(deg < (n - 1))
      if (length(eligible) == 0) break
      
      i <- sample(eligible, 1)
      deg[i] <- deg[i] + 1L
      remaining <- remaining - 1L
    }
    
    deg <- sort(deg, decreasing = TRUE)
    
    if (is_graphical_safe(deg)) {
      return(deg)
    }
  }
  
  stop("Could not generate a graphical random initial sequence.")
}

degree_sequence_target_centralisation <- function(
    size,
    average_degree,
    freeman_centralisation,
    tolerance = 0.01,
    min_degree = 1L,
    max_iter = 50000,
    n_restarts = 20,
    seed = NULL,
    verbose = FALSE,
    allow_sideways = TRUE,
    p_worse = 0.02,
    collect_n = 10,
    extra_iter = 5000,
    return_random = TRUE
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- as.integer(size)
  if (n < 3) stop("size must be at least 3.")
  if (average_degree < 0 || average_degree > (n - 1)) {
    stop("average_degree must be between 0 and n - 1.")
  }
  if (freeman_centralisation < 0 || freeman_centralisation > 1) {
    stop("freeman_centralisation must be between 0 and 1.")
  }
  
  target_total_degree <- round(n * average_degree)
  if ((target_total_degree %% 2) != 0) {
    target_total_degree <- target_total_degree + 1L
  }
  
  if (target_total_degree < n * min_degree) {
    stop("Requested average degree is too small for min_degree.")
  }
  
  best_deg <- NULL
  best_diff <- Inf
  best_cent <- NA_real_
  
  acceptable <- list()
  acceptable_keys <- character(0)
  first_hit_found <- FALSE
  post_hit_counter <- 0L
  
  for (restart in seq_len(n_restarts)) {
    
    deg <- make_random_initial_degseq(
      n = n,
      total_degree = target_total_degree,
      min_degree = min_degree
    )
    
    current_cent <- freeman_from_degree(deg)
    current_diff <- abs(current_cent - freeman_centralisation)
    
    if (current_diff < best_diff) {
      best_deg <- deg
      best_diff <- current_diff
      best_cent <- current_cent
    }
    
    if (current_diff <= tolerance) {
      key <- paste(deg, collapse = "-")
      if (!(key %in% acceptable_keys)) {
        acceptable[[length(acceptable) + 1L]] <- list(
          degree_sequence = deg,
          realised_centralisation = current_cent,
          abs_diff = current_diff
        )
        acceptable_keys <- c(acceptable_keys, key)
      }
      first_hit_found <- TRUE
    }
    
    for (iter in seq_len(max_iter)) {
      proposal <- deg
      
      donors <- which(proposal > min_degree)
      receivers <- which(proposal < (n - 1))
      
      if (length(donors) == 0 || length(receivers) == 0) break
      
      from <- sample(donors, 1)
      to   <- sample(receivers, 1)
      
      if (from == to) next
      
      proposal[from] <- proposal[from] - 1L
      proposal[to]   <- proposal[to] + 1L
      proposal <- sort(proposal, decreasing = TRUE)
      
      if (any(proposal < min_degree)) next
      if (!is_graphical_safe(proposal)) next
      
      prop_cent <- freeman_from_degree(proposal)
      prop_diff <- abs(prop_cent - freeman_centralisation)
      
      accept <- FALSE
      
      if (prop_diff < current_diff) {
        accept <- TRUE
      } else if (allow_sideways && prop_diff == current_diff) {
        accept <- TRUE
      } else if (runif(1) < p_worse) {
        accept <- TRUE
      }
      
      if (!accept) next
      
      deg <- proposal
      current_cent <- prop_cent
      current_diff <- prop_diff
      
      if (current_diff < best_diff) {
        best_deg <- deg
        best_diff <- current_diff
        best_cent <- current_cent
      }
      
      if (current_diff <= tolerance) {
        key <- paste(deg, collapse = "-")
        
        if (!(key %in% acceptable_keys)) {
          acceptable[[length(acceptable) + 1L]] <- list(
            degree_sequence = deg,
            realised_centralisation = current_cent,
            abs_diff = current_diff
          )
          acceptable_keys <- c(acceptable_keys, key)
          
          if (verbose) {
            message(
              "Accepted unique candidate ", length(acceptable),
              " at restart ", restart,
              ", iter ", iter,
              ", cent = ", round(current_cent, 4),
              ", dmax = ", max(deg)
            )
          }
        }
        
        if (!first_hit_found) {
          first_hit_found <- TRUE
          post_hit_counter <- 0L
        }
      }
      
      if (first_hit_found) {
        post_hit_counter <- post_hit_counter + 1L
      }
      
      if (length(acceptable) >= collect_n) break
      if (first_hit_found && post_hit_counter >= extra_iter) break
    }
    
    if (length(acceptable) >= collect_n) break
    if (first_hit_found && post_hit_counter >= extra_iter) break
  }
  
  if (length(acceptable) > 0) {
    chosen <- if (return_random) {
      acceptable[[sample.int(length(acceptable), 1)]]
    } else {
      acceptable[[1]]
    }
    
    return(list(
      degree_sequence = chosen$degree_sequence,
      realised_centralisation = chosen$realised_centralisation,
      abs_diff = chosen$abs_diff,
      converged = TRUE,
      n_acceptable_found = length(acceptable),
      acceptable_centralisations = sapply(acceptable, `[[`, "realised_centralisation"),
      acceptable_dmax = sapply(acceptable, function(x) max(x$degree_sequence)),
      acceptable_sequences = lapply(acceptable, `[[`, "degree_sequence")
    ))
  }
  
  list(
    degree_sequence = best_deg,
    realised_centralisation = best_cent,
    best_abs_diff = best_diff,
    converged = FALSE,
    n_acceptable_found = 0
  )
}

degree_sequence_sample <- function(
    nsim = 20,
    size,
    average_degree,
    freeman_centralisation,
    tolerance = 0.01,
    min_degree = 1L,
    max_iter = 50000,
    n_restarts = 20,
    max_attempts = 500,
    seed = NULL,
    verbose = FALSE,
    collect_n = 10,
    extra_iter = 5000
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  results <- list()
  seen <- character(0)
  n_found <- 0
  attempts <- 0
  
  while (n_found < nsim && attempts < max_attempts) {
    attempts <- attempts + 1
    
    res <- degree_sequence_target_centralisation(
      size = size,
      average_degree = average_degree,
      freeman_centralisation = freeman_centralisation,
      tolerance = tolerance,
      min_degree = min_degree,
      max_iter = max_iter,
      n_restarts = n_restarts,
      seed = sample.int(.Machine$integer.max, 1),
      verbose = FALSE,
      collect_n = collect_n,
      extra_iter = extra_iter,
      return_random = TRUE
    )
    
    if (!isTRUE(res$converged)) next
    
    key <- paste(res$degree_sequence, collapse = "-")
    
    if (!(key %in% seen)) {
      n_found <- n_found + 1
      results[[n_found]] <- res$degree_sequence
      seen <- c(seen, key)
      
      if (verbose) {
        message(
          "Found unique sequence ", n_found, "/", nsim,
          " at attempt ", attempts,
          " (generator found ", res$n_acceptable_found, " acceptable candidates)"
        )
      }
    }
  }
  
  if (n_found < nsim) {
    warning("Only found ", n_found, " unique valid sequences out of requested ", nsim)
  }
  
  list(
    degree_sequences = results,
    n_found = n_found,
    attempts = attempts,
    success_rate = n_found / attempts
  )
}

out2 <- degree_sequence_sample(
  nsim = 5,
  size = 100,
  average_degree = 3,
  freeman_centralisation = 0.7,
  tolerance = 0.05,
  seed = 123,
  collect_n = 20,
  extra_iter = 10000,
  verbose = TRUE
)

out2
