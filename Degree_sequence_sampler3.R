################################################################################
# Degree-sequence MCMC sampler with mixed move sizes
# - fixed size
# - target average degree WITH optional tolerance
# - target Freeman centralisation within tolerance
# - minimum degree enforced
# - graphicality enforced
# - mixed proposal moves, including total-degree-changing moves
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
      
      if (
        abs(prop_c - freeman_centralisation) <= tolerance &&
        within_avg_degree_band(prop, average_degree, average_degree_tolerance) &&
        is_graphical_safe(prop)
      ) {
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
        average_degree = mean(current),
        total_degree = sum(current),
        dmax = max(current),
        accepted = accepted,
        move_type = move_type
      )
    }
    
    if (total_steps > burnin && ((total_steps - burnin) %% thin == 0)) {
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

plot_degseq_trace <- function(out) {
  if (is.null(out$trace)) stop("No trace stored in output.")
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(1, 3))
  
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
    out$trace$average_degree,
    type = "l",
    xlab = "Step",
    ylab = "Average degree",
    main = "Trace: average degree"
  )
  abline(h = out$average_degree_target, lty = 2)
  abline(h = out$average_degree_target + out$average_degree_tolerance, lty = 3)
  abline(h = out$average_degree_target - out$average_degree_tolerance, lty = 3)
  
  plot(
    out$trace$step,
    out$trace$dmax,
    type = "l",
    xlab = "Step",
    ylab = "Maximum degree",
    main = "Trace: dmax"
  )
}

################################################################################
# Test example
################################################################################

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

unique(summarise_degseq_sample(out))
plot_degseq_trace(out)
