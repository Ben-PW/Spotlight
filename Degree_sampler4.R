degree_sequence_profile_sampler <- function(
    nsim,
    size,
    average_degree,
    average_degree_tolerance = 0.3,
    freeman_centralisation,
    tolerance = 0.05,
    min_degree = 1L,
    max_steps = 500000,
    max_start_tries = 10000,
    move_probs = c(
      move1 = 0.40,
      split2 = 0.20,
      merge2 = 0.20,
      add2 = 0.10,
      drop2 = 0.10
    ),
    profile_digits_c = 3,
    profile_digits_iqr = 2,
    seed = NULL,
    verbose = TRUE
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- as.integer(size)
  start_total <- get_total_degree(n, average_degree)
  
  is_valid <- function(deg) {
    c_val <- freeman_from_degree(deg)
    avg_val <- mean(deg)
    
    abs(c_val - freeman_centralisation) <= tolerance &&
      abs(avg_val - average_degree) <= average_degree_tolerance &&
      min(deg) >= min_degree &&
      is_graphical_safe(deg)
  }
  
  make_profile <- function(deg) {
    paste(
      max(deg),
      round(freeman_from_degree(deg), profile_digits_c),
      #round(IQR(deg), profile_digits_iqr),
      round(mean(deg), 3),
      sep = "_"
    )
  }
  
  save_sequence <- function(deg, step, profile) {
    list(
      degree_sequence = deg,
      dmax = max(deg),
      realised_centralisation = freeman_from_degree(deg),
      realised_average_degree = mean(deg),
      realised_total_degree = sum(deg),
      sd_degree = sd(deg),
      degree_iqr = IQR(deg),
      n_degree_values = length(unique(deg)),
      profile = profile,
      step = step
    )
  }
  
  # ---- initial sequence ----
  
  current <- construct_initial_degseq(
    size = n,
    total_degree = start_total,
    min_degree = min_degree,
    max_tries = max_start_tries
  )
  
  if (is.null(current)) {
    stop("Could not construct initial graphical degree sequence.")
  }
  
  current_c <- freeman_from_degree(current)
  
  # If initial sequence is outside centralisation band, walk until valid
  init_steps <- 0L
  
  while (!is_valid(current)) {
    init_steps <- init_steps + 1L
    
    if (init_steps > max_steps) {
      stop("Could not find initial sequence inside target bands.")
    }
    
    move <- propose_degseq_move_mixed(
      current,
      min_degree = min_degree,
      move_probs = move_probs
    )
    
    prop <- move$prop
    if (is.null(prop)) next
    
    prop_c <- freeman_from_degree(prop)
    
    current_dist <-
      abs(current_c - freeman_centralisation) +
      abs(mean(current) - average_degree)
    
    prop_dist <-
      abs(prop_c - freeman_centralisation) +
      abs(mean(prop) - average_degree)
    
    if (prop_dist < current_dist && is_graphical_safe(prop)) {
      current <- prop
      current_c <- prop_c
    }
  }
  
  if (verbose) {
    message(
      "Initial valid state found",
      " | C = ", round(freeman_from_degree(current), 4),
      " | avg degree = ", round(mean(current), 4),
      " | dmax = ", max(current),
      " | total degree = ", sum(current),
      " | init_steps = ", init_steps
    )
  }
  
  # ---- profile collection ----
  
  out <- vector("list", nsim)
  seen_profiles <- character(0)
  seen_sequences <- character(0)
  
  saved <- 0L
  total_steps <- 0L
  attempted_moves <- 0L
  accepted_moves <- 0L
  
  while (saved < nsim && total_steps < max_steps) {
    
    total_steps <- total_steps + 1L
    
    move <- propose_degseq_move_mixed(
      current,
      min_degree = min_degree,
      move_probs = move_probs
    )
    
    prop <- move$prop
    if (is.null(prop)) next
    
    attempted_moves <- attempted_moves + 1L
    
    if (!is_valid(prop)) next
    
    current <- prop
    accepted_moves <- accepted_moves + 1L
    
    profile <- make_profile(current)
    key <- paste(current, collapse = "-")
    
    if (profile %in% seen_profiles) next
    if (key %in% seen_sequences) next
    
    saved <- saved + 1L
    
    out[[saved]] <- save_sequence(
      deg = current,
      step = total_steps,
      profile = profile
    )
    
    seen_profiles <- c(seen_profiles, profile)
    seen_sequences <- c(seen_sequences, key)
    
    if (verbose) {
      message(
        "Saved ", saved, "/", nsim,
        " | step = ", total_steps,
        " | dmax = ", max(current),
        " | C = ", round(freeman_from_degree(current), 4),
        " | avg degree = ", round(mean(current), 4),
        " | total degree = ", sum(current),
        " | sd = ", round(sd(current), 3),
        " | IQR = ", round(IQR(current), 3),
        " | nvals = ", length(unique(current))
      )
    }
  }
  
  out <- out[seq_len(saved)]
  
  list(
    degree_sequences = lapply(out, `[[`, "degree_sequence"),
    dmax = sapply(out, `[[`, "dmax"),
    realised_centralisation = sapply(out, `[[`, "realised_centralisation"),
    realised_average_degree = sapply(out, `[[`, "realised_average_degree"),
    realised_total_degree = sapply(out, `[[`, "realised_total_degree"),
    sd_degree = sapply(out, `[[`, "sd_degree"),
    degree_iqr = sapply(out, `[[`, "degree_iqr"),
    n_degree_values = sapply(out, `[[`, "n_degree_values"),
    profiles = sapply(out, `[[`, "profile"),
    steps_saved = sapply(out, `[[`, "step"),
    n_found = saved,
    nsim_requested = nsim,
    size = n,
    average_degree_target = average_degree,
    average_degree_tolerance = average_degree_tolerance,
    freeman_target = freeman_centralisation,
    tolerance = tolerance,
    total_steps = total_steps,
    attempted_moves = attempted_moves,
    accepted_moves = accepted_moves,
    acceptance_rate = if (attempted_moves > 0) accepted_moves / attempted_moves else NA_real_,
    init_steps = init_steps
  )
}

n30_ad3_c01 <- degree_sequence_profile_sampler(
  nsim = 50,
  size = 30,
  average_degree = 3,
  average_degree_tolerance = 0.3,
  freeman_centralisation = 0.1,
  tolerance = 0.05,
  min_degree = 1,
  seed = 123
)

summarise_degseq_features(n30_ad3_c01)

summarise_degseq_features(n30_ad3_c01) %>%
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