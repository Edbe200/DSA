# Run Markov trace
run_markov_trace <- function(m_TP, v_init_state, n_cycles){
  m_trace <- matrix(0, nrow = n_cycles, ncol = length(v_init_state))
  colnames(m_trace) <- names(v_init_state)
  m_trace[1, ] <- v_init_state
  
  for(i in 2:n_cycles){
    m_trace[i, ] <- m_trace[i-1, ] %*% m_TP
  }
  
  m_trace
}

# Calculate QALYs
calc_qalys <- function(m_trace, v_utilities, disc_rate) {
  n_cycles <- nrow(m_trace)
  v_disc_factors <- 1 / (1 + disc_rate)^(0:(n_cycles - 1))
  sum((m_trace * v_disc_factors) %*% v_utilities)
}

# Calculate costs
calc_costs <- function(m_trace, v_costs, disc_rate) {
  n_cycles <- nrow(m_trace)
  v_disc_factors <- 1 / (1 + disc_rate)^(0:(n_cycles - 1))
  sum((m_trace * v_disc_factors) %*% v_costs)
}

# Apply treatment effects to transition probabilities
apply_treatment_effect <- function(m_TP, transitions) {
  m_TP_new <- m_TP
  
  for(trans_name in names(transitions)) {
    states <- strsplit(trans_name, "_")[[1]]
    m_TP_new[states[1], states[2]] <- m_TP_new[states[1], states[2]] * transitions[[trans_name]]
  }
  
  affected_states <- unique(sapply(strsplit(names(transitions), "_"), `[`, 1))
  
  for(state in affected_states) {
    m_TP_new[state, state] <- 1 - sum(m_TP_new[state, colnames(m_TP_new) != state])
  }
  
  m_TP_new
}

# Run complete model
run_model <- function(m_TP, v_init_state, n_cycles, v_utilities, v_costs, r, treatment = NULL) {
  
  if (!is.null(treatment)) {
    m_TP <- apply_treatment_effect(m_TP, treatment$transitions)
  }
  
  m_trace <- run_markov_trace(m_TP, v_init_state, n_cycles)
  total_qalys <- calc_qalys(m_trace, v_utilities, r)
  total_costs <- calc_costs(m_trace, v_costs, r)
  
  if (!is.null(treatment)) {
    person_time_treated <- sum(m_trace[, treatment$apply_states])
    total_costs <- total_costs + (person_time_treated * treatment$cost)
  }
  
  list(qalys = total_qalys, costs = total_costs, trace = m_trace)
}

# Update parameters for sensitivity analysis
update_parameters <- function(base_params, param_name, new_value) {
  params <- base_params
  
  switch(param_name,
         u_H = { params$v_utilities["H"] <- new_value },
         u_S1 = { params$v_utilities["S1"] <- new_value },
         u_S2 = { params$v_utilities["S2"] <- new_value },
         u_R = { params$v_utilities["R"] <- new_value },
         c_S1 = { params$v_costs["S1"] <- new_value },
         c_S2 = { params$v_costs["S2"] <- new_value },
         c_trt_standard = { params$trt_standard$cost <- new_value },
         c_trt_advanced = { params$trt_advanced$cost <- new_value },
         rr_S1_S2_standard = { params$trt_standard$transitions$S1_S2 <- new_value },
         rr_S1_D_standard = { params$trt_standard$transitions$S1_D <- new_value },
         rr_S2_D_standard = { params$trt_standard$transitions$S2_D <- new_value },
         rr_S1_S2_advanced = { params$trt_advanced$transitions$S1_S2 <- new_value },
         rr_S1_D_advanced = { params$trt_advanced$transitions$S1_D <- new_value },
         rr_S2_D_advanced = { params$trt_advanced$transitions$S2_D <- new_value },
         discount_rate = { params$r <- new_value }
  )
  
  params
}

# Compare two strategies
run_comparison <- function(m_TP, v_init_state, n_cycles, params, 
                           strategy1_name, strategy2_name) {
  
  get_treatment <- function(strategy_name) {
    switch(strategy_name,
           "Baseline" = NULL,
           "Standard Care" = params$trt_standard,
           "Advanced Therapy" = params$trt_advanced
    )
  }
  
  results1 <- run_model(m_TP, v_init_state, n_cycles, 
                        params$v_utilities, params$v_costs, params$r,
                        treatment = get_treatment(strategy1_name))
  
  results2 <- run_model(m_TP, v_init_state, n_cycles,
                        params$v_utilities, params$v_costs, params$r,
                        treatment = get_treatment(strategy2_name))
  
  inc_costs <- results2$costs - results1$costs
  inc_qalys <- results2$qalys - results1$qalys
  icer <- if(inc_qalys > 0) inc_costs / inc_qalys else NA
  
  list(icer = icer, inc_costs = inc_costs, inc_qalys = inc_qalys)
}

# One-way sensitivity analysis
run_owsa <- function(m_TP, v_init_state, n_cycles, base_params, dsa_params,
                     strategy1 = "Baseline", strategy2 = "Standard Care") {
  
  dsa_results <- data.frame()
  
  for(i in 1:nrow(dsa_params)) {
    param_name <- dsa_params$parameter[i]
    cat("Testing parameter:", param_name, "\n")
    
    params_low <- update_parameters(base_params, param_name, dsa_params$lower[i])
    results_low <- run_comparison(m_TP, v_init_state, n_cycles, params_low, 
                                  strategy1, strategy2)
    
    params_high <- update_parameters(base_params, param_name, dsa_params$upper[i])
    results_high <- run_comparison(m_TP, v_init_state, n_cycles, params_high, 
                                   strategy1, strategy2)
    
    dsa_results <- rbind(dsa_results, 
                         data.frame(
                           parameter = param_name,
                           bound = "lower",
                           value = dsa_params$lower[i],
                           icer = results_low$icer,
                           inc_costs = results_low$inc_costs,
                           inc_qalys = results_low$inc_qalys
                         ),
                         data.frame(
                           parameter = param_name,
                           bound = "upper",
                           value = dsa_params$upper[i],
                           icer = results_high$icer,
                           inc_costs = results_high$inc_costs,
                           inc_qalys = results_high$inc_qalys
                         )
    )
  }
  
  dsa_results
}
