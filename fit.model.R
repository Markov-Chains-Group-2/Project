fit.model <- function(n.sims, burn.in, inits = NULL, n_states, observed_data, M = 144) {
  library(MCMCpack)
  library(MASS)  # For mvrnorm
  
  Y <- observed_data
  n <- length(Y)

  
  
  ### FUNCTION TO EVALUATE THE LOG-POSTERIOR ###
  post <- function(log.lhood, emission_params, duration_params, p.init, trans_matrix) {
    if (any(emission_params <= 0) || any(duration_params <= 0) ) {return(-Inf)}

    log_prior <- 0
    
    # Log-prior for the transition matrix (only states 2 to n_states-1, ignoring states 1 and n_states)
    for (state in 2:(n_states - 1)) {
      trans_row <- trans_matrix[state, c(state - 1, state + 1)]
      # Ensure probabilities sum to 1 before computing log-prior
      trans_row <- trans_row / sum(trans_row)
      log_prior <- log_prior + log(ddirichlet(trans_row, c(1, 1)))  # Dirichlet prior for contiguous transitions
    }
    
    return(
      log.lhood + 
        sum(dgamma(emission_params[, 1], 0.001, 0.001, log = TRUE)) +  # Shape for emissions
        sum(dgamma(emission_params[, 2], 0.001, 0.001, log = TRUE)) +  # Scale for emissions
        sum(dgamma(duration_params[, 1], 0.001, 0.001, log = TRUE)) +  # Shape for durations
        sum(dgamma(duration_params[, 2], 0.001, 0.001, log = TRUE)) +  # Scale for durations
        log(ddirichlet(p.init, rep(1, n_states))) 
    )
  }

  ### LOGLIKELIHOOD VERIFY

  loglikelihood <- function(p.init, m.init, params_dist_emission, params_dist_sojourn, n_states, n, y) {
    # Verificar que todos los parámetros sean positivos
    #if (any(params_dist_emission<= 0)) {
    #    print("F")
    #    return(-Inf)
    #}
    
    #if (any(params_dist_sojourn<= 0)) {
    #    print("F")
    #    return(-Inf)
    #}
  
  # Si todos los parámetros son positivos, calcular la verosimilitud
  log_likelihood <- for.alg(p.init = p.init,
                            m.init = m.init,
                            params_dist_emission = params_dist_emission,
                            params_dist_sojourn = params_dist_sojourn,
                            n_states = n_states,
                            n = n,
                            y = y)
  
  return(log_likelihood)
}

  print("===== Iniciando MCMC =====")
  flush.console()
  ### INITIAL VALUES ###
  if (is.null(inits)) {
    emission_params <- matrix(c(runif(n_states, 1, 2), runif(n_states, 1, 2)), ncol = 2)
    duration_params <- matrix(c(runif(n_states, 1, 2), runif(n_states, 1, 2)), ncol = 2)
    
    p.initA <- rep(1 / n_states, n_states)
    trans_matrix <- diag(0, n_states)
    trans_matrix[row(trans_matrix) == col(trans_matrix) - 1] <- 0.5  
    trans_matrix[row(trans_matrix) == col(trans_matrix) + 1] <- 0.5
    trans_matrix[1, 2] <- 1
    trans_matrix[n_states, n_states - 1] <- 1
  } else {
    emission_params <- inits$emission_params
    duration_params <- inits$duration_params
    p.initA <- inits$p.init
    trans_matrix <- inits$trans_matrix
  }

  print("===== Parámetros inicializados ======")
  print("Parámetros de emisión (shape y scale):")
  print(emission_params)
  print("Parámetros de duración (shape y scale):")
  print(duration_params)
  print("Distribución inicial (p.initA):")
  print(p.initA)
  print("Matriz de transición:")
  print(trans_matrix)
  flush.console()

  paramA <- c(emission_params[, 1],
              emission_params[, 2],
              duration_params[, 1],
              duration_params[, 2])

  ### Setup for MCMC Sampling ###
  n_params <- 7 * n_states - 4
  #sims.matrix <- matrix(nrow = n.sims, ncol = n_params + 2)  # parámetros + likelihood + posterior
  sims_list <- vector("list", n.sims)  # parámetros + likelihood + posterior
  
  log.likelihoodA <- loglikelihood(p.init = p.initA,
                             m.init = trans_matrix,
                             params_dist_emission = emission_params,
                             params_dist_sojourn = duration_params,
                             n_states = n_states,
                             n = n,
                             y = Y)
  
  print(paste("Initial Log-likelihood:", log.likelihoodA))
  flush.console()
  log.postA <- post(log.likelihoodA, emission_params, duration_params, p.initA, trans_matrix)   
  
  accept1 <- accept2 <- accept3 <- 0
  
  U1 <- runif(n.sims)
  U2 <- runif(n.sims)
  U3 <- runif(n.sims)

  paramB <- paramA
  
  ### Variances of proposals for each parameter
  ss_shape_emission <- 1e-5  # Variance for shape of emission (n_states)
  ss_scale_emission <- 1e-5   # Variance for scale of emission (n_states)
  ss_shape_duration <- 5e-5    # Variance for shape of duration (n_states)
  ss_scale_duration <- 1e-5    # Variance for scale of duration (n_states)

  # Combine variances for each vector into a single vector
  VAR <- c(rep(ss_shape_emission, n_states),      # Shape for emissions
           rep(ss_scale_emission, n_states),     # Scale for emissions
           rep(ss_shape_duration, n_states),     # Shape for duration
           rep(ss_scale_duration, n_states))     # Scale for duration

  # Adjust the covariance matrix for the proposal distribution
  Sigma <- diag(VAR)

  # Generate proposal samples
  Q <- mvrnorm(n.sims, rep(0, length(VAR)), Sigma)  # Generate proposals
  Z_trans_matrix <- 200	### For adjusting acceptance rate of P
  Z_p.init <- 100	### For adjusting acceptance rate of p.init

  ### MCMC Sampling Process ###
  for (i in 1:n.sims) {
    if (i %% 100 == 0) {
      print("Verosimilitud:")
      print(log.likelihoodA)
      print(paste("Iteration:", i))
      flush.console()
    }
    
    paramB <- exp(log(paramA) + Q[i, ])  # Proposed values
    flush.console()
    ###################################
    # Update for emission params#
    ###################################
    
    emission_paramsB <- matrix(
      c(paramB[1:n_states], paramB[(n_states + 1):(2 * n_states)]),
      ncol = 2
    )
    
    log.likelihoodB <- loglikelihood(p.init = p.initA,
                               m.init = trans_matrix,
                               params_dist_emission = emission_paramsB,
                               params_dist_sojourn = duration_params,
                               n_states = n_states,
                               n = n,
                               y = Y)
    
    log.postB <- post(log.likelihoodB, emission_paramsB, duration_params, p.initA, trans_matrix)
    numer <- log.postB + sum(paramB[1:n_states]) + sum((n_states + 1):(2 * n_states))           
    denom <- log.postA + sum(paramA[1:n_states]) + sum((n_states + 1):(2 * n_states))        
    acceptance <- min(1, exp(numer - denom))
    #print("%%%%%%%%%%%%%%%%%%")
    #print("LOG.POSTB EMISSION")
    #print(log.postB)
    #print("%%%%%%%%%%%%%%%%%%")
    #print("%%%%%%%%%%%%%%%%%%")
    #print("PARAMS")
    #print(paramB[(n_states + 1):(2 * n_states)])
    #print("%%%%%%%%%%%%%%%%%%")
    if (log.postB > log.postA || U1[i] <= acceptance) {
      emission_params[, 1] <- emission_paramsB[, 1]
      paramA[1:n_states] <- paramB[1:n_states]
      paramA[(n_states + 1):(2 * n_states)] <- paramB[(n_states + 1):(2 * n_states)]
      accept1 <- accept1 + 1
      log.likelihoodA <- log.likelihoodB
      log.postA <- log.postB
    }
    
    #############################
    # Update for duration params #
    #############################

    duration_paramsB <- matrix(
      c(paramB[(2 * n_states + 1):(3 * n_states)], paramB[(3 * n_states + 1):(4 * n_states)]),
      ncol = 2
    )
    
    log.likelihoodB <- loglikelihood(p.init = p.initA,
                               m.init = trans_matrix,
                               params_dist_emission = emission_params,
                               params_dist_sojourn = duration_paramsB,
                               n_states = n_states,
                               n = n,
                               y = Y)
    
    log.postB <- post(log.likelihoodB, emission_params, duration_paramsB, p.initA, trans_matrix)
    numer <- log.postB + sum(paramB[(2 * n_states + 1):(3 * n_states)]) + sum(paramB[(3 * n_states + 1):(4 * n_states)])           
    denom <- log.postA + sum(paramA[(2 * n_states + 1):(3 * n_states)]) + sum(paramA[(3 * n_states + 1):(4 * n_states)])     
    
    acceptance <- min(1, exp(numer - denom))
    #print("%%%%%%%%%%%%%%%%%%")
    #print("LOG.POSTB SOJOURN")
    #print(log.postB)
    #print("%%%%%%%%%%%%%%%%%%")
    if (log.postB > log.postA || U2[i] <= acceptance) {
      duration_params[, 1] <- duration_paramsB[, 1]
      paramA[(2 * n_states + 1):(3 * n_states)] <- paramB[(2 * n_states + 1):(3 * n_states)]
      paramA[(3 * n_states + 1):(4 * n_states)] <- paramB[(3 * n_states + 1):(4 * n_states)]
      accept2 <- accept2 + 1
      log.likelihoodA <- log.likelihoodB
      log.postA <- log.postB
    }

    
    #####################
    # UPDATE FOR p.init and transition matrix using Dirichlet#
    #####################
    p.initB <- as.vector(rDirichlet(1, Z_p.init * p.initA))
    trans_matrixB <- trans_matrix
    
    for (state in 2:(n_states - 1)) {
      # Only update for intermediate states
      trans_probs <- rDirichlet(1, c(1, 1))  # Dirichlet for two transitions (to previous and next states)
      trans_matrixB[state, state - 1] <- trans_probs[1]
      trans_matrixB[state, state + 1] <- trans_probs[2]
      trans_matrixB[state, state] <- 0  # No self-transition
    }

    trans_matrixB[1, 2] <- 1
    trans_matrixB[n_states, n_states - 1] <- 1
    log.likelihoodB <- loglikelihood(p.init = p.initB,
                               m.init = trans_matrixB,
                               params_dist_emission = emission_params,
                               params_dist_sojourn = duration_params,
                               n_states = n_states,
                               n = n,
                               y = Y)

    log.postB <- post(log.likelihoodB, emission_params, duration_params, p.initB, trans_matrixB)
    numer <- log.postB           
    denom <- log.postA 

    numer <- log.postB  + log(ddirichlet(p.initA, Z_p.init * p.initB))
    denom <- log.postA  + log(ddirichlet(p.initB, Z_p.init * p.initA))

    
    for (state in 2:(n_states - 1)) {
    
       
        current_trans <- trans_matrix[state, c(state - 1, state + 1)]
        proposed_trans <- trans_matrixB[state, c(state - 1, state + 1)]
        
       
        log_prior_current_numer <- log(ddirichlet(current_trans, Z_trans_matrix * proposed_trans))
        
        log_prior_current_denom <- log(ddirichlet(proposed_trans, Z_trans_matrix * current_trans))
        

        numer <- numer + log_prior_current_numer
        denom <- denom + log_prior_current_denom
    }

    acceptance <- min(1, exp(numer - denom))
    
    #print("%%%%%%%%%%%%%%%%%%")
    #print("LOG.POSTB PROBS")
    #print(log.postB)
    #print("%%%%%%%%%%%%%%%%%%")
    if (!is.na(acceptance) && (log.postB > log.postA || U3[i] <= acceptance)) {
      p.initA <- p.initB
      trans_matrix <- trans_matrixB
      accept3 <- accept3 + 1
      log.likelihoodA <- log.likelihoodB
      log.postA <- log.postB
    }

    
    
    ##########################
    # Store the sampled values#
    ##########################

    sims_list[[i]] <- list(
    params_emission = emission_paramsB,
    params_sojourn = duration_params,
    trans_matrix = trans_matrix,
    p.init = p.initA,
    log.like = log.likelihoodB,
    log.pos = log.postB
  )

  }
  
  print("===== Muestreo Finalizado =====")
  print("Verosimilitud:")
  print(log.likelihoodA)
  return(list(sims_list = sims_list, acceptance_rates = c(accept1, accept2, accept3) / n.sims))
}
