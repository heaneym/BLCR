LCR_Gibbs <- function(X, Y, G, 
                      beta_prior_cov, beta_prior_mean,
                      theta_hyperparam, clust_var_prior,
                      item.sel = FALSE, cov.sel = FALSE, 
                      verbose = FALSE, relabel = TRUE,
                      n_samples = 1000, burnin = 500, 
                      thinby = 1,
                      sfm = FALSE, sparse_prior = 1/G, sparse_sigma_MH = 0.1,
                      init_seed = FALSE, beta_init = NULL,
                      z_init = NULL, theta_init = NULL,
                      nu_init = NULL, gamma_init = NULL,
                      K_specify = NULL){
  args <- match.call()
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (is.data.frame(Y)) {
    Y <- as.matrix(Y)
  }
  init <- LCR_init(X = X, Y = Y, G = G, 
                   beta_prior_mean = beta_prior_mean,
                   beta_prior_cov = beta_prior_cov,
                   theta_hyperparam = theta_hyperparam,
                   clust_var_prior = clust_var_prior,
                   K_specify = K_specify)
  list2env(init, envir = environment())
  if (init_seed){
    beta <- beta_init
    theta <- theta_init
    nu <- nu_init
    gamma <- gamma_init
    z <- z_init
  }
  X_current <- X_current[,which(gamma == 1), drop = FALSE]

  #want to create functions that initialise the sample arrays and select relevant functions
  z_update <- get_z_update_function(item.sel)
  impute_missing_values <- get_imputation_function(item.sel)
  nu_update <- get_nu_update(item.sel)
  S_update <- get_S_update(item.sel)
  theta_update <- get_theta_update(item.sel)
  gamma_update <- get_gamma_update(cov.sel)
  log_post_compute <- get_log_post_function(item.sel)
  n_iter <- thinby*n_samples + burnin
  sample_arrays <- LCR_init_sample_arrays(beta = beta, n_samples = n_samples, omega = omega,
                                          w = w, z = z, theta = theta, gamma = gamma, N_g = N_g,
                                          N_gjk = N_gjk, nu = nu, Y)
  list2env(sample_arrays, envir = environment())
  G_eff_samples <- numeric(n_samples) 
  sparse_weights_accept <- 0   
  if (sfm) relabel <- FALSE
  sample_count <- 0
  start_time <- Sys.time()
  progress_interval <- 500
  for (count in 1:n_iter){
    N_up <- N_updates(z = z, Y_indicator = Y_indicator, Y_missing_indicator = Y_missing_indicator, N_jk = N_jk, nu = nu)
    N_g <- N_up$N_g
    N_gjk <- N_up$N_gjk
    Y_current <- impute_missing_values(Y_current = Y_current, Y_missing_indicator = Y_missing_indicator, z = z, theta = theta, nu = nu, N_gjk = N_up$N_gjk_array, N_g = N_g, theta_hyperparam = theta_hyperparam, M = M, K = K)
    Y_indicator <- 1*array(outer(Y_current, 1:max_K, "=="), dim = dim_Y_indicator)
    S <- S_update(S = S, Y_indicator = Y_indicator, Y_missing_indicator = Y_missing_indicator, z = z, M = M, G = G, K = K)
    theta <- theta_update(S = S, theta_hyperparam = theta_hyperparam, K = K, G = G, M = M, theta = theta)
    mu <- mu_update(X = X, beta = beta)
    exp_mu <- exp(mu)
    nu <- nu_update(nu = nu, M = M, K = K, G = G, theta_hyperparam = theta_hyperparam, N_g = N_g, N_gjk = N_up$N_gjk_array, Y_missing_indicator = Y_missing_indicator, inclusion_sum = inclusion_sum, exclusion_sum = exclusion_sum, count = count, z = z)
    C <- C_update(mu = mu, exp_mu = exp_mu, G = G, C = C)
    eta <- eta_update(mu = mu, C = C)
    logit_probs <- logit_probs_update(eta = eta)
    log_logit_probs <- log(logit_probs)
    log_theta <- log(theta)
    z_params <- list(
      z = z, nu = nu, mu = mu, K = K, theta_hyperparam = theta_hyperparam,
      N_gjk = N_gjk, N_g = N_g, Y_indicator = Y_indicator, n = n, G = G,
      omega = omega, C = C, p = p, beta = beta, beta_cov_inv = beta_cov_inv,
      beta_cov_inv_chol = beta_cov_inv_chol, beta_prior_mean = beta_prior_mean,
      beta_prior_cov_inv = beta_prior_cov_inv, gamma = gamma,
      X_current = X_current,
      
      log_logit_probs = log_logit_probs, log_theta = log_theta,
      Y = Y, M = M, Y_missing_indicator = Y_missing_indicator
    )
    z_up <- z_update(z_params)
    list2env(z_up, envir = environment())
    kappa <- kappa_update(z = z)
    omega <- omega_update(eta = eta, mu = mu, G = G, n = n, omega = omega)
    A <- A_update(kappa = kappa, omega = omega, C = C, G = G)
    gamma_up <- gamma_update(gamma = gamma, p = p, beta_prior_cov_inv = beta_prior_cov_inv, beta_prior_mean = beta_prior_mean, kappa = kappa, G = G, omega = omega, tau = tau, X = X, A = A, beta_mean = beta_mean, beta_cov_inv = beta_cov_inv, X_current = X_current)
    list2env(gamma_up, envir = environment())
    if (sfm) {
      beta_up <- beta_update_sfm(beta = beta, X_current = X_current, G = G, p = p,
                                 beta_prior_cov_inv = beta_prior_cov_inv, A = A,
                                 beta_prior_mean = beta_prior_mean, omega = omega)
      list2env(beta_up, envir = environment())
      mu   <- mu_update(X = X_current, beta = beta)    
      sparse_weights_up <- sparse_weights_update(beta = beta, mu = mu, z = z, G = G, sparse_prior = sparse_prior, sparse_sigma_MH = sparse_sigma_MH, n = n)
      beta <- sparse_weights_up$beta
      sparse_weights_accept <- sparse_weights_accept + sparse_weights_up$accepted
      if (count <= burnin)   #adaptive update for weight tuning                            
        #sparse_sigma_MH <- exp(log(sparse_sigma_MH) + (sparse_weights_up$accepted - 0.30) / sqrt(count))
        sparse_sigma_MH <- exp(log(sparse_sigma_MH) + (sparse_weights_up$accept_prob - 0.30) / sqrt(count))
    } else {
      beta_up <- beta_update(gamma = gamma, X_current = X_current, G = G, p = p,
                             beta_prior_cov_inv = beta_prior_cov_inv, A = A,
                             beta_prior_mean = beta_prior_mean, omega = omega)
      list2env(beta_up, envir = environment())
    }
    log_post_params <- list(
      mu = mu, Y_indicator = Y_indicator, 
      log_theta = log_theta, z = z, 
      theta_hyperparam = theta_hyperparam, 
      beta_prior_mean = beta_prior_mean, 
      beta_prior_cov_inv = beta_prior_cov_inv, 
      beta = beta, nu = nu, M = M, 
      clust_var_prior = clust_var_prior,
      K = K, N_jk = N_jk, N_gjk = N_gjk, 
      N_g = N_g, n = n, G = G
    )
    log_post_like <- log_post_compute(log_post_params)
    list2env(log_post_like, envir = environment())
    if (count > burnin && ((count - burnin) %% thinby == 0)) {
      sample_count <- sample_count + 1
      beta_samples[,,sample_count] <- beta
      omega_samples[,,sample_count] <- omega
      z_samples[,,sample_count] <- z
      w_samples[,,sample_count] <- w
      theta_samples[,,,sample_count] <- theta
      nu_samples[,sample_count] <- nu
      gamma_samples[,sample_count] <- gamma
      N_g_samples[,sample_count] <- N_g
      N_gjk_samples[,,,sample_count] <- N_gjk
      log_post_samples[sample_count] <- log_post
      log_like_samples[sample_count] <- log_like
      Y_imputed_samples[,,sample_count] <- Y_current
      G_eff_samples[sample_count] <- sum(colSums(z) > 0)
    }
    if (verbose && count %% progress_interval == 0 && count > burnin) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      avg_time_per_sample <- elapsed / sample_count
      est_remaining <- round((n_samples - sample_count) * avg_time_per_sample, 1)
      percent <- round(100 * sample_count / n_samples)
      bar_width <- 40
      n_hashes <- floor(bar_width * sample_count / n_samples)
      bar <- paste0("[", paste0(rep("=", n_hashes), collapse = ""), 
                    paste0(rep(" ", bar_width - n_hashes), collapse = ""), "]")
      
      line <- sprintf("\r%s %3d%% | %d of %d samples | ~%.1f sec remaining", 
                      bar, percent, sample_count, n_samples, est_remaining)
      cat(sprintf("%-*s", 80, line))  
      flush.console()
      #print(sum(colSums(z) > 0))
    }
  }
  if (verbose) {
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("\nSampling complete! Total time: %.1f seconds (%.1f minutes)\n", 
                total_time, total_time / 60))
  }
  if (relabel){
    cat("relabelling...")
    relabelled_samples <- relabel_outputs(beta_samples = beta_samples, z_samples = z_samples, w_samples = w_samples, theta_samples = theta_samples, N_g_samples = N_g_samples, N_gjk_samples = N_gjk_samples, n_samples = n_samples, n = n, G = G)
    list2env(relabelled_samples, envir = environment())
  }
  gating <- extract_gating(beta_samples = beta_samples, G = G)
  
  post_process <- LCR_post_process(N_g_samples = N_g_samples,
                                   z_samples = z_samples,
                                   beta_samples = beta_samples, 
                                   theta_samples = theta_samples,
                                   nu_samples = nu_samples,
                                   gamma_samples = gamma_samples,
                                   N_jk = N_jk,
                                   Y_indicator = Y_indicator,
                                   X = X, 
                                   beta_prior_mean = beta_prior_mean, 
                                   beta_prior_cov_inv = beta_prior_cov_inv,
                                   M = M, K = K, n = n,
                                   log_like_samples = log_like_samples, 
                                   p = p, G = G, 
                                   log_post_compute = log_post_compute,
                                   item.sel = item.sel,
                                   n_samples = n_samples,
                                   N_gjk_samples = N_gjk_samples,
                                   theta_hyperparam = theta_hyperparam,
                                   clust_var_prior = clust_var_prior,
                                   Y_imputed_samples = Y_imputed_samples,
                                   Y_missing_indicator = Y_missing_indicator,
                                   Y = Y)
  list2env(post_process, envir = environment())
  samples <- list(beta_samples = beta_samples, 
                  omega_samples = omega_samples,
                  z_samples = z_samples,
                  w_samples = w_samples,
                  theta_samples = theta_samples,
                  nu_samples = nu_samples,
                  gamma_samples = gamma_samples,
                  N_g_samples = N_g_samples,
                  N_gjk_samples = N_gjk_samples,
                  pi_samples = pi_samples,
                  log_post_samples = log_post_samples,
                  log_like_samples = log_like_samples,
                  Y_imputed_samples = Y_imputed_samples)
  result <- list()
  result$samples <- samples
  result$pi <- pi_estimate
  result$pi.sd <- pi_sd
  result$itemprob <- theta_estimate
  result$itemprob.sd <- theta_sd
  result$assignment_prob <- z_estimate
  result$assignment_sd <- z_sd
  result$assignment_entropy <- z_entropy
  result$beta_estimate <- beta_estimate
  result$beta_sd <- beta_sd
  result$item_inclusion_prob <- nu_estimate
  result$item_inclusion_sd <- nu_sd
  result$cov_inclusion_prob <- gamma_estimate
  result$cov_inclusion_sd <- gamma_sd
  result$Z <- Z
  result$logpost <- log_post
  result$loglike <- log_like
  result$DIC <- DIC
  result$AICM <- AIC
  result$BICM <- BIC
  result$iter <- n_iter
  result$burn.in <- burnin
  result$thin.by <- thinby
  result$n_samples <- n_samples
  result$relabel <- relabel
  result$item.sel <- item.sel
  result$cov.sel <- cov.sel
  result$cov.ind <- cov_ind
  result$item.ind <- item_ind
  result$N_gjk_estimate <- N_gjk_estimate
  result$N_g_estimate <- N_g_estimate
  result$Y_imputed_mode <- Y_imputed_mode  
  result$Y_imputed_probs <- Y_imputed_probs  
  result$Y_imputed_max_prob <- Y_imputed_max_prob  
  result$Y_imputed_entropy <- Y_imputed_entropy  
  result$Y_missing_indicator <- Y_missing_indicator
  result$Y <- Y
  result$X <- X
  
  result$G_eff_samples <- G_eff_samples
  result$G_eff_posterior <- summarise_G_eff(N_g_samples, n)
  result$sparse_weights_accept_rate <- sparse_weights_accept/n_iter
  
  
  
  result$gating_intercepts <- gating$intercept_samples   
  result$gating_slopes     <- gating$slope_samples       
  result$weights           <- gating$weight_samples      
  result$weight_profile    <- gating$weight_ordered_mean
  
  return(result)
}
