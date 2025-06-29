LCR_post_process <- function(N_g_samples, z_samples, beta_samples, 
                             theta_samples, nu_samples,
                             gamma_samples, N_jk, Y_indicator,
                             X, beta_prior_mean, beta_prior_cov_inv,
                             M, K, n, log_like_samples, 
                             p, G, log_post_compute, item.sel, 
                             n_samples, N_gjk_samples, 
                             theta_hyperparam, clust_var_prior){
  
  if(item.sel){
    theta_estimate <- array(0, dim = c(G, M, max(K)))
    theta_var <- array(0, dim = c(G, M, max(K)))
    for (g in 1:G){
      for (j in 1:M){
        for (k in 1:K[j]){
          theta_estimate[g,j,k] <- (1/n_samples)*(sum((N_gjk_samples[g,j,k,] + theta_hyperparam)/(N_g_samples[g,] + K[j]*theta_hyperparam)))
          var_term1_summand_num <- (N_gjk_samples[g,j,k,] + theta_hyperparam)*(N_g_samples[g,] + (K[j] - 1)*theta_hyperparam - N_gjk_samples[g,j,k,])
          var_term1_summand_denom <- ((N_g_samples[g,] + K[j]*theta_hyperparam)^2)*(N_g_samples[g,] + K[j]*theta_hyperparam + 1)
          var_term1 <- (1/n_samples)*sum(var_term1_summand_num/var_term1_summand_denom)
          var_term2_summand <- (N_gjk_samples[g,j,k,] + theta_hyperparam)/(N_g_samples[g,] + K[j]*theta_hyperparam) - theta_estimate[g,j,k]
          theta_var[g,j,k] <- var_term1 + (1/n_samples)*sum((var_term2_summand)^2)
        }
      }
    }
    theta_sd <- sqrt(theta_var)
    #finish this!!!
  } else {
    theta_estimate <- apply(theta_samples, c(1,2,3), mean)
    theta_sd <- apply(theta_samples, c(1,2,3), sd)
  }
  pi_samples <- N_g_samples/colSums(N_g_samples)
  pi_estimate <- apply(pi_samples, 1, mean)
  pi_sd <- apply(pi_samples, 1, sd)
  Z <- apply(z_samples, c(1,2), mean)
  z_sd <- apply(z_samples, c(1,2), sd)
  z_estimate_temp <- apply(Z,1,which.max)
  one_hot <- matrix(0, nrow = n, ncol = G)
  one_hot[cbind(1:n,z_estimate_temp)] <- 1
  z_estimate <- one_hot
  z_entropy <- -rowSums(Z * log(Z + 1e-12))
  beta_estimate <- apply(beta_samples, c(1,2), mean)
  beta_sd <- apply(beta_samples, c(1,2), sd)
  log_theta_estimate <- log(theta_estimate)
  nu_estimate <- apply(nu_samples, 1, mean)
  nu_estimate_hard <- as.integer(nu_estimate > 0.5)
  nu_sd <- apply(nu_samples, 1, sd)
  gamma_estimate <- apply(gamma_samples, 1, mean)
  gamma_estimate_hard <- as.integer(gamma_estimate>0.5)
  gamma_sd <- apply(gamma_samples, 1, sd)
  N_up_estimate <- N_updates(z = z_estimate, Y_indicator = Y_indicator, N_jk = N_jk, nu = nu_estimate_hard)
  N_gjk_estimate <- N_up_estimate$N_gjk
  N_g_estimate <- N_up_estimate$N_g
  mu_estimate <- X%*%beta_estimate
  
  #do I need to include the log like/log post function as an argument?
  
  log_post_params_estimate <- list(
    mu = mu_estimate, Y_indicator = Y_indicator, 
    log_theta = log_theta_estimate, z = z_estimate, 
    theta_hyperparam = theta_hyperparam, 
    beta_prior_mean = beta_prior_mean, 
    beta_prior_cov_inv = beta_prior_cov_inv, 
    beta = beta_estimate, nu = nu_estimate_hard, M = M, 
    clust_var_prior = clust_var_prior,
    K = K, N_jk = N_jk, N_gjk = N_gjk_estimate, 
    N_g = N_g_estimate, n = n, G = G
  )
  
  log_post_like_estimate <- log_post_compute(log_post_params_estimate) 
  
  log_post_estimate <- log_post_like_estimate$log_post
  log_like_estimate <- log_post_like_estimate$log_like
  
  num_params <- (G-1)*(p+1) + G*(sum(K - 1))
  AIC <- -2*log_like_estimate + 2*num_params
  BIC <- -2*log_like_estimate + num_params*log(n)
  
  mean_log_like <- mean(log_like_samples)
  mean_deviance <- -2*mean_log_like
  deviance_estimate <- -2*log_like_estimate
  DIC <- 2*mean_deviance - deviance_estimate
  
  
  theta_list <- list()
  theta_sd_list <- list()
  for (j in 1:M){
    theta_list[[j]] <- matrix(0,nrow = G, ncol = K[j])
    theta_sd_list[[j]] <- matrix(0,nrow = G, ncol = K[j])
    for (g in 1:G){
      for (k in 1:K[j]){
        theta_list[[j]][g,k] <- theta_estimate[g,j,k]
        theta_sd_list[[j]][g,k] <- theta_sd[g,j,k]
      }
    }
  }
  
  #Going to permute the indices so that g indices are sorted by the group proportions
  descend_perm <- order(pi_estimate, decreasing = TRUE)
  pi_samples_perm <- pi_samples[descend_perm, , drop = FALSE]
  pi_estimate_perm <- pi_estimate[descend_perm]
  pi_sd_perm <- pi_sd[descend_perm]
  Z_perm <- Z[,descend_perm]
  z_sd_perm <- z_sd[,descend_perm]
  z_estimate_perm <- z_estimate[,descend_perm]
  beta_estimate_perm <- beta_estimate[,descend_perm]
  beta_sd_perm <- beta_sd[,descend_perm]
  theta_list_perm <- lapply(theta_list, function(mat) mat[descend_perm, , drop = FALSE])
  theta_row_names <- paste0("Group ", seq_len(G))
  theta_list_named <- lapply(theta_list_perm, function(mat) {
    rownames(mat) <- theta_row_names
    mat
  })
  theta_sd_list_perm <- lapply(theta_sd_list, function(mat) mat[descend_perm, , drop = FALSE])
  theta_sd_list_named <- lapply(theta_sd_list_perm, function(mat) {
    rownames(mat) <- theta_row_names
    mat
  })
  N_gjk_estimate_perm <- N_gjk_estimate[descend_perm,,]
  N_g_estimate_perm <- N_g_estimate[descend_perm]
  
  
  outputs <- list(
    pi_samples = pi_samples_perm,
    pi_estimate = pi_estimate_perm,
    pi_sd = pi_sd_perm,
    Z = Z_perm,
    z_sd = z_sd_perm,
    z_estimate = z_estimate_perm,
    z_entropy = z_entropy,
    beta_estimate = beta_estimate_perm,
    beta_sd = beta_sd_perm,
    theta_estimate = theta_list_named,
    theta_sd = theta_sd_list_named,
    nu_estimate = nu_estimate,
    item_ind = nu_estimate_hard,
    nu_sd = nu_sd,
    gamma_estimate = gamma_estimate,
    gamma_sd = gamma_sd,
    cov_ind = gamma_estimate_hard,
    N_gjk_estimate = N_gjk_estimate_perm,
    N_g_estimate = N_g_estimate_perm,
    log_post = log_post_estimate,
    log_like = log_like_estimate,
    AIC = AIC,
    BIC = BIC,
    DIC = DIC
  )
  return(outputs)
}




