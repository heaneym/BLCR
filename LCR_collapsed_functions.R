# N_updates <- function(z, Y_indicator, N_jk, nu){
#   N_g <- colSums(z)
#   N_gjk <- apply(Y_indicator, c(2, 3), function(S_jk) t(z) %*% S_jk)
#   return(list(N_g = N_g, N_gjk = N_gjk))
# }

# N_updates <- function(z, Y_indicator, Y_missing_indicator, N_jk, nu){
#   N_g <- colSums(z)
#   N_gjk <- apply(Y_indicator, c(2, 3), function(S_jk) {
#     observed_S_jk <- S_jk * (1 - Y_missing_indicator)
#     t(z) %*% observed_S_jk
#   })
#   
#   return(list(N_g = N_g, N_gjk = N_gjk))
# }

# N_updates <- function(z, Y_indicator, Y_missing_indicator, N_jk, nu){
#   N_g <- colSums(z)
#   
#   # Get dimensions
#   G <- ncol(z)
#   M <- dim(Y_indicator)[2]
#   max_K <- dim(Y_indicator)[3]
#   
#   # Initialize array with correct dimensions: G x M x max_K
#   N_gjk <- array(0, dim = c(G, M, max_K))
#   
#   # Compute counts manually to ensure correct dimensions
#   for (j in 1:M) {
#     for (k in 1:max_K) {
#       # Get indicator for item j, category k
#       S_jk <- Y_indicator[, j, k]
#       
#       # Zero out missing entries
#       observed_S_jk <- S_jk * (1 - Y_missing_indicator[, j])
#       
#       # Compute counts per class: gives vector of length G
#       N_gjk[, j, k] <- as.vector(t(z) %*% observed_S_jk)
#     }
#   }
#   
#   # Convert to vector for C++
#   N_gjk_vec <- as.vector(N_gjk)
#   
#   return(list(N_g = N_g, N_gjk = N_gjk_vec))
# }

N_updates <- function(z, Y_indicator, Y_missing_indicator, N_jk, nu){
  N_g <- colSums(z)
  

  G <- ncol(z)
  M <- dim(Y_indicator)[2]
  max_K <- dim(Y_indicator)[3]
  

  N_gjk_array <- array(0, dim = c(G, M, max_K))
  

  for (j in 1:M) {
    for (k in 1:max_K) {
      S_jk <- Y_indicator[, j, k]
      observed_S_jk <- S_jk * (1 - Y_missing_indicator[, j])
      N_gjk_array[, j, k] <- as.vector(t(z) %*% observed_S_jk)
    }
  }
  

  return(list(
    N_g = N_g, 
    N_gjk = as.vector(N_gjk_array),  # Vector for C++
    N_gjk_array = N_gjk_array        # Array for R functions
  ))
}





# nu_update_collapsed <- function(nu, M, K, G, theta_hyperparam, N_g, N_gjk, inclusion_sum, exclusion_sum, count){
#   nu_prop <- nu
#   j_prop <- sample(1:M,1)
#   if(nu[j_prop] == 1){
#     nu_prop[j_prop] <- 0
#     log_gamma_N_alpha <- lgamma(N_gjk[,j_prop,] + theta_hyperparam)
#     log_sum1 <- sum(log_gamma_N_alpha)
#     diff_log <- sum(lgamma(N_g + K[j_prop]*theta_hyperparam)) - log_sum1
#     log_accept_ratio <- exclusion_sum[j_prop] + diff_log
#   } else {
#     nu_prop[j_prop] <- 1
#     log_gamma_N_alpha <- lgamma(N_gjk[,j_prop,] + theta_hyperparam)
#     log_sum1 <- sum(log_gamma_N_alpha)
#     diff_log <- log_sum1 - sum(lgamma(N_g + K[j_prop]*theta_hyperparam))
#     log_accept_ratio <- inclusion_sum[j_prop] + diff_log
#     
#   }
#   accept_ratio <- exp(log_accept_ratio)
#   accept_ratio <- as.numeric(accept_ratio)
#   if(runif(1) < min(1,accept_ratio)){
#     nu <- nu_prop
#   }
#   return(nu)
# }


nu_update_collapsed <- function(nu, M, K, G, theta_hyperparam, N_g, N_gjk,
                                Y_missing_indicator, inclusion_sum, exclusion_sum, count, z){
  nu_prop <- nu
  j_prop <- sample(1:M, 1)

  # Compute N_g,obs(j) = number of observations in each class with observed value for item j
  N_g_obs_j <- colSums((1 - Y_missing_indicator[, j_prop]) * z)

  if(nu[j_prop] == 1){
    # Proposing to exclude variable j_prop
    nu_prop[j_prop] <- 0

    log_gamma_N_alpha <- lgamma(N_gjk[, j_prop, 1:K[j_prop]] + theta_hyperparam)
    log_sum1 <- sum(log_gamma_N_alpha)

    # Use observed counts only in denominator
    diff_log <- sum(lgamma(N_g_obs_j + K[j_prop]*theta_hyperparam)) - log_sum1
    log_accept_ratio <- exclusion_sum[j_prop] + diff_log

  } else {
    # Proposing to include variable j_prop
    nu_prop[j_prop] <- 1

    log_gamma_N_alpha <- lgamma(N_gjk[, j_prop, 1:K[j_prop]] + theta_hyperparam)
    log_sum1 <- sum(log_gamma_N_alpha)

    # Use observed counts only in denominator
    diff_log <- log_sum1 - sum(lgamma(N_g_obs_j + K[j_prop]*theta_hyperparam))
    log_accept_ratio <- inclusion_sum[j_prop] + diff_log
  }

  accept_ratio <- exp(log_accept_ratio)
  accept_ratio <- as.numeric(accept_ratio)

  if(runif(1) < min(1, accept_ratio)){
    nu <- nu_prop
  }

  return(nu)
}







# impute_missing_values_collapsed <- function(Y_current, Y_missing_indicator, z, theta, nu, N_gjk,N_g,theta_hyperparam, M, K) {
#   n <- nrow(Y_current)
#   Y_imputed <- Y_current
#   
#   for (i in 1:n) {
#     class_i <- which(z[i,] == 1)[1]
#     
#     for (j in 1:M) {
#       if (nu[j] == 1 && Y_missing_indicator[i, j] == 1) {
#         
#         N_g_obs_j <- sum(z[, class_i] * (1 - Y_missing_indicator[, j]))
#         
#         probs <- (N_gjk[class_i, j, 1:K[j]] + theta_hyperparam) / 
#           (N_g_obs_j + K[j] * theta_hyperparam)
#         
#         probs <- probs / sum(probs)
#         
#         Y_imputed[i, j] <- sample(1:K[j], size = 1, prob = probs)
#       }
#     }
#   }
#   
#   return(Y_imputed)
# }


impute_missing_values_collapsed <- function(Y_current, Y_missing_indicator, z, nu, 
                                            N_gjk, N_g, theta_hyperparam, M, K, theta) {
  n <- nrow(Y_current)
  Y_imputed <- Y_current
  
  for (i in 1:n) {
    class_i <- which(z[i,] == 1)[1]
    
    for (j in 1:M) {
      if (nu[j] == 1 && Y_missing_indicator[i, j] == 1) {
        
        # N_g_obs_j EXCLUDING observation i
        N_g_obs_j_minus_i <- sum(z[-i, class_i] * (1 - Y_missing_indicator[-i, j]))
        
        # N_gjk EXCLUDING observation i (but this is already computed without i 
        # since i has missing value for j, it wasn't counted in N_gjk)
        # So just use N_gjk directly
        probs <- (N_gjk[class_i, j, 1:K[j]] + theta_hyperparam) / 
          (N_g_obs_j_minus_i + K[j] * theta_hyperparam)
        
        probs <- probs / sum(probs)
        
        Y_imputed[i, j] <- sample(1:K[j], size = 1, prob = probs)
      }
    }
  }
  
  return(Y_imputed)
}



get_imputation_function <- function(item.sel){
  if (item.sel){
    impute_missing_values <- impute_missing_values_collapsed
  } else {
    impute_missing_values <- impute_missing_values_uncollapsed
  }
  return(impute_missing_values)
}



get_nu_update <- function(item.sel){
  if (item.sel){
    nu_update <- nu_update_collapsed
  } else {
    nu_update <- nu_update_uncollapsed
  }
  return(nu_update)
}

S_update_collapsed <- function(Y_indicator,Y_missing_indicator,z,M,G,K,S){
  return(S)
}

theta_update_collapsed <- function(S,theta_hyperparam,K,G,M,theta){
  return(theta)
}