# S_update_uncollapsed <- function(Y_indicator,z,M,G,K,S){
#   S <- array(0, dim = c(G, M, max(K)))
#   for (k in 1:max(K)) {
#     S[,,k] <- t(z)%*%Y_indicator[,,k]
#   }
#   return(S)
# }

S_update_uncollapsed <- function(Y_indicator, Y_missing_indicator, z, M, G, K, S) {
  S <- array(0, dim = c(G, M, max(K)))
  
  # Only accumulate counts for observed values
  for (k in 1:max(K)) {
    # Multiply indicator by (1 - missing_indicator) to zero out missing values
    observed_indicator <- Y_indicator[,,k] * (1 - Y_missing_indicator)
    S[,,k] <- t(z) %*% observed_indicator
  }
  
  return(S)
}


# theta_update_uncollapsed <- function(S,theta_hyperparam,K,G,M,theta){
#   alpha_S <- array(theta_hyperparam, dim = dim(S)) + S
#   gam_samples <- array(rgamma(prod(dim(S)), shape = as.vector(alpha_S), scale = 1),dim = dim(S))
#   gam_sums <- apply(gam_samples, c(1, 2), sum)
#   theta <- sweep(gam_samples, MARGIN = c(1, 2), gam_sums, FUN = "/")
#   return(theta)
# }


theta_update_uncollapsed <- function(S, theta_hyperparam, K, G, M, theta) {
  alpha_S <- array(theta_hyperparam, dim = dim(S)) + S
  
  # # Diagnostic checks
  # cat("Checking alpha_S:\n")
  # cat("Range of S:", range(S, na.rm = TRUE), "\n")
  # cat("theta_hyperparam:", theta_hyperparam, "\n")
  # cat("Range of alpha_S:", range(alpha_S, na.rm = TRUE), "\n")
  # cat("Any NA in S:", any(is.na(S)), "\n")
  # cat("Any NA in alpha_S:", any(is.na(alpha_S)), "\n")
  # cat("Any Inf in alpha_S:", any(is.infinite(alpha_S)), "\n")
  # cat("Any negative in alpha_S:", any(alpha_S < 0, na.rm = TRUE), "\n")
  # 
  # # Check for specific problematic values
  # if (any(is.na(alpha_S))) {
  #   cat("NA positions:\n")
  #   print(which(is.na(alpha_S), arr.ind = TRUE))
  # }
  # if (any(is.infinite(alpha_S))) {
  #   cat("Inf positions:\n")
  #   print(which(is.infinite(alpha_S), arr.ind = TRUE))
  # }
  
  gam_samples <- array(0, dim = dim(S))
  
  # Only generate gamma samples for valid categories
  for (j in 1:M) {
    for (k in 1:K[j]) {
      # Add safety check before rgamma
      shapes <- alpha_S[, j, k]
      if (any(is.na(shapes)) || any(is.infinite(shapes)) || any(shapes <= 0)) {
        cat("Problem at j=", j, ", k=", k, "\n")
        cat("shapes:", shapes, "\n")
        # Use small positive value as fallback
        shapes[is.na(shapes) | is.infinite(shapes) | shapes <= 0] <- theta_hyperparam
      }
      gam_samples[, j, k] <- rgamma(G, shape = shapes, scale = 1)
    }
  }
  
  # Normalize within each (g, j) only over valid categories
  theta <- array(0, dim = dim(S))
  for (j in 1:M) {
    gam_sums <- rowSums(gam_samples[, j, 1:K[j], drop = FALSE])
    for (k in 1:K[j]) {
      theta[, j, k] <- gam_samples[, j, k] / gam_sums
    }
  }
  
  return(theta)
}




nu_update_uncollapsed <- function(nu, M, K, G, theta_hyperparam, N_g, N_gjk, 
                                  Y_missing_indicator, inclusion_sum, exclusion_sum, count, z){
  return(nu)
}

get_S_update <- function(item.sel){
  if (item.sel){
    S_update <- S_update_collapsed
  } else {
    S_update <- S_update_uncollapsed
  }
  return(S_update)
}

get_theta_update <- function(item.sel){
  if (item.sel){
    theta_update <- theta_update_collapsed
  } else {
    theta_update <- theta_update_uncollapsed
  }
  return(theta_update)
}


impute_missing_values_uncollapsed <- function(Y_current, Y_missing_indicator, z, theta, nu, N_gjk,N_g,theta_hyperparam, M, K) {
  n <- nrow(Y_current)
  Y_imputed <- Y_current
  
  for (i in 1:n) {
    # Get assigned class for this observation
    class_i <- which(z[i,] == 1)[1]  # Take first if multiple (shouldn't happen)
    
    for (j in 1:M) {
      if (Y_missing_indicator[i, j] == 1) {  # If missing
        # Extract probabilities more explicitly
        probs <- as.numeric(theta[class_i, j, 1:K[j]])
        
        # Ensure valid probability vector
        if (any(is.na(probs)) || any(probs < 0)) {
          cat("Invalid probs at i=", i, ", j=", j, ", class=", class_i, "\n")
          probs <- rep(1/K[j], K[j])  # Use uniform as fallback
        }
        
        # Normalize to ensure sum to 1
        probs <- probs / sum(probs)
        
        Y_imputed[i, j] <- sample(1:K[j], size = 1, prob = probs)
      }
    }
  }
  
  return(Y_imputed)
}





