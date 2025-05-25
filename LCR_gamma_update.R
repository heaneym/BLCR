gamma_update_cov_sel <- function(gamma,p,beta_prior_cov_inv, beta_prior_mean, kappa, G, X_current, omega, tau, X, A, beta_mean, beta_cov_inv ){
  if (p == 1){
    gamma_prop_index <- 2
  } else{
    gamma_prop_index <- sample(2:(p+1), 1) 
  }
  gamma_prop <- gamma
  if(gamma[gamma_prop_index] == 1){
    gamma_prop[gamma_prop_index] <- 0 
  } else {
    gamma_prop[gamma_prop_index] <- 1
  }
  X_current <- X[,gamma == 1, drop = FALSE]
  X_prop <- X[,gamma_prop == 1, drop = FALSE]
  beta_prior_cov_inv_prop <- beta_prior_cov_inv[gamma_prop == 1,gamma_prop == 1, drop = FALSE]
  beta_prior_cov_inv_current <- beta_prior_cov_inv[gamma == 1,gamma == 1, drop = FALSE]
  beta_prior_mean_prop <- beta_prior_mean[gamma_prop == 1, drop = FALSE]
  beta_prior_mean_current <- beta_prior_mean[gamma == 1, drop = FALSE]
  beta_cov_inv_prop <- array(0, dim = c(dim(beta_prior_cov_inv_prop),G-1))
  beta_cov_inv_current <- array(0, dim = c(dim(beta_prior_cov_inv_current),G-1))
  beta_cov_inv_chol_prop <- array(0, dim = c(dim(beta_prior_cov_inv_prop),G-1))
  beta_cov_inv_chol_current <- array(0, dim = c(dim(beta_prior_cov_inv_current),G-1))
  beta_mean_prop <- array(0, dim = c(length(beta_prior_mean_prop),G-1))
  beta_mean_current <- array(0, dim = c(length(beta_prior_mean_current),G-1))
  for (g in 1:(G-1)){
    #print(dim(X_prop))
    #print(dim(beta_prior_cov_inv_prop))
    beta_cov_inv_prop[,,g] <- crossprod(X_prop,omega[,g]*X_prop) + beta_prior_cov_inv_prop
    beta_cov_inv_current[,,g] <- crossprod(X_current,omega[,g]*X_current) + beta_prior_cov_inv_current
    beta_cov_inv_chol_prop[,,g] <- chol(beta_cov_inv_prop[,,g])
    beta_cov_inv_chol_current[,,g] <- chol(beta_cov_inv_current[,,g])
    #A <- kappa[,g] + omega[,g]*C[,g]
    fs_prop <- forwardsolve(t(beta_cov_inv_chol_prop[,,g]),crossprod(X_prop,A[,g]) + beta_prior_cov_inv_prop%*%beta_prior_mean_prop)
    fs_current <- forwardsolve(t(beta_cov_inv_chol_current[,,g]), crossprod(X_current,A[,g]) + beta_prior_cov_inv_current%*%beta_prior_mean_current)
    beta_mean_prop[,g] <- backsolve(beta_cov_inv_chol_prop[,,g], fs_prop) 
    beta_mean_current[,g] <- backsolve(beta_cov_inv_chol_current[,,g], fs_current)
  }
  beta_cov_inv_block_prop <- Matrix::bdiag(lapply(1:(G-1), function(g) beta_cov_inv_prop[,,g]))
  beta_cov_inv_block_current <- Matrix::bdiag(lapply(1:(G-1), function(g) beta_cov_inv_current[,,g]))
  beta_cov_inv_chol_block_prop <- Matrix::bdiag(lapply(1:(G-1), function(g) beta_cov_inv_chol_prop[,,g]))
  beta_cov_inv_chol_block_current <- Matrix::bdiag(lapply(1:(G-1), function(g) beta_cov_inv_chol_current[,,g]))
  beta_prior_cov_inv_block_prop <- Matrix::bdiag(replicate(G - 1, beta_prior_cov_inv_prop, simplify = FALSE))
  beta_prior_cov_inv_block_current <- Matrix::bdiag(replicate(G - 1, beta_prior_cov_inv_current, simplify = FALSE))
  beta_mean_prop_vec <- as.vector(beta_mean_prop)
  beta_mean_current_vec <- as.vector(beta_mean_current)
  beta_prior_mean_prop_vec <- rep(beta_prior_mean_prop, G - 1)
  beta_prior_mean_current_vec <- rep(beta_prior_mean_current, G-1)
  #print(beta_cov_inv_chol_block_current)
  #print(beta_cov_inv_chol_block_prop)
  log_det_prop <- -2*sum(log(Matrix::diag(beta_cov_inv_chol_block_prop)))
  log_det_current <- -2*sum(log(Matrix::diag(beta_cov_inv_chol_block_current)))
  log_det_prior_prop <- -2*sum(log(Matrix::diag(beta_prior_cov_inv_block_prop)))
  log_det_prior_current <- -2*sum(log(Matrix::diag(beta_prior_cov_inv_block_current)))
  #print(log_det_prop)
  #print(log_det_current)
  #print(log_det_prior_prop)
  #print(log_det_prior_current)
  if (is.na(log_det_prop)){print('log_det_prop is NA')}
  if (is.na(log_det_current)){print('log_det_current is NA')}
  if (is.na(log_det_prior_prop)){print('log_det_prior_prop is NA')}
  if (is.na(log_det_prior_current)){print('log_det_prior_current is NA')}
  log_det_sum <- log_det_prop - log_det_current - log_det_prior_prop + log_det_prior_current
  mean_var_contribution <- t(beta_mean_prop_vec)%*%(beta_cov_inv_block_prop%*%beta_mean_prop_vec) + t(beta_prior_mean_prop_vec)%*%(beta_prior_cov_inv_block_prop%*%beta_prior_mean_prop_vec) - t(beta_mean_current_vec)%*%(beta_cov_inv_block_current%*%beta_mean_current_vec) - t(beta_prior_mean_current_vec)%*%(beta_prior_cov_inv_block_current%*%beta_prior_mean_current_vec)
  #print(log_det_sum)
  log_accept_ratio <- (1/2)*(log_det_sum + mean_var_contribution)
  sum_gamma_prop <- sum(gamma_prop[2:(p+1)])
  sum_gamma_current <- sum(gamma[2:(p+1)])
  log_accept_ratio <- log_accept_ratio #+ (sum_gamma_prop-sum_gamma_current)*log(tau) + (sum_gamma_current-sum_gamma_prop)*log(1-tau)
  accept_ratio <- exp(log_accept_ratio)
  accept_ratio <- as.numeric(accept_ratio)
  #print(log_accept_ratio)
  #print(accept_ratio)
  if (runif(1) <= min(1,accept_ratio)){
    gamma <- gamma_prop
    X_current <- X_prop
    beta_mean <- beta_mean_prop
    beta_cov_inv <- beta_cov_inv_prop
    #print('accepted')
    #print(accept_ratio)
  }
  return(list(gamma = gamma, X_current = X_current, beta_mean = beta_mean, beta_cov_inv = beta_cov_inv))
}


gamma_update_no_sel <- function(gamma,p,beta_prior_cov_inv, beta_prior_mean, kappa, G, X_current, omega, tau, X, A, beta_mean, beta_cov_inv){
  return(list(gamma = gamma, X_current = X_current, beta_mean = beta_mean, beta_cov_inv = beta_cov_inv))
}

get_gamma_update <- function(cov.sel){
  if (cov.sel){
    gamma_update <- gamma_update_cov_sel
  } else {
    gamma_update <- gamma_update_no_sel
  }
  return(gamma_update)
}