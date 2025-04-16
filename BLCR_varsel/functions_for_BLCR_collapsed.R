initialise_variables_BLCR_collapsed <- function(G, X, Y, beta_prior_var_factor = 10, beta_prior_mean, clust_var_prior = 0.5, alpha  = 1){
  if (is.vector(X)){
    ones <- rep(1,length(X))
    X <- cbind(ones,X)
  } else{
    ones <- rep(1,nrow(X))
    X <- cbind(ones,X)
  }
  
  p <- ncol(X) - 1
  M <- ncol(Y)
  n <- nrow(Y)
  K <- apply(Y, 2, function(x) length(unique(x)))
  
  
  #-----------------------------------------------------------------------------
  #-------------------------beta initialisation---------------------------------
  #-----------------------------------------------------------------------------
  
  beta <- cbind(matrix(rnorm((G-1)*(p+1)), nrow = (p+1)),0)
  mu <- X%*%beta
  logit_probs <- t(apply(matrix(runif(n * G), nrow = n), 1, function(x) x / sum(x)))
  log_logit_probs <- log(logit_probs)
  
  
  #-----------------------------------------------------------------------------
  #-------------------------omega, Omega,...initialisation here-----------------
  #-----------------------------------------------------------------------------
  
  C <- matrix(rnorm((n)*(G)), nrow = n)
  omega <- matrix(rpg(n*(G-1)), nrow = n)
  eta <- matrix(0, nrow = n, ncol = G)
  
  beta_cov <- array(0, dim = c(p+1,p+1,G-1))
  beta_cov_inv <- array(0, dim = c(p+1,p+1,G-1))
  beta_cov_inv_chol <- array(0, dim = c(p+1,p+1,G-1))
  beta_mean <- matrix(rnorm((G-1)*(p+1)), ncol = (G-1), nrow = (p+1) )
  
  
  #-----------------------------------------------------------------------------
  #-------------------------z initialisation------------------------------------
  #-----------------------------------------------------------------------------
  
  w <- array(0, dim = c(n,G)) 
  for (i in 1:n) {
    # Randomly assign initial group
    w[i,] <- rep(1/G,G)
  }
  
  z <- array(0, dim = c(n, G))
  # Use K-means for better initial clustering
  km <- kmeans(X, G)
  for (i in 1:n) {
    z[i, km$cluster[i]] <- 1
  }
  kappa <- z - 1/2
  
  #-----------------------------------------------------------------------------
  #-------------------------gamma, tau initialisation---------------------------
  #-----------------------------------------------------------------------------
  
  
  gamma <- rep(1,p+1)
  tau <- runif(1)
  nu <- rep(1,M)
  
  
  
  
  
  
  Y_indicator <- 1*array(outer(Y, 1:max(K), "=="), dim = c(nrow(Y), ncol(Y), length(1:max(K))))
  N_jk <- colSums(Y_indicator, dims = 1)
  
  
  sum1_inclusion <- (G-1)*(lgamma(K*alpha) - K*lgamma(alpha))
  log_gamma_N_jk_alpha <- lgamma(N_jk + alpha)
  sum2_inclusion <- lgamma(n + K*alpha) - rowSums(log_gamma_N_jk_alpha)
  sum3_inclusion <- log(clust_var_prior) - log(1-clust_var_prior)
  inclusion_sum <- sum1_inclusion + sum2_inclusion + sum3_inclusion
  exclusion_sum <- -inclusion_sum
  
  
  
  
  
  
  beta_prior_mean <- rep(0,ncol(X)) 
  beta_prior_cov <- diag(beta_prior_var_factor^2,ncol(X))
  beta_prior_cov_inv <- solve(beta_prior_cov)
    
  N_g <- colSums(z)
  N_gjk <- apply(Y_indicator, c(2, 3), function(S_jk) t(z) %*% S_jk)
  
  A <- kappa[,1:(G-1)]+omega*C[,1:(G-1)]
  
  return(list(X = X, Y = Y, p = p, M = M, n = n, K = K, beta = beta, mu = mu, C = C, omega = omega, eta = eta, beta_cov = beta_cov, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, beta_mean = beta_mean, w = w, z = z, gamma = gamma, tau = tau, Y_indicator = Y_indicator, N_jk = N_jk, inclusion_sum = inclusion_sum, exclusion_sum = exclusion_sum, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, beta_prior_cov_inv = beta_prior_cov_inv, nu = nu, alpha = alpha, N_g = N_g, N_gjk = N_gjk, X_current = X, logit_probs = logit_probs, log_logit_probs = log_logit_probs, kappa = kappa, A = A))
}


N_updates <- function(z, Y_indicator, N_jk, nu){
  N_g <- colSums(z)
  N_gjk <- apply(Y_indicator, c(2, 3), function(S_jk) t(z) %*% S_jk)
  return(list(N_g = N_g, N_gjk = N_gjk))
}


nu_update <- function(nu, M, K, G, alpha, N_g, N_gjk, inclusion_sum, exclusion_sum, count){
  nu_prop <- nu
  j_prop <- sample(1:M,1)
  if(nu[j_prop] == 1){
    nu_prop[j_prop] <- 0
    log_gamma_N_alpha <- lgamma(N_gjk[,j_prop,] + alpha)
    log_sum1 <- sum(log_gamma_N_alpha)
    diff_log <- sum(lgamma(N_g + K[j_prop]*alpha)) - log_sum1
    log_accept_ratio <- exclusion_sum[j_prop] + diff_log
    #print(paste('exclusion', log_accept_ratio))
  } else {
    nu_prop[j_prop] <- 1
    log_gamma_N_alpha <- lgamma(N_gjk[,j_prop,] + alpha)
    log_sum1 <- sum(log_gamma_N_alpha)
    diff_log <- log_sum1 - sum(lgamma(N_g + K[j_prop]*alpha))
    log_accept_ratio <- inclusion_sum[j_prop] + diff_log
    #print(paste('inclusion', log_accept_ratio))
  }
  accept_ratio <- exp(log_accept_ratio)
  accept_ratio <- as.numeric(accept_ratio)
  if(runif(1) < min(1,accept_ratio)){
    nu <- nu_prop
    #print('accepted')
  }
  #print(accept_ratio)
  return(nu)
}




# THIS z_update FUNCTION NEEDS OPTIMISING


z_update <- function(z, nu, mu, K, alpha, N_gjk, N_g, Y_indicator, n, G,
                     omega, C, p, beta, beta_cov_inv, beta_cov_inv_chol,
                     beta_prior_mean, beta_prior_cov_inv, gamma, X_current) {
  w <- matrix(0, nrow = n, ncol = G)
  which_item_var <- which(nu == 1)
  #K_current <- K[nu == 1]
  K_current <- K[which_item_var]
  for (i in 1:n) {
    curr_z_i <- z[i,]  # Store current assignment

    # Remove contribution of observation i from counts
    N_g_minus_i <- N_g - curr_z_i

    for (g in 1:G) {
      # Calculate new counts if observation i is assigned to group g
      new_z_i <- rep(0, G)
      new_z_i[g] <- 1

      # Update counts
      N_g_temp <- N_g_minus_i + new_z_i

      N_gjk_temp <- N_gjk
      for (j in which_item_var) {
        for (k in 1:dim(Y_indicator)[3]) {
          if (Y_indicator[i,j,k] == 1) {
            for (h in 1:G) {
              # Remove current assignment
              N_gjk_temp[h,j,k] <- N_gjk_temp[h,j,k] - curr_z_i[h]
              # Add new assignment
              N_gjk_temp[h,j,k] <- N_gjk_temp[h,j,k] + (h == g)
            }
          }
        }
      }


      log_gamma_N_gjk_temp <- lgamma(N_gjk_temp[,which_item_var,, drop=FALSE] + alpha)
      term2 <- sum(log_gamma_N_gjk_temp)

      term3 <- sum(sapply(1:G, function(h) {
        sum(lgamma(N_g_temp[h] + K_current * alpha))
      }))



      nonzero_indices <- which(gamma == 1)

      term0 <- mu[i,g]


      w[i,g] <- term0 + term2 - term3
    }


    max_w_i <- max(w[i,])
    exp_w_i <- exp(w[i,] - max_w_i)
    w[i,] <- exp_w_i/sum(exp_w_i)


    z[i,] <- rmultinom(1, 1, w[i,])

    # Update overall counts after this reassignment
    N_g <- N_g_minus_i + z[i,]

    # Update N_gjk
    for (j in 1:dim(Y_indicator)[2]) {
      for (k in 1:dim(Y_indicator)[3]) {
        if (Y_indicator[i,j,k] == 1) {
          N_gjk[,j,k] <- N_gjk[,j,k] - curr_z_i + z[i,]
        }
      }
    }
  }

  # Calculate final kappa
  kappa <- z - 0.5

  return(list(z = z, w = w, kappa = kappa, N_g = N_g, N_gjk = N_gjk))
}






















logit_probs_update <- function(eta){
  logit_probs <- exp(eta)/(1+exp(eta)) 
  return(logit_probs)
}



omega_update <- function(eta,mu,G,n,omega){
  omega <- matrix(rpg(n * (G-1), 1, as.vector(eta[,1:(G-1)])), nrow = n, ncol = (G-1))
  return(omega)
}

mu_update <- function(X,beta){
  return(X%*%beta)
}






C_update <- function(mu, exp_mu, G, C) {
  if (ncol(exp_mu) > 2) {
    epsilon <- 1e-10
    row_sums <- rowSums(exp_mu)
    C <- log(pmax(row_sums - exp_mu, epsilon))
  } else {
    C <- mu[, c(2, 1)]
  }
  return(C)
}








kappa_update <- function(z){
  return(z - 1/2)
}

eta_update <- function(mu,C){
  return(mu - C)
}



A_update <- function(kappa, omega, C){
  return(kappa[,1:(G-1)]+omega*C[,1:(G-1)])
}


beta_update <- function(gamma, X_current, G, p, beta_prior_cov_inv, A, beta_prior_mean, omega){
  beta <- matrix(0, nrow = (p+1), ncol = G)
  nonzero_indices <- which(gamma == 1)
  beta_prior_cov_inv1 <- beta_prior_cov_inv[nonzero_indices,nonzero_indices]
  beta_prior_cov1 <- chol2inv(chol(beta_prior_cov_inv1))
  beta_prior_mean1 <- beta_prior_mean[nonzero_indices]
  n_nonzero <- length(nonzero_indices)
  
  beta_mean <- matrix(0, nrow = (p+1), ncol = (G-1))
  beta_cov_inv <- array(0, dim = c(n_nonzero, n_nonzero, G-1))
  beta_cov_inv_chol <- array(0, dim = c(n_nonzero, n_nonzero, G-1))
  
  for(g in 1:(G-1)){
    beta_cov_inv[,,g] <- crossprod(X_current, omega[,g]*X_current) + beta_prior_cov_inv1
    beta_cov_inv_chol[,,g] <- chol(beta_cov_inv[,,g])
    fs <- forwardsolve(t(beta_cov_inv_chol[,,g]),crossprod(X_current,A[,g]) + beta_prior_cov_inv1%*%beta_prior_mean1)
    beta_mean[nonzero_indices,g] <- backsolve(beta_cov_inv_chol[,,g], fs)
    beta[nonzero_indices,g] <- rMVNormP(1,beta_mean[nonzero_indices,g],beta_cov_inv_chol[,,g])
  }  
  return(list(beta = beta, beta_cov_inv = beta_cov_inv, beta_mean = beta_mean, beta_cov_inv_chol = beta_cov_inv_chol))
}





logsumexp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}


rMVNormP <- function(n, mu, chol_precision) {
  p <- length(mu)
  L <- chol_precision
  Z <- matrix(rnorm(p*n), p, n)
  Y <- backsolve(L, Z, transpose = TRUE)
  X <- sweep(Y, 1, mu, FUN = "+")
  return(X)
}










gamma_update <- function(gamma,p,beta_prior_cov_inv, beta_prior_mean, kappa, G, X_current, omega, tau, X, temp){
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
  return(list(gamma = gamma,X_current = X_current, beta_mean = beta_mean, beta_cov_inv = beta_cov_inv))
}

tau_update <- function(gamma,a,b){
  p <- length(gamma)-1
  sum_gamma <- sum(gamma[2:(p+1)])
  tau <- rbeta(1,sum_gamma + a, p - sum_gamma + b)
  return(tau)
}

BLCR_varsel_burn <- function(burnin, beta_samples, nu_samples, gamma_samples, z_samples, w_samples, N_g_samples, N_gjk_samples, omega_samples){
  n_iter <- dim(beta_samples)[3]
  burned_indices <- (burnin+1):n_iter
  beta_samples_burned <- beta_samples[,,burned_indices]
  nu_samples_burned <- nu_samples[,burned_indices]
  gamma_samples_burned <- gamma_samples[,burned_indices]
  z_samples_burned <- z_samples[,,burned_indices]
  w_samples_burned <- w_samples[,,burned_indices]
  N_g_samples_burned <- N_g_samples[,burned_indices] 
  N_gjk_samples_burned <- N_gjk_samples[,,,burned_indices]
  omega_samples_burned <- omega_samples[,,burned_indices]
  return(list(beta_samples_burned = beta_samples_burned, nu_samples_burned = nu_samples_burned, gamma_samples_burned = gamma_samples_burned, z_samples_burned = z_samples_burned, w_samples_burned = w_samples_burned, N_gjk_samples_burned = N_gjk_samples_burned, N_g_samples_burned = N_g_samples_burned, omega_samples_burned = omega_samples_burned))
}


BLCR_varsel_relabel <- function(perm, beta_samples, z_samples, w_samples, N_gjk_samples, N_g_samples, omega_samples){
  relabelled_beta_samples <- array(0, dim = dim(beta_samples))
  relabelled_z_samples <- array(0, dim = dim(z_samples))
  relabelled_w_samples <- array(0, dim = dim(w_samples))
  relabelled_N_gjk_samples <- array(0, dim = dim(N_gjk_samples))
  relabelled_N_g_samples <- array(0, dim = dim(N_g_samples))
  relabelled_omega_samples <- array(0, dim = dim(omega_samples))
  samples_length <- dim(beta_samples)[3]
  for (sample in 1:samples_length) {
    relabelled_beta_samples[,,sample] <- beta_samples[, perm[sample, ], sample]
    relabelled_z_samples[,,sample] <- z_samples[,perm[sample,], sample]
    relabelled_w_samples[,,sample] <- w_samples[,perm[sample,], sample]
    relabelled_N_gjk_samples[,,,sample] <- N_gjk_samples[perm[sample,],,,sample]
    relabelled_N_g_samples[,sample] <- N_g_samples[perm[sample,],sample]
  }
  #the omega_samples needs to be dealt with for >2 groups, but we'll look at that later. For now I'll just return omega_samples
  relabelled_omega_samples <- omega_samples
  return(list(relabelled_omega_samples = relabelled_omega_samples, relabelled_N_g_samples = relabelled_N_g_samples, relabelled_N_gjk_samples = relabelled_N_gjk_samples, relabelled_w_samples = relabelled_w_samples, relabelled_z_samples = relabelled_z_samples, relabelled_beta_samples = relabelled_beta_samples))
}



BLCR_varsel_thin <- function(thin, beta_samples, nu_samples, gamma_samples, z_samples, w_samples, N_g_samples, N_gjk_samples, omega_samples){
  samples_length <- dim(beta_samples)[3]
  thin_indices <- seq(1, samples_length, by = thin)
  beta_samples_thin <- beta_samples[,,thin_indices]
  nu_samples_thin <- nu_samples[,thin_indices]
  gamma_samples_thin <- gamma_samples[,thin_indices]
  z_samples_thin <- z_samples[,,thin_indices]
  w_samples_thin <- w_samples[,,thin_indices]
  N_g_samples_thin <- N_g_samples[,thin_indices]
  N_gjk_samples_thin <- N_gjk_samples[,,,thin_indices]
  
  #Again, this needs modification for G>2
  omega_samples_thin <- omega_samples[,thin_indices]
  
  return(list(beta_samples_thin = beta_samples_thin, nu_samples_thin = nu_samples_thin, gamma_samples_thin = gamma_samples_thin, z_samples_thin = z_samples_thin, w_samples_thin = w_samples_thin, N_g_samples_thin = N_g_samples_thin, N_gjk_samples_thin = N_gjk_samples_thin, omega_samples_thin = omega_samples_thin))
}


#making a function here that relevels the beta samples AFTER RELABELLING

BLCR_varsel_beta_relevel <- function(beta_samples){
  beta_samples_relevelled <- beta_samples
  G <- dim(beta_samples)[2]
  for (count in 1:(dim(beta_samples)[3])){
    baseline <- beta_samples[,G,count]
    beta_samples_relevelled[,,count] <- beta_samples[,,count] - baseline
  }
  return(beta_samples_relevelled = beta_samples_relevelled)
}









