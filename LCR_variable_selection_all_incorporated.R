library(MCMCpack)
library(BayesLogit)
library(BayesLCA)
library(einsum)

### Functions for initialising variables before starting the sampler

# 1. Function for the sampler without the item selection

initialise_variables_polyagamma_varsel <- function(G,alpha, X,Y, beta_prior_mean = NULL, beta_prior_cov = NULL){
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
  #-------------------------theta initialisation--------------------------------
  #-----------------------------------------------------------------------------
  
  theta <- array(0, dim = c(G, M, max(K)))  #Initialising theta
  for (g in 1:G) { #generating a sample from theta_{gj vectors for each g and j}
    for (j in 1:M) {
      theta[g, j, 1:K[j]] <- rdirichlet(1, alpha)  
    }
  }
  
  
  #-----------------------------------------------------------------------------
  #-------------------------beta initialisation---------------------------------
  #-----------------------------------------------------------------------------
  
  beta <- cbind(matrix(rnorm((G-1)*(p+1)), nrow = (p+1)),0)
  mu <- X%*%beta
  
  
  #-----------------------------------------------------------------------------
  #-------------------------omega, Omega,...initialisation here-----------------
  #-----------------------------------------------------------------------------
  
  C <- matrix(rnorm((n)*(G)), nrow = n)
  omega <- matrix(rpg(n*(G-1)), nrow = n)
  eta <- matrix(0, nrow = n, ncol = G)
  
  beta_cov <- array(0, dim = c(p+1,p+1,G-1))
  beta_cov_inv <- array(0, dim = c(p+1,p+1,G-1))
  beta_cov_inv_chol <- array(0, dim = c(p+1,p+1,G-1))
  beta_mean <- matrix(0, ncol = (G-1), nrow = (p+1) )
  
  
  #-----------------------------------------------------------------------------
  #------------------beta prior parameters here---------------------------------
  #-----------------------------------------------------------------------------
  if (is.null(beta_prior_mean)){
    beta_prior_mean <- rep(0,p+1)
  }
  if (is.null(beta_prior_cov)){
    beta_prior_cov <- diag(10^2, p+1)
  }
  beta_prior_cov_inv <- solve(beta_prior_cov)
  
  #-----------------------------------------------------------------------------
  #-------------------------z initialisation------------------------------------
  #-----------------------------------------------------------------------------
  
  w <- array(0, dim = c(n,G)) 
  z <- array(0, dim = c(n,G)) 
  for(i in 1:n){
    w[i,] <- rep(1/G,G)
    z[i,] <- rmultinom(1,1,w[i,])
  }
  
  #logit probability matrix update
  
  logit_probs <- matrix(0,nrow = n, ncol = G)
  
  
  #-----------------------------------------------------------------------------
  #-------------------------gamma, tau initialisation---------------------------
  #-----------------------------------------------------------------------------
  
  
  gamma <- rep(1,p+1)
  tau <- runif(1)
  kappa <- z-1/2
  A <- kappa[,1:(G-1)]+omega*C[,1:(G-1)]
  
  

  #print(K)
  Y_indicator <- 1*array(outer(Y, 1:max(K), "=="), dim = c(nrow(Y), ncol(Y), length(1:max(K))))
  
  return(list(Y = Y, Y_indicator = Y_indicator, X = X, X_current =X, theta = theta, alpha  = alpha, beta = beta, w = w, z = z, C = C, n = n, M = M,G = G, K = K, p = p, omega = omega,beta_cov = beta_cov, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov ,beta_prior_cov_inv = beta_prior_cov_inv , gamma = gamma, beta_mean = beta_mean, eta = eta, logit_probs = logit_probs, tau = tau, mu = mu, A = A))
  
}





# 2. Function for the sampler with the item selection


initialise_variables_BLCR_collapsed <- function(G, X, Y, beta_prior_cov = NULL, beta_prior_mean = NULL, clust_var_prior = 0.5, alpha  = 1){
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
  
  #need to sort out the alpha (theta hyperparameter) stufff
  sum1_inclusion <- (G-1)*(lgamma(K*alpha) - K*lgamma(alpha))
  log_gamma_N_jk_alpha <- lgamma(N_jk + alpha)
  sum2_inclusion <- lgamma(n + K*alpha) - rowSums(log_gamma_N_jk_alpha)
  sum3_inclusion <- log(clust_var_prior) - log(1-clust_var_prior)
  inclusion_sum <- sum1_inclusion + sum2_inclusion + sum3_inclusion
  exclusion_sum <- -inclusion_sum
  
  
  
  
  
  
  if (is.null(beta_prior_mean)){
    beta_prior_mean <- rep(0,p+1)
  }
  if (is.null(beta_prior_cov)){
    beta_prior_cov <- diag(10^2, p+1)
    beta_prior_cov_inv <- solve(beta_prior_cov)
  } else {
    beta_prior_cov_inv <- solve(beta_prior_cov)
  }
  
  N_g <- colSums(z)
  N_gjk <- apply(Y_indicator, c(2, 3), function(S_jk) t(z) %*% S_jk)
  
  A <- kappa[,1:(G-1)]+omega*C[,1:(G-1)]
  
  return(list(X = X, Y = Y, p = p, M = M, n = n, K = K, beta = beta, mu = mu, C = C, omega = omega, eta = eta, beta_cov = beta_cov, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, beta_mean = beta_mean, w = w, z = z, gamma = gamma, tau = tau, Y_indicator = Y_indicator, N_jk = N_jk, inclusion_sum = inclusion_sum, exclusion_sum = exclusion_sum, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, beta_prior_cov_inv = beta_prior_cov_inv, nu = nu, alpha = alpha, N_g = N_g, N_gjk = N_gjk, X_current = X, logit_probs = logit_probs, log_logit_probs = log_logit_probs, kappa = kappa, A = A))
}


# 3. Function for the regular sampler without variable selection


initialise_variables_polyagamma <- function(G,alpha, X,Y, beta_prior_mean = NULL, beta_prior_cov = NULL){
  if (is.vector(X)){
    ones <- rep(1,length(X))
    X1 <- cbind(ones,X)
  } else{
    ones <- rep(1,nrow(X))
    X1 <- cbind(ones,X)
  }
  p <- ncol(X1) - 1
  M <- ncol(Y)
  n <- nrow(Y)
  K <- apply(Y, 2, function(x) length(unique(x)))
  
  #-----------------------------------------------------------------------------
  #-------------------------theta initialisation--------------------------------
  #-----------------------------------------------------------------------------
  
  theta <- array(0, dim = c(G, M, max(K)))  #Initialising theta
  for (g in 1:G) { #generating a sample from theta_{gj vectors for each g and j}
    for (j in 1:M) {
      theta[g, j, 1:K[j]] <- rdirichlet(1, alpha)  
    }
  }
  
  
  #-----------------------------------------------------------------------------
  #-------------------------beta initialisation---------------------------------
  #-----------------------------------------------------------------------------
  
  beta <- cbind(matrix(rnorm((G-1)*(p+1)), nrow = (p+1)),0)
  logit_probs <- t(apply(matrix(runif(n * G), nrow = n), 1, function(x) x / sum(x)))
  log_logit_probs <- log(logit_probs)
  
  
  #-----------------------------------------------------------------------------
  #-------------------------omega, Omega,...initialisation here------------------------
  #-----------------------------------------------------------------------------
  
  C <- matrix(rnorm((n)*(G)), nrow = n)
  omega <- matrix(rpg(n*(G-1)), nrow = n)
  eta <- matrix(0, nrow = n, ncol = G)
  
  beta_cov <- array(0, dim = c(p+1,p+1,G-1))
  beta_cov_inv <- array(0, dim = c(p+1,p+1,G-1))
  beta_cov_inv_chol <- array(0, dim = c(p+1,p+1,G-1))
  beta_mean <- matrix(0, nrow = (p+1), ncol = (G-1) )
  
  
  #-----------------------------------------------------------------------------
  #------------------beta prior parameters here---------------------------------
  #-----------------------------------------------------------------------------
  if (is.null(beta_prior_mean)){
    beta_prior_mean <- rep(0,p+1)
  }
  if (is.null(beta_prior_cov)){
    beta_prior_cov <- diag(10^2, p+1)
  } 
  beta_prior_cov_inv <- solve(beta_prior_cov)
  
  #-----------------------------------------------------------------------------
  #-------------------------z initialisation------------------------------------
  #-----------------------------------------------------------------------------
  
  w <- array(0, dim = c(n,G)) 
  z <- array(0, dim = c(n,G)) 
  for(i in 1:n){
    w[i,] <- rep(1/G,G)
    z[i,] <- rmultinom(1,1,w[i,])
  }
  
  Y_indicator <- 1*array(outer(Y, 1:max(K), "=="), dim = c(nrow(Y), ncol(Y), length(1:max(K))))
  gamma <- rep(1, p+1)
  
  return(list(Y = Y, X = X1, theta = theta, alpha  = alpha,
              beta = beta, w = w, z = z, C = C, n = n, M = M, G = G,
              K = K, p = p, Y_indicator = Y_indicator, omega = omega,
              beta_cov = beta_cov, beta_cov_inv = beta_cov_inv, 
              beta_cov_inv_chol = beta_cov_inv_chol, 
              beta_prior_mean = beta_prior_mean, 
              beta_prior_cov = beta_prior_cov ,
              beta_prior_cov_inv = beta_prior_cov_inv, 
              beta_mean = beta_mean, eta = eta, logit_probs = logit_probs,
              log_logit_probs = log_logit_probs,
              gamma = gamma))
  
}







##Functions specific to the model with the theta parameters included (uncollapsed) 

S_update <- function(Y_indicator,z,M,G,K){
  S <- array(0, dim = c(G, M, max(K)))
  for (k in 1:max(K)) {
    S[,,k] <- t(z)%*%Y_indicator[,,k]
  }
  return(S)
}

theta_update <- function(S,alpha,K,G,M){
  alpha_S <- array(alpha, dim = dim(S)) + S
  gam_samples <- array(rgamma(prod(dim(S)), shape = as.vector(alpha_S), scale = 1),dim = dim(S))
  gam_sums <- apply(gam_samples, c(1, 2), sum)
  theta <- sweep(gam_samples, MARGIN = c(1, 2), gam_sums, FUN = "/")
  return(theta)
}

# w_update_uncollapsed <- function(log_logit_probs,log_theta,Y,G,n,M,K){
#   log_theta_list <- lapply(1:M, function(j) {
#     log_theta[, j, Y[, j]]
#   })
#   sum_log_theta <- Reduce(`+`, log_theta_list)
#   log_w <- log_logit_probs + t(sum_log_theta)
#   w <- exp(log_w)
#   return(w)
# }
# 
# 
# 
# 
# z_update_uncollapsed <- function(w){
#   z <- t(apply(w, 1, function(row) rmultinom(1, 1, row)))
#   return(z)
# }

z_update_uncollapsed <- function(log_logit_probs,log_theta,Y,G,n,M,K){
  log_theta_list <- lapply(1:M, function(j) {
    log_theta[, j, Y[, j]]
  })
  sum_log_theta <- Reduce(`+`, log_theta_list)
  log_w <- log_logit_probs + t(sum_log_theta)
  w <- exp(log_w)
  z <- t(apply(w, 1, function(row) rmultinom(1, 1, row)))
  return(list(w = w, z = z))
}






##Functions specific to the collapsed model

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
  #print(exclusion_sum)
  return(nu)
}



z_update_collapsed <- function(z, nu, mu, K, alpha, N_gjk, N_g, Y_indicator, n, G,
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
    
    mask_vec <- as.logical(Y_indicator[i, , ] == 1)

    for (g in 1:G) {
      # Calculate new counts if observation i is assigned to group g
      new_z_i <- rep(0, G)
      new_z_i[g] <- 1

      # Update counts
      N_g_temp <- N_g_minus_i + new_z_i

      N_gjk_temp <- N_gjk
      
      
      N_gjk_mat      <- matrix(N_gjk, nrow = G)
      N_gjk_mat_temp <- matrix(N_gjk_temp, nrow = G)
      delta <- -curr_z_i + as.numeric(seq_len(G) == g)
      N_gjk_mat_temp[ , mask_vec] <- N_gjk_mat[ , mask_vec] + delta
      dim(N_gjk_mat_temp) <- c(G,M,max(K))
      N_gjk_temp <- N_gjk_mat_temp


      log_gamma_N_gjk_temp <- lgamma(N_gjk_temp[,which_item_var,, drop=FALSE] + alpha)
      term2 <- sum(log_gamma_N_gjk_temp)

      temp_mat <- outer(N_g_temp, K_current * alpha, "+")
      term3 <- sum(lgamma(temp_mat))
      
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
    delta <- z[i, ] - curr_z_i
    N_gjk_mat[ , mask_vec] <- N_gjk_mat[ , mask_vec] + delta
    dim(N_gjk_mat) <- c(G,M,max(K))
    N_gjk <- N_gjk_mat
  }

  # Calculate final kappa
  kappa <- z - 0.5

  return(list(z = z, w = w, kappa = kappa, N_g = N_g, N_gjk = N_gjk))
}
















## Functions common to both collapsed and uncollapsed

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



A_update <- function(kappa, omega, C, G){
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


gamma_update <- function(gamma,p,beta_prior_cov_inv, beta_prior_mean, kappa, G, X_current, omega, tau, X, A, beta_mean, beta_cov_inv ){
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







## General functions


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

log_softmax <- function(mu) {
  max_per_row <- apply(mu, 1, max)
  mu_centered <- mu - max_per_row
  logsumexp <- log(rowSums(exp(mu_centered)))
  result <- mu_centered - logsumexp
  return(result)
}









## Post processing functions - need to make these for the uncollapsed ones as well

BLCR_varsel_burn_collapsed <- function(burnin, beta_samples, nu_samples, gamma_samples, z_samples, w_samples, N_g_samples, N_gjk_samples, omega_samples){
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




BLCR_varsel_relabel_collapsed <- function(perm, beta_samples, z_samples, w_samples, N_gjk_samples, N_g_samples, omega_samples){
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



BLCR_varsel_thin_collapsed <- function(thin, beta_samples, nu_samples, gamma_samples, z_samples, w_samples, N_g_samples, N_gjk_samples, omega_samples){
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


BLCR_varsel_beta_relevel_collapsed <- function(beta_samples){
  beta_samples_relevelled <- beta_samples
  G <- dim(beta_samples)[2]
  for (count in 1:(dim(beta_samples)[3])){
    baseline <- beta_samples[,G,count]
    beta_samples_relevelled[,,count] <- beta_samples[,,count] - baseline
  }
  return(beta_samples_relevelled = beta_samples_relevelled)
}




LCR_log_post_compute <- function(mu, Y_indicator, log_theta, z, theta_prior_param, beta_prior_mean, beta_prior_cov_inv, beta){
  term1_mat <- log_softmax(mu)
  term2_array <- einsum("ijk,gjk -> igj", Y_indicator, log_theta)
  term2_mat <- apply(term2_array, c(1, 2), sum)
  combined_terms <- z*(term1_mat + term2_mat)
  sum_combined <- sum(combined_terms)
  theta_prior_term1 <- theta_prior_param - 1 
  #log_theta_scaled <- log_theta * array(theta_prior_term1, dim = c(1, 1, length(theta_prior_term1)))
  log_theta_scaled <- log_theta * array(theta_prior_term1, dim = dim(log_theta))
  prior_term1 <- sum(log_theta_scaled)
  beta_prior_term1 <- (-1/2)*sum(colSums(beta * (beta_prior_cov_inv %*% beta)))
  beta_prior_term2 <- sum(crossprod(beta_prior_mean, beta_prior_cov_inv %*% beta))
  log_like <- sum_combined
  log_post <- sum_combined + prior_term1 + beta_prior_term1 + beta_prior_term2
  return(list(log_like = log_like, log_post = log_post))
}







LCR_collapsed_log_post_compute <- function(beta, beta_prior_cov_inv, beta_prior_mean, nu, mu, M, clust_var_prior, z, K, N_jk, N_gjk, N_g, n, theta_hyperparam){
  beta_prior_term1 <- (-1/2)*sum(colSums(beta * (beta_prior_cov_inv %*% beta)))
  beta_prior_term2 <- sum(crossprod(beta_prior_mean, beta_prior_cov_inv %*% beta))
  beta_prior_term <- beta_prior_term1 + beta_prior_term2
  current_indices <- which(nu == 1)
  excl_indices <- which(nu == 0)
  sum_current_indices <- sum(nu)
  nu_prior_term <- sum_current_indices*log(clust_var_prior) + (M-sum_current_indices)*log(1-clust_var_prior)
  logit_term <- sum(z*log_softmax(mu))
  K_current <- K[current_indices]
  K_excl <- K[excl_indices]
  N_jk_excl <- N_jk[excl_indices,]
  N_jk_current <- N_jk[current_indices,]
  N_gjk_current <- N_gjk[,current_indices,]
  term1 <- sum(lgamma(K_excl*theta_hyperparam))
  term2 <- sum(K_excl*lgamma(theta_hyperparam))
  term3 <- sum(lgamma(N_jk_excl + theta_hyperparam))
  term4 <- sum(lgamma(n + K_excl*theta_hyperparam))
  log_like_term1 <- term1 - term2 + term3 - term4
  term5 <- G*sum(lgamma(K_current*theta_hyperparam))
  term6 <- G*sum(K_current*lgamma(theta_hyperparam))
  term7 <- sum(lgamma(N_gjk_current + theta_hyperparam))
  term8 <- sum_current_indices*sum(lgamma(N_g+theta_hyperparam))
  log_like_term2 <- term5 - term6 + term7 - term8
  log_like <- log_like_term1 + log_like_term2
  log_post <- log_like + beta_prior_term + nu_prior_term + logit_term
  return(list(log_like = log_like, log_post = log_post))
}











## Creating 4 different functions, with the aim of combining them later into a single function
## 1. LCR sampler using Polya-Gamma augmentation.
## 2. LCR sampler with item selection.
## 3. LCR sampler with covariate selection.
## 4. LCR sampler with item and covariate selection.



## 1. LCR sampler using Polya-Gamma augmentation. Few bits to be added in this
# needs to return stuff, and relabelling needs to be included somehow. 

LCR_Gibbs_no_sel <- function(X, Y, G = 2, theta_prior_param, beta_prior_mean, beta_prior_cov, burnin = 500, n_samples = 1000, thinby = 1, verbose = FALSE, relabel = TRUE){
  init <- initialise_variables_polyagamma(G = G, a = theta_prior_param, X = X, Y = Y, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov)
  list2env(init, envir = environment())
  n_iter <- thinby*n_samples + burnin
  beta_samples <- array(0, dim = c(dim(beta), n_samples))
  omega_samples <- array(0, dim = c(dim(omega), n_samples))
  w_samples <- array(0, dim = c(dim(w), n_samples))
  z_samples <- array(0, dim = c(dim(z), n_samples))
  theta_samples <- array(0, dim = c(dim(theta), n_samples))
  log_post_samples <- numeric(n_samples)
  log_like_samples <- numeric(n_samples)
  sample_count <- 0
  start_time <- Sys.time()
  progress_interval <- 500
  for (count in 1:n_iter){
    S <- S_update(Y_indicator = Y_indicator, z = z, M = M, G = G, K = K)
    theta <- theta_update(S = S, alpha = theta_prior_param, K = K, G = G, M = M)
    mu <- mu_update(X = X, beta = beta)
    exp_mu <- exp(mu)
    C <- C_update(mu = mu, exp_mu = exp_mu, G = G, C = C)
    eta <- eta_update(mu = mu, C = C)
    logit_probs <- logit_probs_update(eta = eta)
    log_logit_probs <- log(logit_probs)
    log_theta <- log(theta)
    #w <- w_update_uncollapsed(log_logit_probs = log_logit_probs,log_theta = log_theta, Y = Y, G = G, n = n, M = M, K = K)
    #z <- z_update_uncollapsed(w = w)
    z_up <- z_update_uncollapsed(log_logit_probs = log_logit_probs, log_theta = log_theta, Y = Y, G = G, n = n, M = M, K = K)
    list2env(z_up, envir = environment())
    kappa <- kappa_update(z = z)
    omega <- omega_update(eta = eta, mu = mu, G = G, n = n, omega = omega)
    A <- A_update(kappa = kappa, omega = omega, C = C, G = G)
    beta_up <- beta_update(gamma = gamma , X_current = X, G = G, p = p, beta_prior_cov_inv = beta_prior_cov_inv , A = A, beta_prior_mean = beta_prior_mean , omega = omega)
    list2env(beta_up, envir = environment())
    log_post_compute <- LCR_log_post_compute(mu = mu, Y_indicator = Y_indicator, log_theta = log_theta, z = z, theta_prior_param = theta_prior_param, beta_prior_mean = beta_prior_mean, beta_prior_cov_inv = beta_prior_cov_inv, beta = beta)
    list2env(log_post_compute, envir = environment())
    if (count > burnin && ((count - burnin) %% thinby == 0)) {
      sample_count <- sample_count + 1
      beta_samples[,,sample_count] <- beta
      omega_samples[,,sample_count] <- omega
      z_samples[,,sample_count] <- z
      w_samples[,,sample_count] <- w
      theta_samples[,,,sample_count] <- theta
      log_post_samples[sample_count] <- log_post
      log_like_samples[sample_count] <- log_like
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
      cat(sprintf("\r%s %3d%% | %d of %d samples | ~%.1f sec remaining", 
                  bar, percent, sample_count, n_samples, est_remaining))
      flush.console()
    }
  }
  if (verbose) {
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("\nSampling complete! Total time: %.1f seconds (%.1f minutes)\n", 
                total_time, total_time / 60))
  }
  if (relabel){
    print("relabelling")
    z_mat <- matrix(0, nrow = n_samples, ncol = n)
    for (iter in 1:n_samples){
      temp <- z_samples[,,iter]
      temp_row <- apply(temp,1,which.max)
      z_mat[iter,] <- temp_row 
    }
    w_samples_reshape <- aperm(w_samples,c(3,1,2))
    ls_perm <- ls$permutations$STEPHENS
    reordered_beta_samples <- array(0, dim = dim(beta_samples))
    reordered_z_samples <- array(0, dim = dim(z_samples))
    reordered_w_samples <- array(0, dim = dim(w_samples))
    reordered_theta_samples <- array(0, dim = dim(theta_samples))
    for (i in 1:n_samples) {
      reordered_beta_samples[,,i] <- beta_samples[,ls_perm[i,],i]
      reordered_z_samples[,,i] <- z_samples[,ls_perm[i,],i]
      reordered_w_samples[,,i] <- w_samples[,ls_perm[i,],i]
      reordered_theta_samples[,,,i] <- theta_samples[ls_perm[i,], , ,i]
    }
    beta_samples <- reordered_beta_samples
    z_samples <- reordered_z_samples
    w_samples <- reordered_w_samples
    theta_samples <- reordered_theta_samples
  }
  N_g_samples <- apply(z_samples, 3, colSums)
  pi_samples <- N_g_samples/colSums(N_g_samples)
  #Estimates for each parameter
  beta_estimate <- apply(beta_samples, c(1,2), mean)
  beta_sd <- apply(beta_samples, c(1,2), sd)
  mu_estimate <-  mu_update(X = X, beta = beta_estimate)
  #need to reformat the theta matrix into the list as before
  theta_estimate <- apply(theta_samples, c(1,2,3), mean)
  log_theta_estimate <- log(theta_estimate)
  theta_sd <- apply(theta_samples, c(1,2,3), sd)
  
  Z <- apply(z_samples, c(1,2), mean)
  z_estimate <- apply(Z,1,which.max)
  
  log_post_estimate_compute <- LCR_log_post_compute(mu = mu_estimate, Y_indicator = Y_indicator, log_theta = log_theta_estimate, z = z_estimate, theta_prior_param = theta_prior_param, beta_prior_mean = beta_prior_mean, beta_prior_cov_inv = beta_prior_cov_inv, beta = beta_estimate)
  log_post_estimate <- log_post_estimate_compute$log_post
  log_like_estimate <- log_post_estimate_compute$log_like
  deviance_estimated_params <- -2*log_like_estimate
  deviance_samples <- -2*log_post_samples
  mean_deviance <- mean(deviance_samples)
  DIC <- 2*mean_deviance - deviance_estimated_params
  samples <- list(logpost = log_post_samples, loglik = log_like_samples) #FINISH THIS OFF - need to reformat the item probabilities
  return(list(beta_samples = beta_samples, theta_samples = theta_samples, z_samples = z_samples, w_samples = w_samples, pi_samples = pi_samples,log_post_samples = log_post_samples, log_like_samples = log_like_samples))
}



## 2. LCR sampler using Polya-Gamma augmentation, with covariate selection. Few bits to be added in this
# needs to return stuff, and relabelling needs to be included somehow. 


LCR_Gibbs_cov_sel <- function(X, Y, G = 2, theta_prior_param = NULL, beta_prior_mean, beta_prior_cov, burnin = 500, n_samples = 1000, thinby = 1, verbose = FALSE, tau_prior_a = 1, tau_prior_b = 1, relabel = TRUE){
  init <- initialise_variables_polyagamma_varsel(G = G, a = theta_prior_param, X = X, Y = Y, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov)
  list2env(init, envir = environment())
  n_iter <- thinby*n_samples + burnin
  beta_samples <- array(0, dim = c(dim(beta), n_samples))
  omega_samples <- array(0, dim = c(dim(omega), n_samples))
  w_samples <- array(0, dim = c(dim(w), n_samples))
  z_samples <- array(0, dim = c(dim(z), n_samples))
  theta_samples <- array(0, dim = c(dim(theta), n_samples))
  gamma_samples <- matrix(0, nrow = length(gamma), ncol = n_samples)
  log_post_samples <- numeric(n_samples)
  log_like_samples <- numeric(n_samples)
  tau_samples <- numeric(n_samples)
  sample_count <- 0
  start_time <- Sys.time()
  progress_interval <- 500
  for (count in 1:n_iter){
    S <- S_update(Y_indicator = Y_indicator, z = z, M = M, G = G, K = K)
    theta <- theta_update(S = S, alpha = theta_prior_param, K = K, G = G, M = M)
    mu <- mu_update(X = X, beta = beta)
    exp_mu <- exp(mu)
    C <- C_update(mu = mu, exp_mu = exp_mu, G = G, C = C)
    eta <- eta_update(mu = mu, C = C)
    logit_probs <- logit_probs_update(eta = eta)
    log_logit_probs <- log(logit_probs)
    log_theta <- log(theta)
    #w <- w_update_uncollapsed(log_logit_probs = log_logit_probs,log_theta = log_theta, Y = Y, G = G, n = n, M = M, K = K)
    #z <- z_update_uncollapsed(w = w)
    z_up <- z_update_uncollapsed(log_logit_probs = log_logit_probs, log_theta = log_theta, Y = Y, G = G, n = n, M = M, K = K)
    list2env(z_up, envir = environment())
    kappa <- kappa_update(z = z)
    omega <- omega_update(eta = eta, mu = mu, G = G, n = n, omega = omega)
    A <- A_update(kappa = kappa, omega = omega, C = C, G = G)
    beta_up <- beta_update(gamma = gamma , X_current = X_current, G = G, p = p, beta_prior_cov_inv = beta_prior_cov_inv , A = A, beta_prior_mean = beta_prior_mean , omega = omega)
    list2env(beta_up, envir = environment())
    gamma_up <- gamma_update(gamma = gamma, p = p, beta_prior_cov_inv = beta_prior_cov_inv, beta_prior_mean = beta_prior_mean, kappa = kappa, G = G, omega = omega, tau = tau, X = X, A = A, beta_mean = beta_mean, beta_cov_inv = beta_cov_inv)
    list2env(gamma_up, envir = environment())
    log_post_compute <- LCR_log_post_compute(mu = mu, Y_indicator = Y_indicator, log_theta = log_theta, z = z, theta_prior_param = theta_prior_param, beta_prior_mean = beta_prior_mean, beta_prior_cov_inv = beta_prior_cov_inv, beta = beta)
    list2env(log_post_compute, envir = environment())
    if (count > burnin && ((count - burnin) %% thinby == 0)) {
      sample_count <- sample_count + 1
      beta_samples[,,sample_count] <- beta
      omega_samples[,,sample_count] <- omega
      z_samples[,,sample_count] <- z
      w_samples[,,sample_count] <- w
      theta_samples[,,,sample_count] <- theta
      gamma_samples[,sample_count] <- gamma
      tau_samples[sample_count] <- tau
      log_post_samples[sample_count] <- log_post
      log_like_samples[sample_count] <- log_like
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
      cat(sprintf("\r%s %3d%% | %d of %d samples | ~%.1f sec remaining", 
                  bar, percent, sample_count, n_samples, est_remaining))
      flush.console()
    }
  }
  if (verbose) {
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("\nSampling complete! Total time: %.1f seconds (%.1f minutes)\n", 
                total_time, total_time / 60))
  }
  if (relabel){
    print("relabelling")
    z_mat <- matrix(0, nrow = n_samples, ncol = n)
    for (iter in 1:n_samples){
      temp <- z_samples[,,iter]
      temp_row <- apply(temp,1,which.max)
      z_mat[iter,] <- temp_row 
    }
    w_samples_reshape <- aperm(w_samples,c(3,1,2))
    ls <- label.switching(method = "STEPHENS", z = z_mat, K = G, p = w_samples_reshape)
    ls_perm <- ls$permutations$STEPHENS
    reordered_beta_samples <- array(0, dim = dim(beta_samples))
    reordered_z_samples <- array(0, dim = dim(z_samples))
    reordered_w_samples <- array(0, dim = dim(w_samples))
    reordered_theta_samples <- array(0, dim = dim(theta_samples))
    for (i in 1:n_samples) {
      reordered_beta_samples[,,i] <- beta_samples[,ls_perm[i,],i]
      reordered_z_samples[,,i] <- z_samples[,ls_perm[i,],i]
      reordered_w_samples[,,i] <- w_samples[,ls_perm[i,],i]
      reordered_theta_samples[,,,i] <- theta_samples[ls_perm[i,], , ,i]
    }
    beta_samples <- reordered_beta_samples
    z_samples <- reordered_z_samples
    w_samples <- reordered_w_samples
    theta_samples <- reordered_theta_samples
  }
  N_g_samples <- apply(z_samples, 3, colSums)
  pi_samples <- N_g_samples/colSums(N_g_samples)
  deviance_samples <- -2*log_post_samples
  mean_deviance <- mean(deviance_samples)
  #need to calculate the mean of parameter values to use for the deviance calculation
  
  #DIC <- 
  return(list(beta_samples = beta_samples, gamma_samples = gamma_samples, theta_samples = theta_samples, z_samples = z_samples, w_samples = w_samples, log_post_samples = log_post_samples, log_like_samples = log_like_samples, pi_samples = pi_samples))
}


## 3. LCR sampler using Polya-Gamma augmentation, with item selection. Few bits to be added in this
# needs to return stuff, and relabelling needs to be included somehow.


LCR_Gibbs_item_sel <- function(X, Y, G = 2, theta_prior_param, beta_prior_mean, beta_prior_cov, burnin = 500, n_samples = 1000, thinby = 1, verbose = FALSE, tau_prior_a = 1, tau_prior_b = 1, clust_var_prior = 0.5, relabel = TRUE){
  init <- initialise_variables_BLCR_collapsed(G = G, X = X, Y = Y, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, clust_var_prior = clust_var_prior, alpha = theta_prior_param)
  list2env(init, envir = environment())
  n_iter <- thinby*n_samples + burnin
  beta_samples <- array(0, dim = c(dim(beta), n_samples))
  omega_samples <- array(0, dim = c(dim(omega), n_samples))
  nu_samples <- array(0, dim = c(length(nu), n_samples))
  z_samples <- array(0, dim = c(dim(z), n_samples))
  w_samples <- array(0, dim = c(dim(w), n_samples))
  N_g_samples <- array(0, dim = c(length(N_g), n_samples))
  N_gjk_samples <- array(0, dim = c(dim(N_gjk), n_samples))
  log_post_samples <- numeric(n_samples)
  log_like_samples <- numeric(n_samples)
  sample_count <- 0
  start_time <- Sys.time()
  progress_interval <- 500
  for (count in 1:n_iter){
    N_up <- N_updates(z = z, Y_indicator = Y_indicator, N_jk = N_jk, nu = nu)
    N_g <- N_up$N_g
    N_gjk <- N_up$N_gjk
    mu <- mu_update(X = X, beta = beta)
    exp_mu <- exp(mu)
    nu <- nu_update(nu = nu, M = M, K = K, G = G, alpha = theta_prior_param, N_g, N_gjk = N_gjk, inclusion_sum = inclusion_sum, exclusion_sum = exclusion_sum, count = count)
    C <- C_update(mu = mu, exp_mu = exp_mu, G = G, C = C)
    eta <- eta_update(mu = mu, C = C)
    logit_probs <- logit_probs_update(eta = eta)
    log_logit_probs <- log(logit_probs)
    z_up <- z_update_collapsed(z = z, nu = nu, mu = mu, K = K, alpha = theta_prior_param, N_gjk = N_gjk, N_g = N_g, n = n, G = G, Y_indicator = Y_indicator, gamma = gamma, X_current = X_current,  beta = beta, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, omega = omega, C = C)
    list2env(z_up, envir = environment())
    kappa <- kappa_update(z = z)
    omega <- omega_update(eta = eta, mu = mu, G = G, n = n, omega = omega)
    A <- A_update(kappa = kappa, omega = omega, C = C, G = G)
    beta_up <- beta_update(gamma = gamma , X_current = X, G = G, p = p, beta_prior_cov_inv = beta_prior_cov_inv , A = A, beta_prior_mean = beta_prior_mean , omega = omega)
    list2env(beta_up, envir = environment())
    log_post_compute <- LCR_collapsed_log_post_compute(beta = beta, beta_prior_cov_inv = beta_prior_cov_inv, beta_prior_mean = beta_prior_mean, nu = nu, mu = mu, M = M, clust_var_prior = clust_var_prior, z = z, K = K, N_jk = N_jk, N_gjk = N_gjk, N_g = N_g, n = n, theta_hyperparam = theta_prior_param)
    list2env(log_post_compute, envir = environment())
    if (count > burnin && ((count - burnin) %% thinby == 0)) {
      sample_count <- sample_count + 1
      beta_samples[,,sample_count] <- beta
      omega_samples[,,sample_count] <- omega
      z_samples[,,sample_count] <- z
      w_samples[,,sample_count] <- w
      nu_samples[,sample_count] <- nu
      N_gjk_samples[,,,sample_count] <- N_gjk
      N_g_samples[,sample_count] <- N_g
      log_post_samples[sample_count] <- log_post
      log_like_samples[sample_count] <- log_like
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
      cat(sprintf("\r%s %3d%% | %d of %d samples | ~%.1f sec remaining", 
                  bar, percent, sample_count, n_samples, est_remaining))
      flush.console()
    }
  }
  if (verbose) {
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("\nSampling complete! Total time: %.1f seconds (%.1f minutes)\n", 
                total_time, total_time / 60))
  }
  if (relabel){
    print("relabelling")
    z_mat <- matrix(0, nrow = n_samples, ncol = n)
    for (iter in 1:n_samples){
      temp <- z_samples[,,iter]
      temp_row <- apply(temp,1,which.max)
      z_mat[iter,] <- temp_row 
    }
    w_samples_reshape <- aperm(w_samples,c(3,1,2))
    ls <- label.switching(method = "STEPHENS", z = z_mat, K = G, p = w_samples_reshape)
    ls_perm <- ls$permutations$STEPHENS
    reordered_beta_samples <- array(0, dim = dim(beta_samples))
    reordered_z_samples <- array(0, dim = dim(z_samples))
    reordered_w_samples <- array(0, dim = dim(w_samples))
    reordered_N_gjk_samples <- array(0, dim = dim(N_gjk_samples))
    reordered_N_g_samples <- array(0, dim = dim(N_g_samples))
    for (i in 1:n_samples) {
      reordered_beta_samples[,,i] <- beta_samples[,ls_perm[i,],i]
      reordered_z_samples[,,i] <- z_samples[,ls_perm[i,],i]
      reordered_w_samples[,,i] <- w_samples[,ls_perm[i,],i]
      reordered_N_gjk_samples[,,,i] <- N_gjk_samples[ls_perm[i,],,,i]
      reordered_N_g_samples[,i] <- N_g_samples[ls_perm[i,],i]
    }
    beta_samples <- reordered_beta_samples
    z_samples <- reordered_z_samples
    w_samples <- reordered_w_samples
    N_gjk_samples <- reordered_N_gjk_samples
    N_g_samples <- reordered_N_g_samples
  }
  pi_samples <- N_g_samples/colSums(N_g_samples)
  deviance_samples <- -2*log_post_samples
  mean_deviance <- mean(deviance_samples)
  #need to calculate the mean of parameter values to use for the deviance calculation
  
  DIC <- 
}


## 4. LCR sampler using Polya-Gamma augmentation, with covariate and item selection. Few bits to be added in this
# needs to return stuff, and relabelling needs to be included somehow. 


#need to test the performance of the variable selection method as well.
LCR_Gibbs_both_sel <- function(X, Y, G = 2, theta_prior_param, beta_prior_mean, beta_prior_cov, burnin = 500, n_samples = 1000, thinby = 1, verbose = FALSE, tau_prior_a = 1, tau_prior_b = 1, clust_var_prior = 0.5, relabel = TRUE){
  init <- initialise_variables_BLCR_collapsed(G = G, X = X, Y = Y, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, clust_var_prior = clust_var_prior, alpha = theta_prior_param)
  list2env(init, envir = environment())
  n_iter <- thinby*n_samples + burnin
  beta_samples <- array(0, dim = c(dim(beta), n_samples))
  omega_samples <- array(0, dim = c(dim(omega), n_samples))
  nu_samples <- array(0, dim = c(length(nu), n_samples))
  z_samples <- array(0, dim = c(dim(z), n_samples))
  w_samples <- array(0, dim = c(dim(w), n_samples))
  N_g_samples <- array(0, dim = c(length(N_g), n_samples))
  N_gjk_samples <- array(0, dim = c(dim(N_gjk), n_samples))
  gamma_samples <- matrix(0, nrow = length(gamma), ncol = n_samples)
  log_post_samples <- numeric(n_samples)
  log_like_samples <- numeric(n_samples)
  sample_count <- 0
  start_time <- Sys.time()
  progress_interval <- 500
  for (count in 1:n_iter){
    N_up <- N_updates(z = z, Y_indicator = Y_indicator, N_jk = N_jk, nu = nu)
    N_g <- N_up$N_g
    N_gjk <- N_up$N_gjk
    mu <- mu_update(X = X, beta = beta)
    exp_mu <- exp(mu)
    nu <- nu_update(nu = nu, M = M, K = K, G = G, alpha = theta_prior_param, N_g, N_gjk = N_gjk, inclusion_sum = inclusion_sum, exclusion_sum = exclusion_sum, count = count)
    C <- C_update(mu = mu, exp_mu = exp_mu, G = G, C = C)
    eta <- eta_update(mu = mu, C = C)
    logit_probs <- logit_probs_update(eta = eta)
    log_logit_probs <- log(logit_probs)
    z_up <- z_update_collapsed(z = z, nu = nu, mu = mu, K = K, alpha = theta_prior_param, N_gjk = N_gjk, N_g = N_g, n = n, G = G, Y_indicator = Y_indicator, gamma = gamma, X_current = X_current,  beta = beta, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, omega = omega, C = C)
    list2env(z_up, envir = environment())
    kappa <- kappa_update(z = z)
    omega <- omega_update(eta = eta, mu = mu, G = G, n = n, omega = omega)
    A <- A_update(kappa = kappa, omega = omega, C = C, G = G)
    beta_up <- beta_update(gamma = gamma , X_current = X_current, G = G, p = p, beta_prior_cov_inv = beta_prior_cov_inv , A = A, beta_prior_mean = beta_prior_mean , omega = omega)
    list2env(beta_up, envir = environment())
    gamma_up <- gamma_update(gamma = gamma, p = p, beta_prior_cov_inv = beta_prior_cov_inv, beta_prior_mean = beta_prior_mean, kappa = kappa, G = G, omega = omega, tau = tau, X = X, A = A, beta_mean = beta_mean, beta_cov_inv = beta_cov_inv)
    list2env(gamma_up, envir = environment())
    log_post_compute <- LCR_collapsed_log_post_compute(beta = beta, beta_prior_cov_inv = beta_prior_cov_inv, beta_prior_mean = beta_prior_mean, nu = nu, mu = mu, M = M, clust_var_prior = clust_var_prior, z = z, K = K, N_jk = N_jk, N_gjk = N_gjk, N_g = N_g, n = n, theta_hyperparam = theta_prior_param)
    list2env(log_post_compute, envir = environment())
    if (count > burnin && ((count - burnin) %% thinby == 0)) {
      sample_count <- sample_count + 1
      beta_samples[,,sample_count] <- beta
      omega_samples[,,sample_count] <- omega
      z_samples[,,sample_count] <- z
      w_samples[,,sample_count] <- w
      nu_samples[,sample_count] <- nu
      N_gjk_samples[,,,sample_count] <- N_gjk
      N_g_samples[,sample_count] <- N_g
      gamma_samples[,sample_count] <- gamma
      log_post_samples[sample_count] <- log_post
      log_like_samples[sample_count] <- log_like
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
      cat(sprintf("\r%s %3d%% | %d of %d samples | ~%.1f sec remaining", 
                  bar, percent, sample_count, n_samples, est_remaining))
      flush.console()
    }
  }
  if (verbose) {
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("\nSampling complete! Total time: %.1f seconds (%.1f minutes)\n", 
                total_time, total_time / 60))
  }
  if (relabel){
    print("relabelling")
    z_mat <- matrix(0, nrow = n_samples, ncol = n)
    for (iter in 1:n_samples){
      temp <- z_samples[,,iter]
      temp_row <- apply(temp,1,which.max)
      z_mat[iter,] <- temp_row 
    }
    w_samples_reshape <- aperm(w_samples,c(3,1,2))
    ls <- label.switching(method = "STEPHENS", z = z_mat, K = G, p = w_samples_reshape)
    ls_perm <- ls$permutations$STEPHENS
    reordered_beta_samples <- array(0, dim = dim(beta_samples))
    reordered_z_samples <- array(0, dim = dim(z_samples))
    reordered_w_samples <- array(0, dim = dim(w_samples))
    reordered_N_gjk_samples <- array(0, dim = dim(N_gjk_samples))
    reordered_N_g_samples <- array(0, dim = dim(N_g_samples))
    for (i in 1:n_samples) {
      reordered_beta_samples[,,i] <- beta_samples[,ls_perm[i,],i]
      reordered_z_samples[,,i] <- z_samples[,ls_perm[i,],i]
      reordered_w_samples[,,i] <- w_samples[,ls_perm[i,],i]
      reordered_N_gjk_samples[,,,i] <- N_gjk_samples[ls_perm[i,],,,i]
      reordered_N_g_samples[,i] <- N_g_samples[ls_perm[i,],i]
    }
    beta_samples <- reordered_beta_samples
    z_samples <- reordered_z_samples
    w_samples <- reordered_w_samples
    N_gjk_samples <- reordered_N_gjk_samples
    N_g_samples <- reordered_N_g_samples
  }
  pi_samples <- N_g_samples/colSums(N_g_samples)
  deviance_samples <- -2*log_post_samples
  mean_deviance <- mean(deviance_samples)
  #need to calculate the mean of parameter values to use for the deviance calculation
  
  DIC <- 
  return(list(beta_samples = beta_samples, gamma_samples = gamma_samples, nu_samples = nu_samples, z_samples = z_samples, w_samples = w_samples, N_gjk_samples = N_gjk_samples, N_g_samples = N_g_samples, pi_samples = pi_samples))
}










#Simulating some data - testing the covariate selection



p <- 4
G <- 2
n <- 500

#item probabilities 
theta1 <- matrix(c(0.6,0.4,0.2,0.8), nrow = G, byrow = TRUE)
theta2 <- matrix(c(0.8,0.2,0.5,0.5), nrow = G, byrow = TRUE)
theta3 <- matrix(c(0.7,0.3,0.4,0.6), nrow = G, byrow = TRUE)
theta4 <- matrix(c(0.6,0.4,0.9,0.1), nrow = G, byrow = TRUE)
theta5 <- matrix(c(0.5,0.5,0.5,0.5), nrow = G, byrow = TRUE)
theta6 <- matrix(c(0.4,0.6,0.4,0.6), nrow = G, byrow = TRUE)
theta7 <- matrix(c(0.3,0.7,0.3,0.7), nrow = G, byrow = TRUE)
theta8 <- matrix(c(0.2,0.8,0.2,0.8), nrow = G, byrow = TRUE)
theta9 <- matrix(c(0.9,0.1,0.9,0.1), nrow = G, byrow = TRUE)
theta10 <- matrix(c(0.6,0.4,0.6,0.4), nrow = G, byrow = TRUE)
theta11 <- matrix(c(0.7,0.3,0.7,0.3), nrow = G, byrow = TRUE)
theta12 <- matrix(c(0.8,0.2,0.8,0.2), nrow = G, byrow = TRUE)
theta13 <- matrix(c(0.1,0.9,0.1,0.9), nrow = G, byrow = TRUE)

realtheta <- list(theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10, theta11, theta12, theta13)
M <- length(realtheta)
K <- sapply(realtheta, ncol)

X_raw <- matrix(rnorm(n * (p)), nrow = n)

X <- scale(X_raw)

#X[,2] <- 0.7*X[,1] + rnorm(n,0,0.1)

#realbeta <- c(0,0.7,1,-0.8,0.5,0)
#realbeta <- c(0,0.5,1,-0.5,0)
#realbeta <- c(0, 0.3, 0.5, -0.4, 0.2)
#realbeta <- c(0,0.7,1,-0.8,0.5,0)
realbeta <- c(0,0.7,1,-0.8,0.5) 
#realbeta <- c(0,0.3,0.7,-0.4,0.2)
#realbeta <- c(0,1.5,2,-3,-1)
#realbeta <- c(0,-0.8,1,0.7,0.5)

mu <- cbind(1,X)%*%realbeta

logit_prob <- exp(mu)/(1+exp(mu))

logit_prob <- cbind(logit_prob, 1-logit_prob)

class <- apply(logit_prob, 1, function(x) sample(1:G, 1, prob = x))

Y <- matrix(0,nrow = n,ncol = M)

for(i in 1:n){
  for (j in 1:M){
    theta_mat_row <- realtheta[[j]][class[i],]
    Y[i,j] <- sample(1:K[j], 1, prob = theta_mat_row)
  }
}


#X <- X[,2:6]


#testing the effect that the last covariate has on the logit probabilities
mu_excl <- cbind(1,X)[,2:4]%*%realbeta[2:4]
logit_prob_excl <- exp(mu_excl)/(1+exp(mu_excl))
logit_prob_excl <- cbind(logit_prob_excl, 1-logit_prob_excl)











beta_prior_mean <- rep(0, p+1)
beta_prior_cov <- diag(10^2, p+1)


#Y <- Y[,1:4]

#an interesting thing is happening when I include a noise variable along with the other ones (in particular beta = (0, 0.7, 1, -0.8, 0.5) vs beta = (0, 0.7, 1, -0.8, 0.5, 0))
#We get that the 2nd variable (corresponding to beta = 0.7) is excluded almost all of the time in the model where we have the extra noise covariate. Needs to be looked into further.


#LCR_fit1 <- LCR_Gibbs_cov_sel(X = X, Y = Y, G = 2, theta_prior_param = c(1,1), beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, burnin = 2000, thinby = 10, n_samples = 10000, verbose = TRUE)
LCR_fit1 <- LCR_Gibbs_both_sel(X = X, Y = Y, G = 2, theta_prior_param = c(1,1), beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, burnin = 2000, thinby = 10, n_samples = 5000, verbose = TRUE)
cov_incl_prop1 <- apply(LCR_fit1$gamma_samples, 1, mean)
cov_incl_prop1
coincidence <- ((LCR_fit1$gamma_samples)%*%t(LCR_fit1$gamma_samples))/(dim(LCR_fit1$gamma_samples)[2])

#Gonna see how much these variables actually impact the probabilities.
mu_excl1 <- cbind(1,X)[,-c(2)]%*%realbeta[-c(2)]
logit_prob_excl1 <- exp(mu_excl1)/(1+exp(mu_excl1))
logit_prob_excl1 <- cbind(logit_prob_excl1, 1-logit_prob_excl1)
logit_prob[1:10,]
logit_prob_excl1[1:10,]

mu_excl2 <- cbind(1,X)[,-c(3)]%*%realbeta[-c(3)]
logit_prob_excl2 <- exp(mu_excl2)/(1+exp(mu_excl2))
logit_prob_excl2 <- cbind(logit_prob_excl2, 1-logit_prob_excl2)
logit_prob[1:10,]
logit_prob_excl2[1:10,]

mu_excl3 <- cbind(1,X)[,-c(4)]%*%realbeta[-c(4)]
logit_prob_excl3 <- exp(mu_excl3)/(1+exp(mu_excl3))
logit_prob_excl3 <- cbind(logit_prob_excl3, 1-logit_prob_excl3)
logit_prob[1:10,]
logit_prob_excl3[1:10,]

mu_excl4 <- cbind(1,X)[,-c(5)]%*%realbeta[-c(5)]
logit_prob_excl4 <- exp(mu_excl4)/(1+exp(mu_excl4))
logit_prob_excl4 <- cbind(logit_prob_excl4, 1-logit_prob_excl4)
logit_prob[1:10,]
logit_prob_excl4[1:10,]

js_div_1 <- js_divergence(logit_prob,logit_prob_excl1)
js_div_2 <- js_divergence(logit_prob,logit_prob_excl2)
js_div_3 <- js_divergence(logit_prob,logit_prob_excl3)
js_div_4 <- js_divergence(logit_prob,logit_prob_excl4)

kl_div_1 <- kl_divergence(logit_prob,logit_prob_excl1)
kl_div_2 <- kl_divergence(logit_prob,logit_prob_excl2)
kl_div_3 <- kl_divergence(logit_prob,logit_prob_excl3)
kl_div_4 <- kl_divergence(logit_prob,logit_prob_excl4)

dens1 <- density(js_div_1)
dens2 <- density(js_div_2)
dens3 <- density(js_div_3)
dens4 <- density(js_div_4)

dens1_kl <- density(kl_div_1)
dens2_kl <- density(kl_div_2)
dens3_kl <- density(kl_div_3)
dens4_kl <- density(kl_div_4)

plot(dens1, col = "blue", lwd = 2, main = "JS Divergence Densities", xlab = "Value", ylim = c(0, max(dens1$y, dens2$y, dens3$y, dens4$y)))
lines(dens2, col = "red", lwd = 2)
lines(dens3, col = "green", lwd = 2)
lines(dens4, col = "purple", lwd = 2)

plot(dens1_kl, col = "blue", lwd = 2, main = "KL Divergence Densities", xlab = "Value", ylim = c(0, max(dens1_kl$y, dens2_kl$y, dens3_kl$y, dens4_kl$y)))
lines(dens2_kl, col = "red", lwd = 2)
lines(dens3_kl, col = "green", lwd = 2)
lines(dens4_kl, col = "purple", lwd = 2)

# LCR_fit2 <- LCR_Gibbs_cov_sel(X = X, Y = Y, G = 2, theta_prior_param = c(1,1), beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, burnin = 1000, thinby = 10, n_samples = 5000, verbose = TRUE)
# cov_incl_prop2 <- apply(LCR_fit2$gamma_samples, 1, mean)
# LCR_fit3 <- LCR_Gibbs_cov_sel(X = X, Y = Y, G = 2, theta_prior_param = c(1,1), beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, burnin = 1000, thinby = 10, n_samples = 5000, verbose = TRUE)
# cov_incl_prop3 <- apply(LCR_fit3$gamma_samples, 1, mean)
# LCR_fit4 <- LCR_Gibbs_cov_sel(X = X, Y = Y, G = 2, theta_prior_param = c(1,1), beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, burnin = 1000, thinby = 10, n_samples = 5000, verbose = TRUE)
# cov_incl_prop4 <- apply(LCR_fit4$gamma_samples, 1, mean)
# LCR_fit5 <- LCR_Gibbs_cov_sel(X = X, Y = Y, G = 2, theta_prior_param = c(1,1), beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, burnin = 1000, thinby = 10, n_samples = 5000, verbose = TRUE)
# cov_incl_prop5 <- apply(LCR_fit5$gamma_samples, 1, mean)


js_divergence <- function(P, Q) {
  eps <- 1e-12
  P <- pmax(P, eps)
  Q <- pmax(Q, eps)
  M <- 0.5 * (P + Q)
  0.5 * rowSums(P * log(P / M)) + 0.5 * rowSums(Q * log(Q / M))
}

kl_divergence <- function(P, Q) {
  eps <- 1e-12
  P <- pmax(P, eps)
  Q <- pmax(Q, eps)
  rowSums(P * log(P / Q))
}




LCR_Gibbs <- function(X, Y, ){
  args <- match.call()
  init <- 
}













