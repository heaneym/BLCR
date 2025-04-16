initialise_variables_polyagamma_varsel <- function(G,a = NULL, X,Y, beta_prior_mean = NULL, beta_prior_cov = NULL){
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
  if (is.null(a)){
    alpha <- rep(1, max(K))
  } else {
    if(length(a) != max(K)){
      print("The parameter a has the wrong length")
    } else{
      alpha <- a
    }
  }
  for (g in 1:G) { #generating a sample from theta_{gj vectors for each g and j}
    for (j in 1:M) {
      theta[g, j, 1:K[j]] <- rdirichlet(1, alpha[1:K[j]])  
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
  omega <- matrix(rpg(n*G), nrow = n)
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
    beta_prior_cov <- diag(p+1)
    beta_prior_cov_inv <- solve(beta_prior_cov)
  } else {
    beta_prior_cov_inv <- solve(beta_prior_cov)
  
  }
  
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
  
  
  
  
  #Precalculating various quantities needed in the sampler
  
  Y_flat <- as.vector(t(Y))
  I <- rep(1:n, each = M * G)
  J <- rep(rep(1:M, each = G), n)
  G_0 <- rep(1:G, times = n * M)
  K_0 <- rep(Y_flat, each = G)
  
  return(list(Y = Y, X = X, X_current =X, theta = theta, alpha  = alpha, beta = beta, w = w, z = z, C = C, n = n, M = M,G = G, K = K, p = p, Y_flat = Y_flat, I = I, J = J, G_0 = G_0, K_0 = K_0, omega = omega,beta_cov = beta_cov, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov ,beta_prior_cov_inv = beta_prior_cov_inv , gamma = gamma, beta_mean = beta_mean, eta = eta, logit_probs = logit_probs, tau = tau, mu = mu))
  
}



S_update <- function(Y,z,M,G,K){
  Y_indicator <- array(0, dim = c(n, M, max(K)))
  Y_indicator[cbind(rep(1:nrow(Y), ncol(Y)), rep(1:ncol(Y), each = nrow(Y)), as.vector(Y))] <- 1
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




logit_probs_update <- function(eta){
  logit_probs <- exp(eta)/(1+exp(eta)) 
  return(logit_probs)
}

w_update <- function(log_logit_probs,log_theta,Y,G,n,M,K){
  log_theta_list <- lapply(1:M, function(j) {
    log_theta[, j, Y[, j]]
  })
  sum_log_theta <- Reduce(`+`, log_theta_list)
  log_w <- log_logit_probs + t(sum_log_theta)
  w <- exp(log_w)
  return(w)
}

z_update <- function(w){
  z <- t(apply(w, 1, function(row) rmultinom(1, 1, row)))
  return(z)
}

mu_update <- function(X,beta){
  return(X%*%beta)
}

eta_update <- function(mu,C){
  return(mu - C)
}

kappa_update <- function(z){
  return(z - 1/2)
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

#C_update <- function(mu,exp_mu,G,C){
  #if (ncol(exp_mu) > 2) {
    # Use apply to vectorize the operation across columns
    #epsilon <- 1e-10
    #C <- t(apply(exp_mu, 1, function(row) log(pmax(sum(row) - row, epsilon))))
  #} else {
    # For the case when ncol(exp_mu) <= 2, we can directly assign
    #C <- mu[,c(2,1)]
  #}
 #return(C)
#}

omega_update <- function(eta,mu,G,n,omega){
  omega <- matrix(rpg(n * G, 1, as.vector(eta)), nrow = n, ncol = G)
  return(omega)
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










beta_update <- function(gamma, X_current, G, p, beta_prior_cov_inv, A){
  beta <- matrix(0, nrow = (p+1), ncol = G)
  nonzero_indices <- which(gamma == 1)
  beta_prior_cov_inv1 <- beta_prior_cov_inv[nonzero_indices,nonzero_indices]
  beta_prior_cov1 <- chol2inv(chol(beta_prior_cov_inv1))
  beta_prior_mean1 <- beta_prior_mean[nonzero_indices]
  for(g in 1:(G-1)){
    #can probably remove this calculation here
    beta_cov_inv <- crossprod(X_current, omega[,g]*X_current) + beta_prior_cov_inv1
    beta_cov_inv_chol <- chol(beta_cov_inv)
    #A <- kappa[,g] + omega[,g]*C[,g]
    fs <- forwardsolve(t(beta_cov_inv_chol),crossprod(X_current,A[,g]) + beta_prior_cov_inv1%*%beta_prior_mean1)
    beta_mean <- backsolve(beta_cov_inv_chol, fs)
    beta[nonzero_indices,g] <- rMVNormP(1,beta_mean,beta_cov_inv_chol)
  }  
  return(beta)
}

tau_update <- function(gamma,a,b){
  p <- length(gamma)-1
  sum_gamma <- sum(gamma[2:(p+1)])
  tau <- rbeta(1,sum_gamma + a, p - sum_gamma + b)
  return(tau)
}

#function for drawing a mvnorm sample using the Cholesky decomp of precision rather than sigma

rMVNormP <- function(n, mu, chol_precision) {
  p <- length(mu)
  L <- chol_precision
  Z <- matrix(rnorm(p*n), p, n)
  Y <- backsolve(L, Z, transpose = TRUE)
  X <- sweep(Y, 1, mu, FUN = "+")
  return(X)
}


log_beta_function <- function(x){
  sum_x <- sum(x)
  gamma_sum_x <- gamma(sum_x)
  log_gamma_sum_x <- log(gamma_sum_x)
  gamma_x <- gamma(x)
  log_gamma_x <- log(gamma_x)
  sum_log_gamma_x <- sum(log_gamma_x)
  return(sum_log_gamma_x - log_gamma_sum_x)
}



  






