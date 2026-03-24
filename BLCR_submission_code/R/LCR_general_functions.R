
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





