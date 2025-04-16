PG_varsel_mult_init <- function(X,Y, beta_prior_mean = NULL, beta_prior_cov = NULL){
  if (is.vector(X)){
    ones <- rep(1,length(X))
    X <- cbind(ones,X)
  } else{
    ones <- rep(1,nrow(X))
    X <- cbind(ones,X)
  }
  p <- ncol(X) - 1
  G <- length(unique(Y))
  n <- length(Y)
  
  
  #-----------------------------------------------------------------------------
  #-------------------------beta initialisation---------------------------------
  #-----------------------------------------------------------------------------
  
  beta <- cbind(matrix(rnorm((G-1)*(p+1)), ncol = (G-1)),0)
  
  
  #-----------------------------------------------------------------------------
  #-------------------------omega, Omega,...initialisation here-----------------
  #-----------------------------------------------------------------------------
  
  C <- matrix(rnorm((n)*(G)), nrow = n)
  omega <- matrix(rpg(n*G), nrow = n)
  Omega <- array(0, dim = c(n,n,(G)))
  for(g in 1:(G-1)){
    Omega[,,g] <- diag(omega[,g])
  }
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
    beta_prior_cov_inv <- beta_prior_cov
  }
  
  
  #-----------------------------------------------------------------------------
  #------------------gamma initialisation here----------------------------------
  #-----------------------------------------------------------------------------
  
  gamma <- rep(1, p+1)
  
  
  tau <- runif(1)
  
  #need to create one-hot z matrix as well
  
  z <- model.matrix(~ factor(Y) - 1)
  
  
  return(list(Y = Y, X = X, p = p, G = G, n = n, beta = beta, C = C, omega = omega, Omega = Omega, eta = eta, beta_cov = beta_cov, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, beta_mean = beta_mean, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov, beta_prior_cov_inv = beta_prior_cov_inv, gamma = gamma, z = z, X_current = X, tau = tau))
  
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


#beta_cov_mean_update <- function(Omega, beta_prior_cov, beta_prior_mean, kappa,C,X, G, beta_prior_cov_inv, beta_cov_inv, beta_cov_inv_chol, beta_cov, beta_mean){
  #for (g in 1:(G-1)){
    #beta_cov_inv[,,g] <- t(X)%*%(Omega[,,g]%*%X) + beta_prior_cov_inv
    #beta_cov_inv_chol[,,g] <- chol(beta_cov_inv[,,g])
    #beta_cov[,,g] <- chol2inv(beta_cov_inv_chol[,,g])
    #A <- kappa[,g] + Omega[,,g]%*%C[,g]
    #beta_mean[,g] <- beta_cov[,,g]%*%(t(X)%*%A - beta_prior_cov_inv%*%beta_prior_mean)
  #}
  #return(list(beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, beta_cov = beta_cov, beta_mean = beta_mean))
#}




#beta_update <- function(beta_mean, beta_cov, G = G, beta){
  #for (g in 1:(G-1)){
    #beta[,g] <- mvrnorm(1,mu = beta_mean[,g], Sigma = beta_cov[,,g])
  #}
  #return(beta)
#}


C_update <- function(mu,exp_mu,G,C){
  if (ncol(exp_mu)>2){
    for (g in 1:G){
      C[, g] <- log(rowSums(exp_mu[,-g, drop = FALSE]))
    }
  } else {
    for (g in 1:G){
      C[, g] <- mu[,-g, drop = FALSE]
    }
  }
  return(C)
}

omega_update <- function(eta,mu,G,n,omega,Omega){
  for(g in 1:G){
    for (i in 1:n){
      omega[i,g] <- rpg(1,1,eta[i,g])
    }
    Omega[,,g] <- diag(omega[,g])
  }
  return(list(omega = omega, Omega = Omega))
}


#gamma_update <- function(){
  #if (p == 1){
    #gamma_prop_index <- 2
  #} else {
    #gamma_prop_index <- sample((2:(p+1)),1)
  #}
  #gamma_prop <- gamma
  #if(gamma[gamma_prop_index] == 1){
    #gamma[gamma_prop_index] <- 0
  #} else {
    #gamma[gamma_prop_index] <- 1
  #}
  #X_current <- X[gamma == 1, drop = FALSE]
  #X_prop <- X[gamma_prop == 1, drop = FALSE]
  #beta_prior_precision_prop <- beta_prior_precision[gamma_prop == 1, gamma_prop ==1, drop = FALSE]
  #beta_prior_precision_current <- beta_prior_precision[gamma == 1, gamma == 1, drop = FALSE]
  #beta_prior_mean_prop <- beta_prior_mean[gamma_prop == 1, drop = FALSE]
  #beta_prior_mean_current <- beta_prior_precision[gamma == 1, drop = FALSE]
#}

gamma_update <- function(gamma,p,beta_prior_cov_inv, beta_prior_mean, kappa, G, X_current, tau){
  gamma_prop_index <- sample(2:(p+1), 1)
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
  beta_cov_prop <- array(0, dim = c(dim(beta_prior_cov_inv_prop),G-1))
  beta_cov_current <- array(0, dim = c(dim(beta_prior_cov_inv_current),G-1))
  beta_mean_prop <- array(0, dim = c(G-1,length(beta_prior_mean_prop)))
  beta_mean_current <- array(0, dim = c(G-1,length(beta_prior_mean_current)))
  
  for (g in 1:(G-1)){
    beta_cov_inv_prop[,,g] <- t(X_prop)%*%(Omega[,,g]%*%X_prop) + beta_prior_cov_inv_prop
    beta_cov_inv_current[,,g] <- t(X_current)%*%(Omega[,,g]%*%X_current) + beta_prior_cov_inv_current
    beta_cov_inv_chol_prop[,,g] <- chol(beta_cov_inv_prop[,,g])
    beta_cov_inv_chol_current[,,g] <- chol(beta_cov_inv_current[,,g])
    beta_cov_prop[,,g] <- chol2inv(beta_cov_inv_chol_prop[,,g])
    beta_cov_current[,,g] <- chol2inv(beta_cov_inv_chol_current[,,g])
    A <- kappa[,g] + Omega[,,g]%*%C[,g]
    beta_mean_prop[g,] <- beta_cov_prop[,,g]%*%(t(X_prop)%*%A + beta_prior_cov_inv_prop%*%beta_prior_mean_prop)
    beta_mean_current[g,] <- beta_cov_current[,,g]%*%(t(X_current)%*%A + beta_prior_cov_inv_current%*%beta_prior_mean_current)
  }
  beta_cov_inv_list_prop <- lapply(1:dim(beta_cov_inv_prop)[3], function(g) beta_cov_inv_prop[,,g])
  beta_cov_inv_block_prop <- do.call(bdiag, beta_cov_inv_list_prop)
  beta_cov_inv_list_current <- lapply(1:dim(beta_cov_inv_current)[3], function(g) beta_cov_inv_current[,,g])
  beta_cov_inv_block_current <- do.call(bdiag, beta_cov_inv_list_current)
  
  beta_cov_list_prop <- lapply(1:dim(beta_cov_prop)[3], function(g) beta_cov_prop[,,g])
  beta_cov_block_prop <- do.call(bdiag, beta_cov_list_prop)
  beta_cov_list_current <- lapply(1:dim(beta_cov_current)[3], function(g) beta_cov_current[,,g])
  beta_cov_block_current <- do.call(bdiag, beta_cov_list_current)
  
  beta_mean_prop_vec <- as.vector(beta_mean_prop)
  beta_mean_current_vec <- as.vector(beta_mean_current)
  det_cov_prop <- det(beta_cov_block_prop)
  det_cov_current <- det(beta_cov_block_current)
  log_det_current <- log(det_cov_current)
  log_det_prop <- log(det_cov_prop)
  
  log_accept_ratio1 <- (1/2)*(t(beta_mean_prop_vec)%*%(beta_cov_inv_block_prop%*%beta_mean_prop_vec) - t(beta_mean_current_vec)%*%(beta_cov_inv_block_current%*%beta_mean_current_vec) + log_det_prop - log_det_current)
  
  sum_gamma_prop <- sum(gamma_prop[2:(p+1)])
  sum_gamma_current <- sum(gamma[2:(p+1)])
  
  #print((sum_gamma_prop-sum_gamma_current)*log(tau) + (sum_gamma_current-sum_gamma_prop)*log(1-tau))
  
  log_accept_ratio <- log_accept_ratio1 + (sum_gamma_prop-sum_gamma_current)*log(tau) + (sum_gamma_current-sum_gamma_prop)*log(1-tau)
  
  accept_ratio <- exp(log_accept_ratio)
  accept_ratio <- as.numeric(accept_ratio)
  #print(accept_ratio)
  #print(gamma_prop_index)
  if (runif(1) <= min(1,accept_ratio)){
    gamma <- gamma_prop
    X_current <- X_prop
    beta_mean <- beta_mean_prop
    beta_cov <- beta_cov_prop
    beta_cov_inv <- beta_cov_inv_prop
  }
  return(list(gamma = gamma,X_current = X_current, beta_mean = beta_mean, beta_cov = beta_cov, beta_cov_inv = beta_cov_inv))
}


beta_update <- function(gamma, X_current, G, p){
  beta <- matrix(0, nrow = (p+1), ncol = G)
  for(g in 1:(G-1)){
    nonzero_indices <- which(gamma == 1)
    #print(nonzero_indices)
    beta_prior_cov_inv1 <- beta_prior_cov_inv[nonzero_indices,nonzero_indices]
    #print(beta_prior_cov_inv1)
    beta_prior_cov1 <- chol2inv(chol(beta_prior_cov_inv1))
    #print(beta_prior_cov1)
    beta_prior_mean1 <- beta_prior_mean[nonzero_indices]
    #print(beta_prior_mean)
    #X_current <- X[,nonzero_indices, drop = FALSE]
    #print(ncol(X_current))
    beta_cov_inv <- t(X_current)%*%Omega[,,g]%*%X_current + beta_prior_cov_inv1
    #print(beta_cov_inv)
    beta_cov <- chol2inv(chol(beta_cov_inv))
    #print(beta_cov)
    beta_mean <- beta_cov%*%(t(X_current)%*%(kappa[,g] + Omega[,,g]%*%C[,g]) + beta_prior_cov_inv1%*%beta_prior_mean1)
    #print(beta_mean)
    beta[nonzero_indices,g] <- mvrnorm(n = 1, mu = beta_mean, Sigma = beta_cov)
  }  
  return(beta)
}

tau_update <- function(gamma,a,b){
  p <- length(gamma)-1
  sum_gamma <- sum(gamma[2:(p+1)])
  tau <- rbeta(1,sum_gamma + a, p - sum_gamma + b)
  return(tau)
}






