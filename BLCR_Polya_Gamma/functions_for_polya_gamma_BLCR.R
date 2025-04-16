initialise_variables_polyagamma <- function(G,a = NULL, X,Y, beta_prior_mean = NULL, beta_prior_cov = NULL){
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
  
  
  #-----------------------------------------------------------------------------
  #-------------------------omega, Omega,...initialisation here------------------------
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
  beta_mean <- matrix(0, nrow = (p+1), ncol = (G-1) )
  
  
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
  #-------------------------z initialisation------------------------------------
  #-----------------------------------------------------------------------------
  
  w <- array(0, dim = c(n,G)) 
  z <- array(0, dim = c(n,G)) 
  for(i in 1:n){
    w[i,] <- rep(1/G,G)
    z[i,] <- rmultinom(1,1,w[i,])
  }
  
  #gamma matrix update
  
  gamma <- matrix(0,nrow = n, ncol = G)
  
  
  
  #Precalculating various quantities needed in the sampler
  
  Y_flat <- as.vector(t(Y))
  I <- rep(1:n, each = M * G)
  J <- rep(rep(1:M, each = G), n)
  G_0 <- rep(1:G, times = n * M)
  K_0 <- rep(Y_flat, each = G)
  
  return(list(Y = Y, X = X1, theta = theta, alpha  = alpha, beta = beta, w = w, z = z, C = C, n = n, M = M,G = G, K = K, p = p, Y_flat = Y_flat, I = I, J = J, G_0 = G_0, K_0 = K_0, Omega = Omega, omega = omega,beta_cov = beta_cov, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, beta_prior_mean = beta_prior_mean, beta_prior_cov = beta_prior_cov ,beta_prior_cov_inv = beta_prior_cov_inv , gamma = gamma, beta_mean = beta_mean, eta = eta))
  
}




prepare_data_poLCA_covariates <- function(seed, n, M, p, beta,theta = NULL, G){
  #seed is the seed we want to simulate the data with, n is the number of observations
  #M is the number of clustering variables, p is the number of covariates, and beta is the vector of these covariates
  set.seed(seed)
  if(is.null(beta)){
    b <- NULL
  } else {
    b <- matrix(beta,nrow = p+1, ncol = G-1)
  }
  dat <- poLCA.simdata(N = n, probs = theta, nclass = G, ndv = M, nresp = NULL, x = NULL, niv = p, b = b, P = NULL, missval = FALSE, pctmiss = NULL)
  Y_columns <- lapply(1:M, function(i) dat[["dat"]][[paste0("Y", i)]])
  Y <- do.call(cbind, Y_columns)
  X_columns <- lapply(1:p, function(i) dat[["dat"]][[paste0("X", i)]])
  X <- do.call(cbind, X_columns)
  tc <- dat[["trueclass"]]
  return(list(tc = tc, Y = Y, X = X, data = dat))
}



S_update <- function(Y_flat,I,J,G_0,K_0,M,z,G,K){
  z_binded_M <- do.call(cbind, replicate(M, z, simplify = FALSE))
  z_flat <- as.vector(t(z_binded_M))
  S_df <- data.frame(
    i = I,         
    j = J,      
    g = G_0,         
    k = K_0,       
    z_value = z_flat 
  )
  S_df <- S_df[S_df$z_value == 1, ]
  S <- with(S_df, table(g, j, k))
  if(dim(S)[1] < G){
    newrows <- array(0,dim = c((G-dim(S)[1]),M,max(K)))
    S <- abind(S,newrows,along = 1)
  }
  return(S)
}

theta_update <- function(alpha,G,M,K,S,theta){
  sapply(1:G, function(g) {
    sapply(1:M, function(j) {
      # Generate theta values for each g and j
      theta[g, j, 1:K[j]] <<- rdirichlet(1, alpha[1:K[j]] + S[g, j, 1:K[j]])
    })
  })
  return(theta)
}

#gamma_update <- function(exp_mu){
  #gamma <- exp_mu/rowSums(exp_mu)
  #return(gamma)
#}

gamma_update <- function(eta){
  gamma <- exp(eta)/(1+exp(eta)) 
  return(gamma)
}

w_update <- function(log_gamma,log_theta,Y,G,n,M,K){
  #optimise this bit
  log_sum <- matrix(0, nrow = n, ncol = G)
  for (i in 1:n){
    for(g in 1:G){
      for (j in 1:M){
        y_ij = Y[i,j]
        #remember here that the y_ij can be zero in some binary datasets and similar, so make sure to account for this by addng 1
        log_sum[i,g] <- log_sum[i,g] + (y_ij <= K[j])*log_theta[g,j,y_ij]
      }
    }
  }
  log_w <- log_gamma + log_sum
  #log_w <- log_gamma + A
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

beta_cov_mean_update <- function(Omega, beta_prior_cov, beta_prior_mean, kappa,C,X, G, beta_prior_cov_inv, beta_cov_inv, beta_cov_inv_chol, beta_cov, beta_mean){
  beta_cov_inv <- array(0,dim=c(dim(beta_prior_cov),G-1))
  beta_cov_inv_chol <- array(0,dim = c(dim(beta_prior_cov),G-1))
  beta_cov <- array(0, dim = c(dim(beta_prior_cov),G-1))
  beta_mean <- matrix(0, nrow = ncol(X), ncol = G-1)
  for (g in 1:(G-1)){
    beta_cov_inv[,,g] <- t(X)%*%(Omega[,,g]%*%X) + beta_prior_cov_inv
    beta_cov_inv_chol[,,g] <- chol(beta_cov_inv[,,g])
    beta_cov[,,g] <- chol2inv(beta_cov_inv_chol[,,g])
    A <- kappa[,g] + Omega[,,g]%*%C[,g]
    beta_mean[,g] <- beta_cov[,,g]%*%(t(X)%*%A - beta_prior_cov_inv%*%beta_prior_mean)
  }
  return(list(beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, beta_cov = beta_cov, beta_mean = beta_mean))
}

beta_update <- function(beta_mean, beta_cov, G = G, beta){
  for (g in 1:(G-1)){
    beta[,g] <- mvrnorm(1,mu = beta_mean[,g], Sigma = beta_cov[,,g])
  }
  return(beta)
}


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




