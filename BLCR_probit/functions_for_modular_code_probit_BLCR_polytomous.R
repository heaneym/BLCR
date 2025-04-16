library(MASS)
library(numDeriv)
#library(Rcpp)
#library(RcppArmadillo)
#sourceCpp("./BLCR/BLCR_probit/BLCR_probit_cpp_functions.cpp", verbose = TRUE)

#This file contains the functions for the multinomial probit LCR model

#Update for S
S_update <- function(Y_flat,I,J,G_0,K_0,M,z){
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
  return(S)
}


#Update function for theta
theta_update <- function(alpha,G,M,K,S,theta){
  sapply(1:G, function(g) {
    sapply(1:M, function(j) {
      # Generate theta values for each g and j
      theta[g, j, 1:K[j]] <<- rdirichlet(1, alpha[1:K[j]] + S[g, j, 1:K[j]])
    })
  })
  return(theta)
}


#This is to update C - it is very slow
C_update <- function(z, mu, Sigma, n, G, C) {
  C <- C
  max_indices <- max.col(z)
  for (i in 1:n) {
    repeat {
      C_prop <- mvrnorm(1, mu = mu[i,], Sigma = Sigma)
      if (which.max(C_prop) == max_indices[i]) {
        C[i,] <- C_prop
        print(i)
        break
      }
    }
  }
  return(C)
}


# Function to calculate proportions for each row
calculate_proportions <- function(row) {
  table_counts <- table(factor(row, levels = 1:3)) 
  proportions <- table_counts / length(row)      
  return(as.numeric(proportions))                 
}



#Function to get the multinomial probit probabilities - This needs to be optimised later on

#gamma_update <- function(mu, Sigma, n,G,n_samples){
  #gamma <- matrix(0,nrow = n,ncol = G)
  #vec_epsilon_samples <- mvrnorm(n_samples, mu = rep(0, n*G), Sigma = kronecker(diag(n),Sigma))
  #epsilon_samples <- array(0, dim = c(n,G,n_samples))
  #linear_predictors <- array(0, dim = c(n,G,n_samples))
  #chosen_group <- matrix(0,nrow = n,ncol = n_samples)
  #Optimise this!
  #for (sample in 1:n_samples){
    #epsilon_samples[,,sample] <- matrix(vec_epsilon_samples[sample,], nrow = n, byrow = FALSE)
    #linear_predictors[,,sample] <- mu + epsilon_samples[,,sample]
    #chosen_group[,sample] <- apply(linear_predictors[,,sample], 1, which.max)
  #}
  #gamma_mat <- t(apply(matrix_data, 1, calculate_proportions))
  #return(gamma_mat)
#}



#The gamma_update function below is the one that currently works but is quite slow

gamma_update <- function(mu, Sigma, n, G, n_samples) {
  vec_epsilon_samples <- mvrnorm(n_samples, mu = rep(0, n*G), Sigma = kronecker(diag(n), Sigma))
  epsilon_samples <- array(vec_epsilon_samples, dim = c(n, G, n_samples))
  linear_predictors <- sweep(epsilon_samples, c(1,2), mu, '+')
  chosen_group <- t(apply(linear_predictors, 3, function(x) apply(x, 1, which.max)))
  gamma_mat <- t(apply(chosen_group, 1, calculate_proportions))
  return(gamma_mat)
}


#This calls the gamma_update from a cpp implementation - cannot get this to work - issues with cpp file paths and stuff

#gamma_update <- function(mu, Sigma, n, G, n_samples) {
#  gamma_update_cpp(mu, Sigma, n, G, n_samples)
#}





#Function for updating w - COME BACK TO THIS - could do with optimising
w_update <- function(gamma,theta){
  log_gamma <- log(gamma)
  log_theta <- log(theta)
  log_sum <- matrix(0, nrow = n, ncol = G)
  for (i in 1:n){
    for(g in 1:G){
      for (j in 1:M){
        y_ij = Y[i,j]
        log_sum[i,g] <- log_sum[i,g] + (y_ij <= K[j])*log_theta[g,j,y_ij]
      }
    }
  }
  log_w <- log_gamma + log_sum
  w <- exp(log_w)
  return(w)
}




z_update <- function(w){
  z <- t(apply(w, 1, function(row) rmultinom(1, 1, row)))
  return(z)
}



# should the last row of beta be 0? Also should maybe get the Omega and Omega inverse matrix outside of the function since we need it for the eta update as well. Same for C_vec
beta_update <- function(U,Omega,Omega_inv,C_vec){
  mat1 <- t(U)%*%(Omega_inv%*%U)
  mat1_chol <- chol(mat1)
  beta_cov <- chol2inv(mat1_chol)
  mean_factor <- beta_cov%*%t(U)%*%Omega_inv
  beta_mean <- mean_factor%*%C_vec
  vec_beta <- MASS::mvrnorm(1, mu = beta_mean, Sigma = beta_cov)
  beta <- matrix(vec_beta, nrow = G, byrow = FALSE) 
  return(beta)
}


#Need to have precomputed Omega_inv and Omega_det,C_vec,mu
eta_density_unnormalised <- function(eta,C_vec,mu,n){
  Sig_up <- Sigma_matrix_update(eta, G)
  Sigma <- Sig_up$Sigma
  Omega <- kronecker(diag(n), Sigma)
  Omega_chol <- chol(Omega)
  Omega_inv <- chol2inv(Omega_chol)
  log_Omega_det <- -2 * sum(log(diag(Omega_chol)))  # Log determinant
  mu_vec <- as.vector(mu)
  exp_arg <- (-1/2) * t(C_vec - mu_vec) %*% (Omega_inv %*% (C_vec - mu_vec))
  log_result <- 0.5 * log_Omega_det + exp_arg
  return(log_result)
}

#eta_density_wrapper <- function(eta) {
  #return(-eta_density_unnormalised(eta, C_vec, mu, n))
#}

eta_density_wrapper <- function(eta) {
  return(-eta_density_unnormalised(eta, C_vec, mu, n))
}

find_mode <- function(f, init_eta) {
  result <- optimParallel(
    par = init_eta,
    fn = f,
    method = "BFGS",
    control = list(fnscale = 1, trace = 2),
    parallel = list(cl = cl, forward = TRUE),
    hessian = TRUE
  )
  return(list(eta = result$par, hessian = result$hessian)) 
}


#compute_hessian <- function(f, eta_mode) {
  #hessian_matrix <- hessian(f, eta_mode, method = "forward")
  #return(hessian_matrix)
#}



eta_update <- function(init_eta){
  mode_calc <- find_mode(eta_density_wrapper,init_eta)
  #eta_mode <- find_mode(eta_density_wrapper,init_eta)
  eta_mode <- mode_calc$eta
  hessian_matrix <- mode_calc$hessian
  print(hessian_matrix)
  #hessian_matrix <- compute_hessian(eta_density_wrapper,eta_mode)
  #hessian_matrix_chol <- chol(hessian_matrix)
  hessian_matrix_inv <- solve(hessian_matrix)
  #hessian_matrix_inv <- chol2inv(hessian_matrix_chol)
  result <- mvrnorm(1,mu = eta_mode, Sigma = hessian_matrix_inv)
  return(list(eta = result, eta_mode = eta_mode))
  return(result)
  #mode_Sigma_chol <- t(chol(mode_Sigma))
  #lower_indices <- which(row(mode_Sigma_chol) > col(mode_Sigma_chol))
  #mode_eta <- mode_Sigma_chol[lower_indices]
  
  #want to get the hessian at that mode and then get the corresponding normal and draw from that
}



Sigma_matrix_initialise <- function(G){
  eta_length <- (G*(G-1))/2
  eta <- rnorm(eta_length)
  L_components <- exp(eta)
  strictly_lower_mat <- matrix(0,nrow = G,ncol = G)
  lower_indices <- which(row(strictly_lower_mat) > col(strictly_lower_mat))
  strictly_lower_mat[lower_indices] <- L_components
  L <- diag(G) + strictly_lower_mat
  Sigma <- L%*%t(L)
  result <- list(eta = eta, Sigma = Sigma)
}


#Im not entirely sure how the Sigma matrix should be parametrised by the eta vector
Sigma_matrix_update <- function(eta,G){
  #L_components <- exp(eta)
  L_components <- eta
  strictly_lower_mat <- matrix(0,nrow = G,ncol = G)
  lower_indices <- which(row(strictly_lower_mat) > col(strictly_lower_mat))
  strictly_lower_mat[lower_indices] <- L_components
  L <- diag(G) + strictly_lower_mat
  Sigma <- L%*%t(L)
  result <- list(eta = eta, Sigma = Sigma)
  return(result)
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
















initialise_variables <- function(G,a = NULL, X,Y){
  p <- ncol(X)
  if (is.vector(X)){
    ones <- rep(1,length(X))
    X <- cbind(ones,X)
  } else{
    ones <- rep(1,nrow(X))
    X <- cbind(ones,X)
  }
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
  
  beta <- rbind(matrix(rnorm((G-1)*(p+1),2,0.5), nrow = (G-1)),0)
  
  
  #-----------------------------------------------------------------------------
  #-------------------------eta and Sigma initialisation------------------------
  #-----------------------------------------------------------------------------
  
  A <- Sigma_matrix_initialise(G)
  eta <- A$eta
  Sigma <- A$Sigma
  
  #-----------------------------------------------------------------------------
  #-------------------------z initialisation------------------------------------
  #-----------------------------------------------------------------------------
  
  w <- array(0, dim = c(n,G)) 
  z <- array(0, dim = c(n,G)) 
  for(i in 1:n){
    w[i,] <- rep(1/G,G)
    z[i,] <- rmultinom(1,1,w[i,])
  }
  
  
  
  
  #-----------------------------------------------------------------------------
  #-------------------------C initialisation------------------------------------
  #-----------------------------------------------------------------------------
  
  C <- array(rnorm(n*G), dim = c(n,G))
  
  #Precalculating various quantities needed in the sampler
  
  Y_flat <- as.vector(t(Y))
  I <- rep(1:n, each = M * G)
  J <- rep(rep(1:M, each = G), n)
  G_0 <- rep(1:G, times = n * M)
  K_0 <- rep(Y_flat, each = G)
  
  return(list(Y = Y, X = X, theta = theta, alpha  = alpha, beta = beta, w = w, z = z, C = C, n = n, M = M, K = K, p = p, Y_flat = Y_flat, I = I, J = J, G_0 = G_0, K_0 = K_0, eta = eta, Sigma = Sigma))
  
}












