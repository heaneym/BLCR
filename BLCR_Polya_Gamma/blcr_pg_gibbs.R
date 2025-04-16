library(BayesLogit)
library(MCMCpack)
library(nnet)
library(abind)
source("C:/Users/matth/Documents/AIM CP Project/LCA/LCA/BLCR/BLCR_Polya_Gamma/functions_for_polya_gamma_BLCR.R")

blcr_pg_gibbs <- function(X, Y, G = 2, burnin = 500, n_samples = 1000, thin = 1, verbose = TRUE){
  #X is an (nxp)-matrix of covariates.
  #Y is an (nxM)-matrix of categorical responses.
  #G is the number of latent classes.
  
  #getting the number of iterations we need
  
  n_iter <- thin*n_samples + burnin
  
  #Initialising variables
  init <- initialise_variables_polyagamma(G = G,X = X, Y = Y)
  Y_flat <- init$Y_flat
  I <- init$I
  J <- init$J
  G_0 <- init$G_0
  K_0 <- init$K_0
  M <- init$M
  n <- init$n
  z <- init$z
  w <- init$w
  alpha <- init$alpha
  K <- init$K
  theta <- init$theta
  gamma <- init$gamma
  X <- init$X
  Y <- init$Y
  beta <- init$beta
  C <- init$C
  p <- init$p
  beta_prior_cov <- init$beta_prior_cov
  beta_prior_cov_inv <- init$beta_prior_cov_inv
  beta_prior_mean <- init$beta_prior_mean
  beta_cov <- init$beta_cov
  beta_cov_inv <- init$beta_cov_inv
  beta_cov_inv_chol <- init$beta_cov_inv_chol
  beta_mean <- init$beta_mean
  eta <- init$eta
  omega <- init$omega
  Omega <- init$Omega
  
  beta_samples <- array(0, dim = c(dim(beta),n_iter))
  omega_samples <- array(0, dim = c(dim(omega),n_iter))
  w_samples <- array(0, dim = c(dim(w),n_iter))
  z_samples <- array(0, dim = c(dim(z),n_iter))
  theta_samples <- array(0, dim = c(dim(theta),n_iter))
  
  beta_samples_unfixed <- array(0, dim = c(dim(beta),n_iter))
  omega_samples_unfixed <- array(0, dim = c(dim(omega),n_iter))
  w_samples_unfixed <- array(0, dim = c(dim(w),n_iter))
  z_samples_unfixed <- array(0, dim = c(dim(z),n_iter))
  theta_samples_unfixed <- array(0, dim = c(dim(theta),n_iter))
  
  cost_matrices <- array(0, dim = c(G,G,n_iter))
  perms <- matrix(0,nrow = G, ncol = n_iter)
  permuted_cost_matrices <-array(0, dim = c(G,G,n_iter))
  
  #Update steps
  
  for (count in 1:n_iter){
    
    if(count == (burnin+1)){
      z_init <- z
    }
    
    S <- S_update(Y_flat = Y_flat,I = I, J = J, G_0 = G_0, K_0 = K_0, M = M, z = z, G = G, K = K)
    
    theta <- theta_update(alpha = alpha, G = G, M = M, K = K, S = S, theta = theta)
    
    mu <- mu_update(X,beta)
    
    exp_mu <- exp(mu)
    
    C <- C_update(mu,exp_mu, G = G,C)
    
    eta <- eta_update(mu,C)
    
    gamma <- gamma_update(eta)
    
    log_gamma <- log(gamma)
    
    log_theta <- log(theta)
    
    w <- w_update(log_gamma,log_theta, Y, G,n,M,K)
    
    z <- z_update(w)
    
    kappa <- kappa_update(z)
    
    omega_up <- omega_update(eta,mu, G = G, n = n, omega = omega, Omega = Omega)
    
    omega <- omega_up$omega
    
    Omega <- omega_up$Omega
    
    beta_up <- beta_cov_mean_update(Omega,beta_prior_cov, beta_prior_mean, kappa, C,X, G = G, beta_prior_cov_inv = beta_prior_cov_inv, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, beta_cov = beta_cov, beta_mean = beta_mean)
    
    beta_cov <- beta_up$beta_cov
    
    beta_mean <- beta_up$beta_mean
    
    beta <- beta_update(beta_mean,beta_cov, G, beta = beta)
    
    if (count>=burnin+1){
      A_mat <- t(z_init)%*%z
      optimal_permutation <- matchClasses(A_mat, method = "exact", verbose = FALSE)
      perms[,count] <- optimal_permutation
      cost_matrices[,,count] <- A_mat
      permuted_A_mat <- A_mat[, optimal_permutation]
      permuted_cost_matrices[,,count] <- permuted_A_mat 
      
      beta_samples_unfixed[,,count] <- beta
      beta_samples[,,count] <- beta[optimal_permutation,]
      
      omega_samples_unfixed[,,count] <- omega
      omega_samples[,,count] <- omega[,optimal_permutation]
      
      z_samples_unfixed[,,count] <- z
      z_samples[,,count] <- z[,optimal_permutation]
      
      w_samples_unfixed[,,count] <- w
      w_samples[,,count] <- w[,optimal_permutation]
      
      theta_samples_unfixed[,,,count] <- theta
      theta_samples[,,,count] <- theta[optimal_permutation,,]
    } else {
      beta_samples[,,count] <- beta
      
      omega_samples[,,count] <- omega
      
      z_samples[,,count] <- z
      
      w_samples[,,count] <- w
      
      theta_samples[,,,count] <- theta
    }
    
    if (verbose == TRUE && count %% 500 == 0){
      print(paste(count,' of ',n_iter, ' completed.'))
    }
  }
  
  #post-processing steps
  #removing burnin
  beta_samples_burned <- beta_samples[,,-(1:burnin)]
  omega_samples_burned <- omega_samples[,,-(1:burnin)]
  z_samples_burned <- z_samples[,,-(1:burnin)]
  w_samples_burned <- w_samples[,,-(1:burnin)]
  theta_samples_burned <- theta_samples[,,,-(1:burnin)]
  #thinning
  beta_samples_clean <- beta_samples_burned[,,seq(1,dim(beta_samples_burned)[3], by = thin)]
  omega_samples_clean <- omega_samples_burned[,,seq(1,dim(omega_samples_burned)[3], by = thin)]
  z_samples_clean <- z_samples_burned[,,seq(1,dim(z_samples_burned)[3], by = thin)]
  w_samples_clean <- w_samples_burned[,,seq(1,dim(w_samples_burned)[3], by = thin)]
  theta_samples_clean <- theta_samples_burned[,,,seq(1,dim(theta_samples_burned)[4], by = thin)]
  
  #parameter means and standard deviations
  beta_estimate <- apply(beta_samples_clean,c(1,2), mean)
  omega_estimate <- apply(omega_samples_clean,c(1,2), mean)
  theta_estimate <- apply(theta_samples_clean,c(1,2,3), mean)
  beta_sd <- apply(beta_samples_clean,c(1,2), sd)
  omega_sd <- apply(omega_samples_clean,c(1,2), sd)
  theta_sd <- apply(theta_samples_clean,c(1,2,3), sd)
  
  mu_est <- X%*%t(beta_estimate)
  exp_mu_est <- exp(mu_est)
  gamma_est <- exp_mu_est/rowSums(exp_mu_est)
  class_membership <- apply(gamma_est,1,which.max)
  
  pi <- numeric(G)
  for (g in 1:G){
    pi[g] <- length(class_membership[class_membership == g])/length(class_membership)
  }
  
  return(list(beta_estimate = beta_estimate, omega_estimate = omega_estimate, theta_estimate = theta_estimate, beta_sd = beta_sd, omega_sd = omega_sd, theta_sd = theta_sd, class_membership = class_membership, pi = pi, beta_samples = beta_samples_clean, omega_samples = omega_samples_clean, theta_samples = theta_samples_clean, perms = perms, cost_matrices = cost_matrices, permuted_cost_matrices = permuted_cost_matrices, beta_samples_unfixed, theta_samples_unfixed, omega_samples_unfixed, z_samples_unfixed, w_samples_unfixed))
  
  #probably a few other outputs that I could include in here as well
  
  #Note that the estimated beta is actually the transpose, since the sim_data function that I made takes the columns to be beta_g and the sampler takes the rows to be beta_g
  
  #This can also do with some optimisation since its not particularly fast at the moment
}