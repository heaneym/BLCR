LCR_init <- function(X, Y, G, beta_prior_mean, beta_prior_cov, theta_hyperparam, clust_var_prior){
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
  K <- apply(Y, 2, function(x) length(unique(x[!is.na(x)])))
  max_K <- max(K)
  theta <- array(0, dim = c(G, M, max(K)))  #Initialising theta
  for (g in 1:G) { 
    for (j in 1:M) {
        theta[g, j, 1:K[j]] <- rdirichlet(1,rep(theta_hyperparam,K[j]))
    }
  }
  beta <- cbind(matrix(rnorm((G-1)*(p+1)), nrow = (p+1)),0)
  mu <- X1%*%beta
  logit_probs <- t(apply(matrix(runif(n * G), nrow = n), 1, function(x) x / sum(x)))
  log_logit_probs <- log(logit_probs)
  C <- matrix(rnorm((n)*(G)), nrow = n)
  omega <- matrix(rpg(n*(G-1)), nrow = n)
  eta <- matrix(0, nrow = n, ncol = G)
  beta_cov <- array(0, dim = c(p+1,p+1,G-1))
  beta_cov_inv <- array(0, dim = c(p+1,p+1,G-1))
  beta_cov_inv_chol <- array(0, dim = c(p+1,p+1,G-1))
  beta_mean <- matrix(rnorm((G-1)*(p+1)), ncol = (G-1), nrow = (p+1))
  beta_prior_cov_inv <- solve(beta_prior_cov)
  w <- array(0, dim = c(n,G)) 
  for (i in 1:n) {
    w[i,] <- rep(1/G,G)
  }
  z <- array(0, dim = c(n, G))
  if (nrow(unique(X)) >= G){
    km <- kmeans(X, G)
    z_init_assign <- km$cluster
    z[cbind(1:n, z_init_assign)] <- 1
  } else {
    z_init_assign <- sample(1:G, size = n, replace = TRUE)
    z[cbind(1:n, z_init_assign)] <- 1
  }
  kappa <- z - 1/2
  gamma <- rep(1,p+1)
  tau <- runif(1)
  nu <- rep(1,M)
  Y_missing_indicator <- is.na(Y) * 1
  Y_for_indicator <- Y
  Y_for_indicator[is.na(Y)] <- 1
  dim_Y_indicator <- c(nrow(Y), ncol(Y), length(1:max(K)))
  Y_indicator <- 1 * array(outer(Y_for_indicator, 1:max_K, "=="), dim = dim_Y_indicator)
  for (j in 1:M) {
    Y_indicator[Y_missing_indicator[, j] == 1, j, ] <- 0
  }
  S <- array(0, dim = c(G, M, max(K)))
  for (k in 1:max(K)) {
    S[,,k] <- t(z)%*%Y_indicator[,,k]
  }
  Y_current <- LCR_initialize_missing_values(Y = Y, Y_missing_indicator = Y_missing_indicator, K = K, n = n, M = M)

  N_jk <- array(0, dim = c(M, max_K))
  for (j in 1:M) {
    for (k in 1:K[j]) {
      N_jk[j, k] <- sum(Y_indicator[, j, k])  
    }
  }
  

  n_obs_j <- colSums(1 - Y_missing_indicator)  
  
  sum1_inclusion <- (G-1) * (lgamma(K*theta_hyperparam) - K*lgamma(theta_hyperparam))
  
  log_gamma_N_jk_alpha <- lgamma(N_jk + theta_hyperparam)
  
  sum2_inclusion <- lgamma(n_obs_j + K*theta_hyperparam) - rowSums(log_gamma_N_jk_alpha)
  
  sum3_inclusion <- log(clust_var_prior) - log(1-clust_var_prior)
  
  inclusion_sum <- sum1_inclusion + sum2_inclusion + sum3_inclusion
  exclusion_sum <- -inclusion_sum
  
  N_g <- colSums(z)
  N_gjk <- array(0, dim = c(G, M, max_K))
  

  for (j in 1:M) {
    for (k in 1:max_K) {

      S_jk <- Y_indicator[, j, k]
      

      observed_S_jk <- S_jk * (1 - Y_missing_indicator[, j])
      

      N_gjk[, j, k] <- as.vector(t(z) %*% observed_S_jk)
    }
  }
  
  A <- kappa[,1:(G-1)]+omega*C[,1:(G-1)]
  init <- list(Y = Y ,X = X1, X_current = X1, p = p, 
               M = M, n = n, K = K, 
               theta = theta, beta = beta, mu = mu, 
               logit_probs = logit_probs,
               log_logit_probs = log_logit_probs, 
               C = C, omega = omega, eta = eta, beta_cov = beta_cov,
               beta_cov_inv = beta_cov_inv, 
               beta_cov_inv_chol = beta_cov_inv_chol,
               beta_mean = beta_mean,
               w = w,
               z  = z,
               kappa = kappa,
               gamma = gamma,
               tau = tau,
               nu = nu,
               Y_indicator = Y_indicator,
               N_jk = N_jk,
               inclusion_sum = inclusion_sum,
               exclusion_sum = exclusion_sum,
               N_g = N_g,
               N_gjk = N_gjk,
               A = A,
               theta_hyperparam = theta_hyperparam,
               clust_var_prior = clust_var_prior,
               beta_prior_cov_inv = beta_prior_cov_inv,
               S = S,
               Y_missing_indicator = Y_missing_indicator,
               Y_current = Y_current,
               dim_Y_indicator = dim_Y_indicator,
               max_K = max_K
               )
  return(init)
}


LCR_init_sample_arrays <- function(beta, n_samples, omega,
                                   w, z, theta, gamma, N_g,
                                   N_gjk, nu, Y){
  beta_samples <- array(0, dim = c(dim(beta), n_samples))
  omega_samples <- array(0, dim = c(dim(omega), n_samples))
  w_samples <- array(0, dim = c(dim(w), n_samples))
  z_samples <- array(0, dim = c(dim(z), n_samples))
  theta_samples <- array(0, dim = c(dim(theta), n_samples))
  gamma_samples <- matrix(0, nrow = length(gamma), ncol = n_samples)
  N_g_samples <- array(0, dim = c(length(N_g), n_samples))
  N_gjk_samples <- array(0, dim = c(dim(N_gjk), n_samples))
  nu_samples <- matrix(0, nrow = length(nu), ncol = n_samples)
  Y_imputed_samples <- array(0, dim = c(nrow(Y),ncol(Y),n_samples))
  log_post_samples <- numeric(n_samples)
  log_like_samples <- numeric(n_samples)
  return(list(beta_samples = beta_samples, omega_samples = omega_samples,
         w_samples = w_samples, z_samples = z_samples, 
         theta_samples = theta_samples, gamma_samples = gamma_samples,
         N_g_samples = N_g_samples, N_gjk_samples = N_gjk_samples,
         nu_samples = nu_samples, log_post_samples = log_post_samples,
         log_like_samples = log_like_samples, Y_imputed_samples = Y_imputed_samples))
}



LCR_initialize_missing_values <- function(Y, Y_missing_indicator, K, n, M) {
  n <- nrow(Y)
  M <- ncol(Y)
  Y_initialized <- Y
  
  for (j in 1:M) {
    observed_vals <- Y[Y_missing_indicator[, j] == 0, j]
    
    if (length(observed_vals) > 0) {
      obs_table <- table(factor(observed_vals, levels = 1:K[j]))
      probs <- as.numeric(obs_table / sum(obs_table))
      missing_rows <- which(Y_missing_indicator[, j] == 1)
      if (length(missing_rows) > 0) {
        Y_initialized[missing_rows, j] <- sample(1:K[j], 
                                                 size = length(missing_rows), 
                                                 replace = TRUE, 
                                                 prob = probs)
      }
    } else {
      missing_rows <- which(Y_missing_indicator[, j] == 1)
      Y_initialized[missing_rows, j] <- sample(1:K[j], 
                                               size = length(missing_rows), 
                                               replace = TRUE)
    }
  }
  
  return(Y_initialized)
}




LCR_select_functions <- function(item.sel, cov.sel){
  
}

