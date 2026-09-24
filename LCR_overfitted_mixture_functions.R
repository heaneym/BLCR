
beta_update_sfm <- function(beta, X_current, G, p, beta_prior_cov_inv,
                            A, beta_prior_mean, omega){
  X_cov <- X_current[, -1, drop = FALSE]              
  q <- ncol(X_cov)                                
  cov_inv <- beta_prior_cov_inv[-1, -1, drop = FALSE]   
  m0 <- beta_prior_mean[-1]
  intercept <- beta[1, ]                                  
  
  beta_mean <- matrix(0, nrow = q, ncol = (G - 1))
  beta_cov_inv <- array(0, dim = c(q, q, G - 1))
  beta_cov_inv_chol <- array(0, dim = c(q, q, G - 1))
  for (g in 1:(G - 1)){
    A_g <- A[, g] - intercept[g] * omega[, g]               
    prec <- crossprod(X_cov, omega[, g] * X_cov) + cov_inv
    L <- chol(prec)
    fs <- forwardsolve(t(L), crossprod(X_cov, A_g) + cov_inv %*% m0)
    mg <- backsolve(L, fs)
    beta_mean[, g] <- mg
    beta_cov_inv[, , g] <- prec
    beta_cov_inv_chol[, , g] <- L
    beta[2:(p + 1), g] <- rMVNormP(1, mg, L)
  }
  beta[2:(p + 1), G] <- 0                             
  list(beta = beta, beta_cov_inv = beta_cov_inv,
       beta_mean = beta_mean, beta_cov_inv_chol = beta_cov_inv_chol)
}

#Function for sampling the auxiliary weights (for sparsity) - uses Metropolis with reparametrised weights and a random walk - updates the parameter in the form of the offset to the intercept in beta

sparse_weights_update <- function(beta, mu, z, G, sparse_prior, sparse_sigma_MH, has_covariates = TRUE, n){
  N_g <- colSums(z)
  if (!has_covariates){
    sparse_weights_current <- as.numeric(rdirichlet(1, sparse_prior + N_g))
    accepted <- TRUE
  } else{
    #Removing the intercept contribution from mu for use in the acceptance criterion
    a <- exp(mu - matrix(beta[1, ], n, G, byrow = TRUE)) 
    #Forming the log-target (we'll take the difference of this function evaluated at current and proposed values)
    log_target <- function(weights) sum((sparse_prior + N_g) * log(weights)) - sum(log(as.vector(a %*% weights)))
    #Extracting the current values for the weights - including a shift by the max column for stability
    intercept  <- beta[1, ] - max(beta[1, ]) 
    sparse_weights_current <- exp(intercept)/sum(exp(intercept))
    #proposing a new value for the parameter v (v_g = log(eta_g/eta_G), however this denominator of eta_G is a common shift across groups so doesnt impact the expression)
    v_prop <- log(sparse_weights_current) + rnorm(G, 0, sparse_sigma_MH)
    #Shifting for stability
    v_prop <- v_prop - max(v_prop)
    #transforming back
    sparse_weights_prop <- exp(v_prop) / sum(exp(v_prop))
    log_alpha <- log_target(sparse_weights_prop) - log_target(sparse_weights_current)
    accepted <- (log(runif(1)) < log_alpha)
    accept_prob <- exp(min(0, log_alpha))
    if (accepted) sparse_weights_current <- sparse_weights_prop
  }
  #We want it in terms of a zero-baseline, so we transform back to this
  beta[1, ] <- log(sparse_weights_current) - log(sparse_weights_current[G]) 
  list(beta = beta, sparse_weights = sparse_weights_current, accepted = accepted, accept_prob = accept_prob)
}


summarise_G_eff <- function(N_g_samples, n, min_frac = 0.01){
  thr  <- max(1, ceiling(min_frac * n))
  occ  <- apply(N_g_samples, 2, function(ng) sum(ng > 0))      
  mass <- apply(N_g_samples, 2, function(ng) sum(ng >= thr))   
  pmf  <- function(x){ t <- table(factor(x, levels = sort(unique(x)))); t / sum(t) }
  list(threshold = thr, occupancy_posterior = pmf(occ),
       cluster_posterior = pmf(mass),
       cluster_mode = as.integer(names(which.max(table(mass)))),
       cluster_mean = mean(mass))
}


extract_gating <- function(beta_samples, G){
  intercept_samples <- matrix(beta_samples[1, , ], nrow = G)        
  slope_samples <- beta_samples[-1, , , drop = FALSE]           
  weight_samples <- apply(intercept_samples, 2, function(b){     
    b <- b - max(b)
    e <- exp(b)
    e / sum(e)
  })
  weight_ordered_mean <- rowMeans(apply(weight_samples, 2, sort, decreasing = TRUE))
  list(intercept_samples = intercept_samples,   
       slope_samples = slope_samples,      
       weight_samples = weight_samples,      
       weight_ordered_mean = weight_ordered_mean) 
}


