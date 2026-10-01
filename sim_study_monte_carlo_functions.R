# This is a script defining the functions for running the Monte Carlo Simulation for LCR

# Functions for creating arrays to store quantities of interest


make_beta_arrays <- function(p, G, n_replicates) {
  list(
    coverage = array(NA, dim = c(p + 1, G - 1, n_replicates)),
    bias = array(NA, dim = c(p + 1, G - 1, n_replicates)),
    mse = array(NA, dim = c(p + 1, G - 1, n_replicates)),
    coverage_cond = array(NA, dim = c(p + 1, G - 1, n_replicates)),
    bias_cond = array(NA, dim = c(p + 1, G - 1, n_replicates)),
    mse_cond = array(NA, dim = c(p + 1, G - 1, n_replicates)),
    n_included = array(NA, dim = c(p + 1, n_replicates))  
  )
}

make_theta_arrays <- function(K, G, n_replicates) {
  list(
    coverage = lapply(K, function(Kj) array(NA, dim = c(G, Kj, n_replicates))),
    bias = lapply(K, function(Kj) array(NA, dim = c(G, Kj, n_replicates))),
    mse = lapply(K, function(Kj) array(NA, dim = c(G, Kj, n_replicates))),
    coverage_cond = lapply(K, function(Kj) array(NA, dim = c(G, Kj, n_replicates))),
    bias_cond = lapply(K, function(Kj) array(NA, dim = c(G, Kj, n_replicates))),
    mse_cond = lapply(K, function(Kj) array(NA, dim = c(G, Kj, n_replicates))),
    mean = lapply(K, function(Kj) array(NA, dim = c(G, Kj, n_replicates))),
    var = lapply(K, function(Kj) array(NA, dim = c(G, Kj, n_replicates))),
    sd = lapply(K, function(Kj) array(NA, dim = c(G, Kj, n_replicates))),
    n_included = array(NA, dim = c(length(K), n_replicates))  
  )
}

make_selection_arrays <- function(n_vars, n_replicates) {
  list(
    pip = array(NA, dim = c(n_vars, n_replicates)),
    indicator = array(NA, dim = c(n_vars, n_replicates))  
  )
}

#Function for computing the coverage, bias, MSE
summarise_draws <- function(draws, truth, min_draws = 10) {
  n <- length(draws)
  if (n < min_draws) {
    return(list(coverage = NA, bias = NA, mse = NA, n = n))
  }
  ci <- hdi(draws, credMass = 0.95)
  list(
    coverage = (truth >= ci["lower"]) & (truth <= ci["upper"]),
    bias = mean(draws) - truth,
    mse  = (mean(draws) - truth)^2,
    n = n
  )
}




# function for getting correspondence between true labels and the estimated partition 
# so that we have beta in terms of the correct baseline

get_relabel_perm <- function(true_labels, est_labels, G) {
  conf_mat <- table(factor(est_labels, levels = 1:G),
                    factor(true_labels, levels = 1:G))
  perm <- solve_LSAP(conf_mat, maximum = TRUE)
  as.integer(perm)  
}

# function to apply a class permutation to theta matrix

relabel_theta_mat <- function(theta_mat, perm) {
  out <- theta_mat
  out[perm, ] <- theta_mat      
  out
}

# function to apply a class permutation to beta
relabel_beta_mat <- function(beta_mat_G, perm) {
  beta_mat_G[, perm, drop = FALSE]
}


#function for taking the desired quantities from the MCMC output
extract_replicate_summary <- function(fit, true_z, true_beta, true_theta, G, K, p, M) {
  
  perm <- get_relabel_perm(true_labels = true_z, est_labels = max.col(fit$assignment_prob), G = G)
  
  n_iter <- dim(fit$samples$beta_samples)[3]
  
  beta_draws_relab <- array(NA, dim = dim(fit$samples$beta_samples))
  for (i in 1:n_iter) {
    beta_draws_relab[, , i] <- relabel_beta_mat(fit$samples$beta_samples[, , i], perm)
  }
  
  n_iter_theta <- dim(fit$samples$theta_samples)[4]
  
  theta_draws_relab <- vector("list", M)
  for (j in 1:M) {
    Kj <- K[j]
    arr_j <- fit$samples$theta_samples[, j, 1:Kj, , drop = TRUE]
    if (Kj == 1) dim(arr_j) <- c(G, 1, n_iter_theta)
    arr_j_relab <- array(NA, dim = dim(arr_j))
    for (i in 1:n_iter_theta) {
      arr_j_relab[, , i] <- relabel_theta_mat(arr_j[, , i], perm)
    }
    theta_draws_relab[[j]] <- arr_j_relab
  }
  
  gamma_samples <- fit$samples$gamma_samples
  
  
  nu_samples <- fit$samples$nu_samples
  
  beta_cov  <- array(NA, dim = c(p + 1, G - 1))
  beta_bias <- array(NA, dim = c(p + 1, G - 1))
  beta_mse  <- array(NA, dim = c(p + 1, G - 1))
  
  beta_cov_cond  <- array(NA, dim = c(p + 1, G - 1))
  beta_bias_cond <- array(NA, dim = c(p + 1, G - 1))
  beta_mse_cond  <- array(NA, dim = c(p + 1, G - 1))
  
  beta_n_included <- numeric(p + 1)
  
  for (j in 1:(p + 1)) {
    
    mask <- gamma_samples[j, ] == 1
    beta_n_included[j] <- sum(mask)
    
    for (g in 1:(G - 1)) {
      
      draws <- beta_draws_relab[j, g, ] - beta_draws_relab[j, G, ]
      truth <- true_beta[j, g] - true_beta[j, G]
      
      uncond <- summarise_draws(draws, truth)
      beta_cov[j, g]  <- uncond$coverage
      beta_bias[j, g] <- uncond$bias
      beta_mse[j, g]  <- uncond$mse
      
      cond <- summarise_draws(draws[mask], truth)
      beta_cov_cond[j, g]  <- cond$coverage
      beta_bias_cond[j, g] <- cond$bias
      beta_mse_cond[j, g]  <- cond$mse
    }
  }
  
  
  theta_cov <- vector("list", M) 
  theta_bias <- vector("list", M)
  theta_mse <- vector("list", M)
  theta_cov_cond <- vector("list", M)
  theta_bias_cond <- vector("list", M)
  theta_mse_cond <- vector("list", M)
  theta_n_included <- numeric(M)
  theta_mean <- lapply(fit$itemprob, relabel_theta_mat, perm = perm)
  theta_sd <- lapply(fit$itemprob.sd, relabel_theta_mat, perm = perm)
  
  for (j in 1:M) {
    Kj <- K[j]
    
    mask <- nu_samples[j, ] == 1
    theta_n_included[j] <- sum(mask)
    
    cov_j  <- matrix(NA, G, Kj)
    bias_j  <- matrix(NA, G, Kj)
    mse_j  <- matrix(NA, G, Kj)
    cov_j_cond <- matrix(NA, G, Kj)
    bias_j_cond <- matrix(NA, G, Kj)
    mse_j_cond <- matrix(NA, G, Kj)
    
    for (g in 1:G) {
      for (k in 1:Kj) {
        
        draws <- theta_draws_relab[[j]][g, k, ]
        truth <- true_theta[[j]][g, k]
        
        uncond <- summarise_draws(draws, truth)
        cov_j[g, k] <- uncond$coverage
        bias_j[g, k] <- uncond$bias
        mse_j[g, k] <- uncond$mse
        
        cond <- summarise_draws(draws[mask], truth)
        cov_j_cond[g, k] <- cond$coverage
        bias_j_cond[g, k] <- cond$bias
        mse_j_cond[g, k] <- cond$mse
      }
    }
    
    theta_cov[[j]] <- cov_j 
    theta_bias[[j]] <- bias_j
    theta_mse[[j]] <- mse_j
    theta_cov_cond[[j]] <- cov_j_cond
    theta_bias_cond[[j]] <- bias_j_cond
    theta_mse_cond[[j]] <- mse_j_cond
  }
  
  item_pip <- fit$item_inclusion_prob
  pred_pip <- fit$cov_inclusion_prob[-1]
  
  list(
    beta_cov = beta_cov, beta_bias = beta_bias, beta_mse = beta_mse,
    beta_cov_cond = beta_cov_cond, beta_bias_cond = beta_bias_cond, beta_mse_cond = beta_mse_cond,
    beta_n_included = beta_n_included,
    
    theta_cov = theta_cov, theta_bias = theta_bias, theta_mse = theta_mse,
    theta_cov_cond = theta_cov_cond, theta_bias_cond = theta_bias_cond, theta_mse_cond = theta_mse_cond,
    theta_n_included = theta_n_included, theta_mean = theta_mean, theta_sd = theta_sd,
    
    item_pip = item_pip, pred_pip = pred_pip
  )
}


# # Summary objects
# # beta: average coverage/bias/mse across active vs inactive coefficients
# summarise_beta <- function(res, group) {
#   
#   cov_mean_intercept <- mean(res$beta$coverage[group == "intercept", , ], na.rm = TRUE)
#   cov_mean_active <- mean(res$beta$coverage[group == "active", , ], na.rm = TRUE)
#   cov_mean_inactive <- mean(res$beta$coverage[group == "inactive", , ], na.rm = TRUE)
#   
#   bias_mean_intercept <- mean(res$beta$bias[group == "intercept", , ], na.rm = TRUE)
#   bias_mean_active <- mean(res$beta$bias[group == "active", , ], na.rm = TRUE)
#   bias_mean_inactive <- mean(res$beta$bias[group == "inactive", , ], na.rm = TRUE)
#   
#   mse_mean_intercept <- mean(res$beta$mse[group == "intercept", , ], na.rm = TRUE)
#   mse_mean_active <- mean(res$beta$mse[group == "active", , ], na.rm = TRUE)
#   mse_mean_inactive <- mean(res$beta$mse[group == "inactive", , ], na.rm = TRUE)
#   
#   data.frame(
#     group = c("intercept", "active", "inactive"),
#     coverage = c(cov_mean_intercept, cov_mean_active, cov_mean_inactive),
#     bias = c(bias_mean_intercept, bias_mean_active, bias_mean_inactive),
#     mse = c(mse_mean_intercept, mse_mean_active, mse_mean_inactive)
#   )
# }
# 
# # item/predictor selection: TPR / FPR across replications 
# summarise_selection <- function(sel_arrays, true_active) {
#   ind <- sel_arrays$indicator  
#   tpr <- mean(ind[true_active, , drop = FALSE] == 1)
#   fpr <- mean(ind[!true_active, , drop = FALSE] == 1)
#   data.frame(TPR = tpr, FPR = fpr)
# }












