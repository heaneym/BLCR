# This script runs both simulation 1 and 2 for 300 different replicates, N = 500, 300, 150 
# (TBC whether this is enough) 


# We don't want to store the MCMC objects due to their size, so we should compute what we need and store those
# before moving to the following iteration


# What we need to compute (recommended by the reviewer):
# 1. The proportion of MC replications in which the 95% HDI contains the true parameter value
# 2. Bias and MSE of posterior point estimates
# 3. Variable Selection Accuracy: aggregated metrics such as 
# the average true positive rate and false positive rate for 
# both item and predictor selection computed across all MC replications

library(e1071)      
library(clue)        
library(HDInterval)  
library(future.apply)
library(progressr)

#Simulation study 1

#Defining model parameters

G <- 2
p <- 6          
M <- 8
K <- rep(3, 8)  

sim1_beta <- matrix(c(-0.5, 0.8, 1.2, 1, 0.4, 0, 0), nrow = p + 1)

sim1_theta1 <- matrix(c(0.15,0.25,0.6,0.7,0.2,0.1), nrow = G, ncol = K[1], byrow = TRUE)
sim1_theta2 <- matrix(c(0.2,0.35,0.45,0.55,0.3,0.15), nrow = G, ncol = K[2], byrow = TRUE)
sim1_theta3 <- matrix(c(0.1,0.15,0.75,0.8,0.15,0.05), nrow = G, ncol = K[3], byrow = TRUE)
sim1_theta4 <- matrix(c(0.25,0.4,0.35,0.45,0.35,0.2), nrow = G, ncol = K[4], byrow = TRUE)
sim1_theta5 <- matrix(c(0.4, 0.5, 0.1, 0.4, 0.5, 0.1), nrow = G, ncol = K[5], byrow = TRUE)
sim1_theta6 <- matrix(c(0.7, 0.1, 0.2, 0.7, 0.1, 0.2), nrow = G, ncol = K[6], byrow = TRUE)
sim1_theta7 <- matrix(c(0.1, 0.5, 0.4, 0.1, 0.5, 0.4), nrow = G, ncol = K[7], byrow = TRUE)
sim1_theta8 <- matrix(1/K[8], nrow = G, ncol = K[8])

sim1_theta <- list(sim1_theta1, sim1_theta2, sim1_theta3, sim1_theta4,
                   sim1_theta5, sim1_theta6, sim1_theta7, sim1_theta8)


true_item_active <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)  
true_pred_active <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)               

n_replicates <- 300

#Simulating datasets

sim_study1_dataset_list_N500 <- vector("list", n_replicates)
sim_study1_dataset_list_N300 <- vector("list", n_replicates)
sim_study1_dataset_list_N150 <- vector("list", n_replicates)

set.seed(123)
for (t in 1:n_replicates) {
  sim_study1_dataset_list_N500[[t]] <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 500)
  sim_study1_dataset_list_N300[[t]] <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 300)
  sim_study1_dataset_list_N150[[t]] <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 150)
}


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
    n_included = array(NA, dim = c(length(K), n_replicates))  
  )
}

make_selection_arrays <- function(n_vars, n_replicates) {
  list(
    pip = array(NA, dim = c(n_vars, n_replicates)),
    indicator = array(NA, dim = c(n_vars, n_replicates))  
  )
}


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


sample_sizes <- c("N500" = 500, "N300" = 300, "N150" = 150)

# Creating a list for all results from simulation 1
sim1_results <- lapply(names(sample_sizes), function(nm) {
  list(
    beta = make_beta_arrays(p, G, n_replicates),
    theta = make_theta_arrays(K, G, n_replicates),
    item_sel = make_selection_arrays(M, n_replicates),
    pred_sel = make_selection_arrays(p, n_replicates)
  )
})
names(sim1_results) <- names(sample_sizes)

dataset_lists <- list(
  N500 = sim_study1_dataset_list_N500,
  N300 = sim_study1_dataset_list_N300,
  N150 = sim_study1_dataset_list_N150
)

# function for getting correspondence between true labels and the estimated partition 
# so that we have beta in terms of the correct baseline

get_relabel_perm <- function(true_labels, est_labels, G) {
  conf_mat <- table(factor(est_labels, levels = 1:G),
                    factor(true_labels, levels = 1:G))
  perm <- solve_LSAP(conf_mat, maximum = TRUE)
  as.integer(perm)  
}

# function to apply a class permutation to theta

relabel_theta_mat <- function(theta_mat, perm) {
  theta_mat[perm, , drop = FALSE]
}

# function to apply a class permutation to beta
relabel_beta_mat <- function(beta_mat_G, perm) {
  beta_mat_G[, perm, drop = FALSE]
}



extract_replicate_summary <- function(fit, true_z, true_beta, true_theta, G, K, p) {
  
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

  
  theta_cov  <- vector("list", M); theta_bias  <- vector("list", M); theta_mse  <- vector("list", M)
  theta_cov_cond <- vector("list", M); theta_bias_cond <- vector("list", M); theta_mse_cond <- vector("list", M)
  theta_n_included <- numeric(M)
  theta_mean <- fit$itemprob
  theta_sd <- fit$itemprob.sd
  
  for (j in 1:M) {
    Kj <- K[j]
    
    mask <- nu_samples[j, ] == 1
    theta_n_included[j] <- sum(mask)
    
    cov_j  <- matrix(NA, G, Kj); bias_j  <- matrix(NA, G, Kj); mse_j  <- matrix(NA, G, Kj)
    cov_j_cond <- matrix(NA, G, Kj); bias_j_cond <- matrix(NA, G, Kj); mse_j_cond <- matrix(NA, G, Kj)
    
    for (g in 1:G) {
      for (k in 1:Kj) {
        
        draws <- theta_draws_relab[[j]][g, k, ]
        truth <- true_theta[[j]][g, k]
        
        uncond <- summarise_draws(draws, truth)
        cov_j[g, k]  <- uncond$coverage
        bias_j[g, k] <- uncond$bias
        mse_j[g, k]  <- uncond$mse
        
        cond <- summarise_draws(draws[mask], truth)
        cov_j_cond[g, k]  <- cond$coverage
        bias_j_cond[g, k] <- cond$bias
        mse_j_cond[g, k]  <- cond$mse
      }
    }
    
    theta_cov[[j]] <- cov_j;  theta_bias[[j]]  <- bias_j;  theta_mse[[j]]  <- mse_j
    theta_cov_cond[[j]] <- cov_j_cond; theta_bias_cond[[j]] <- bias_j_cond; theta_mse_cond[[j]] <- mse_j_cond
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


#Main loop - fiting model, extracting quantities we want, discard other less important stuff


# for (nm in names(sample_sizes)) {
#   
#   cat("=== Sample size:", nm, "===\n")
#   
#   for (t in 1:n_replicates) {
#     cat("Replicate:", t, "\n")
#     
#     dat <- dataset_lists[[nm]][[t]]
#     
#     fit <- LCR_Gibbs(X = dat$X, Y = dat$Y,
#                      G = G, beta_prior_cov = diag(10^2, p + 1), beta_prior_mean = rep(0, p + 1),
#                      theta_hyperparam = 1, clust_var_prior = 0.5,
#                      item.sel = TRUE, cov.sel = TRUE, verbose = TRUE,
#                      relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)
#     
# 
#     summ <- extract_replicate_summary(
#       fit = fit, true_z = dat$class,
#       true_beta = cbind(sim1_beta, 0),   
#       true_theta = sim1_theta, G = G, K = K
#     )
#     
#     sim1_results[[nm]]$beta$coverage[, , t] <- summ$beta_cov
#     sim1_results[[nm]]$beta$bias[, , t] <- summ$beta_bias
#     sim1_results[[nm]]$beta$mse[, , t] <- summ$beta_mse
# 
#     for (m in 1:M) {
#       sim1_results[[nm]]$theta$coverage[[m]][, , t] <- summ$theta_cov[[m]]
#       sim1_results[[nm]]$theta$bias[[m]][, , t] <- summ$theta_bias[[m]]
#       sim1_results[[nm]]$theta$mse[[m]][, , t] <- summ$theta_mse[[m]]
#     }
# 
#     sim1_results[[nm]]$item_sel$pip[, t] <- summ$item_pip
#     sim1_results[[nm]]$item_sel$indicator[, t] <- summ$item_pip > 0.5
#     
#     sim1_results[[nm]]$pred_sel$pip[, t] <- summ$pred_pip
#     sim1_results[[nm]]$pred_sel$indicator[, t] <- summ$pred_pip > 0.5
#     
#     
#     print(sim1_results[[nm]]$beta$coverage[, , t])
#     print(sim1_results[[nm]]$beta$bias[, , t])
#     print(sim1_results[[nm]]$beta$mse[, , t])
#     
#     rm(fit, dat, summ)
#     gc()
#   }
# 
#   saveRDS(sim1_results[[nm]], paste0("sim1_results_", nm, ".rds"))
# }


plan(multisession, workers = 8)

handlers(global = TRUE)
handlers("progress")

for (nm in names(sample_sizes)) {
  
  cat("=== Sample size:", nm, "===\n")
  
  with_progress({
    
    pr <- progressor(n_replicates)
    
    results <- future_lapply(
      1:n_replicates,
      function(t) {
        
        Rcpp::sourceCpp("./z_update_collapsed.cpp")
        Rcpp::sourceCpp("./z_update_uncollapsed.cpp")
        
        dat <- dataset_lists[[nm]][[t]]
        
        fit <- LCR_Gibbs(
          X = dat$X,
          Y = dat$Y,
          G = G,
          beta_prior_cov = diag(10^2, p + 1),
          beta_prior_mean = rep(0, p + 1),
          theta_hyperparam = 1,
          clust_var_prior = 0.5,
          item.sel = TRUE,
          cov.sel = TRUE,
          verbose = FALSE,
          relabel = TRUE,
          n_samples = 5000,
          burnin = 1000,
          thinby = 10
        )
        
        summ <- extract_replicate_summary(
          fit = fit,
          true_z = dat$class,
          true_beta = cbind(sim1_beta, 0),
          true_theta = sim1_theta,
          G = G,
          K = K
        )
        
        pr(sprintf("Replicate %d finished", t))
        
        summ
      },
      future.seed = TRUE
    )
  })
  
  # Put the results from each replicate into the existing arrays
  
  # for (t in 1:n_replicates) {
  #   
  #   summ <- results[[t]]
  #   
  #   sim1_results[[nm]]$beta$coverage[, , t] <- summ$beta_cov
  #   sim1_results[[nm]]$beta$bias[, , t] <- summ$beta_bias
  #   sim1_results[[nm]]$beta$mse[, , t] <- summ$beta_mse
  #   
  #   for (m in 1:M) {
  #     sim1_results[[nm]]$theta$coverage[[m]][, , t] <-
  #       summ$theta_cov[[m]]
  #     
  #     sim1_results[[nm]]$theta$bias[[m]][, , t] <-
  #       summ$theta_bias[[m]]
  #     
  #     sim1_results[[nm]]$theta$mse[[m]][, , t] <-
  #       summ$theta_mse[[m]]
  #   }
  #   
  #   sim1_results[[nm]]$item_sel$pip[, t] <- summ$item_pip
  #   sim1_results[[nm]]$item_sel$indicator[, t] <-
  #     summ$item_pip > 0.5
  #   
  #   sim1_results[[nm]]$pred_sel$pip[, t] <- summ$pred_pip
  #   sim1_results[[nm]]$pred_sel$indicator[, t] <-
  #     summ$pred_pip > 0.5
  # }
  
  #Updated code computes coverage, bias, MSE for the iterations for all iterations 
  #as well as those where a variable is included
  for (t in 1:n_replicates) {
    summ <- results[[t]]

    sim1_results[[nm]]$beta$coverage[, , t] <- summ$beta_cov
    sim1_results[[nm]]$beta$bias[, , t] <- summ$beta_bias
    sim1_results[[nm]]$beta$mse[, , t] <- summ$beta_mse
    sim1_results[[nm]]$beta$coverage_cond[, , t] <- summ$beta_cov_cond
    sim1_results[[nm]]$beta$bias_cond[, , t] <- summ$beta_bias_cond
    sim1_results[[nm]]$beta$mse_cond[, , t] <- summ$beta_mse_cond
    sim1_results[[nm]]$beta$n_included[, t] <- summ$beta_n_included

    for (m in 1:M) {
      sim1_results[[nm]]$theta$coverage[[m]][, , t] <- summ$theta_cov[[m]]
      
      sim1_results[[nm]]$theta$mean[[m]][, , t] <- summ$theta_mean[[m]]
      sim1_results[[nm]]$theta$sd[[m]][, , t]  <- summ$theta_sd[[m]]
      
      sim1_results[[nm]]$theta$bias[[m]][, , t] <- summ$theta_bias[[m]]
      sim1_results[[nm]]$theta$mse[[m]][, , t] <- summ$theta_mse[[m]]
      sim1_results[[nm]]$theta$coverage_cond[[m]][, , t] <- summ$theta_cov_cond[[m]]
      sim1_results[[nm]]$theta$bias_cond[[m]][, , t] <- summ$theta_bias_cond[[m]]
      sim1_results[[nm]]$theta$mse_cond[[m]][, , t] <- summ$theta_mse_cond[[m]]
    }
    sim1_results[[nm]]$theta$n_included[, t] <- summ$theta_n_included

    sim1_results[[nm]]$item_sel$pip[, t] <- summ$item_pip
    sim1_results[[nm]]$item_sel$indicator[, t] <- summ$item_pip > 0.5
    sim1_results[[nm]]$pred_sel$pip[, t] <- summ$pred_pip
    sim1_results[[nm]]$pred_sel$indicator[, t] <- summ$pred_pip > 0.5
  }
  
  
  
  
  
  saveRDS(
    sim1_results[[nm]],
    paste0("sim1_results_", nm, ".rds")
  )
}

plan(sequential)













# Summary objects
# beta: average coverage/bias/mse across active vs inactive coefficients
summarise_beta <- function(res, group) {
  
  cov_mean_intercept <- mean(res$beta$coverage[group == "intercept", , ], na.rm = TRUE)
  cov_mean_active <- mean(res$beta$coverage[group == "active", , ], na.rm = TRUE)
  cov_mean_inactive <- mean(res$beta$coverage[group == "inactive", , ], na.rm = TRUE)
  
  bias_mean_intercept <- mean(res$beta$bias[group == "intercept", , ], na.rm = TRUE)
  bias_mean_active <- mean(res$beta$bias[group == "active", , ], na.rm = TRUE)
  bias_mean_inactive <- mean(res$beta$bias[group == "inactive", , ], na.rm = TRUE)
  
  mse_mean_intercept <- mean(res$beta$mse[group == "intercept", , ], na.rm = TRUE)
  mse_mean_active <- mean(res$beta$mse[group == "active", , ], na.rm = TRUE)
  mse_mean_inactive <- mean(res$beta$mse[group == "inactive", , ], na.rm = TRUE)
  
  data.frame(
    group = c("intercept", "active", "inactive"),
    coverage = c(cov_mean_intercept, cov_mean_active, cov_mean_inactive),
    bias = c(bias_mean_intercept, bias_mean_active, bias_mean_inactive),
    mse = c(mse_mean_intercept, mse_mean_active, mse_mean_inactive)
  )
}

# item/predictor selection: TPR / FPR across replications 
summarise_selection <- function(sel_arrays, true_active) {
  ind <- sel_arrays$indicator  
  tpr <- mean(ind[true_active, , drop = FALSE] == 1)
  fpr <- mean(ind[!true_active, , drop = FALSE] == 1)
  data.frame(TPR = tpr, FPR = fpr)
}

# Example usage once all N's are done:
# active_idx_beta <- c(TRUE, true_pred_active)  # TRUE for intercept, then predictors
# beta_summary_N500 <- summarise_beta(sim1_results$N500, active_idx_beta)
# pred_sel_summary_N500 <- summarise_selection(sim1_results$N500$pred_sel, true_pred_active)
# item_sel_summary_N500 <- summarise_selection(sim1_results$N500$item_sel, true_item_active)