# This is a script for runnning some summary output from the Monte Carlo Simulation for LCR

# We'll make a function for each of the things we want to report, which are

# 1. beta : coverage, bias, MSE

# 2. theta: coverage, bias, MSE

# 3. predictor variable selection: aggregated metrics, average true positive, false positive (anything else?)

# 3. item variable selection: aggregated metrics, average true positive, false positive (anything else?)



beta_summarise <- function(fit){
  
  beta_list <- fit$beta
  
  p <- dim(beta_list$coverage)[1] - 1
  
  G <- dim(beta_list$coverage)[2] + 1
  
  n_replicates <- dim(beta_list$coverage)[3]
  
  #Computing the mean of beta (ignoring NAs)
  mean_coverage_cond <- matrix(NA, nrow = p+1, ncol = G-1)
  mean_bias_cond <- matrix(NA, nrow = p+1, ncol = G-1)
  mean_mse_cond <- matrix(NA, nrow = p+1, ncol = G-1)
  
  for (l in 1:(p+1)){
    for (g in 1:(G-1)){
      mean_coverage_cond[l,g] <- mean(beta_list$coverage_cond[l,g,], na.rm = TRUE)
      mean_bias_cond[l,g] <- mean(beta_list$bias_cond[l,g,], na.rm = TRUE)
      mean_mse_cond[l,g] <- mean(beta_list$mse_cond[l,g,], na.rm = TRUE)
    }
  }
  
  beta_cond <- list(
    coverage = beta_list$coverage_cond,
    bias = beta_list$bias_cond,
    mse = beta_list$mse_cond,
    mean_coverage = mean_coverage_cond,
    mean_bias = mean_bias_cond,
    mean_mse = mean_mse_cond,
    n_included = beta_list$n_included
  )
  beta_uncond <- list(
    coverage = beta_list$coverage,
    bias = beta_list$bias,
    mse = beta_list$mse,
    mean_coverage = apply(beta_list$coverage, c(1,2), mean),
    mean_bias = apply(beta_list$bias, c(1,2), mean),
    mean_mse = apply(beta_list$mse, c(1,2), mean)
  )
  return(list(
    conditional = beta_cond,
    unconditional = beta_uncond
  ))
}


# Function for computing the quantile based credible interval for theta and then checking whether it contains the true value 
quantile_95_cred_interval_coverage_theta <- function(theta_means, theta_sd, true_theta){
  
  M <- length(theta_means)
  
  n_replicate <- dim(theta_means[[1]])[3]
  
  G <- dim(theta_means[[1]])[1]
  
  K <- unlist(lapply(theta_means, function(x) dim(x)[2]))
  
  coverage <- list()
  
  for (j in 1:M){
    coverage[[j]] <- array(NA, dim = c(G, K[j], n_replicate))
    for (t in 1:n_replicate){
      for (g in 1:G){
        for (k in 1:K[j]){
          mu <- theta_means[[j]][g,k,t]
          v <- (theta_sd[[j]][g,k,t])^2
          S <- (mu*(1-mu))/(v) - 1
          q_025 <- qbeta(0.025, S*mu, S - S*mu)
          q_975 <- qbeta(0.975, S*mu, S - S*mu)
          
          coverage[[j]][g, k, t] <- (true_theta[[j]][g,k] >= q_025) & (true_theta[[j]][g,k] <= q_975)
        }
      }
    }
  }
  return(coverage)
}






# I think we may need to modify this in order to split into conditional and unconditional on inclusion 
# (if this is possible - may not be due to only having a point estimate and variance for theta)
theta_summarise <- function(fit, true_theta){
  theta_list <- fit$theta
  
  theta_means <- theta_list$mean
  theta_sd <- theta_list$sd
  coverage <- quantile_95_cred_interval_coverage_theta(theta_means = theta_means, theta_sd = theta_sd, true_theta = true_theta)
  
  
  M <- length(coverage)
  
  mean_coverage <- list()
  mean_bias <- list()
  mean_mse <- list()
  for (j in 1:M){
    mean_coverage[[j]] <- apply(coverage[[j]], c(1,2), mean)
    mean_bias[[j]] <- apply(theta_list$bias[[j]], c(1,2), mean)
    mean_mse[[j]] <- apply(theta_list$mse[[j]], c(1,2), mean)
  }
  
  
  result <- list(
    coverage = coverage,
    bias = theta_list$bias,
    mse = theta_list$mse,
    mean_coverage = mean_coverage,
    mean_bias = mean_bias,
    mean_mse = mean_mse,
    n_included = theta_list$n_included
  )
  return(result)
}

pred_var_sel_summarise <- function(fit, true_pred_active, threshold = 0.5){
  
  pip <- fit$pred_sel$pip
  
  p <- nrow(pip)
  
  n_replicates <- ncol(pip)
  
  ind <- pip>threshold
  
  TP <- colSums(true_pred_active & ind)
  FP <- colSums(!true_pred_active & ind)
  FN <- colSums(true_pred_active & !ind)
  TN <- colSums(!true_pred_active & !ind)
  
  TPR <- TP/(TP + FN)
  FPR <- FP/(FP + TN)
  
  per_replicate <- data.frame(
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    TPR = TPR,
    FPR = FPR
  )
  
  mean_pip <- rowMeans(pip)
  mean_ind <- rowMeans(ind)
  
  per_pred <- data.frame(
    variable = seq_len(p),
    true_active = true_pred_active,
    mean_pip = mean_pip,
    mean_ind = mean_ind
  )
  
  result <- list(
    per_replicate = per_replicate,
    per_predictor = per_pred
  ) 
  
  return(result)
}

item_var_sel_summarise <- function(fit, true_item_active, threshold = 0.5){
  
  pip <- fit$item_sel$pip
  
  M <- nrow(pip)
  
  n_replicates <- ncol(pip)
  
  ind <- pip>threshold
  
  TP <- colSums(true_item_active & ind)
  FP <- colSums(!true_item_active & ind)
  FN <- colSums(true_item_active & !ind)
  TN <- colSums(!true_item_active & !ind)
  
  TPR <- TP/(TP + FN)
  FPR <- FP/(FP + TN)
  
  per_replicate <- data.frame(
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    TPR = TPR,
    FPR = FPR
  )
  
  mean_pip <- rowMeans(pip)
  mean_ind <- rowMeans(ind)
  
  per_item <- data.frame(
    variable = seq_len(M),
    true_active = true_item_active,
    mean_pip = mean_pip,
    mean_ind = mean_ind
  )
  
  result <- list(
    per_replicate = per_replicate,
    per_item = per_item
  ) 
  
  return(result)
}




