# This script runs both simulation 1 and 2 for 300 different replicates, N = 500, 300, 150 

library(here)
setwd(here())
source('./sim_study_monte_carlo_functions.R')
source('./LCR_sim_data.R')
source('./LCR_Gibbs.R')

library(parallelly)
library(e1071)      
library(clue)        
library(HDInterval)  
library(future.apply)
library(progressr)
library(MCMCpack)

cpp_files <- normalizePath(c("./z_update_collapsed.cpp",
                             "./z_update_uncollapsed.cpp"))
cpp_cache <- file.path(getwd(), "rcpp_cache")
dir.create(cpp_cache, showWarnings = FALSE)

load_cpp <- function(files, cache) {
  for (f in files) {
    Rcpp::sourceCpp(f, cacheDir = cache, env = globalenv())
  }
  invisible(TRUE)
}

load_cpp(cpp_files, cpp_cache)


ensure_cpp_loaded <- function(files, cache) {
  if (!isTRUE(get0(".cpp_loaded", envir = globalenv(), inherits = FALSE))) {
    load_cpp(files, cache)
    assign(".cpp_loaded", TRUE, envir = globalenv())
  }
  invisible(TRUE)
}

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


true_item_active <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)  
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







#Main loop - fitting model, extracting quantities we want, discard other less important stuff

#Simulation 1


plan(multisession, workers = 64)

handlers(global = TRUE)
handlers("progress")

for (nm in names(sample_sizes)) {
  
  cat("=== Sample size:", nm, "===\n")
  
  with_progress({
    
    pr <- progressor(n_replicates)
    
    results <- future_lapply(
      1:n_replicates,
      function(t) {
        
        ensure_cpp_loaded(cpp_files, cpp_cache)
        # Rcpp::sourceCpp("./z_update_collapsed.cpp")
        # Rcpp::sourceCpp("./z_update_uncollapsed.cpp")
        
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
          K = K,
          p = p,
          M = M
        )
        
        pr(sprintf("Replicate %d finished", t))
        
        summ
      },
      future.seed = TRUE
    )
  })
  
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







#Running the same thing now for simulation 2



G <- 3
p <- 6
M <- 13
K <- c(2,2,2,3,3,3,4,4,3,3,5,5,5)

# covariate effects matrix
sim2_beta <- matrix(c(0, 0, 1, -1, -1, 1, 0.5, -0.5, -0.4, 0.4, 0, 0, 0, 0), nrow = p+1, byrow = TRUE)



sim2_theta1 <- matrix(c(0.15, 0.6, 0.8, 0.85, 0.4, 0.2), nrow = G, ncol = K[1])
sim2_theta2 <- matrix(c(0.25, 0.45, 0.7, 0.75, 0.55, 0.3), nrow = G, ncol = K[2])
sim2_theta3 <- matrix(c(0.7, 0.2, 0.65, 0.3, 0.8, 0.35), nrow = G, ncol = K[3])
sim2_theta4 <- matrix(c(0.1, 0.35, 0.7, 0.25, 0.4, 0.2, 0.65, 0.25, 0.1), nrow = G, ncol = K[4])
sim2_theta5 <- matrix(c(0.2, 0.25, 0.65, 0.15, 0.6, 0.25, 0.65, 0.15, 0.1), nrow = G, ncol = K[5])
sim2_theta6 <- matrix(c(0.15, 0.5, 0.75, 0.2, 0.35, 0.15, 0.65, 0.15, 0.1), nrow = G, ncol = K[6])
sim2_theta7 <- matrix(c(0.1, 0.25, 0.6, 0.15, 0.35, 0.25, 0.25, 0.25, 0.1, 0.5, 0.15, 0.05), nrow = G, ncol = K[7])
sim2_theta8 <- matrix(c(0.15, 0.2, 0.55, 0.2, 0.45, 0.2, 0.2, 0.25, 0.15, 0.45, 0.1, 0.1), nrow = G, ncol = K[8])
sim2_theta9 <- matrix(c(0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1), nrow = G, ncol = K[9])
sim2_theta10 <- matrix(c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2), nrow = G, ncol = K[10])
sim2_theta11 <- matrix(1/K[11], nrow = G, ncol = K[11])
sim2_theta12 <- matrix(c(0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.2, 0.2, 0.2, 0.25, 0.25, 0.25, 0.3, 0.3, 0.3), nrow = G, ncol = K[12])
sim2_theta13 <- matrix(c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = G, ncol = K[13])


sim2_theta <- list(sim2_theta1, sim2_theta2, sim2_theta3, sim2_theta4, sim2_theta5, sim2_theta6, sim2_theta7, sim2_theta8, sim2_theta9, sim2_theta10, sim2_theta11, sim2_theta12, sim2_theta13)



true_item_active <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)  
true_pred_active <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)               

n_replicates <- 300

#Simulating datasets

sim_study2_dataset_list_N500 <- vector("list", n_replicates)
sim_study2_dataset_list_N300 <- vector("list", n_replicates)
sim_study2_dataset_list_N150 <- vector("list", n_replicates)

set.seed(123)
for (t in 1:n_replicates) {
  sim_study2_dataset_list_N500[[t]] <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 500)
  sim_study2_dataset_list_N300[[t]] <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 300)
  sim_study2_dataset_list_N150[[t]] <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 150)
}


sample_sizes <- c("N500" = 500, "N300" = 300, "N150" = 150)


sim2_results <- lapply(names(sample_sizes), function(nm) {
  list(
    beta = make_beta_arrays(p, G, n_replicates),
    theta = make_theta_arrays(K, G, n_replicates),
    item_sel = make_selection_arrays(M, n_replicates),
    pred_sel = make_selection_arrays(p, n_replicates)
  )
})
names(sim2_results) <- names(sample_sizes)

dataset_lists <- list(
  N500 = sim_study2_dataset_list_N500,
  N300 = sim_study2_dataset_list_N300,
  N150 = sim_study2_dataset_list_N150
)



plan(multisession, workers = 32)

handlers(global = TRUE)
handlers("progress")

for (nm in names(sample_sizes)) {
  
  cat("=== Sample size:", nm, "===\n")
  
  with_progress({
    
    pr <- progressor(n_replicates)
    
    results <- future_lapply(
      1:n_replicates,
      function(t) {
        
        ensure_cpp_loaded(cpp_files, cpp_cache)
        # Rcpp::sourceCpp("./z_update_collapsed.cpp")
        # Rcpp::sourceCpp("./z_update_uncollapsed.cpp")
        
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
          true_beta = cbind(sim2_beta, 0),
          true_theta = sim2_theta,
          G = G,
          K = K,
          p = p,
          M = M
        )
        
        pr(sprintf("Replicate %d finished", t))
        
        summ
      },
      future.seed = TRUE
    )
  })
  for (t in 1:n_replicates) {
    summ <- results[[t]]
    
    sim2_results[[nm]]$beta$coverage[, , t] <- summ$beta_cov
    sim2_results[[nm]]$beta$bias[, , t] <- summ$beta_bias
    sim2_results[[nm]]$beta$mse[, , t] <- summ$beta_mse
    sim2_results[[nm]]$beta$coverage_cond[, , t] <- summ$beta_cov_cond
    sim2_results[[nm]]$beta$bias_cond[, , t] <- summ$beta_bias_cond
    sim2_results[[nm]]$beta$mse_cond[, , t] <- summ$beta_mse_cond
    sim2_results[[nm]]$beta$n_included[, t] <- summ$beta_n_included
    
    for (m in 1:M) {
      sim2_results[[nm]]$theta$coverage[[m]][, , t] <- summ$theta_cov[[m]]
      
      sim2_results[[nm]]$theta$mean[[m]][, , t] <- summ$theta_mean[[m]]
      sim2_results[[nm]]$theta$sd[[m]][, , t]  <- summ$theta_sd[[m]]
      
      sim2_results[[nm]]$theta$bias[[m]][, , t] <- summ$theta_bias[[m]]
      sim2_results[[nm]]$theta$mse[[m]][, , t] <- summ$theta_mse[[m]]
      sim2_results[[nm]]$theta$coverage_cond[[m]][, , t] <- summ$theta_cov_cond[[m]]
      sim2_results[[nm]]$theta$bias_cond[[m]][, , t] <- summ$theta_bias_cond[[m]]
      sim2_results[[nm]]$theta$mse_cond[[m]][, , t] <- summ$theta_mse_cond[[m]]
    }
    sim2_results[[nm]]$theta$n_included[, t] <- summ$theta_n_included
    
    sim2_results[[nm]]$item_sel$pip[, t] <- summ$item_pip
    sim2_results[[nm]]$item_sel$indicator[, t] <- summ$item_pip > 0.5
    sim2_results[[nm]]$pred_sel$pip[, t] <- summ$pred_pip
    sim2_results[[nm]]$pred_sel$indicator[, t] <- summ$pred_pip > 0.5
  }
  
  
  
  
  
  saveRDS(
    sim2_results[[nm]],
    paste0("sim2_results_", nm, ".rds")
  )
}

plan(sequential)













#Same thing for simulation 3 







G <- 4
n <- 500

#We consider 5 predictors, one representing age (N(0,1)), and 4 representing clinical dianoses 
# (binary, proportions 0.2, 0,15, 0.1, 0.5 - replicating ASD, other, ID diagnoses, gender respectively)
p <- 5
sim3_X <- matrix(NA, ncol = p, nrow = n)

# Continuous, representing age
sim3_X[,1] <- rnorm(n)
# Binary (0.5 proportion) representing gender
sim3_X[,2] <- sample(c(0,1), n, replace = TRUE, prob = c(0.5,0.5))
# Binary (0.2 proportion) representing ASD
sim3_X[,3] <- sample(c(0,1), n, replace = TRUE, prob = c(0.8,0.2))
# Binary (0.15 proportion) representing other diagnoses
sim3_X[,4] <- sample(c(0,1), n, replace = TRUE, prob = c(0.85,0.15))
# Binary (0.06 proportion) representing ID diagnosis
sim3_X[,5] <- sample(c(0,1), n, replace = TRUE, prob = c(0.94,0.06))

M <- 40
K <- rep(3, M)

# covariate effects matrix
sim3_beta <- matrix(c(
  -1.6, -1.2, -1.6,   # intercept
  1.0, -0.8,  0.6,   # standard normal
  1.2,  0.0, -1.0,   # binary 
  -1.5,  1.2,  0.8,   # binary 
  1.0, -1.2,  1.5,   # binary 
  0.0,  0.0,  0.0    # binary - non-informative
), nrow = p + 1, byrow = TRUE)

#We let items 1-20 be informative with distinct parameters, and 21-40 non-informative

# sim3_theta1 <- matrix(c(
#   0.4, 0.3, 0.3,
#   0.2, 0.3, 0.5,
#   0.1, 0.2, 0.7,
#   0.6, 0.2, 0.2
# ), 
# nrow = G, ncol = 3, byrow = TRUE)
# sim3_theta2 <- matrix(c(
#   0.2, 0.4, 0.4,
#   0.5, 0.3, 0.2,
#   0.1, 0.3, 0.6,
#   0.7, 0.2, 0.1
# ), 
# nrow = G, ncol = 3, byrow = TRUE)
# sim3_theta3 <- matrix(c(
#   0.5, 0.2, 0.3,
#   0.1, 0.4, 0.5,
#   0.2, 0.2, 0.6,
#   0.4, 0.3, 0.3
# ), 
# nrow = G, ncol = 3, byrow = TRUE)
# sim3_theta4 <- matrix(c(
#   0.15, 0.3, 0.55,
#   0.6, 0.25, 0.15,
#   0.1, 0.3, 0.6,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta5 <- matrix(c(
#   0.2, 0.3, 0.5,
#   0.55, 0.3, 0.15,
#   0.15, 0.25, 0.6,
#   0.6, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta6 <- matrix(c(
#   0.15, 0.35, 0.5,
#   0.65, 0.2, 0.15,
#   0.1, 0.3, 0.6,
#   0.45, 0.35, 0.2
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta7 <- matrix(c(
#   0.2, 0.25, 0.55,
#   0.6, 0.25, 0.15,
#   0.4, 0.35, 0.25,
#   0.7, 0.2, 0.1
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta8 <- matrix(c(
#   0.15, 0.30, 0.55,
#   0.55, 0.25, 0.20,
#   0.15, 0.30, 0.55,
#   0.60, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta9 <- matrix(c(
#   0.2, 0.3, 0.5,
#   0.65, 0.2, 0.15,
#   0.1, 0.25, 0.65,
#   0.55, 0.3, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta10 <- matrix(c(
#   0.15, 0.35, 0.5,
#   0.6, 0.2, 0.2,
#   0.15, 0.3, 0.55,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta11 <- matrix(c(
#   0.2, 0.3, 0.5,
#   0.55, 0.3, 0.15,
#   0.45, 0.3, 0.25,
#   0.6, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta12 <- matrix(c(
#   0.6, 0.25, 0.15,
#   0.15, 0.3, 0.55,
#   0.1, 0.3, 0.6,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta13 <- matrix(c(
#   0.55, 0.3, 0.15,
#   0.2, 0.3, 0.5,
#   0.15, 0.25, 0.6,
#   0.45, 0.35, 0.2
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta14 <- matrix(c(
#   0.65, 0.2, 0.15,
#   0.15, 0.35, 0.5,
#   0.1, 0.3, 0.6,
#   0.6, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta15 <- matrix(c(
#   0.6, 0.25, 0.15,
#   0.2, 0.3, 0.5,
#   0.4, 0.35, 0.25,
#   0.7, 0.2, 0.1
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta16 <- matrix(c(
#   0.55, 0.25, 0.2,
#   0.15, 0.3, 0.55,
#   0.15, 0.3, 0.55,
#   0.55, 0.3, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta17 <- matrix(c(
#   0.6, 0.2, 0.2,
#   0.2, 0.3, 0.5,
#   0.1, 0.25, 0.65,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta18 <- matrix(c(
#   0.6, 0.2, 0.2,
#   0.15, 0.35, 0.5,
#   0.15, 0.3, 0.55,
#   0.6, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta19 <- matrix(c(
#   0.55, 0.3, 0.15,
#   0.2, 0.3, 0.5,
#   0.45, 0.3, 0.25,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta20 <- matrix(c(
#   0.5, 0.3, 0.2,
#   0.5, 0.3, 0.2,
#   0.15, 0.30, 0.55,
#   0.45, 0.35, 0.2
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# 

sim3_theta1 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta2 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.05, 0.25, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta3 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.20, 0.10, 0.70,
  0.55, 0.30, 0.15     # G4 slightly elevated
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta4 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,    # G3 only moderate
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta5 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.10, 0.20, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta6 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.10, 0.15, 0.75,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

# Items 7-12: G1 low, G2 high, G3 high, G4 low
sim3_theta7 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta8 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta9 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta10 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.15, 0.35, 0.50,
  0.55, 0.30, 0.15     # G4 slightly elevated
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta11 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.05, 0.05, 0.90,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta12 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.15, 0.15, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

# Items 13-16: G1 moderate, G2 low, G3 high, G4 low
sim3_theta13 <- matrix(c(
  0.35, 0.40, 0.25,
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta14 <- matrix(c(
  0.35, 0.40, 0.25,
  0.80, 0.15, 0.05,
  0.10, 0.10, 0.80,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta15 <- matrix(c(
  0.35, 0.40, 0.25,
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,    # G3 only moderate
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta16 <- matrix(c(
  0.35, 0.40, 0.25,
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

# Items 17-20: G1 low, G2 moderate, G3 high, G4 low
sim3_theta17 <- matrix(c(
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,
  0.10, 0.20, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta18 <- matrix(c(
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,
  0.10, 0.30, 0.60,
  0.55, 0.30, 0.15     # G4 slightly elevated
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta19 <- matrix(c(
  0.55, 0.30, 0.15,    # G1 slightly elevated
  0.35, 0.40, 0.25,
  0.10, 0.10, 0.80,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta20 <- matrix(c(
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,
  0.10, 0.20, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

# ---------------- Non-informative items 21-40 ----------------
sim3_theta21 <- matrix(rep(c(0.5, 0.4, 0.1), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta22 <- matrix(rep(c(0.5, 0.4, 0.1), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta23 <- matrix(rep(c(0.3, 0.1, 0.6), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta24 <- matrix(rep(c(0.5, 0.2, 0.3), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta25 <- matrix(rep(c(0.2, 0.05, 0.75), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta26 <- matrix(rep(c(0.7, 0.05, 0.25), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta27 <- matrix(rep(c(0.25, 0.5, 0.25), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta28 <- matrix(rep(c(0.2, 0.3, 0.5), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta29 <- matrix(rep(c(0.6, 0.2, 0.2), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta30 <- matrix(rep(c(0.7, 0.1, 0.2), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta31 <- matrix(rep(c(0.65, 0.15, 0.2), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta32 <- matrix(rep(c(0.15, 0.2, 0.65), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta33 <- matrix(rep(c(0.45, 0.25, 0.3), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta34 <- matrix(rep(c(0.05, 0.35, 0.6), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta35 <- matrix(rep(c(0.85, 0.1, 0.05), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta36 <- matrix(rep(c(0.1, 0.35, 0.55), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta37 <- matrix(rep(c(0.3, 0.67, 0.03), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta38 <- matrix(rep(c(0.8, 0.1, 0.1), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta39 <- matrix(rep(c(1/3, 1/3, 1/3), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta40 <- matrix(rep(c(0.1, 0.5, 0.4), G), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta <- list(
  sim3_theta1,  sim3_theta2,  sim3_theta3,  sim3_theta4,  sim3_theta5,
  sim3_theta6,  sim3_theta7,  sim3_theta8,  sim3_theta9,  sim3_theta10,
  sim3_theta11, sim3_theta12, sim3_theta13, sim3_theta14, sim3_theta15,
  sim3_theta16, sim3_theta17, sim3_theta18, sim3_theta19, sim3_theta20,
  sim3_theta21, sim3_theta22, sim3_theta23, sim3_theta24, sim3_theta25,
  sim3_theta26, sim3_theta27, sim3_theta28, sim3_theta29, sim3_theta30,
  sim3_theta31, sim3_theta32, sim3_theta33, sim3_theta34, sim3_theta35,
  sim3_theta36, sim3_theta37, sim3_theta38, sim3_theta39, sim3_theta40
)



true_item_active <- c(rep(TRUE, 20), rep(FALSE, 20))  
true_pred_active <- rep(TRUE, 5)               

#We'll try with 100 replicates for the moment and will see about scaling up later on
n_replicates <- 100

#Simulating datasets

sim_study3_dataset_list_N500 <- vector("list", n_replicates)
sim_study3_dataset_list_N300 <- vector("list", n_replicates)
sim_study3_dataset_list_N150 <- vector("list", n_replicates)

set.seed(123)
for (t in 1:n_replicates) {
  sim_study3_dataset_list_N500[[t]] <- LCR_sim_data(theta = sim3_theta, beta = sim3_beta, n_samples = 500)
  sim_study3_dataset_list_N300[[t]] <- LCR_sim_data(theta = sim3_theta, beta = sim3_beta, n_samples = 300)
  sim_study3_dataset_list_N150[[t]] <- LCR_sim_data(theta = sim3_theta, beta = sim3_beta, n_samples = 150)
}


sample_sizes <- c("N500" = 500, "N300" = 300, "N150" = 150)

# Creating a list for all results from simulation 1
sim3_results <- lapply(names(sample_sizes), function(nm) {
  list(
    beta = make_beta_arrays(p, G, n_replicates),
    theta = make_theta_arrays(K, G, n_replicates),
    item_sel = make_selection_arrays(M, n_replicates),
    pred_sel = make_selection_arrays(p, n_replicates)
  )
})
names(sim3_results) <- names(sample_sizes)

dataset_lists <- list(
  N500 = sim_study3_dataset_list_N500,
  N300 = sim_study3_dataset_list_N300,
  N150 = sim_study3_dataset_list_N150
)



plan(multisession, workers = 32)

handlers(global = TRUE)
handlers("progress")

for (nm in names(sample_sizes)) {
  
  cat("=== Sample size:", nm, "===\n")
  
  with_progress({
    
    pr <- progressor(n_replicates)
    
    results <- future_lapply(
      1:n_replicates,
      function(t) {
        
        # Rcpp::sourceCpp("./z_update_collapsed.cpp")
        # Rcpp::sourceCpp("./z_update_uncollapsed.cpp")
        ensure_cpp_loaded(cpp_files, cpp_cache)
        
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
          true_beta = cbind(sim3_beta, 0),
          true_theta = sim3_theta,
          G = G,
          K = K,
          p = p,
          M = M
        )
        
        pr(sprintf("Replicate %d finished", t))
        
        summ
      },
      future.seed = TRUE
    )
  })

  for (t in 1:n_replicates) {
    summ <- results[[t]]
    
    sim3_results[[nm]]$beta$coverage[, , t] <- summ$beta_cov
    sim3_results[[nm]]$beta$bias[, , t] <- summ$beta_bias
    sim3_results[[nm]]$beta$mse[, , t] <- summ$beta_mse
    sim3_results[[nm]]$beta$coverage_cond[, , t] <- summ$beta_cov_cond
    sim3_results[[nm]]$beta$bias_cond[, , t] <- summ$beta_bias_cond
    sim3_results[[nm]]$beta$mse_cond[, , t] <- summ$beta_mse_cond
    sim3_results[[nm]]$beta$n_included[, t] <- summ$beta_n_included
    
    for (m in 1:M) {
      sim3_results[[nm]]$theta$coverage[[m]][, , t] <- summ$theta_cov[[m]]
      
      sim3_results[[nm]]$theta$mean[[m]][, , t] <- summ$theta_mean[[m]]
      sim3_results[[nm]]$theta$sd[[m]][, , t]  <- summ$theta_sd[[m]]
      
      sim3_results[[nm]]$theta$bias[[m]][, , t] <- summ$theta_bias[[m]]
      sim3_results[[nm]]$theta$mse[[m]][, , t] <- summ$theta_mse[[m]]
      sim3_results[[nm]]$theta$coverage_cond[[m]][, , t] <- summ$theta_cov_cond[[m]]
      sim3_results[[nm]]$theta$bias_cond[[m]][, , t] <- summ$theta_bias_cond[[m]]
      sim3_results[[nm]]$theta$mse_cond[[m]][, , t] <- summ$theta_mse_cond[[m]]
    }
    sim3_results[[nm]]$theta$n_included[, t] <- summ$theta_n_included
    
    sim3_results[[nm]]$item_sel$pip[, t] <- summ$item_pip
    sim3_results[[nm]]$item_sel$indicator[, t] <- summ$item_pip > 0.5
    sim3_results[[nm]]$pred_sel$pip[, t] <- summ$pred_pip
    sim3_results[[nm]]$pred_sel$indicator[, t] <- summ$pred_pip > 0.5
  }
  
  
  
  
  
  saveRDS(
    sim3_results[[nm]],
    paste0("sim3_results_", nm, ".rds")
  )
}

plan(sequential)



















