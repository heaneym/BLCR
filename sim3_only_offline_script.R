# This script runs both simulation 1 and 2 for 300 different replicates, N = 500, 300, 150 

library(here)
setwd(here())
source('./sim_study_monte_carlo_functions.R')
source('./LCR_sim_data.R')
source('./LCR_Gibbs.R')

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





#Same thing for simulation 3 

G <- 4
n <- 500
M <- 40
K <- rep(3,M)

p <- 5

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
true_pred_active <- c(rep(TRUE, 4), FALSE)               

#We'll try with 100 replicates for the moment and will see about scaling up later on
n_replicates <- 300

#Simulating datasets

sim_study3_dataset_list_N500 <- vector("list", n_replicates)
sim_study3_dataset_list_N300 <- vector("list", n_replicates)
sim_study3_dataset_list_N150 <- vector("list", n_replicates)

sim_study3_X_list_N500 <- vector("list", n_replicates)
sim_study3_X_list_N300 <- vector("list", n_replicates)
sim_study3_X_list_N150 <- vector("list", n_replicates)

p <- 5
sim3_X <- matrix(NA, ncol = p, nrow = n)



set.seed(123)
for (t in 1:n_replicates) {
  sim_study3_X_list_N500[[t]] <- matrix(NA, ncol = p, nrow = 500)
  sim_study3_X_list_N500[[t]][,1] <- rnorm(500)
  sim_study3_X_list_N500[[t]][,2] <- sample(c(0,1), 500, replace = TRUE, prob = c(0.5,0.5))
  sim_study3_X_list_N500[[t]][,3] <- sample(c(0,1), 500, replace = TRUE, prob = c(0.8,0.2))
  sim_study3_X_list_N500[[t]][,4] <- sample(c(0,1), 500, replace = TRUE, prob = c(0.85,0.15))
  sim_study3_X_list_N500[[t]][,5] <- sample(c(0,1), 500, replace = TRUE, prob = c(0.94,0.06))
  
  sim_study3_X_list_N300[[t]] <- matrix(NA, ncol = p, nrow = 300)
  sim_study3_X_list_N300[[t]][,1] <- rnorm(300)
  sim_study3_X_list_N300[[t]][,2] <- sample(c(0,1), 300, replace = TRUE, prob = c(0.5,0.5))
  sim_study3_X_list_N300[[t]][,3] <- sample(c(0,1), 300, replace = TRUE, prob = c(0.8,0.2))
  sim_study3_X_list_N300[[t]][,4] <- sample(c(0,1), 300, replace = TRUE, prob = c(0.85,0.15))
  sim_study3_X_list_N300[[t]][,5] <- sample(c(0,1), 300, replace = TRUE, prob = c(0.94,0.06))
  
  sim_study3_X_list_N150[[t]] <- matrix(NA, ncol = p, nrow = 150)
  sim_study3_X_list_N150[[t]][,1] <- rnorm(150)
  sim_study3_X_list_N150[[t]][,2] <- sample(c(0,1), 150, replace = TRUE, prob = c(0.5,0.5))
  sim_study3_X_list_N150[[t]][,3] <- sample(c(0,1), 150, replace = TRUE, prob = c(0.8,0.2))
  sim_study3_X_list_N150[[t]][,4] <- sample(c(0,1), 150, replace = TRUE, prob = c(0.85,0.15))
  sim_study3_X_list_N150[[t]][,5] <- sample(c(0,1), 150, replace = TRUE, prob = c(0.94,0.06))
  
  sim_study3_dataset_list_N500[[t]] <- LCR_sim_data_given_X(theta = sim3_theta, beta = sim3_beta, n_samples = 500, X = sim_study3_X_list_N500[[t]])
  sim_study3_dataset_list_N300[[t]] <- LCR_sim_data_given_X(theta = sim3_theta, beta = sim3_beta, n_samples = 300, X = sim_study3_X_list_N300[[t]])
  sim_study3_dataset_list_N150[[t]] <- LCR_sim_data_given_X(theta = sim3_theta, beta = sim3_beta, n_samples = 150, X = sim_study3_X_list_N150[[t]])
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
handlers("txtprogressbar")

for (nm in names(sample_sizes)) {
  
  cat("=== Sample size:", nm, "===\n")
  
  progress_file <- file.path(getwd(), paste0("sim3_progress_", nm, ".log"))
  writeLines(character(0), progress_file)   # reset/create the file
  start_time <- Sys.time()
  
  results <- future_lapply(
    1:n_replicates,
    function(t) {
      
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
        thinby = 10,
        K_specify = K
      )
      
      summ <- extract_replicate_summary(
        fit = fit,
        true_z = dat$class,
        true_beta = cbind(sim3_beta, 0),
        true_theta = sim3_theta,
        G = G, K = K, p = p, M = M
      )
      
      # one line per finished replicate; append = TRUE is safe enough for this use
      cat(sprintf("%s | %s | replicate %d finished\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), nm, t),
          file = progress_file, append = TRUE)
      
      summ
    },
    future.seed = TRUE,
    future.scheduling = Inf     # one task per replicate, so workers pick up jobs dynamically
  )
  
  cat(sprintf("%s | %s done in %.1f min\n", format(Sys.time()), nm,
              as.numeric(difftime(Sys.time(), start_time, units = "mins"))))

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



















