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


sim_study1_dataset_list_N500 <- list()
sim_study2_dataset_list_N500 <- list()
sim_study3_dataset_list_N500 <- list()

sim_study1_dataset_list_N300 <- list()
sim_study2_dataset_list_N300 <- list()
sim_study3_dataset_list_N300 <- list()

sim_study1_dataset_list_N150 <- list()
sim_study2_dataset_list_N150 <- list()
sim_study3_dataset_list_N150 <- list()

# Defining model parameters


G <- 2
p <- 6
M <- 8
K <- rep(3,8)

# covariate effects vector
sim1_beta <- matrix(c(-0.5, 0.8, 1.2, 1, 0.4, 0, 0), nrow = p+1)

#item probability parameters
sim1_theta1 <- matrix(c(0.15,0.25,0.6,0.7,0.2,0.1), nrow = G, ncol = K[1], byrow = TRUE)
sim1_theta2 <- matrix(c(0.2,0.35,0.45,0.55,0.3,0.15), nrow = G, ncol = K[2], byrow = TRUE)
sim1_theta3 <- matrix(c(0.1,0.15,0.75,0.8,0.15,0.05), nrow = G, ncol = K[3], byrow = TRUE)
sim1_theta4 <- matrix(c(0.25,0.4,0.35,0.45,0.35,0.2), nrow = G, ncol = K[4], byrow = TRUE)
sim1_theta5 <- matrix(c(0.4, 0.5, 0.1, 0.4, 0.5, 0.1), nrow = G, ncol = K[5], byrow = TRUE)
sim1_theta6 <- matrix(c(0.7, 0.1, 0.2, 0.7, 0.1, 0.2), nrow = G, ncol = K[6], byrow = TRUE)
sim1_theta7 <- matrix(c(0.1, 0.5, 0.4, 0.1, 0.5, 0.4), nrow = G, ncol = K[7], byrow = TRUE)
sim1_theta8 <- matrix(1/K[8], nrow = G, ncol = K[8])

sim1_theta <- list(sim1_theta1, sim1_theta2, sim1_theta3, sim1_theta4, sim1_theta5, sim1_theta6, sim1_theta7, sim1_theta8)


# Simulating replicates for simulation 1
n_replicates <- 300

set.seed(123)
for (t in 1:n_replicates){
  sim_study1_dataset_list_N500[[t]] <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 500)
  sim_study1_dataset_list_N300[[t]] <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 300)
  sim_study1_dataset_list_N150[[t]] <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 150)
}


# Creating storage arrays

# beta coverage arrays
sim1_beta_coverage_array_N500 <- array(0, dim = c(P+1, G-1, n_replicates))
sim1_beta_coverage_array_N300 <- array(0, dim = c(P+1, G-1, n_replicates))
sim1_beta_coverage_array_N150 <- array(0, dim = c(P+1, G-1, n_replicates))

# theta coverage arrays
sim1_theta_coverage_array_N500 <- array()
sim1_theta_coverage_array_N300 <- array()
sim1_theta_coverage_array_N150 <- array()
  
# beta bias arrays
sim1_beta_bias_array_N500 <- array(0, dim = c(P+1, G-1, n_replicates))
sim1_beta_bias_array_N300 <- array(0, dim = c(P+1, G-1, n_replicates))
sim1_beta_bias_array_N150 <- array(0, dim = c(P+1, G-1, n_replicates))

# theta bias arrays
sim1_theta_bias_array_N500 <- array()
sim1_theta_bias_array_N300 <- array()
sim1_theta_bias_array_N150 <- array()
  
# beta MSE arrays
sim1_beta_mse_array_N500 <- array(0, dim = c(P+1, G-1, n_replicates))
sim1_beta_mse_array_N300 <- array(0, dim = c(P+1, G-1, n_replicates)) 
sim1_beta_mse_array_N150 <- array(0, dim = c(P+1, G-1, n_replicates)) 

# theta MSE arrays
sim1_theta_mse_array_N500 <- array()  
sim1_theta_mse_array_N300 <- array()
sim1_theta_mse_array_N150 <- array()

# item selection indicator arrays
sim1_item_sel_ind_array_N500 <- array(0, dim = )
sim1_item_sel_ind_array_N300 <- array(0, dim = )
sim1_item_sel_ind_array_N150 <- array(0, dim = )

# item selection probabiliity arrays
sim1_item_sel_ind_array_N500 <- array(0, dim = )
sim1_item_sel_ind_array_N300 <- array(0, dim = )
sim1_item_sel_ind_array_N150 <- array(0, dim = )
  
# predictor selection indicator arrays
sim1_pred_sel_array_N500 <- array(0, dim = )
sim1_pred_sel_array_N300 <- array(0, dim = )
sim1_pred_sel_array_N150 <- array(0, dim = )

# predictor selection probability arrays
sim1_pred_sel_array_N500 <- array(0, dim = )
sim1_pred_sel_array_N300 <- array(0, dim = )
sim1_pred_sel_array_N150 <- array(0, dim = )


# Running the sampler
# We want to double check convergernce metrics etc. before running all of this...
for (t in 1:n_replicates){
  print(paste('Replicate: ', t))
  fit <- LCR_Gibbs(X = sim_study1_dataset_list_N500[[t]]$X, Y = sim_study1_dataset_list_N500[[t]]$Y, 
                   G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                   theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                   relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)
  
  
}










