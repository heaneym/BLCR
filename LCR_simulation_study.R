source('LCR_Gibbs.R')
source('LCR_sim_data.R')

#Simulation Study 1

# 2 groups
# 8 item variables, 4 of which are useful for clustering (all item variables have 3 possible outcomes)
# 6 covariates, of varying effect size
# 300 observations

# structure parameters
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
sim1_theta5 <- matrix(1/K[5], nrow = G, ncol = K[5])
sim1_theta6 <- matrix(1/K[6], nrow = G, ncol = K[6])
sim1_theta7 <- matrix(1/K[7], nrow = G, ncol = K[7])
sim1_theta8 <- matrix(1/K[8], nrow = G, ncol = K[8])

sim1_theta <- list(sim1_theta1, sim1_theta2, sim1_theta3, sim1_theta4, sim1_theta5, sim1_theta6, sim1_theta7, sim1_theta8)

sim1_data <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 300, n_noise_var = 0)

sim1_LCR_fit <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                          theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                          relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


#Simulation Study 2

# 3 groups
# 13 item variables, 8 of which are useful for clustering 
# variables 1,2,3 have 2 possible responses, 4,5,6 have 3 possible responses, and 7 and 8 have 4 possible responses
# variables 9,10 have 3 possible responses, 11,12,13 have 5 possible responses
# 6 covariates, of varying effect size for different groups
# 300 observations

# structure parameters
G <- 3
p <- 6
M <- 13
K <- c(2,2,2,3,3,3,4,4,3,3,5,5,5)

# covariate effects matrix
sim2_beta <- matrix(c(-1, 0.8, 1.1, 1.2, 0.3, 0, 0, -0.3, 0.5, 0, 0.4, 0, 0, 0), nrow = p+1)


sim2_theta1 <- matrix(c(0.15, 0.6, 0.8, 0.85, 0.4, 0.2), nrow = G, ncol = K[1])
sim2_theta2 <- matrix(c(0.25, 0.45, 0.7, 0.75, 0.55, 0.3), nrow = G, ncol = K[2])
sim2_theta3 <- matrix(c(0.7, 0.2, 0.65, 0.3, 0.8, 0.35), nrow = G, ncol = K[3])
sim2_theta4 <- matrix(c(0.1, 0.35, 0.7, 0.25, 0.4, 0.2, 0.65, 0.25, 0.1), nrow = G, ncol = K[4])
sim2_theta5 <- matrix(c(0.2, 0.25, 0.65, 0.15, 0.6, 0.25, 0.65, 0.15, 0.1), nrow = G, ncol = K[5])
sim2_theta6 <- matrix(c(0.15, 0.5, 0.75, 0.2, 0.35, 0.15, 0.65, 0.15, 0.1), nrow = G, ncol = K[6])
sim2_theta7 <- matrix(c(0.1, 0.25, 0.6, 0.15, 0.35, 0.25, 0.25, 0.25, 0.1, 0.5, 0.15, 0.05), nrow = G, ncol = K[7])
sim2_theta8 <- matrix(c(0.15, 0.2, 0.55, 0.2, 0.45, 0.2, 0.2, 0.25, 0.15, 0.45, 0.1, 0.1), nrow = G, ncol = K[8])
sim2_theta9 <- matrix(1/K[9], nrow = G, ncol = K[9])
sim2_theta10 <- matrix(1/K[10], nrow = G, ncol = K[10])
sim2_theta11 <- matrix(1/K[11], nrow = G, ncol = K[11])
sim2_theta12 <- matrix(1/K[12], nrow = G, ncol = K[12])
sim2_theta13 <- matrix(1/K[13], nrow = G, ncol = K[13])


sim2_theta <- list(sim2_theta1, sim2_theta2, sim2_theta3, sim2_theta4, sim2_theta5, sim2_theta6, sim2_theta7, sim2_theta8, sim2_theta9, sim2_theta10, sim2_theta11, sim2_theta12, sim2_theta13)


sim2_data <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 300, n_noise_var = 0)


sim2_LCR_fit <- LCR_Gibbs(X = sim2_data$X, Y = sim2_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                          theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                          relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)



