#Missing data test simulation script

#I'll take the data from simulation 1 and inject some missingness, and see how it performs

G <- 2
p <- 6
M <- 8
K <- rep(3,8)

# covariate effects vector
sim1_MD_beta <- matrix(c(-0.5, 0.8, 1.2, 1, 0.4, 0, 0), nrow = p+1)

#item probability parameters
sim1_MD_theta1 <- matrix(c(0.15,0.25,0.6,0.7,0.2,0.1), nrow = G, ncol = K[1], byrow = TRUE)
sim1_MD_theta2 <- matrix(c(0.2,0.35,0.45,0.55,0.3,0.15), nrow = G, ncol = K[2], byrow = TRUE)
sim1_MD_theta3 <- matrix(c(0.1,0.15,0.75,0.8,0.15,0.05), nrow = G, ncol = K[3], byrow = TRUE)
sim1_MD_theta4 <- matrix(c(0.25,0.4,0.35,0.45,0.35,0.2), nrow = G, ncol = K[4], byrow = TRUE)
sim1_MD_theta5 <- matrix(c(0.4, 0.5, 0.1, 0.4, 0.5, 0.1), nrow = G, ncol = K[5], byrow = TRUE)
sim1_MD_theta6 <- matrix(c(0.7, 0.1, 0.2, 0.7, 0.1, 0.2), nrow = G, ncol = K[6], byrow = TRUE)
sim1_MD_theta7 <- matrix(c(0.1, 0.5, 0.4, 0.1, 0.5, 0.4), nrow = G, ncol = K[7], byrow = TRUE)
sim1_MD_theta8 <- matrix(1/K[8], nrow = G, ncol = K[8])

sim1_MD_theta <- list(sim1_theta1, sim1_theta2, sim1_theta3, sim1_theta4, sim1_theta5, sim1_theta6, sim1_theta7, sim1_theta8)

set.seed(123)

#Simulating dataset from using the above parameters
sim1_MD_data <- LCR_sim_data(theta = sim1_MD_theta, beta = sim1_MD_beta, n_samples = 500)
sim1_MD_data_full <- sim1_MD_data

set.seed(123)
sim1_MD_data$Y <- LCR_inject_missingness(Y = sim1_MD_data$Y, prop_missing = 0.05, max_missing_per_row = 3)

set.seed(123)

sim1_MD_LCR_fit_predsel <- LCR_Gibbs(X = sim1_MD_data$X, Y = sim1_MD_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(123)

sim1_MD_LCR_fit_varsel <- LCR_Gibbs(X = sim1_MD_data$X, Y = sim1_MD_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                     theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                     relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)
