Y_CSHQ <- as.matrix(read.csv('./CSHQ_response_Y.csv'))
X_CSHQ <- as.matrix(read.csv('./CSHQ_predictor_X.csv'))[,-1]




#Testing the overfitted mixture approach for LCR for simulation study 1 - for the moment excluding the noise covariates to make things simpler

G <- 2
p <- 4
M <- 8
K <- rep(3,8)

# covariate effects vector
sim1_beta <- matrix(c(-0.5, 0.8, 1.2, 1, 0.4), nrow = p+1)

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

set.seed(123)

#Simulatiing dataset from using the above parameters
sim1_data <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 500)

set.seed(123)

sim1_data_test_fit_overfitted <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 12,                 
                                           beta_prior_cov = diag(10^2,p+1), beta_prior_mean = rep(0,p+1), 
                                           theta_hyperparam = 1, clust_var_prior = 0.5,
                                           item.sel = TRUE, cov.sel = FALSE,   
                                           sfm = TRUE, sparse_prior = 1/12,
                                           relabel = TRUE,
                                           n_samples = 5000, burnin = 5000,
                                           verbose = TRUE)
sim1_data_test_fit_overfitted$G_eff_posterior        
sim1_data_test_fit_overfitted$sparse_weights_accept_rate    
adjustedRandIndex(max.col(sim1_data_test_fit_overfitted$Z), sim1_data$class)

sim1_relevelled_beta_samples <- sim1_data_test_fit_overfitted$samples$beta_samples

for (t in 1:dim(sim1_relevelled_beta_samples)[3]){
  sim1_relevelled_beta_samples[,,t] <- sim1_relevelled_beta_samples[,,t] - sim1_relevelled_beta_samples[,1,t]
}



#Testing for simulation study 2 - for the moment excluding the noise covariates to make things simpler



G <- 3
p <- 4
M <- 13
K <- c(2,2,2,3,3,3,4,4,3,3,5,5,5)

# covariate effects matrix
sim2_beta <- matrix(c(0, 0, 1, -1, -1, 1, 0.5, -0.5, -0.4, 0.4), nrow = p+1, byrow = TRUE)



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


set.seed(129)


sim2_data <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 500)

sim2_data_test_fit_overfitted <- LCR_Gibbs(X = sim2_data$X, Y = sim2_data$Y, G = 12,                 
                                           beta_prior_cov = diag(10^2,p+1), beta_prior_mean = rep(0,p+1), 
                                           theta_hyperparam = 1, clust_var_prior = 0.5,
                                           item.sel = TRUE, cov.sel = FALSE,   
                                           sfm = TRUE, sparse_prior = 1/12,
                                           relabel = TRUE,
                                           n_samples = 5000, burnin = 5000,
                                           verbose = TRUE)


CSHQ_data_test_fit_overfitted <- LCR_Gibbs(X = X_CSHQ, Y = Y_CSHQ, G = 12,                 
                                           beta_prior_cov = diag(10^2,ncol(X_CSHQ)+1), beta_prior_mean = rep(0,ncol(X_CSHQ)+1), 
                                           theta_hyperparam = 1, clust_var_prior = 0.5,
                                           item.sel = TRUE, cov.sel = FALSE,   
                                           sfm = TRUE, sparse_prior = 0.005,
                                           relabel = FALSE,
                                           n_samples = 5000, burnin = 5000,
                                           verbose = TRUE)

#Running Wade and Ghahramani clustering on this


g_vec_samples <- apply(CSHQ_data_test_fit_overfitted$samples$z_samples, c(1,3), which.max)

psm_mat_overfitted <- comp.psm(t(g_vec_samples))

start_cl <- as.integer(factor(max.col(CSHQ_data_test_fit_overfitted$Z)))

minVI_est_overfitted_wade <- minVI(psm = psm_mat_overfitted, method = 'greedy', start.cl = start_cl)

minVI_est_overfitted_salso <- salso(t(g_vec_samples))













