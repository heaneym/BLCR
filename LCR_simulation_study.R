working_dir <- getwd()
source('LCR_Gibbs.R')
source('LCR_sim_data.R')
library(ggplot2)
library(RColorBrewer)
library(reshape2) 
library(BayesLCA)

kl_divergence_rowwise <- function(p, q) {
  if (!identical(dim(p), dim(q))) {
    stop("Matrices must have the same dimensions")
  }
  epsilon <- 1e-10
  p <- pmax(p, epsilon)
  q <- pmax(q, epsilon)
  kl_values <- numeric(nrow(p))
  for (i in 1:nrow(p)) {
    kl_values[i] <- sum(p[i, ] * log(p[i, ] / q[i, ]))
  }
  kl_values[kl_values < 0 & abs(kl_values) < 1e-10] <- 0
  return(kl_values)
}



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
sim1_theta5 <- matrix(c(0.4, 0.5, 0.1, 0.4, 0.5, 0.1), nrow = G, ncol = K[5], byrow = TRUE)
sim1_theta6 <- matrix(c(0.7, 0.1, 0.2, 0.7, 0.1, 0.2), nrow = G, ncol = K[6], byrow = TRUE)
sim1_theta7 <- matrix(c(0.1, 0.5, 0.4, 0.1, 0.5, 0.4), nrow = G, ncol = K[7], byrow = TRUE)
sim1_theta8 <- matrix(1/K[8], nrow = G, ncol = K[8])

sim1_theta <- list(sim1_theta1, sim1_theta2, sim1_theta3, sim1_theta4, sim1_theta5, sim1_theta6, sim1_theta7, sim1_theta8)

sim1_data <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 500, n_noise_var = 0)

write.csv(data.frame(cbind(sim1_data$X, sim1_data$Y)), "sim_data1.csv", row.names = FALSE)

sim1_prelim_collapsed_fit <- blca.collapsed(X = sim1_data$Y, G = 1, iter = 50000, burn.in = 1000, thin = 0.1, G.sel = TRUE, verbose = TRUE)

sim1_LCR_fit_no_varsel <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                          theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                          relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim1_LCR_fit_varsel <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                    theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                    relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim1_LCR_fit_covsel <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim1_LCR_fit_itemsel <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim1_true_labels <- sim1_data$class
sim1_LCR_fit_no_varsel_crossclass_table <- table(sim1_true_labels, max.col(sim1_LCR_fit_no_varsel$assignment_prob))
sim1_LCR_fit_varsel_crossclass_table <- table(sim1_true_labels, max.col(sim1_LCR_fit_varsel$assignment_prob))
sim1_LCR_fit_covsel_crossclass_table <- table(sim1_true_labels, max.col(sim1_LCR_fit_covsel$assignment_prob))
sim1_LCR_fit_itemsel_crossclass_table <- table(sim1_true_labels, max.col(sim1_LCR_fit_covsel$assignment_prob))

sim1_LCR_fit_varsel_vs_no_varsel_table <- table(max.col(sim1_LCR_fit_no_varsel$assignment_prob), max.col(sim1_LCR_fit_varsel$assignment_prob))

sim1_LCR_fit_no_varsel_ARI <- adj.rand.index(sim1_true_labels, max.col(sim1_LCR_fit_no_varsel$assignment_prob))
sim1_LCR_fit_varsel_ARI <- adj.rand.index(sim1_true_labels, max.col(sim1_LCR_fit_varsel$assignment_prob))
sim1_LCR_fit_covsel_ARI <- adj.rand.index(sim1_true_labels, max.col(sim1_LCR_fit_covsel$assignment_prob))
sim1_LCR_fit_itemsel_ARI <- adj.rand.index(sim1_true_labels, max.col(sim1_LCR_fit_covsel$assignment_prob))

sim1_LCR_fit_varsel_vs_no_varsel_ARI <- adj.rand.index(max.col(sim1_LCR_fit_no_varsel$assignment_prob), max.col(sim1_LCR_fit_varsel$assignment_prob))

png("./sim_study_plots/sim_study1_plots/sim1_no_varsel_log_post.png", width = 800, height = 600)
plot(sim1_LCR_fit_no_varsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 1', xlab = 'Iteration Number', ylab = 'log posterior')
dev.off()

png("./sim_study_plots/sim_study1_plots/sim1_varsel_log_post.png", width = 800, height = 600)
plot(sim1_LCR_fit_varsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 1', xlab = 'Iteration Number', ylab = 'log posterior')
dev.off()

png("./sim_study_plots/sim_study1_plots/sim1_covsel_log_post.png", width = 800, height = 600)
plot(sim1_LCR_fit_covsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 1', xlab = 'Iteration Number', ylab = 'log posterior')
dev.off()

png("./sim_study_plots/sim_study1_plots/sim1_itemsel_log_post.png", width = 800, height = 600)
plot(sim1_LCR_fit_itemsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 1', xlab = 'Iteration Number', ylab = 'log posterior')
dev.off()

png("./sim_study_plots/sim_study1_plots/sim1_no_varsel_log_post_acf.png", width = 800, height = 600)
acf(sim1_LCR_fit_no_varsel$samples$log_post_samples, main = 'Simulated Data 1')
dev.off()

png("./sim_study_plots/sim_study1_plots/sim1_varsel_log_post_acf.png", width = 800, height = 600)
acf(sim1_LCR_fit_varsel$samples$log_post_samples, main = 'Simulated Data 1')
dev.off()

png("./sim_study_plots/sim_study1_plots/sim1_covsel_log_post_acf.png", width = 800, height = 600)
acf(sim1_LCR_fit_covsel$samples$log_post_samples, main = 'Simulated Data 1')
dev.off()

png("./sim_study_plots/sim_study1_plots/sim1_itemsel_log_post_acf.png", width = 800, height = 600)
acf(sim1_LCR_fit_itemsel$samples$log_post_samples, main = 'Simulated Data 1')
dev.off()


png("./sim_study_plots/sim_study1_plots/beta_plots/beta0_density_plot.png", width = 800, height = 600)
plot(density(sim1_LCR_fit_no_varsel$samples$beta_samples[1,1,]), main = expression(paste('Density plot of ', beta[0])), xlab = expression(paste('Value of ', beta[0])))
rug(sim1_LCR_fit_no_varsel$samples$beta_samples[1,1,])
dev.off()

png("./sim_study_plots/sim_study1_plots/beta_plots/beta1_density_plot.png", width = 800, height = 600)
plot(density(sim1_LCR_fit_no_varsel$samples$beta_samples[2,1,]), main = expression(paste('Density plot of ', beta[1])), xlab = expression(paste('Value of ', beta[1])))
rug(sim1_LCR_fit_no_varsel$samples$beta_samples[2,1,])
dev.off()

png("./sim_study_plots/sim_study1_plots/beta_plots/beta2_density_plot.png", width = 800, height = 600)
plot(density(sim1_LCR_fit_no_varsel$samples$beta_samples[3,1,]), main = expression(paste('Density plot of ', beta[2])), xlab = expression(paste('Value of ', beta[2])))
rug(sim1_LCR_fit_no_varsel$samples$beta_samples[3,1,])
dev.off()

png("./sim_study_plots/sim_study1_plots/beta_plots/beta3_density_plot.png", width = 800, height = 600)
plot(density(sim1_LCR_fit_no_varsel$samples$beta_samples[4,1,]), main = expression(paste('Density plot of ', beta[3])), xlab = expression(paste('Value of ', beta[3])))
rug(sim1_LCR_fit_no_varsel$samples$beta_samples[4,1,])
dev.off()

png("./sim_study_plots/sim_study1_plots/beta_plots/beta4_density_plot.png", width = 800, height = 600)
plot(density(sim1_LCR_fit_no_varsel$samples$beta_samples[5,1,]), main = expression(paste('Density plot of ', beta[4])), xlab = expression(paste('Value of ', beta[4])))
rug(sim1_LCR_fit_no_varsel$samples$beta_samples[5,1,])
dev.off()

png("./sim_study_plots/sim_study1_plots/beta_plots/beta5_density_plot.png", width = 800, height = 600)
plot(density(sim1_LCR_fit_no_varsel$samples$beta_samples[6,1,]), main = expression(paste('Density plot of ', beta[5])), xlab = expression(paste('Value of ', beta[5])))
rug(sim1_LCR_fit_no_varsel$samples$beta_samples[6,1,])
dev.off()

png("./sim_study_plots/sim_study1_plots/beta_plots/beta6_density_plot.png", width = 800, height = 600)
plot(density(sim1_LCR_fit_no_varsel$samples$beta_samples[7,1,]), main = expression(paste('Density plot of ', beta[6])), xlab = expression(paste('Value of ', beta[6])))
rug(sim1_LCR_fit_no_varsel$samples$beta_samples[7,1,])
dev.off()


#Checking variable importance using KL divergence


assignment_probs_full_sim1 <- sim1_data$class_prob


assignment_probs_X1_sim1_mu <- cbind(1,sim1_data$X[,-1])%*%sim1_beta[-2,]
assignment_probs_X2_sim1_mu <- cbind(1,sim1_data$X[,-2])%*%sim1_beta[-3,]
assignment_probs_X3_sim1_mu <- cbind(1,sim1_data$X[,-3])%*%sim1_beta[-4,]
assignment_probs_X4_sim1_mu <- cbind(1,sim1_data$X[,-4])%*%sim1_beta[-5,]
assignment_probs_X5_sim1_mu <- cbind(1,sim1_data$X[,-5])%*%sim1_beta[-6,]
assignment_probs_X6_sim1_mu <- cbind(1,sim1_data$X[,-6])%*%sim1_beta[-7,]

assignment_matrix_X1_sim1 <- cbind(exp(assignment_probs_X1_sim1_mu)/(1+exp(assignment_probs_X1_sim1_mu)),1-exp(assignment_probs_X1_sim1_mu)/(1+exp(assignment_probs_X1_sim1_mu)))
assignment_matrix_X2_sim1 <- cbind(exp(assignment_probs_X2_sim1_mu)/(1+exp(assignment_probs_X2_sim1_mu)),1-exp(assignment_probs_X2_sim1_mu)/(1+exp(assignment_probs_X2_sim1_mu)))
assignment_matrix_X3_sim1 <- cbind(exp(assignment_probs_X3_sim1_mu)/(1+exp(assignment_probs_X3_sim1_mu)),1-exp(assignment_probs_X3_sim1_mu)/(1+exp(assignment_probs_X3_sim1_mu)))
assignment_matrix_X4_sim1 <- cbind(exp(assignment_probs_X4_sim1_mu)/(1+exp(assignment_probs_X4_sim1_mu)),1-exp(assignment_probs_X4_sim1_mu)/(1+exp(assignment_probs_X4_sim1_mu)))
assignment_matrix_X5_sim1 <- cbind(exp(assignment_probs_X5_sim1_mu)/(1+exp(assignment_probs_X5_sim1_mu)),1-exp(assignment_probs_X5_sim1_mu)/(1+exp(assignment_probs_X5_sim1_mu)))
assignment_matrix_X6_sim1 <- cbind(exp(assignment_probs_X6_sim1_mu)/(1+exp(assignment_probs_X6_sim1_mu)),1-exp(assignment_probs_X6_sim1_mu)/(1+exp(assignment_probs_X6_sim1_mu)))

X1_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_probs_full_sim1, assignment_matrix_X1_sim1)
X2_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_probs_full_sim1, assignment_matrix_X2_sim1)
X3_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_probs_full_sim1, assignment_matrix_X3_sim1)
X4_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_probs_full_sim1, assignment_matrix_X4_sim1)
X5_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_probs_full_sim1, assignment_matrix_X5_sim1)
X6_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_probs_full_sim1, assignment_matrix_X6_sim1)

kl_df_sim1 <- data.frame(X1 = X1_kl_divergence_sim1,
                         X2 = X2_kl_divergence_sim1,
                         X3 = X3_kl_divergence_sim1,
                         X4 = X4_kl_divergence_sim1) 
data_long_kl_div_sim1 <- melt(kl_df_sim1)


ggplot(data_long_kl_div_sim1, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  coord_cartesian(xlim = c(0, 0.2)) +  
  labs(
    title = "Density Plots of KL Divergence between Full Model and Reduced Model",
    x = "KL divergence value",
    y = "Density",
    fill = "variable excluded"   
  ) +
  theme_classic()

ggsave(filename = './sim_study_plots/sim_study1_plots/KL_divergence_plots_sim_study1/KL_diverg_sim_study1.png', dpi = 600, width = 7, height = 5)

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
sim2_beta <- matrix(c(0, 0, 1, -1, -1, 1, 0.5, -0.5, -0.4, 0.4, 0, 0, 0, 0), nrow = p+1, byrow = TRUE)
#sim2_beta <- matrix(c(-1, 0.8, 1.1, 1.2, 0.3, 0, 0, -0.3, 0.5, 0, 0.4, 0, 0, 0), nrow = p+1)


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


sim2_data1 <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 500)
sim2_data2 <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 500)
sim2_data3 <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 500)
sim2_data4 <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 500)
sim2_data5 <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 500)

sim2_data_frame_binded1 <- data.frame(Y = sim2_data1$Y, X = sim2_data1$X)
sim2_data_frame_binded2 <- data.frame(Y = sim2_data2$Y, X = sim2_data2$X)
sim2_data_frame_binded3 <- data.frame(Y = sim2_data3$Y, X = sim2_data3$X)
sim2_data_frame_binded4 <- data.frame(Y = sim2_data4$Y, X = sim2_data4$X)
sim2_data_frame_binded5 <- data.frame(Y = sim2_data5$Y, X = sim2_data5$X)

write.csv(data.frame(cbind(sim2_data1$X, sim2_data1$Y)), "./sim_study_plots/sim_study2_plots/sim2_data1.csv", row.names = FALSE)
write.csv(data.frame(cbind(sim2_data2$X, sim2_data2$Y)), "./sim_study_plots/sim_study2_plots/sim2_data2.csv", row.names = FALSE)
write.csv(data.frame(cbind(sim2_data3$X, sim2_data3$Y)), "./sim_study_plots/sim_study2_plots/sim2_data3.csv", row.names = FALSE)
write.csv(data.frame(cbind(sim2_data4$X, sim2_data4$Y)), "./sim_study_plots/sim_study2_plots/sim2_data4.csv", row.names = FALSE)
write.csv(data.frame(cbind(sim2_data5$X, sim2_data5$Y)), "./sim_study_plots/sim_study2_plots/sim2_data5.csv", row.names = FALSE)

sim2_prelim_collapsed_fit1 <- blca.collapsed(X = sim2_data1$Y, G = 1, iter = 50000, burn.in = 1000, thin = 0.1, G.sel = TRUE, verbose = TRUE)
sim2_prelim_collapsed_fit2 <- blca.collapsed(X = sim2_data2$Y, G = 1, iter = 50000, burn.in = 1000, thin = 0.1, G.sel = TRUE, verbose = TRUE)
sim2_prelim_collapsed_fit3 <- blca.collapsed(X = sim2_data3$Y, G = 1, iter = 50000, burn.in = 1000, thin = 0.1, G.sel = TRUE, verbose = TRUE)
sim2_prelim_collapsed_fit4 <- blca.collapsed(X = sim2_data4$Y, G = 1, iter = 50000, burn.in = 1000, thin = 0.1, G.sel = TRUE, verbose = TRUE)
sim2_prelim_collapsed_fit5 <- blca.collapsed(X = sim2_data5$Y, G = 1, iter = 50000, burn.in = 1000, thin = 0.1, G.sel = TRUE, verbose = TRUE)

#Running with no variable selection

sim2_LCR_fit_no_varsel1 <- LCR_Gibbs(X = sim2_data1$X, Y = sim2_data1$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                          theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                          relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_no_varsel2 <- LCR_Gibbs(X = sim2_data2$X, Y = sim2_data2$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                     theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                                     relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_no_varsel3 <- LCR_Gibbs(X = sim2_data3$X, Y = sim2_data3$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                     theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                                     relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_no_varsel4 <- LCR_Gibbs(X = sim2_data4$X, Y = sim2_data4$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                     theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                                     relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_no_varsel5 <- LCR_Gibbs(X = sim2_data5$X, Y = sim2_data5$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                     theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                                     relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

#Running with full variable selection


sim2_LCR_fit_varsel1 <- LCR_Gibbs(X = sim2_data1$X, Y = sim2_data1$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                    theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                    relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)
sim2_LCR_fit_varsel2 <- LCR_Gibbs(X = sim2_data2$X, Y = sim2_data2$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)
sim2_LCR_fit_varsel3 <- LCR_Gibbs(X = sim2_data3$X, Y = sim2_data3$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)
sim2_LCR_fit_varsel4 <- LCR_Gibbs(X = sim2_data4$X, Y = sim2_data4$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)
sim2_LCR_fit_varsel5 <- LCR_Gibbs(X = sim2_data5$X, Y = sim2_data5$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

#Running with only predictor selection

sim2_LCR_fit_covsel1 <- LCR_Gibbs(X = sim2_data1$X, Y = sim2_data1$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_covsel2 <- LCR_Gibbs(X = sim2_data2$X, Y = sim2_data2$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_covsel3 <- LCR_Gibbs(X = sim2_data3$X, Y = sim2_data3$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_covsel4 <- LCR_Gibbs(X = sim2_data4$X, Y = sim2_data4$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_covsel5 <- LCR_Gibbs(X = sim2_data5$X, Y = sim2_data5$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

#Running with only item selection


sim2_LCR_fit_itemsel1 <- LCR_Gibbs(X = sim2_data1$X, Y = sim2_data1$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_itemsel2 <- LCR_Gibbs(X = sim2_data2$X, Y = sim2_data2$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                  theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                  relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_itemsel3 <- LCR_Gibbs(X = sim2_data3$X, Y = sim2_data3$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                  theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                  relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_itemsel4 <- LCR_Gibbs(X = sim2_data4$X, Y = sim2_data4$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                  theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                  relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_fit_itemsel5 <- LCR_Gibbs(X = sim2_data5$X, Y = sim2_data5$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                  theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                  relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


#Checking variable importance using KL divergence - need to look at this across all of the datasets I think


assignment_probs_full_sim2 <- sim2_data$class_prob


assignment_probs_X1_sim2_mu <- cbind(1,sim2_data$X[,-1])%*%sim2_beta[-2,]
assignment_probs_X2_sim2_mu <- cbind(1,sim2_data$X[,-2])%*%sim2_beta[-3,]
assignment_probs_X3_sim2_mu <- cbind(1,sim2_data$X[,-3])%*%sim2_beta[-4,]
assignment_probs_X4_sim2_mu <- cbind(1,sim2_data$X[,-4])%*%sim2_beta[-5,]
assignment_probs_X5_sim2_mu <- cbind(1,sim2_data$X[,-5])%*%sim2_beta[-6,]
assignment_probs_X6_sim2_mu <- cbind(1,sim2_data$X[,-6])%*%sim2_beta[-7,]

exp_mu_sim2_X1 <- exp(cbind(assignment_probs_X1_sim2_mu, 0))
exp_mu_sim2_X2 <- exp(cbind(assignment_probs_X2_sim2_mu, 0))
exp_mu_sim2_X3 <- exp(cbind(assignment_probs_X3_sim2_mu, 0))
exp_mu_sim2_X4 <- exp(cbind(assignment_probs_X4_sim2_mu, 0))
exp_mu_sim2_X5 <- exp(cbind(assignment_probs_X5_sim2_mu, 0))
exp_mu_sim2_X6 <- exp(cbind(assignment_probs_X6_sim2_mu, 0))

assignment_matrix_X1_sim2 <- exp_mu_sim2_X1/rowSums(exp_mu_sim2_X1)
assignment_matrix_X2_sim2 <- exp_mu_sim2_X2/rowSums(exp_mu_sim2_X2)
assignment_matrix_X3_sim2 <- exp_mu_sim2_X3/rowSums(exp_mu_sim2_X3)
assignment_matrix_X4_sim2 <- exp_mu_sim2_X4/rowSums(exp_mu_sim2_X4)
assignment_matrix_X5_sim2 <- exp_mu_sim2_X5/rowSums(exp_mu_sim2_X5)
assignment_matrix_X6_sim2 <- exp_mu_sim2_X6/rowSums(exp_mu_sim2_X6)



X1_kl_divergence_sim2 <- kl_divergence_rowwise(assignment_probs_full_sim2, assignment_matrix_X1_sim2)
X2_kl_divergence_sim2 <- kl_divergence_rowwise(assignment_probs_full_sim2, assignment_matrix_X2_sim2)
X3_kl_divergence_sim2 <- kl_divergence_rowwise(assignment_probs_full_sim2, assignment_matrix_X3_sim2)
X4_kl_divergence_sim2 <- kl_divergence_rowwise(assignment_probs_full_sim2, assignment_matrix_X4_sim2)
X5_kl_divergence_sim2 <- kl_divergence_rowwise(assignment_probs_full_sim2, assignment_matrix_X5_sim2)
X6_kl_divergence_sim2 <- kl_divergence_rowwise(assignment_probs_full_sim2, assignment_matrix_X6_sim2)

kl_df_sim2 <- data.frame(X1 = X1_kl_divergence_sim2,
                         X2 = X2_kl_divergence_sim2,
                         X3 = X3_kl_divergence_sim2,
                         X4 = X4_kl_divergence_sim2) 
data_long_kl_div_sim2 <- melt(kl_df_sim2)


ggplot(data_long_kl_div_sim2, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  coord_cartesian(xlim = c(0, 0.2)) +  # Focus on where your data actually is
  labs(title = "Density Plots of KL Divergence between Full Model and Reduced Model",
       x = "Value",
       y = "Density") +
  theme_minimal()



#Looking at the Dean and Raftery non-binary simulated dataset


G <- 3
p <- 6
M <- 10
K <- c(3,2,4,3,3,4,5,2,3,4)

# covariate effects matrix
dr_beta <- matrix(c(-1, 0.8, 1.1, 1.2, 0.3, 0, 0, -0.3, 0.5, 0, 0.4, 0, 0, 0), nrow = p+1)

dr_theta1 <- matrix(c(0.1,0.3,0.6,0.1,0.5,0.2,0.8,0.2,0.2), nrow = G, ncol = K[1])
dr_theta2 <- matrix(c(0.5, 0.1, 0.7, 0.5, 0.9, 0.3), nrow = G, ncol = K[2])
dr_theta3 <- matrix(c(0.2, 0.7, 0.2, 0.2, 0.1, 0.6, 0.3, 0.1, 0.1, 0.3, 0.1, 0.1), nrow = G, ncol = K[3])
dr_theta4 <- matrix(c(0.1, 0.6, 0.4, 0.5, 0.1, 0.4, 0.4, 0.3, 0.2), nrow = G, ncol = K[4])
dr_theta5 <- matrix(c(0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1), nrow = G, ncol = K[5])
dr_theta6 <- matrix(c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3), nrow = G, ncol = K[6])
dr_theta7 <- matrix(c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = G, ncol = K[7])
dr_theta8 <- matrix(c(0.2, 0.2, 0.2, 0.8, 0.8, 0.8), nrow = G, ncol = K[8])
dr_theta9 <- matrix(c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2), nrow = G, ncol = K[9])
dr_theta10 <- matrix(c(0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.6, 0.6, 0.6), nrow = G, ncol = K[10])


dr_theta <- list(dr_theta1, dr_theta2, dr_theta3, dr_theta4, dr_theta5, dr_theta6, dr_theta7)#, dr_theta8, dr_theta9, dr_theta10)


dr_sim_data <- LCR_sim_data(theta = dr_theta, beta = dr_beta, n_samples = 500, n_noise_var = 0)


dr_LCR_fit_no_varsel <- LCR_Gibbs(X = dr_sim_data$X, Y = dr_sim_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                   theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                                   relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 1)

dr_LCR_fit_itemsel <- LCR_Gibbs(X = dr_sim_data$X, Y = dr_sim_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                  theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                  relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 1)





















#Running with the CSHQ sleep survey data, we consider as covariates the variables 
setwd(working_dir)
setwd("..") 
CSHQ_data <- read.spss('./CSHQ_analysis/LynnWalsh thesis data files/Thesis Cleaned File.sav')
setwd(working_dir) 

data <- as.data.frame(data)
sleep_data_responses <- data[,65:99]

mapping <- c(
  "Rarely (0-1 times per week)" = 1,
  "Sometimes (2-4 times per week)" = 2,
  "Usually (5-7 times per week)" = 3
)

sleep_data_responses_numeric <- sleep_data_responses %>% mutate_all(~ mapping[as.character(.)])
#removing bedtime resistance 4 and 6 since they are double counted from sleep anxiety
sleep_data_responses_numeric <- sleep_data_responses_numeric[,-c(4,6)]


df <- cbind(data[,c(35:48,56)],sleep_data_responses_numeric)
df$Diagnosis_ID <- ifelse(is.na(df$Diagnosis_ID),0,1)
df$Diagnosis_CerebralPalsy <- ifelse(is.na(df$Diagnosis_CerebralPalsy),0,1)
df$Diagnosis_ASD <- ifelse(is.na(df$Diagnosis_ASD),0,1)
df$Diagnosis_ADHD <- ifelse(is.na(df$Diagnosis_ADHD),0,1)
df$Diagnosis_DS <- ifelse(is.na(df$Diagnosis_DS),0,1)
df$Diagnosis_Other <- ifelse(is.na(df$Diagnosis_Other),0,1)
df$Diagnosis_None <- ifelse(is.na(df$Diagnosis_None),0,1)


#I'm gonna remove rows that are almost all NA (all indicator variables are missing)

df <- df[complete.cases(df[,16]),]
colnames(df)[16:48] <-c('BR1','BR2','BR3','BR5','SOD1','SD1','SD2','SD3','SA1','SA2','SA3','SA4','NW1','NW2','NW3','P1','P2','P3','P4','P5','P6','P7','SDB1', 'SDB2', 'SDB3', 'DS1','DS2','DS3','DS4','DS5','DS6','DS7','DS8')

#reverse scoring the negative items

df$BR1 <- 4 - df$BR1
df$BR2 <- 4 - df$BR2
df$SOD1 <- 4 - df$SOD1
df$SD2 <- 4 - df$SD2
df$SD3 <- 4 - df$SD3
df$DS1 <- 4 - df$DS1

sleep_data_total_score <- rowSums(df[,16:48])
sleep_disorders <- ifelse(sleep_data_total_score >= 41,1,0)
df <- cbind(df,sleep_data_total_score,sleep_disorders)
df <- df[,-c(8,11)]

# We want to look at a variable for ASD (diagnosis_ASD), one for Intellectual Disability (diagnosis_ID), and then other, which includes anything else.
df$other_diagnoses <- df$Diagnosis_Other + df$Diagnosis_ADHD
df$other_diagnoses <- ifelse(df$other_diagnoses>0, 1, 0)

#Need to check but I'm not sure how relevant the asthma one is here
df$other_diagnoses[107] <- 0

df$interaction_ASD_ID <- (df$Diagnosis_ASD)*(df$Diagnosis_ID)
df$interaction_ASD_other <- (df$Diagnosis_ASD)*(df$other_diagnoses)
df$interaction_ID_other <- (df$Diagnosis_ID)*(df$other_diagnoses)

nd_vec <- rowSums(df[,8:10])
nd_vec[c(83,107,143)] <- 0
nd_vec[which(nd_vec==2 | nd_vec == 3) ] <- 1
nd_vec <- as.numeric(nd_vec)
df <- cbind(df, nd_vec)

CSHQ_Y <- as.matrix(df[,14:46])
#should I include interaction terms here?
df$scaled_age <- scale(df$ChildAge)
CSHQ_X <- as.matrix(cbind(df$Diagnosis_ID, df$Diagnosis_ASD, df$other_diagnoses, df$scaled_age))



CSHQ_LCR_fit_6group <- LCR_Gibbs(X = CSHQ_X, Y = CSHQ_Y, G = 6, beta_prior_cov = diag(10^2,5), beta_prior_mean = rep(0,5), 
                          theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                          relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

CSHQ_LCR_fit_5group <- LCR_Gibbs(X = CSHQ_X, Y = CSHQ_Y, G = 5, beta_prior_cov = diag(10^2,5), beta_prior_mean = rep(0,5), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

CSHQ_LCR_fit_4group <- LCR_Gibbs(X = CSHQ_X, Y = CSHQ_Y, G = 4, beta_prior_cov = diag(10^2,5), beta_prior_mean = rep(0,5), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)



