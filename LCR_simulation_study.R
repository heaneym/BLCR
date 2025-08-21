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
# 500 observations

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

sim1_data <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 500)

write.csv(data.frame(cbind(sim1_data$X, sim1_data$Y, sim1_data$class)), "sim_data1.csv", row.names = FALSE)

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


assignment_probs_full_sim1_mu <- cbind(1,sim1_data$X)%*%sim1_beta


assignment_probs_X1_sim1_mu <- cbind(1,sim1_data$X[,-1])%*%sim1_beta[-2,]
assignment_probs_X2_sim1_mu <- cbind(1,sim1_data$X[,-2])%*%sim1_beta[-3,]
assignment_probs_X3_sim1_mu <- cbind(1,sim1_data$X[,-3])%*%sim1_beta[-4,]
assignment_probs_X4_sim1_mu <- cbind(1,sim1_data$X[,-4])%*%sim1_beta[-5,]
assignment_probs_X5_sim1_mu <- cbind(1,sim1_data$X[,-5])%*%sim1_beta[-6,]
assignment_probs_X6_sim1_mu <- cbind(1,sim1_data$X[,-6])%*%sim1_beta[-7,]

assignment_matrix_full_sim1 <- cbind(exp(assignment_probs_full_sim1_mu)/(1+exp(assignment_probs_full_sim1_mu)),1-exp(assignment_probs_full_sim1_mu)/(1+exp(assignment_probs_full_sim1_mu)))
assignment_matrix_X1_sim1 <- cbind(exp(assignment_probs_X1_sim1_mu)/(1+exp(assignment_probs_X1_sim1_mu)),1-exp(assignment_probs_X1_sim1_mu)/(1+exp(assignment_probs_X1_sim1_mu)))
assignment_matrix_X2_sim1 <- cbind(exp(assignment_probs_X2_sim1_mu)/(1+exp(assignment_probs_X2_sim1_mu)),1-exp(assignment_probs_X2_sim1_mu)/(1+exp(assignment_probs_X2_sim1_mu)))
assignment_matrix_X3_sim1 <- cbind(exp(assignment_probs_X3_sim1_mu)/(1+exp(assignment_probs_X3_sim1_mu)),1-exp(assignment_probs_X3_sim1_mu)/(1+exp(assignment_probs_X3_sim1_mu)))
assignment_matrix_X4_sim1 <- cbind(exp(assignment_probs_X4_sim1_mu)/(1+exp(assignment_probs_X4_sim1_mu)),1-exp(assignment_probs_X4_sim1_mu)/(1+exp(assignment_probs_X4_sim1_mu)))
assignment_matrix_X5_sim1 <- cbind(exp(assignment_probs_X5_sim1_mu)/(1+exp(assignment_probs_X5_sim1_mu)),1-exp(assignment_probs_X5_sim1_mu)/(1+exp(assignment_probs_X5_sim1_mu)))
assignment_matrix_X6_sim1 <- cbind(exp(assignment_probs_X6_sim1_mu)/(1+exp(assignment_probs_X6_sim1_mu)),1-exp(assignment_probs_X6_sim1_mu)/(1+exp(assignment_probs_X6_sim1_mu)))


#Doing this again but with the posterior probabilities

# posterior_assignment_mat_full_sim1 <- LCR_posterior_membership_prob(theta = sim1_theta, beta = cbind(sim1_beta[,,drop = FALSE], 0), X = sim1_data$X, Y = sim1_data$Y)
# posterior_assignment_mat_X1_sim1 <- LCR_posterior_membership_prob(theta = sim1_theta, beta = cbind(sim1_beta[-2,, drop = FALSE],0), X = sim1_data$X[,-1], Y = sim1_data$Y)
# posterior_assignment_mat_X2_sim1 <- LCR_posterior_membership_prob(theta = sim1_theta, beta = cbind(sim1_beta[-3,, drop = FALSE],0), X = sim1_data$X[,-2], Y = sim1_data$Y)
# posterior_assignment_mat_X3_sim1 <- LCR_posterior_membership_prob(theta = sim1_theta, beta = cbind(sim1_beta[-4,, drop = FALSE],0), X = sim1_data$X[,-3], Y = sim1_data$Y)
# posterior_assignment_mat_X4_sim1 <- LCR_posterior_membership_prob(theta = sim1_theta, beta = cbind(sim1_beta[-5,, drop = FALSE],0), X = sim1_data$X[,-4], Y = sim1_data$Y)
# posterior_assignment_mat_X5_sim1 <- LCR_posterior_membership_prob(theta = sim1_theta, beta = cbind(sim1_beta[-6,, drop = FALSE],0), X = sim1_data$X[,-5], Y = sim1_data$Y)
# posterior_assignment_mat_X6_sim1 <- LCR_posterior_membership_prob(theta = sim1_theta, beta = cbind(sim1_beta[-7,, drop = FALSE],0), X = sim1_data$X[,-6], Y = sim1_data$Y)
# 
# posterior_kl_div_X1_sim1 <- kl_divergence_rowwise(posterior_assignment_mat_full_sim1, posterior_assignment_mat_X1_sim1)
# posterior_kl_div_X2_sim1 <- kl_divergence_rowwise(posterior_assignment_mat_full_sim1, posterior_assignment_mat_X2_sim1)
# posterior_kl_div_X3_sim1 <- kl_divergence_rowwise(posterior_assignment_mat_full_sim1, posterior_assignment_mat_X3_sim1)
# posterior_kl_div_X4_sim1 <- kl_divergence_rowwise(posterior_assignment_mat_full_sim1, posterior_assignment_mat_X4_sim1)
# posterior_kl_div_X5_sim1 <- kl_divergence_rowwise(posterior_assignment_mat_full_sim1, posterior_assignment_mat_X5_sim1)
# posterior_kl_div_X6_sim1 <- kl_divergence_rowwise(posterior_assignment_mat_full_sim1, posterior_assignment_mat_X6_sim1)


X1_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_matrix_full_sim1, assignment_matrix_X1_sim1)
X2_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_matrix_full_sim1, assignment_matrix_X2_sim1)
X3_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_matrix_full_sim1, assignment_matrix_X3_sim1)
X4_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_matrix_full_sim1, assignment_matrix_X4_sim1)
X5_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_matrix_full_sim1, assignment_matrix_X5_sim1)
X6_kl_divergence_sim1 <- kl_divergence_rowwise(assignment_matrix_full_sim1, assignment_matrix_X6_sim1)

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



#Plotting densities of kl divergences for posterior
kl_df_sim1_posterior <- data.frame(X1 = posterior_kl_div_X1_sim1,
                         X2 = posterior_kl_div_X2_sim1,
                         X3 = posterior_kl_div_X3_sim1,
                         X4 = posterior_kl_div_X4_sim1) 
data_long_kl_div_sim1_posterior <- melt(kl_df_sim1_posterior)


ggplot(data_long_kl_div_sim1_posterior, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  coord_cartesian(xlim = c(0, 0.10)) +  
  labs(
    title = "Density Plots of KL Divergence between Full Model and Reduced Model",
    x = "KL divergence value",
    y = "Density",
    fill = "variable excluded"   
  ) +
  theme_classic()

ggsave(filename = './sim_study_plots/sim_study1_plots/KL_divergence_plots_sim_study1/KL_diverg_sim_study1_posterior.png', dpi = 600, width = 7, height = 5)




#creating a forest plot to assess HDI overall.

sim1_beta_names <- paste0("beta", 0:6)
true_beta_vals_sim1 <- sim1_beta[,1]
beta_means_sim1 <- apply(-sim1_LCR_fit_no_varsel$samples$beta_samples[,1,], 1, mean)
beta_sd_sim1 <- apply(-sim1_LCR_fit_no_varsel$samples$beta_samples[,1,], 1, sd)
hdi_list_sim1 <- apply(-sim1_LCR_fit_no_varsel$samples$beta_samples[,1,], 1, function(x) hdi(x, ci = 0.95))
hdi_lower_sim1 <- sapply(hdi_list_sim1, function(x) x$CI_low)
hdi_upper_sim1 <- sapply(hdi_list_sim1, function(x) x$CI_high)

plot_data <- data.frame(
  Parameter = sim1_beta_names,
  Mean = beta_means_sim1,
  SD = beta_sd_sim1,
  HDI_low = hdi_lower_sim1,
  HDI_high = hdi_upper_sim1,
  TrueValue = true_beta_vals_sim1
)


plot_data$Parameter <- factor(plot_data$Parameter, levels = rev(sim1_beta_names))

latex_labels <- expression(
  beta[6], beta[5], beta[4], beta[3], beta[2], beta[1], beta[0]
)

cb_blue <- "#4477AA"
cb_orange <- "#EE7733"

ggplot(plot_data, aes(y = Parameter, x = Mean)) +
  geom_point(size = 3, color = cb_blue) +
  geom_errorbarh(aes(xmin = HDI_low, xmax = HDI_high), height = 0.2, color = cb_blue) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(aes(x = TrueValue), shape = 18, size = 3, color = cb_orange) +
  geom_text(aes(x = TrueValue, label = round(TrueValue, 2)), 
            color = cb_orange, vjust = -0.7, hjust = 0.5, size = 3) +
  scale_y_discrete(labels = latex_labels) +
  labs(x = "Logit Coefficient", y = "", 
       title = "Forest Plot of Logit Coefficients with 95% HDI",
       subtitle = "Blue: Posterior mean & 95% HDI, Orange: True value") +
  theme_minimal()

ggsave('./sim_study_plots/sim_study1_plots/beta_plots/sim1_forest_plot.png', height = 8, width = 8, dpi = 600)




#Attempting to create a ridgeline plot to combine density plots with forest plots somehow

posterior_samples <- sim1_LCR_fit_no_varsel$samples$beta_samples[,1,]


sim1_beta_names <- paste0("beta", 0:6)
rownames(posterior_samples) <- sim1_beta_names

long_posterior_data <- as.data.frame(t(posterior_samples)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )

long_posterior_data$Parameter <- factor(long_posterior_data$Parameter, levels = rev(sim1_beta_names))

true_beta_vals_sim1 <- sim1_beta[,1]
true_values_df <- data.frame(
  Parameter = sim1_beta_names,
  TrueValue = true_beta_vals_sim1
)

cb_blue <- "#4477AA"
cb_orange <- "#EE7733"
latex_labels <- TeX(rev(c(
  "$\\beta_0$", "$\\beta_1$", "$\\beta_2$", "$\\beta_3$", 
  "$\\beta_4$", "$\\beta_5$", "$\\beta_6$"
)))

math_labels <- expression(
  beta[6], beta[5], beta[4], beta[3], beta[2], beta[1], beta[0]
)


final_plot <- ggplot(long_posterior_data, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  
  geom_point(
    data = true_values_df,
    aes(x = TrueValue, y = Parameter),
    color = cb_orange,
    shape = 18,
    size = 3
  ) +
  
  geom_text(
    data = true_values_df,
    aes(x = TrueValue, y = Parameter, label = round(TrueValue, 2)),
    color = cb_orange,
    vjust = 2,
    size = 3
  ) +
  
  scale_y_discrete(labels = math_labels) + 
  
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions",
    subtitle = "Posterior density with median & 95% HDI, Orange: True value."
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 

  )

ggsave(
  "./sim_study_plots/sim_study1_plots/beta_plots/ridgeline_plot_sim1.png", 
  plot = final_plot, 
  width = 10, 
  height = 6, 
  dpi = 600 
)



#Clustering point estimates using minimum expected variation of information

classification_mat_sim1_no_varsel <- apply(sim1_LCR_fit_no_varsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
})
classification_mat_sim1_varsel <- apply(sim1_LCR_fit_varsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
})
classification_mat_sim1_covsel <- apply(sim1_LCR_fit_covsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
})
classification_mat_sim1_itemsel <- apply(sim1_LCR_fit_itemsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
})
psm_mat_sim1_no_varsel <- comp.psm(t(classification_mat_sim1_no_varsel))
psm_mat_sim1_varsel <- comp.psm(t(classification_mat_sim1_varsel))
psm_mat_sim1_covsel <- comp.psm(t(classification_mat_sim1_covsel))
psm_mat_sim1_itemsel <- comp.psm(t(classification_mat_sim1_itemsel))
  
minVI_cluster_sim_study1_no_varsel <- minVI(psm_mat_sim1_no_varsel, method = 'greedy', suppress.comment = FALSE, cls.draw = t(classification_mat_sim1_no_varsel), start.cl = max.col(sim1_LCR_fit_no_varsel$Z))
minVI_cluster_sim_study1_varsel <- minVI(psm_mat_sim1_varsel, method = 'greedy', suppress.comment = FALSE, cls.draw = t(classification_mat_sim1_varsel), start.cl = max.col(sim1_LCR_fit_varsel$Z))
minVI_cluster_sim_study1_covsel <- minVI(psm_mat_sim1_covsel, method = 'greedy', suppress.comment = FALSE, cls.draw = t(classification_mat_sim1_covsel), start.cl = max.col(sim1_LCR_fit_covsel$Z))
minVI_cluster_sim_study1_itemsel <- minVI(psm_mat_sim1_itemsel, method = 'greedy', suppress.comment = FALSE, cls.draw = t(classification_mat_sim1_itemsel), start.cl = max.col(sim1_LCR_fit_itemsel$Z))

minVI_cred_ball_sim_study1_no_varsel <- credibleball(minVI_cluster_sim_study1_no_varsel$cl, cls.draw = t(classification_mat_sim1_no_varsel), c.dist = 'VI', alpha = 0.05)
minVI_cred_ball_sim_study1_varsel <- credibleball(minVI_cluster_sim_study1_varsel$cl, cls.draw = t(classification_mat_sim1_varsel), c.dist = 'VI', alpha = 0.05)
minVI_cred_ball_sim_study1_covsel <- credibleball(minVI_cluster_sim_study1_covsel$cl, cls.draw = t(classification_mat_sim1_covsel), c.dist = 'VI', alpha = 0.05)
minVI_cred_ball_sim_study1_itemsel <- credibleball(minVI_cluster_sim_study1_itemsel$cl, cls.draw = t(classification_mat_sim1_itemsel), c.dist = 'VI', alpha = 0.05)





#Getting an estimate of the value of beta coefficients for the predictor and simultaneous selection models (by excluding any samples that have zeros in a given row

X1_inclusion_indices_covsel <- which(sim1_LCR_fit_covsel$samples$gamma_samples[2,] == 1)
X2_inclusion_indices_covsel <- which(sim1_LCR_fit_covsel$samples$gamma_samples[3,] == 1)
X3_inclusion_indices_covsel <- which(sim1_LCR_fit_covsel$samples$gamma_samples[4,] == 1)
X4_inclusion_indices_covsel <- which(sim1_LCR_fit_covsel$samples$gamma_samples[5,] == 1)
X5_inclusion_indices_covsel <- which(sim1_LCR_fit_covsel$samples$gamma_samples[6,] == 1)
X6_inclusion_indices_covsel <- which(sim1_LCR_fit_covsel$samples$gamma_samples[7,] == 1)

sim1_beta1_samples_nonzero_covsel <- sim1_LCR_fit_covsel$samples$beta_samples[2,1,X1_inclusion_indices_covsel] 
sim1_beta2_samples_nonzero_covsel <- sim1_LCR_fit_covsel$samples$beta_samples[3,1,X2_inclusion_indices_covsel] 
sim1_beta3_samples_nonzero_covsel <- sim1_LCR_fit_covsel$samples$beta_samples[4,1,X3_inclusion_indices_covsel] 
sim1_beta4_samples_nonzero_covsel <- sim1_LCR_fit_covsel$samples$beta_samples[5,1,X4_inclusion_indices_covsel] 
sim1_beta5_samples_nonzero_covsel <- sim1_LCR_fit_covsel$samples$beta_samples[6,1,X5_inclusion_indices_covsel] 
sim1_beta6_samples_nonzero_covsel <- sim1_LCR_fit_covsel$samples$beta_samples[7,1,X6_inclusion_indices_covsel] 

sim1_beta1_covsel_estimate <- mean(sim1_beta1_samples_nonzero_covsel)
sim1_beta2_covsel_estimate <- mean(sim1_beta2_samples_nonzero_covsel)
sim1_beta3_covsel_estimate <- mean(sim1_beta3_samples_nonzero_covsel)
sim1_beta4_covsel_estimate <- mean(sim1_beta4_samples_nonzero_covsel)
sim1_beta5_covsel_estimate <- mean(sim1_beta5_samples_nonzero_covsel)
sim1_beta6_covsel_estimate <- mean(sim1_beta6_samples_nonzero_covsel)

beta_covsel_estimate_sim1 <- cbind(c(sim1_LCR_fit_covsel$beta_estimate[1,1], sim1_beta1_covsel_estimate, sim1_beta2_covsel_estimate, sim1_beta3_covsel_estimate, sim1_beta4_covsel_estimate, sim1_beta5_covsel_estimate, sim1_beta6_covsel_estimate),0)
#Should only consider the predictors that had enough posterior probability



X1_inclusion_indices_varsel <- which(sim1_LCR_fit_varsel$samples$gamma_samples[2,] == 1)
X2_inclusion_indices_varsel <- which(sim1_LCR_fit_varsel$samples$gamma_samples[3,] == 1)
X3_inclusion_indices_varsel <- which(sim1_LCR_fit_varsel$samples$gamma_samples[4,] == 1)
X4_inclusion_indices_varsel <- which(sim1_LCR_fit_varsel$samples$gamma_samples[5,] == 1)
X5_inclusion_indices_varsel <- which(sim1_LCR_fit_varsel$samples$gamma_samples[6,] == 1)
X6_inclusion_indices_varsel <- which(sim1_LCR_fit_varsel$samples$gamma_samples[7,] == 1)

sim1_beta1_samples_nonzero_varsel <- sim1_LCR_fit_varsel$samples$beta_samples[2,1,X1_inclusion_indices_covsel] 
sim1_beta2_samples_nonzero_varsel <- sim1_LCR_fit_varsel$samples$beta_samples[3,1,X2_inclusion_indices_covsel] 
sim1_beta3_samples_nonzero_varsel <- sim1_LCR_fit_varsel$samples$beta_samples[4,1,X3_inclusion_indices_covsel] 
sim1_beta4_samples_nonzero_varsel <- sim1_LCR_fit_varsel$samples$beta_samples[5,1,X4_inclusion_indices_covsel] 
sim1_beta5_samples_nonzero_varsel <- sim1_LCR_fit_varsel$samples$beta_samples[6,1,X5_inclusion_indices_covsel] 
sim1_beta6_samples_nonzero_varsel <- sim1_LCR_fit_varsel$samples$beta_samples[7,1,X6_inclusion_indices_covsel] 

sim1_beta1_varsel_estimate <- mean(sim1_beta1_samples_nonzero_varsel)
sim1_beta2_varsel_estimate <- mean(sim1_beta2_samples_nonzero_varsel)
sim1_beta3_varsel_estimate <- mean(sim1_beta3_samples_nonzero_varsel)
sim1_beta4_varsel_estimate <- mean(sim1_beta4_samples_nonzero_varsel)
sim1_beta5_varsel_estimate <- mean(sim1_beta5_samples_nonzero_varsel)
sim1_beta6_varsel_estimate <- mean(sim1_beta6_samples_nonzero_varsel)

beta_varsel_estimate_sim1 <- cbind(c(sim1_LCR_fit_varsel$beta_estimate[1,1], sim1_beta1_varsel_estimate, sim1_beta2_varsel_estimate, sim1_beta3_varsel_estimate, sim1_beta4_varsel_estimate, sim1_beta5_varsel_estimate, sim1_beta6_varsel_estimate),0)












#Simulation Study 2

# 3 groups
# 13 item variables, 8 of which are useful for clustering 
# variables 1,2,3 have 2 possible responses, 4,5,6 have 3 possible responses, and 7 and 8 have 4 possible responses
# variables 9,10 have 3 possible responses, 11,12,13 have 5 possible responses
# 6 covariates, of varying effect size for different groups
# 500 observations

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


png("./sim_study_plots/sim_study2_plots/sim2_itemsel_log_post.png", width = 800, height = 600)
plot(sim2_LCR_fit_itemsel1$samples$log_post_samples, type = 'l', main = 'Simulated Data 2 (Item Selection Model)', xlab = 'Iteration Number', ylab = 'log posterior')
dev.off()

png("./sim_study_plots/sim_study2_plots/sim2_itemsel_log_post_acf.png", width = 800, height = 600)
acf(sim2_LCR_fit_itemsel1$samples$log_post_samples, main = 'Simulated Data 2 (Item Selection Model)')
dev.off()



#Checking variable importance using KL divergence - need to look at this across all of the datasets I think


assignment_probs_full_sim2 <- sim2_data1$class_prob


assignment_probs_X1_sim2_mu <- cbind(1,sim2_data1$X[,-1])%*%sim2_beta[-2,]
assignment_probs_X2_sim2_mu <- cbind(1,sim2_data1$X[,-2])%*%sim2_beta[-3,]
assignment_probs_X3_sim2_mu <- cbind(1,sim2_data1$X[,-3])%*%sim2_beta[-4,]
assignment_probs_X4_sim2_mu <- cbind(1,sim2_data1$X[,-4])%*%sim2_beta[-5,]
assignment_probs_X5_sim2_mu <- cbind(1,sim2_data1$X[,-5])%*%sim2_beta[-6,]
assignment_probs_X6_sim2_mu <- cbind(1,sim2_data1$X[,-6])%*%sim2_beta[-7,]

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

ggsave('sim_study_plots/sim_study2_plots/KL_divergence_plots_sim_study2/KL_diverg_sim_study2.png', dpi = 600, width = 7, height = 5)




#creating a forest plot to assess HDI and estimate quality.

relevelled_beta_samples_sim2_itemsel1 <- sim2_LCR_fit_itemsel1$samples$beta_samples
for (i in 1:dim(relevelled_beta_samples_sim2_itemsel1)[3]){
  relevelled_beta_samples_sim2_itemsel1[,,i] <- relevelled_beta_samples_sim2_itemsel1[,c(3,1,2),i]
  relevelled_beta_samples_sim2_itemsel1[,,i] <- relevelled_beta_samples_sim2_itemsel1[,,i] - relevelled_beta_samples_sim2_itemsel1[,,i][,3]
}

sim2_beta_names <- paste0("beta", 0:6)
true_beta_vals_sim2_group1 <- sim2_beta[,1]
true_beta_vals_sim2_group2 <- sim2_beta[,2]

beta_means_sim2_group1 <- apply(relevelled_beta_samples_sim2_itemsel1[,1,], 1, mean)
beta_means_sim2_group2 <- apply(relevelled_beta_samples_sim2_itemsel1[,2,], 1, mean)

beta_sd_sim2_group1 <- apply(relevelled_beta_samples_sim2_itemsel1[,1,], 1, sd)
beta_sd_sim2_group2 <- apply(relevelled_beta_samples_sim2_itemsel1[,2,], 1, sd)

hdi_list_sim2_group1 <- apply(relevelled_beta_samples_sim2_itemsel1[,1,], 1, function(x) hdi(x, ci = 0.95))
hdi_list_sim2_group2 <- apply(relevelled_beta_samples_sim2_itemsel1[,2,], 1, function(x) hdi(x, ci = 0.95))

hdi_lower_sim2_group1 <- sapply(hdi_list_sim2_group1, function(x) x$CI_low)
hdi_lower_sim2_group2 <- sapply(hdi_list_sim2_group2, function(x) x$CI_low)

hdi_upper_sim2_group1 <- sapply(hdi_list_sim2_group1, function(x) x$CI_high)
hdi_upper_sim2_group2 <- sapply(hdi_list_sim2_group2, function(x) x$CI_high)

plot_data_group1 <- data.frame(
  Parameter = sim2_beta_names,
  Mean = beta_means_sim2_group1,
  SD = beta_sd_sim2_group1,
  HDI_low = hdi_lower_sim2_group1,
  HDI_high = hdi_upper_sim2_group1,
  TrueValue = true_beta_vals_sim2_group1
)

plot_data_group2 <- data.frame(
  Parameter = sim2_beta_names,
  Mean = beta_means_sim2_group2,
  SD = beta_sd_sim2_group2,
  HDI_low = hdi_lower_sim2_group2,
  HDI_high = hdi_upper_sim2_group2,
  TrueValue = true_beta_vals_sim2_group2
)


plot_data_group1$Parameter <- factor(plot_data_group1$Parameter, levels = rev(sim2_beta_names))

latex_labels <- c(
  TeX("$\\beta_6$"),
  TeX("$\\beta_5$"),
  TeX("$\\beta_4$"),
  TeX("$\\beta_3$"),
  TeX("$\\beta_2$"),
  TeX("$\\beta_1$"),
  TeX("$\\beta_0$")
)

cb_blue <- "#4477AA"
cb_orange <- "#EE7733"

ggplot(plot_data_group1, aes(y = Parameter, x = Mean)) +
  geom_point(size = 3, color = cb_blue) +
  geom_errorbarh(aes(xmin = HDI_low, xmax = HDI_high), height = 0.2, color = cb_blue) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(aes(x = TrueValue), shape = 18, size = 3, color = cb_orange) +
  geom_text(aes(x = TrueValue, label = round(TrueValue, 2)), 
            color = cb_orange, vjust = -0.7, hjust = 0.5, size = 3) +
  scale_y_discrete(labels = latex_labels) +
  labs(x = "Logit Coefficient", y = "", 
       title = "Forest Plot of Logit Coefficients with 95% HDI",
       subtitle = "Blue: Posterior mean & 95% HDI, Orange: True value") +
  theme_minimal()

ggsave('./sim_study_plots/sim_study2_plots/beta_plots/sim2_group1_forest_plot.png', height = 8, width = 8, dpi = 600)

ggplot(plot_data_group2, aes(y = Parameter, x = Mean)) +
  geom_point(size = 3, color = cb_blue) +
  geom_errorbarh(aes(xmin = HDI_low, xmax = HDI_high), height = 0.2, color = cb_blue) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(aes(x = TrueValue), shape = 18, size = 3, color = cb_orange) +
  geom_text(aes(x = TrueValue, label = round(TrueValue, 2)), 
            color = cb_orange, vjust = -0.7, hjust = 0.5, size = 3) +
  scale_y_discrete(labels = latex_labels) +
  labs(x = "Logit Coefficient", y = "", 
       title = "Forest Plot of Logit Coefficients with 95% HDI",
       subtitle = "Blue: Posterior mean & 95% HDI, Orange: True value") +
  theme_minimal()

ggsave('./sim_study_plots/sim_study2_plots/beta_plots/sim2_group2_forest_plot.png', height = 8, width = 8, dpi = 600)







#creating ridgeline plots for simulation 2


#group 1

posterior_samples_group1 <- relevelled_beta_samples_sim2_itemsel1[,1,]


sim2_beta_names <- paste0("beta", 0:6)
rownames(posterior_samples_group1) <- sim2_beta_names


long_posterior_data_group1 <- as.data.frame(t(posterior_samples_group1)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "Value"
  )


long_posterior_data_group1$Parameter <- factor(long_posterior_data_group1$Parameter, levels = rev(sim2_beta_names))


true_beta_vals_sim2_group1 <- sim2_beta[,1]
true_values_df_group1 <- data.frame(
  Parameter = sim2_beta_names,
  TrueValue = true_beta_vals_sim2_group1
)


cb_blue <- "#4477AA"
cb_orange <- "#EE7733"
math_labels <- expression(
  beta[6], beta[5], beta[4], beta[3], beta[2], beta[1], beta[0]
)


ridgeline_plot_sim2_group1 <- ggplot(long_posterior_data_group1, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  geom_point(
    data = true_values_df_group1,
    aes(x = TrueValue, y = Parameter),
    color = cb_orange,
    shape = 18,
    size = 3
  ) +
  geom_text(
    data = true_values_df_group1,
    aes(x = TrueValue, y = Parameter, label = round(TrueValue, 2)),
    color = cb_orange,
    vjust = -1.5,
    size = 3
  ) +
  scale_y_discrete(labels = math_labels) +
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions (Group 1)",
    subtitle = "Posterior density with median & 95% credible interval. Orange: true value."
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  "./sim_study_plots/sim_study2_plots/beta_plots/ridgeline_plot_sim2_group1.png", 
  plot = ridgeline_plot_sim2_group1, 
  width = 10, 
  height = 6, 
  dpi = 600 
)


#group 2
 

posterior_samples_group2 <- relevelled_beta_samples_sim2_itemsel1[,2,]


sim2_beta_names <- paste0("beta", 0:6)
rownames(posterior_samples_group2) <- sim2_beta_names


long_posterior_data_group2 <- as.data.frame(t(posterior_samples_group2)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "Value"
  )


long_posterior_data_group2$Parameter <- factor(long_posterior_data_group2$Parameter, levels = rev(sim2_beta_names))


true_beta_vals_sim2_group2 <- sim2_beta[,2]
true_values_df_group2 <- data.frame(
  Parameter = sim2_beta_names,
  TrueValue = true_beta_vals_sim2_group2
)


cb_blue <- "#4477AA"
cb_orange <- "#EE7733"
math_labels <- expression(
  beta[6], beta[5], beta[4], beta[3], beta[2], beta[1], beta[0]
)


ridgeline_plot_sim2_group2 <- ggplot(long_posterior_data_group2, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  geom_point(
    data = true_values_df_group2,
    aes(x = TrueValue, y = Parameter),
    color = cb_orange,
    shape = 18,
    size = 3
  ) +
  geom_text(
    data = true_values_df_group2,
    aes(x = TrueValue, y = Parameter, label = round(TrueValue, 2)),
    color = cb_orange,
    vjust = -1.5,
    size = 3
  ) +
  scale_y_discrete(labels = math_labels) +
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions (Group 2)",
    subtitle = "Posterior density with median & 95% credible interval. Orange: true value."
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  "./sim_study_plots/sim_study2_plots/beta_plots/ridgeline_plot_sim2_group2.png", 
  plot = ridgeline_plot_sim2_group2, 
  width = 10, 
  height = 6, 
  dpi = 600 
)


classification_mat_sim2_no_varsel <- apply(sim2_LCR_fit_no_varsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
})
classification_mat_sim2_varsel <- apply(sim2_LCR_fit_varsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
})
classification_mat_sim2_covsel <- apply(sim2_LCR_fit_covsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
})
classification_mat_sim2_itemsel <- apply(sim2_LCR_fit_itemsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
})
psm_mat_sim2_no_varsel <- comp.psm(t(classification_mat_sim2_no_varsel))
psm_mat_sim2_varsel <- comp.psm(t(classification_mat_sim2_varsel))
psm_mat_sim2_covsel <- comp.psm(t(classification_mat_sim2_covsel))
psm_mat_sim2_itemsel <- comp.psm(t(classification_mat_sim2_itemsel))

minVI_cluster_sim_study2_no_varsel <- minVI(psm_mat_sim2_no_varsel, method = 'draws', suppress.comment = FALSE, cls.draw = t(classification_mat_sim2_no_varsel), start.cl = max.col(sim2_LCR_fit_no_varsel$Z))
minVI_cluster_sim_study2_varsel <- minVI(psm_mat_sim2_varsel, method = 'draws', suppress.comment = FALSE, cls.draw = t(classification_mat_sim2_varsel), start.cl = max.col(sim2_LCR_fit_varsel$Z))
minVI_cluster_sim_study2_covsel <- minVI(psm_mat_sim2_covsel, method = 'draws', suppress.comment = FALSE, cls.draw = t(classification_mat_sim2_covsel), start.cl = max.col(sim2_LCR_fit_covsel$Z))
minVI_cluster_sim_study2_itemsel <- minVI(psm_mat_sim2_itemsel, method = 'draws', suppress.comment = FALSE, cls.draw = t(classification_mat_sim2_itemsel), start.cl = max.col(sim2_LCR_fit_itemsel$Z))


minVI_cred_ball_sim_study2_no_varsel <- credibleball(minVI_cluster_sim_study2_no_varsel$cl, cls.draw = t(classification_mat_sim2_no_varsel), c.dist = 'VI', alpha = 0.05)
minVI_cred_ball_sim_study2_varsel <- credibleball(minVI_cluster_sim_study2_varsel$cl, cls.draw = t(classification_mat_sim2_varsel), c.dist = 'VI', alpha = 0.05)
minVI_cred_ball_sim_study2_covsel <- credibleball(minVI_cluster_sim_study2_covsel$cl, cls.draw = t(classification_mat_sim2_covsel), c.dist = 'VI', alpha = 0.05)
minVI_cred_ball_sim_study2_itemsel <- credibleball(minVI_cluster_sim_study2_itemsel$cl, cls.draw = t(classification_mat_sim2_itemsel), c.dist = 'VI', alpha = 0.05)






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










#simulating one more dataset, similar to dataset 2 but with a different beta value
# 3 groups
# 13 item variables, 8 of which are useful for clustering 
# variables 1,2,3 have 2 possible responses, 4,5,6 have 3 possible responses, and 7 and 8 have 4 possible responses
# variables 9,10 have 3 possible responses, 11,12,13 have 5 possible responses
# 4 covariates, of varying effect size for different groups
# 500 observations

# structure parameters
G <- 3
p <- 4
M <- 13
K <- c(2,2,2,3,3,3,4,4,3,3,5,5,5)

# covariate effects matrix
sim3_beta <- matrix(c(0.5, -0.5, 1, -1, -0.5, 0.5, 0,0,0,0), nrow = p+1, byrow = TRUE)
#sim2_beta <- matrix(c(-1, 0.8, 1.1, 1.2, 0.3, 0, 0, -0.3, 0.5, 0, 0.4, 0, 0, 0), nrow = p+1)


sim3_theta1 <- matrix(c(0.15, 0.6, 0.8, 0.85, 0.4, 0.2), nrow = G, ncol = K[1])
sim3_theta2 <- matrix(c(0.25, 0.45, 0.7, 0.75, 0.55, 0.3), nrow = G, ncol = K[2])
sim3_theta3 <- matrix(c(0.7, 0.2, 0.65, 0.3, 0.8, 0.35), nrow = G, ncol = K[3])
sim3_theta4 <- matrix(c(0.1, 0.35, 0.7, 0.25, 0.4, 0.2, 0.65, 0.25, 0.1), nrow = G, ncol = K[4])
sim3_theta5 <- matrix(c(0.2, 0.25, 0.65, 0.15, 0.6, 0.25, 0.65, 0.15, 0.1), nrow = G, ncol = K[5])
sim3_theta6 <- matrix(c(0.15, 0.5, 0.75, 0.2, 0.35, 0.15, 0.65, 0.15, 0.1), nrow = G, ncol = K[6])
sim3_theta7 <- matrix(c(0.1, 0.25, 0.6, 0.15, 0.35, 0.25, 0.25, 0.25, 0.1, 0.5, 0.15, 0.05), nrow = G, ncol = K[7])
sim3_theta8 <- matrix(c(0.15, 0.2, 0.55, 0.2, 0.45, 0.2, 0.2, 0.25, 0.15, 0.45, 0.1, 0.1), nrow = G, ncol = K[8])
sim3_theta9 <- matrix(c(0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1), nrow = G, ncol = K[9])
sim3_theta10 <- matrix(c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2), nrow = G, ncol = K[10])
sim3_theta11 <- matrix(1/K[11], nrow = G, ncol = K[11])
sim3_theta12 <- matrix(c(0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.2, 0.2, 0.2, 0.25, 0.25, 0.25, 0.3, 0.3, 0.3), nrow = G, ncol = K[12])
sim3_theta13 <- matrix(c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = G, ncol = K[13])


sim3_theta <- list(sim3_theta1, sim3_theta2, sim3_theta3, sim3_theta4, sim3_theta5, sim3_theta6, sim3_theta7, sim3_theta8, sim3_theta9, sim3_theta10, sim3_theta11, sim3_theta12, sim3_theta13)


sim3_data1 <- LCR_sim_data(theta = sim3_theta, beta = sim3_beta, n_samples = 500)


sim3_data_frame_binded1 <- data.frame(Y = sim3_data1$Y, X = sim3_data1$X)


write.csv(data.frame(cbind(sim2_data1$X, sim2_data1$Y)), "./sim_study_plots/sim_study3_plots/sim3_data1.csv", row.names = FALSE)


sim3_prelim_collapsed_fit1 <- blca.collapsed(X = sim3_data1$Y, G = 1, iter = 50000, burn.in = 1000, thin = 0.1, G.sel = TRUE, verbose = TRUE)


#Running with no variable selection

sim3_LCR_fit_no_varsel1 <- LCR_Gibbs(X = sim3_data1$X, Y = sim3_data1$Y, G = 3, beta_prior_cov = diag(10^2,p+1), beta_prior_mean = rep(0,p+1), 
                                     theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                                     relabel = TRUE, n_samples = 5000, burnin = 10, thinby = 10)



#Running with full variable selection


sim3_LCR_fit_varsel1 <- LCR_Gibbs(X = sim3_data1$X, Y = sim3_data1$Y, G = 3, beta_prior_cov = diag(10^2,p+1), beta_prior_mean = rep(0,p+1), 
                                  theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                  relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


#Running with only predictor selection

sim3_LCR_fit_covsel1 <- LCR_Gibbs(X = sim3_data1$X, Y = sim3_data1$Y, G = 3, beta_prior_cov = diag(10^2,p+1), beta_prior_mean = rep(0,p+1), 
                                  theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                  relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)



#Running with only item selection


sim3_LCR_fit_itemsel1 <- LCR_Gibbs(X = sim3_data1$X, Y = sim3_data1$Y, G = 3, beta_prior_cov = diag(10^2,p+1), beta_prior_mean = rep(0,p+1), 
                                   theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                   relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)




































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


