###-------Bayesian Latent Class Regression Script-------###

#Getting functions, packages
working_dir <- getwd()
source('LCR_Gibbs.R')
source('LCR_sim_data.R')
source('beta_summary_table.R')
source('LCR_posterior_membership_prob.R')
library(ggplot2)
library(RColorBrewer)
library(reshape2) 
library(BayesLCA)
library(bayestestR)
library(mcclust.ext)
library(dplyr)
library(tidyr)
library(ggridges)
library(foreign)
library(fossil)
cb_blue <- "#4477AA"
cb_orange <- "#EE7733"


##-------Simulation Studies-------##

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

set.seed(123)

#Simulatiing dataset from using the above parameters
sim1_data <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 500)
sim1_truelabels <- sim1_data$class

#Computing the ARI for the clustering assigned using the true parameters vs the true clustering
true_posterior_cluster_probs_sim1 <- LCR_posterior_membership_prob(theta = sim1_theta, beta = cbind(sim1_beta, 0), X = sim1_data$X, Y = sim1_data$Y)
cat("Adjusted Rand Index for true model parameters:", adj.rand.index(sim1_truelabels, max.col(true_posterior_cluster_probs_sim1)), "\n")



set.seed(124)

#Preliminary collapsed fit of the sampler for simulation 1

sim1_prelim_collapsed_fit <- blca.collapsed(X = sim1_data$Y, G = 1, iter = 50000, burn.in = 1000, thin = 0.1, G.sel = TRUE, var.sel = TRUE, verbose = TRUE)


set.seed(125)

#LCR fit (no variable selection) for simulation 1

sim1_LCR_fit_no_varsel <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                    theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                                    relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

#Creating summary table for beta coefficients (up to permutation of group labels)
sim1_beta_summary_table_no_varsel <- multinomial_logit_coefficient_summary_table_MCMC(sim1_LCR_fit_no_varsel$samples$beta_samples)
print(sim1_beta_summary_table_no_varsel)

#Creating a ridgeline plot for beta coefficients

sim1_beta_samples <- sim1_LCR_fit_no_varsel$samples$beta_samples[,2,]
sim1_beta_names <- paste0("beta", 0:6)
rownames(sim1_beta_samples) <- sim1_beta_names
long_beta_samples_data_sim1 <- as.data.frame(t(sim1_beta_samples)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )
long_beta_samples_data_sim1$Parameter <- factor(long_beta_samples_data_sim1$Parameter, levels = rev(sim1_beta_names))
true_beta_vals_sim1 <- sim1_beta[,1]
true_values_df_sim1 <- data.frame(
  Parameter = sim1_beta_names,
  TrueValue = true_beta_vals_sim1
)
math_labels_sim1 <- expression(
  beta[6], beta[5], beta[4], beta[3], beta[2], beta[1], beta[0]
)
final_ridgeline_plot_sim1 <- ggplot(long_beta_samples_data_sim1, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  geom_point(
    data = true_values_df_sim1,
    aes(x = TrueValue, y = Parameter),
    color = cb_orange,
    shape = 18,
    size = 3
  ) +
  geom_text(
    data = true_values_df_sim1,
    aes(x = TrueValue, y = Parameter, label = round(TrueValue, 2)),
    color = cb_orange,
    vjust = 2,
    size = 3
  ) +
  scale_y_discrete(labels = math_labels_sim1) + 
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

# ggsave(
#   "./sim_study_plots/sim_study1_plots/beta_plots/ridgeline_plot_sim1.png", 
#   plot = final_plot, 
#   width = 10, 
#   height = 6, 
#   dpi = 600 
# )

#Assigning observations to clusters for comparison, by aiming to minimise the posterior expected variation of information
sim1_cluster_mat_no_var_sel <- apply(sim1_LCR_fit_no_varsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
}) 
sim1_psm_mat_no_varsel <- comp.psm(t(sim1_cluster_mat_no_var_sel))
sim1_cluster_minVI_no_var_sel <- minVI(psm = sim1_psm_mat_no_varsel, method = 'greedy', start.cl = max.col(sim1_LCR_fit_no_varsel$Z))

#Cross-classification table and ARI
cat("Cross-classification table for true labels (rows) vs labels assigned from LCR model (columns)\n", 
    paste(capture.output(print(table(sim1_truelabels, sim1_cluster_minVI_no_var_sel$cl))), collapse = "\n"), "\n")
cat("Adjusted Rand Index for LCR model:", adj.rand.index(sim1_truelabels, sim1_cluster_minVI_no_var_sel$cl), "\n")

#Diagnostic plots

plot(sim1_LCR_fit_no_varsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 1', xlab = 'Iteration Number', ylab = 'log posterior')
acf(sim1_LCR_fit_no_varsel$samples$log_post_samples, main = 'Simulated Data 1')






set.seed(126)

#LCR fit with item variable selection for simulation 1

sim1_LCR_fit_itemsel <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                  theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                  relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

#Posterior probability of inclusion for each item variable
sim1_item_inclusion_probs_itemsel <- sim1_LCR_fit_itemsel$item_inclusion_prob
names(sim1_item_inclusion_probs_itemsel) <- paste0('Y', 1:length(sim1_item_inclusion_probs_itemsel))

cat("Posterior inclusion probability for item variables \n",
    paste(names(sim1_item_inclusion_probs_itemsel), collapse = "    "), "\n",
    paste(sprintf("%.4f", sim1_item_inclusion_probs_itemsel), collapse = "  "), "\n")

sim1_cluster_mat_item_sel <- apply(sim1_LCR_fit_itemsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
}) 
sim1_psm_mat_item_sel <- comp.psm(t(sim1_cluster_mat_item_sel))
sim1_cluster_minVI_item_sel <- minVI(psm = sim1_psm_mat_item_sel, method = 'greedy', start.cl = max.col(sim1_LCR_fit_itemsel$Z))

#Cross-classification table and ARI
cat("Cross-classification table for true labels (rows) vs labels assigned from LCR model with item selection (columns)\n", 
    paste(capture.output(print(table(sim1_truelabels, sim1_cluster_minVI_item_sel$cl))), collapse = "\n"), "\n")
cat("Adjusted Rand Index for LCR model with item selection:", adj.rand.index(sim1_truelabels, sim1_cluster_minVI_item_sel$cl), "\n")

#Diagnostic plots

plot(sim1_LCR_fit_itemsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 1 (item selection)', xlab = 'Iteration Number', ylab = 'log posterior')
acf(sim1_LCR_fit_itemsel$samples$log_post_samples, main = 'Simulated Data 1 (item selection)')




set.seed(127)

#LCR fit with predictor variable selection for simulation 1

sim1_LCR_fit_covsel <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

#posterior probability of inclusion for each predictor variable

sim1_predictor_inclusion_probs_covsel <- sim1_LCR_fit_covsel$cov_inclusion_prob
names(sim1_predictor_inclusion_probs_covsel) <- paste0('X', 1:length(sim1_predictor_inclusion_probs_covsel))

cat("Posterior inclusion probability for predictor variables \n",
    paste(names(sim1_predictor_inclusion_probs_covsel), collapse = "    "), "\n",
    paste(sprintf("%.4f", sim1_predictor_inclusion_probs_covsel), collapse = "  "), "\n")

plot(sim1_LCR_fit_covsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 1 (predictor selection)', xlab = 'Iteration Number', ylab = 'log posterior')
acf(sim1_LCR_fit_covsel$samples$log_post_samples, main = 'Simulated Data 1 (predictor selection)')




set.seed(128)

#LCR fit with simultaneous variable selection for simulation 1

sim1_LCR_fit_varsel <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


#Posterior probability of inclusion for each item variable
sim1_item_inclusion_probs_varsel <- sim1_LCR_fit_varsel$item_inclusion_prob
names(sim1_item_inclusion_probs_varsel) <- paste0('Y', 1:length(sim1_item_inclusion_probs_varsel))

cat("Posterior inclusion probability for item variables \n",
    paste(names(sim1_item_inclusion_probs_varsel), collapse = "    "), "\n",
    paste(sprintf("%.4f", sim1_item_inclusion_probs_varsel), collapse = "  "), "\n")


#posterior probability of inclusion for each predictor variable

sim1_predictor_inclusion_probs_varsel <- sim1_LCR_fit_varsel$cov_inclusion_prob
names(sim1_predictor_inclusion_probs_varsel) <- paste0('X', 1:length(sim1_predictor_inclusion_probs_varsel))

cat("Posterior inclusion probability for predictor variables \n",
    paste(names(sim1_predictor_inclusion_probs_varsel), collapse = "    "), "\n",
    paste(sprintf("%.4f", sim1_predictor_inclusion_probs_varsel), collapse = "  "), "\n")

sim1_cluster_mat_varsel <- apply(sim1_LCR_fit_varsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
}) 
sim1_psm_mat_varsel <- comp.psm(t(sim1_cluster_mat_varsel))
sim1_cluster_minVI_varsel <- minVI(psm = sim1_psm_mat_varsel, method = 'greedy', start.cl = max.col(sim1_LCR_fit_varsel$Z))


#Cross classification table and ARI for true labels vs full variable selection clustering
cat("Cross-classification table for true labels (rows) vs labels assigned from LCR model with simultaneous variable selection (columns)\n", 
    paste(capture.output(print(table(sim1_truelabels, sim1_cluster_minVI_varsel$cl))), collapse = "\n"), "\n")
cat("Adjusted Rand Index for true labels vs. LCR model with simultaneous variable selection", adj.rand.index(sim1_truelabels, sim1_cluster_minVI_item_sel$cl), "\n")


#Cross classification table and ARI for full LCR model vs variable selection clustering

cat("Cross-classification table for labels from full LCR model (rows) vs labels assigned from LCR model with simultaneous variable selection (columns)\n", 
    paste(capture.output(print(table(sim1_cluster_minVI_no_var_sel$cl, sim1_cluster_minVI_varsel$cl))), collapse = "\n"), "\n")
cat("Adjusted Rand Index for LCR model (no variable selection) vs. LCR model with simultaneous variable selection:", adj.rand.index(sim1_cluster_minVI_no_var_sel$cl, sim1_cluster_minVI_varsel$cl), "\n")




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

#Simulatiing dataset from using the above parameters
sim2_data <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 500)
sim2_truelabels <- sim2_data$class

#Computing the ARI for the clustering assigned using the true parameters vs the true clustering
true_posterior_cluster_probs_sim2 <- LCR_posterior_membership_prob(theta = sim2_theta, beta = cbind(sim2_beta, 0), X = sim2_data$X, Y = sim2_data$Y)
cat("Adjusted Rand Index for true model parameters:", adj.rand.index(sim2_truelabels, max.col(true_posterior_cluster_probs_sim2)), "\n")


set.seed(130)

#Preliminary collapsed LCA fit of the sampler for simulation 2

sim2_prelim_collapsed_fit <- blca.collapsed(X = sim2_data$Y, G = 1, iter = 50000, burn.in = 1000, thin = 0.1, G.sel = TRUE, var.sel = TRUE, verbose = TRUE)


set.seed(131)

#LCR fit (no variable selection) for simulation 2

sim2_LCR_fit_no_varsel <- LCR_Gibbs(X = sim2_data$X, Y = sim2_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                    theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                                    relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


#We see that the 3 group structure collapses into 2 groups when estimating without variable selection, hence an identifiability issue arises for the beta coefficients
sim2_beta_summary_table_no_varsel <- multinomial_logit_coefficient_summary_table_MCMC(sim2_LCR_fit_no_varsel$samples$beta_samples)
print(sim2_beta_summary_table_no_varsel)



#We can look at the cross-classification table to compare the true labels to the fitted clusters

sim2_cluster_mat_no_var_sel <- apply(sim2_LCR_fit_no_varsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
}) 

sim2_psm_mat_no_varsel <- comp.psm(t(sim2_cluster_mat_no_var_sel))
sim2_cluster_minVI_no_var_sel <- minVI(psm = sim2_psm_mat_no_varsel, method = 'greedy', start.cl = max.col(sim2_LCR_fit_no_varsel$Z))
cat("Cross-classification table for true labels vs labels assigned from LCR model (no variable selection) (columns)\n", 
    paste(capture.output(print(table(sim2_truelabels, sim2_cluster_minVI_no_var_sel$cl))), collapse = "\n"), "\n")


set.seed(132)


#LCR fit with item variable selection for simulation 2

sim2_LCR_fit_itemsel <- LCR_Gibbs(X = sim2_data$X, Y = sim2_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                  theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE, verbose = TRUE, 
                                  relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


#Posterior inclusion probability for each of the item variable in the item selection model
sim2_item_inclusion_probs_itemsel <- sim2_LCR_fit_itemsel$item_inclusion_prob
names(sim2_item_inclusion_probs_itemsel) <- paste0('Y', 1:length(sim2_item_inclusion_probs_itemsel))
cat("Posterior inclusion probability for item variables \n",
    paste(names(sim2_item_inclusion_probs_itemsel), collapse = "    "), "\n",
    paste(sprintf("%.4f", sim2_item_inclusion_probs_itemsel), collapse = "  "), "\n")

#Creating summary table for beta coefficients (up to permutation of group labels)

#Relevelling beta samples for consistency with the simulated values 
relevelled_beta_samples <- sim2_LCR_fit_itemsel$samples$beta_samples
for (sample in 1:dim(relevelled_beta_samples)[3]){
  relevelled_beta_samples[,,sample] <- relevelled_beta_samples[,,sample] - relevelled_beta_samples[,3,sample]
}


sim2_beta_summary_table_itemsel <- multinomial_logit_coefficient_summary_table_MCMC(relevelled_beta_samples)
print(sim2_beta_summary_table_itemsel)



#Creating a ridgeline plot for beta coefficients

sim2_beta1_samples <- relevelled_beta_samples[,2,]
sim2_beta2_samples <- relevelled_beta_samples[,1,]


sim2_beta_names <- paste0("beta", 0:6)
rownames(sim2_beta1_samples) <- sim2_beta_names
rownames(sim2_beta2_samples) <- sim2_beta_names

long_beta1_samples_data_sim2 <- as.data.frame(t(sim2_beta1_samples)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )
long_beta2_samples_data_sim2 <- as.data.frame(t(sim2_beta2_samples)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )

long_beta1_samples_data_sim2$Parameter <- factor(long_beta1_samples_data_sim2$Parameter, levels = rev(sim2_beta_names))
long_beta2_samples_data_sim2$Parameter <- factor(long_beta2_samples_data_sim2$Parameter, levels = rev(sim2_beta_names))

true_beta1_vals_sim2 <- sim2_beta[,1]
true_beta2_vals_sim2 <- sim2_beta[,2]

true_values_beta1_df_sim2 <- data.frame(
  Parameter = sim2_beta_names,
  TrueValue = true_beta1_vals_sim2
)

true_values_beta2_df_sim2 <- data.frame(
  Parameter = sim2_beta_names,
  TrueValue = true_beta2_vals_sim2
)


math_labels_sim2 <- expression(
  beta[6], beta[5], beta[4], beta[3], beta[2], beta[1], beta[0]
)

final_ridgeline_plot_beta_1_sim2 <- ggplot(long_beta1_samples_data_sim2, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  geom_point(
    data = true_values_beta1_df_sim2,
    aes(x = TrueValue, y = Parameter),
    color = cb_orange,
    shape = 18,
    size = 3
  ) +
  geom_text(
    data = true_values_beta1_df_sim2,
    aes(x = TrueValue, y = Parameter, label = round(TrueValue, 2)),
    color = cb_orange,
    vjust = 2,
    size = 3
  ) +
  scale_y_discrete(labels = math_labels_sim2) + 
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions (Group 1)",
    subtitle = "Posterior density with median & 95% HDI, Orange: True value."
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  ) + 
  coord_cartesian(xlim = c(-3, 3))

# ggsave(
#   "./sim_study_plots/sim_study2_plots/beta_plots/ridgeline_plot_sim2_group1.png", 
#   plot = ridgeline_plot_sim2_group1, 
#   width = 10, 
#   height = 6, 
#   dpi = 600 
# )

final_ridgeline_plot_beta_2_sim2 <- ggplot(long_beta2_samples_data_sim2, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  geom_point(
    data = true_values_beta2_df_sim2,
    aes(x = TrueValue, y = Parameter),
    color = cb_orange,
    shape = 18,
    size = 3
  ) +
  geom_text(
    data = true_values_beta2_df_sim2,
    aes(x = TrueValue, y = Parameter, label = round(TrueValue, 2)),
    color = cb_orange,
    vjust = 2,
    size = 3
  ) +
  scale_y_discrete(labels = math_labels_sim2) + 
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions (Group 2)",
    subtitle = "Posterior density with median & 95% HDI, Orange: True value."
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  ) +
coord_cartesian(xlim = c(-3, 3))

# ggsave(
#   "./sim_study_plots/sim_study2_plots/beta_plots/ridgeline_plot_sim2_group2.png", 
#   plot = final_ridgeline_plot_beta_2_sim2 , 
#   width = 10, 
#   height = 6, 
#   dpi = 600 
# )

#Diagnostic plots

plot(sim2_LCR_fit_itemsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 2', xlab = 'Iteration Number', ylab = 'log posterior')
acf(sim2_LCR_fit_itemsel$samples$log_post_samples, main = 'Simulated Data 2')




#Assigning observations to clusters for comparison, by aiming to minimise the posterior expected variation of information
sim2_cluster_mat_itemsel <- apply(sim2_LCR_fit_itemsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
}) 
sim2_psm_mat_itemsel <- comp.psm(t(sim2_cluster_mat_itemsel))
sim2_cluster_minVI_itemsel <- minVI(psm = sim2_psm_mat_itemsel, method = 'greedy', start.cl = max.col(sim2_LCR_fit_itemsel$Z))

#Cross-classification table and ARI
cat("Cross-classification table for true labels (rows) vs labels assigned from LCR model with item selection (columns)\n", 
    paste(capture.output(print(table(sim2_truelabels, sim2_cluster_minVI_itemsel$cl))), collapse = "\n"), "\n")
cat("Adjusted Rand Index for LCR model with item selection:", adj.rand.index(sim2_truelabels, sim2_cluster_minVI_itemsel$cl), "\n")



set.seed(133)



#LCR fit with predictor variable selection for simulation 2

sim2_LCR_fit_covsel <- LCR_Gibbs(X = sim2_data$X, Y = sim2_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

#posterior probability of inclusion for each predictor variable

sim2_predictor_inclusion_probs_covsel <- sim2_LCR_fit_covsel$cov_inclusion_prob
names(sim2_predictor_inclusion_probs_covsel) <- paste0('X', 1:length(sim2_predictor_inclusion_probs_covsel))

cat("Posterior inclusion probability for predictor variables \n",
    paste(names(sim2_predictor_inclusion_probs_covsel), collapse = "    "), "\n",
    paste(sprintf("%.4f", sim2_predictor_inclusion_probs_covsel), collapse = "  "), "\n")

#diagnostic plots
plot(sim2_LCR_fit_covsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 2 (predictor selection)', xlab = 'Iteration Number', ylab = 'log posterior')
acf(sim1_LCR_fit_covsel$samples$log_post_samples, main = 'Simulated Data 2 (predictor selection)')


#LCR fit with simultaneous variable selection for simulation 2

sim2_LCR_fit_varsel <- LCR_Gibbs(X = sim2_data$X, Y = sim2_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


#Posterior probability of inclusion for each item variable
sim2_item_inclusion_probs_varsel <- sim2_LCR_fit_varsel$item_inclusion_prob
names(sim2_item_inclusion_probs_varsel) <- paste0('Y', 1:length(sim2_item_inclusion_probs_varsel))

cat("Posterior inclusion probability for item variables \n",
    paste(names(sim2_item_inclusion_probs_varsel), collapse = "    "), "\n",
    paste(sprintf("%.4f", sim2_item_inclusion_probs_varsel), collapse = "  "), "\n")


#posterior probability of inclusion for each predictor variable

sim2_predictor_inclusion_probs_varsel <- sim2_LCR_fit_varsel$cov_inclusion_prob
names(sim2_predictor_inclusion_probs_varsel) <- paste0('X', 1:length(sim2_predictor_inclusion_probs_varsel))

cat("Posterior inclusion probability for predictor variables \n",
    paste(names(sim2_predictor_inclusion_probs_varsel), collapse = "    "), "\n",
    paste(sprintf("%.4f", sim2_predictor_inclusion_probs_varsel), collapse = "  "), "\n")

sim2_cluster_mat_varsel <- apply(sim2_LCR_fit_varsel$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
}) 
sim2_psm_mat_varsel <- comp.psm(t(sim2_cluster_mat_varsel))
sim2_cluster_minVI_varsel <- minVI(psm = sim2_psm_mat_varsel, method = 'greedy', start.cl = max.col(sim2_LCR_fit_varsel$Z))


#Cross classification table and ARI for true labels vs full variable selection clustering
cat("Cross-classification table for true labels (rows) vs labels assigned from LCR model with simultaneous variable selection (columns)\n", 
    paste(capture.output(print(table(sim2_truelabels, sim2_cluster_minVI_varsel$cl))), collapse = "\n"), "\n")
cat("Adjusted Rand Index for true labels vs. LCR model with simultaneous variable selection", adj.rand.index(sim2_truelabels, sim2_cluster_minVI_itemsel$cl), "\n")


#Cross classification table and ARI for item selection LCR model vs variable selection clustering

cat("Cross-classification table for labels from LCR model with item selection (rows) vs. labels assigned from LCR model with simultaneous variable selection (columns)\n", 
    paste(capture.output(print(table(sim2_cluster_minVI_itemsel$cl, sim2_cluster_minVI_varsel$cl))), collapse = "\n"), "\n")
cat("Adjusted Rand Index for LCR model with item selection vs. LCR model with simultaneous variable selection:", adj.rand.index(sim2_cluster_minVI_itemsel$cl, sim2_cluster_minVI_varsel$cl), "\n")







##-------CSHQ Sleep Patterns Data-------##



#CSV files for response variables and predictors
Y_CSHQ <- as.matrix(read.csv('./CSHQ_response_Y.csv'))
X_CSHQ <- as.matrix(read.csv('./CSHQ_predictor_X.csv'))[,-1]



set.seed(134)


#Model Selection - Running the collapsed sampler and analysing clustering results using the minVI method of Wade and Ghahramani


CSHQ_preliminary_collapsed_run <- blca.collapsed(X = Y_CSHQ, G = 1, iter = 50000, burn.in = 1000, thin = 1/10, G.sel = TRUE, var.sel = TRUE, post.hoc.run = TRUE, verbose = TRUE, relabel = TRUE)

CSHQ_collapsed_group_number_post_prob <- table(CSHQ_preliminary_collapsed_run$samples$G)/sum(table(CSHQ_preliminary_collapsed_run$samples$G)) 

#saving the posterior inclusion probabilities for the item variables for later comparison

LCA_item_inclusion_ind_CSHQ <- do.call(rbind,CSHQ_preliminary_collapsed_run$samples$var.ind) 

LCA_item_inclusion_prob_CSHQ <- apply(LCA_item_inclusion_ind_CSHQ, 2, mean)

#Posterior probability is highest for 6 groups

CSHQ_collapsed_cluster_mat <- do.call(rbind,CSHQ_preliminary_collapsed_run$samples$labels)

CSHQ_collapsed_psm_mat <- comp.psm(CSHQ_collapsed_cluster_mat)

CSHQ_minVI_cluster <- minVI(psm = CSHQ_collapsed_psm_mat, method = 'greedy', start.cl = max.col(CSHQ_preliminary_collapsed_run$Z), suppress.comment = FALSE)

#We can see that the optimal cluster partition found using minimum posterior expected variation of information is actually a 4 group clustering solution

CSHQ_cluster_95_cred_ball <- credibleball(c.star = CSHQ_minVI_cluster$cl, cls.draw = CSHQ_collapsed_cluster_mat, c.dist = 'VI', alpha = 0.05)

#Looking at the 95% credible ball, we can see that the horizontal bound of the ball is a (VI) distance of 1.44 from the point estimate. 
#We can look at the 4 group LCA model clustering solution and check whether it is contained within the 95% credible ball, as well as comparing the correspondence between this and the solution found from the collapsed sampler
#For item variable consistency, we run the collapsed sampler, fixing G = 4.

CSHQ_preliminary_4group_LCA <- blca.collapsed(X = Y_CSHQ, G = 4, iter = 50000, burn.in = 1000, thin = 1/10, G.sel = FALSE, var.sel = TRUE, post.hoc.run = TRUE, verbose = TRUE) 

CSHQ_4group_cluster_mat <- CSHQ_preliminary_4group_LCA$samples$labels

CSHQ_4group_psm_mat <- comp.psm(CSHQ_4group_cluster_mat$`G = 4`)

CSHQ_4group_minVI_cluster <- minVI(psm = CSHQ_4group_psm_mat, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(CSHQ_preliminary_4group_LCA$Z))



CSHQ_4group_minVI_cluster$cl

#VI distance from the 4 group cluster solution to the collapsed model cluster solution

CSHQ_dist_4group_to_point_estimate <- vi.dist(CSHQ_minVI_cluster$cl, CSHQ_4group_minVI_cluster$cl)

#So the 4 group LCA model is contained within the 95% credible ball, hence is a plausible clustering solution
#Checking the cross-classification table and ARI to compare clusterings

table_collapsed_with_4_group <- table(CSHQ_minVI_cluster$cl, CSHQ_4group_minVI_cluster$cl)
ari_collapsed_with_4_group <- adj.rand.index(CSHQ_minVI_cluster$cl, CSHQ_4group_minVI_cluster$cl)

#ARI of 0.86 indicates strong clustering correspondence


set.seed(135)

#Running the LCR model without variable selection first

CSHQ_LCR_no_varsel <- LCR_Gibbs(X = X_CSHQ, Y = Y_CSHQ, G = 4,
                                           beta_prior_cov = diag(c(rep(10^2,4),5^2,5^2)),
                                           beta_prior_mean = rep(0,6), theta_hyperparam = 1,
                                           clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE,
                                           verbose = TRUE, relabel = TRUE, n_samples = 5000,
                                           burnin = 1000, thinby = 10)

#We can see some signs of non-identifiability here, with the log posterior switching between modes
plot(CSHQ_LCR_no_varsel$samples$log_post_samples, type = 'l', main = 'CSHQ Data LCR log posterior (no variable selection)', xlab = 'Iteration Number', ylab = 'log posterior')
acf(CSHQ_LCR_no_varsel$samples$log_post_samples, main = 'CSHQ Data LCR log posterior (no variable selection)')


set.seed(136)

#LCR model with item selection only

CSHQ_LCR_itemsel <- LCR_Gibbs(X = X_CSHQ, Y = Y_CSHQ, G = 4,
                              beta_prior_cov = diag(c(rep(10^2,4),5^2,5^2)),
                              beta_prior_mean = rep(0,6), theta_hyperparam = 1,
                              clust_var_prior = 0.5, item.sel = TRUE, cov.sel = FALSE,
                              verbose = TRUE, relabel = TRUE, n_samples = 5000,
                              burnin = 1000, thinby = 10)

#Saving the item inclusion probabilities
itemsel_item_inclusion_prob_CSHQ <- CSHQ_LCR_itemsel$item_inclusion_prob

#summary table for posterior of beta coefficient parameters
CSHQ_beta_summary_table_itemsel <- multinomial_logit_coefficient_summary_table_MCMC(CSHQ_LCR_itemsel$samples$beta_samples)

CSHQ_beta_summary_table_itemsel



#Creating ridgeline plots


CSHQ_betaB_samples <- CSHQ_LCR_itemsel$samples$beta_samples[,2,]
CSHQ_betaC_samples <- CSHQ_LCR_itemsel$samples$beta_samples[,3,]
CSHQ_betaD_samples <- CSHQ_LCR_itemsel$samples$beta_samples[,4,]


CSHQ_beta_names <- paste0("beta", 0:5)
rownames(CSHQ_betaB_samples) <- CSHQ_beta_names
rownames(CSHQ_betaC_samples) <- CSHQ_beta_names
rownames(CSHQ_betaD_samples) <- CSHQ_beta_names

long_betaB_samples_data_CSHQ <- as.data.frame(t(CSHQ_betaB_samples)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )
long_betaC_samples_data_CSHQ <- as.data.frame(t(CSHQ_betaC_samples)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )
long_betaD_samples_data_CSHQ <- as.data.frame(t(CSHQ_betaD_samples)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )

long_betaB_samples_data_CSHQ$Parameter <- factor(long_betaB_samples_data_CSHQ$Parameter, levels = rev(CSHQ_beta_names))
long_betaC_samples_data_CSHQ$Parameter <- factor(long_betaC_samples_data_CSHQ$Parameter, levels = rev(CSHQ_beta_names))
long_betaD_samples_data_CSHQ$Parameter <- factor(long_betaD_samples_data_CSHQ$Parameter, levels = rev(CSHQ_beta_names))


math_labels_CSHQ <- expression(
  beta[5], beta[4], beta[3], beta[2], beta[1], beta[0]
)

final_ridgeline_plot_betaB_CSHQ <- ggplot(long_betaB_samples_data_CSHQ, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_y_discrete(labels = math_labels_CSHQ) + 
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions (Group B)",
    subtitle = "Posterior density with median & 95% HDI"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  ) + 
  coord_cartesian(xlim = c(-5, 5))


final_ridgeline_plot_betaC_CSHQ <- ggplot(long_betaC_samples_data_CSHQ, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_y_discrete(labels = math_labels_CSHQ) + 
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions (Group C)",
    subtitle = "Posterior density with median & 95% HDI"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  ) + 
  coord_cartesian(xlim = c(-5, 5))

final_ridgeline_plot_betaD_CSHQ <- ggplot(long_betaD_samples_data_CSHQ, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_y_discrete(labels = math_labels_CSHQ) + 
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions (Group D)",
    subtitle = "Posterior density with median & 95% HDI"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  ) + 
  coord_cartesian(xlim = c(-5, 5))
# ggsave(
#   "./sim_study_plots/sim_study2_plots/beta_plots/ridgeline_plot_sim2_group1.png", 
#   plot = ridgeline_plot_sim2_group1, 
#   width = 10, 
#   height = 6, 
#   dpi = 600 
# )

final_ridgeline_plot_beta_2_sim2 <- ggplot(long_beta2_samples_data_sim2, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  geom_point(
    data = true_values_beta2_df_sim2,
    aes(x = TrueValue, y = Parameter),
    color = cb_orange,
    shape = 18,
    size = 3
  ) +
  geom_text(
    data = true_values_beta2_df_sim2,
    aes(x = TrueValue, y = Parameter, label = round(TrueValue, 2)),
    color = cb_orange,
    vjust = 2,
    size = 3
  ) +
  scale_y_discrete(labels = math_labels_sim2) + 
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions (Group 2)",
    subtitle = "Posterior density with median & 95% HDI, Orange: True value."
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  ) +
  coord_cartesian(xlim = c(-3, 3))


final_ridgeline_plot_betaB_CSHQ
final_ridgeline_plot_betaC_CSHQ
final_ridgeline_plot_betaD_CSHQ

#Diagnostic plots

plot(CSHQ_LCR_itemsel$samples$log_post_samples, type = 'l', main = 'CSHQ Data LCR log posterior (item selection)', xlab = 'Iteration Number', ylab = 'log posterior')
acf(CSHQ_LCR_itemsel$samples$log_post_samples, main = 'CSHQ Data LCR log posterior (item selection)')


set.seed(137)

#Running the predictor selection model

CSHQ_LCR_predsel <- LCR_Gibbs(X = X_CSHQ, Y = Y_CSHQ, G = 4,
                              beta_prior_cov = diag(c(rep(10^2,4),5^2,5^2)),
                              beta_prior_mean = rep(0,6), theta_hyperparam = 1,
                              clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE,
                              verbose = TRUE, relabel = TRUE, n_samples = 5000,
                              burnin = 1000, thinby = 10)



#Saving the predictor inclusion probabilities
predsel_predictor_inclusion_prob_CSHQ <- CSHQ_LCR_predsel$cov_inclusion_prob







set.seed(138)

#Running simultaneous variable selection

CSHQ_LCR_varsel <- LCR_Gibbs(X = X_CSHQ, Y = Y_CSHQ, G = 4,
                             beta_prior_cov = diag(c(rep(10^2,4),5^2,5^2)),
                             beta_prior_mean = rep(0,6), theta_hyperparam = 1,
                             clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE,
                             verbose = TRUE, relabel = TRUE, n_samples = 5000,
                             burnin = 1000, thinby = 10)



varsel_predictor_inclusion_prob_CSHQ <- CSHQ_LCR_varsel$cov_inclusion_prob
varsel_item_inclusion_prob_CSHQ <- CSHQ_LCR_varsel$item_inclusion_prob

CSHQ_LCR_varsel_cluster_mat <- matrix(0, nrow = dim(CSHQ_LCR_varsel$samples$z_samples)[3], ncol = nrow(Y_CSHQ))
for (slice in 1:dim(CSHQ_LCR_varsel_cluster_mat)[1]) {
  CSHQ_LCR_varsel_cluster_mat[slice, ] <- max.col(array[, , slice])
}


#Classifying observations using minimum posterior expected variation of information for the variable selection model for comparison with the item selection model

itemsel_vs_varsel_cross_classification_table <- table()



















