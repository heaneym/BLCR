###-------Bayesian Latent Class Regression Script-------###

#Getting functions, packages
working_dir <- getwd()
source('LCR_Gibbs.R')
source('LCR_sim_data.R')
source('beta_summary_table.R')
source('LCR_posterior_membership_prob.R')
source('BLCR_plot_functions/LCA_mosaic_plot.R')
source('BLCR_plot_functions/plot_combined_predictive_T.R')
source('BLCR_plot_functions/compute_expected_subscale_total_across_iterations.R')
source('BLCR_plot_functions/predictive_pmf_T_by_ASD_iter.R')
source('BLCR_plot_functions/compute_pi_asd.R')
source('BLCR_plot_functions/pmf_T_given_class.R')
source('BLCR_plot_functions/predictive_pmf_T_by_ASD.R')
source('BLCR_plot_functions/build_plot_df.R')
source('BLCR_plot_functions/plot_predictive_T.R')
source('BLCR_plot_functions/compute_pi_asd_across_iterations.R')
source('BLCR_plot_functions/compute_theta_samples_from_counts.R')
source('BLCR_plot_functions/plot_subscale_diff_boxplot.R')
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
library(patchwork)
cb_blue <- "#4477AA"
cb_orange <- "#EE7733"
cb_cols <- c(cb_blue, cb_orange)

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
  coord_cartesian(clip = "off") +  # Add this line
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
 #   plot = final_ridgeline_plot_sim1, 
 #   width = 10, 
 #   height = 6, 
 #   dpi = 600 
 # )

# ggsave(
#   "./sim_study_plots/sim_study1_plots/beta_plots/ridgeline_plot_sim1.pdf",
#   plot = final_ridgeline_plot_sim1,
#   width = 8, height = 3.5
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

#png("./sim_study_plots/sim_study1_plots/sim1_no_varsel_log_post.png", width = 800, height = 600)
plot(sim1_LCR_fit_no_varsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 1', xlab = 'Iteration Number', ylab = 'log posterior')
#dev.off()

#png("./sim_study_plots/sim_study1_plots/sim1_no_varsel_log_post_acf.png", width = 800, height = 600)
acf(sim1_LCR_fit_no_varsel$samples$log_post_samples, main = 'Simulated Data 1')
#dev.off()





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

sim1_LCR_no_varsel_vs_itemsel_mean_abs_diff <-  mean(abs(sim1_LCR_fit_no_varsel$beta_estimate - sim1_LCR_fit_itemsel$beta_estimate))

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

#We want to compare the results for different sample sizes, so we'll subsample 150, 200, 250, 300, 350, 400, 450, for comparison with the whole sample for the variable selection results  

sample_indices_150_1 <- sample(1:500, 150)
sample_indices_150_2 <- sample(1:500, 150)
sample_indices_150_3 <- sample(1:500, 150)
sample_indices_150_4 <- sample(1:500, 150)
sample_indices_150_5 <- sample(1:500, 150)
sample_indices_150_6 <- sample(1:500, 150)
sample_indices_150_7 <- sample(1:500, 150)
sample_indices_150_8 <- sample(1:500, 150)
sample_indices_150_9 <- sample(1:500, 150)
sample_indices_150_10 <- sample(1:500, 150)


set.seed(128)

sample_indices_200_1 <- sample(1:500, 200)
sample_indices_200_2 <- sample(1:500, 200)
sample_indices_200_3 <- sample(1:500, 200)
sample_indices_200_4 <- sample(1:500, 200)
sample_indices_200_5 <- sample(1:500, 200)
sample_indices_200_6 <- sample(1:500, 200)
sample_indices_200_7 <- sample(1:500, 200)
sample_indices_200_8 <- sample(1:500, 200)
sample_indices_200_9 <- sample(1:500, 200)
sample_indices_200_10 <- sample(1:500, 200)

set.seed(128)

sample_indices_300_1 <- sample(1:500, 300)
sample_indices_300_2 <- sample(1:500, 300)
sample_indices_300_3 <- sample(1:500, 300)
sample_indices_300_4 <- sample(1:500, 300)
sample_indices_300_5 <- sample(1:500, 300)
sample_indices_300_6 <- sample(1:500, 300)
sample_indices_300_7 <- sample(1:500, 300)
sample_indices_300_8 <- sample(1:500, 300)
sample_indices_300_9 <- sample(1:500, 300)
sample_indices_300_10 <- sample(1:500, 300)

set.seed(128)

sample_indices_400_1 <- sample(1:500, 400)
sample_indices_400_2 <- sample(1:500, 400)
sample_indices_400_3 <- sample(1:500, 400)
sample_indices_400_4 <- sample(1:500, 400)
sample_indices_400_5 <- sample(1:500, 400)
sample_indices_400_6 <- sample(1:500, 400)
sample_indices_400_7 <- sample(1:500, 400)
sample_indices_400_8 <- sample(1:500, 400)
sample_indices_400_9 <- sample(1:500, 400)
sample_indices_400_10 <- sample(1:500, 400)









sim1_X_150_1 <- sim1_data$X[sample_indices_150_1,]
sim1_X_150_2 <- sim1_data$X[sample_indices_150_2,]
sim1_X_150_3 <- sim1_data$X[sample_indices_150_3,]
sim1_X_150_4 <- sim1_data$X[sample_indices_150_4,]
sim1_X_150_5 <- sim1_data$X[sample_indices_150_5,]
sim1_X_150_6 <- sim1_data$X[sample_indices_150_6,]
sim1_X_150_7 <- sim1_data$X[sample_indices_150_7,]
sim1_X_150_8 <- sim1_data$X[sample_indices_150_8,]
sim1_X_150_9 <- sim1_data$X[sample_indices_150_9,]
sim1_X_150_10 <- sim1_data$X[sample_indices_150_10,]


sim1_X_200_1 <- sim1_data$X[sample_indices_200_1,]
sim1_X_200_2 <- sim1_data$X[sample_indices_200_2,]
sim1_X_200_3 <- sim1_data$X[sample_indices_200_3,]
sim1_X_200_4 <- sim1_data$X[sample_indices_200_4,]
sim1_X_200_5 <- sim1_data$X[sample_indices_200_5,]
sim1_X_200_6 <- sim1_data$X[sample_indices_200_6,]
sim1_X_200_7 <- sim1_data$X[sample_indices_200_7,]
sim1_X_200_8 <- sim1_data$X[sample_indices_200_8,]
sim1_X_200_9 <- sim1_data$X[sample_indices_200_9,]
sim1_X_200_10 <- sim1_data$X[sample_indices_200_10,]


sim1_X_300_1 <- sim1_data$X[sample_indices_300_1,]
sim1_X_300_2 <- sim1_data$X[sample_indices_300_2,]
sim1_X_300_3 <- sim1_data$X[sample_indices_300_3,]
sim1_X_300_4 <- sim1_data$X[sample_indices_300_4,]
sim1_X_300_5 <- sim1_data$X[sample_indices_300_5,]
sim1_X_300_6 <- sim1_data$X[sample_indices_300_6,]
sim1_X_300_7 <- sim1_data$X[sample_indices_300_7,]
sim1_X_300_8 <- sim1_data$X[sample_indices_300_8,]
sim1_X_300_9 <- sim1_data$X[sample_indices_300_9,]
sim1_X_300_10 <- sim1_data$X[sample_indices_300_10,]


sim1_X_400_1 <- sim1_data$X[sample_indices_400_1,]
sim1_X_400_2 <- sim1_data$X[sample_indices_400_2,]
sim1_X_400_3 <- sim1_data$X[sample_indices_400_3,]
sim1_X_400_4 <- sim1_data$X[sample_indices_400_4,]
sim1_X_400_5 <- sim1_data$X[sample_indices_400_5,]
sim1_X_400_6 <- sim1_data$X[sample_indices_400_6,]
sim1_X_400_7 <- sim1_data$X[sample_indices_400_7,]
sim1_X_400_8 <- sim1_data$X[sample_indices_400_8,]
sim1_X_400_9 <- sim1_data$X[sample_indices_400_9,]
sim1_X_400_10 <- sim1_data$X[sample_indices_400_10,]












sim1_Y_150_1 <- sim1_data$Y[sample_indices_150_1,]
sim1_Y_150_2 <- sim1_data$Y[sample_indices_150_2,]
sim1_Y_150_3 <- sim1_data$Y[sample_indices_150_3,]
sim1_Y_150_4 <- sim1_data$Y[sample_indices_150_4,]
sim1_Y_150_5 <- sim1_data$Y[sample_indices_150_5,]
sim1_Y_150_6 <- sim1_data$Y[sample_indices_150_6,]
sim1_Y_150_7 <- sim1_data$Y[sample_indices_150_7,]
sim1_Y_150_8 <- sim1_data$Y[sample_indices_150_8,]
sim1_Y_150_9 <- sim1_data$Y[sample_indices_150_9,]
sim1_Y_150_10 <- sim1_data$Y[sample_indices_150_10,]


sim1_Y_200_1 <- sim1_data$Y[sample_indices_200_1,]
sim1_Y_200_2 <- sim1_data$Y[sample_indices_200_2,]
sim1_Y_200_3 <- sim1_data$Y[sample_indices_200_3,]
sim1_Y_200_4 <- sim1_data$Y[sample_indices_200_4,]
sim1_Y_200_5 <- sim1_data$Y[sample_indices_200_5,]
sim1_Y_200_6 <- sim1_data$Y[sample_indices_200_6,]
sim1_Y_200_7 <- sim1_data$Y[sample_indices_200_7,]
sim1_Y_200_8 <- sim1_data$Y[sample_indices_200_8,]
sim1_Y_200_9 <- sim1_data$Y[sample_indices_200_9,]
sim1_Y_200_10 <- sim1_data$Y[sample_indices_200_10,]


sim1_Y_300_1 <- sim1_data$Y[sample_indices_300_1,]
sim1_Y_300_2 <- sim1_data$Y[sample_indices_300_2,]
sim1_Y_300_3 <- sim1_data$Y[sample_indices_300_3,]
sim1_Y_300_4 <- sim1_data$Y[sample_indices_300_4,]
sim1_Y_300_5 <- sim1_data$Y[sample_indices_300_5,]
sim1_Y_300_6 <- sim1_data$Y[sample_indices_300_6,]
sim1_Y_300_7 <- sim1_data$Y[sample_indices_300_7,]
sim1_Y_300_8 <- sim1_data$Y[sample_indices_300_8,]
sim1_Y_300_9 <- sim1_data$Y[sample_indices_300_9,]
sim1_Y_300_10 <- sim1_data$Y[sample_indices_300_10,]


sim1_Y_400_1 <- sim1_data$Y[sample_indices_400_1,]
sim1_Y_400_2 <- sim1_data$Y[sample_indices_400_2,]
sim1_Y_400_3 <- sim1_data$Y[sample_indices_400_3,]
sim1_Y_400_4 <- sim1_data$Y[sample_indices_400_4,]
sim1_Y_400_5 <- sim1_data$Y[sample_indices_400_5,]
sim1_Y_400_6 <- sim1_data$Y[sample_indices_400_6,]
sim1_Y_400_7 <- sim1_data$Y[sample_indices_400_7,]
sim1_Y_400_8 <- sim1_data$Y[sample_indices_400_8,]
sim1_Y_400_9 <- sim1_data$Y[sample_indices_400_9,]
sim1_Y_400_10 <- sim1_data$Y[sample_indices_400_10,]






  

set.seed(128)

sim1_LCR_fit_varsel_full <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_150_1 <- LCR_Gibbs(X = sim1_X_150_1, Y = sim1_Y_150_1, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                      theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                      relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


set.seed(128)
sim1_LCR_fit_varsel_150_2 <- LCR_Gibbs(X = sim1_X_150_2, Y = sim1_Y_150_2, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_150_3 <- LCR_Gibbs(X = sim1_X_150_3, Y = sim1_Y_150_3, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_150_4 <- LCR_Gibbs(X = sim1_X_150_4, Y = sim1_Y_150_4, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_150_5 <- LCR_Gibbs(X = sim1_X_150_5, Y = sim1_Y_150_5, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_150_6 <- LCR_Gibbs(X = sim1_X_150_6, Y = sim1_Y_150_6, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_150_7 <- LCR_Gibbs(X = sim1_X_150_7, Y = sim1_Y_150_7, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_150_8 <- LCR_Gibbs(X = sim1_X_150_8, Y = sim1_Y_150_8, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_150_9 <- LCR_Gibbs(X = sim1_X_150_9, Y = sim1_Y_150_9, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_150_10 <- LCR_Gibbs(X = sim1_X_150_10, Y = sim1_Y_150_10, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


sim1_LCR_pred_inclusion_mean_150 <- apply(cbind(sim1_LCR_fit_varsel_150_1$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_150_2$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_150_3$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_150_4$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_150_5$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_150_6$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_150_7$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_150_8$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_150_9$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_150_10$cov_inclusion_prob), 1, mean)





set.seed(128)
sim1_LCR_fit_varsel_200_1 <- LCR_Gibbs(X = sim1_X_200_1, Y = sim1_Y_200_1, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


set.seed(128)
sim1_LCR_fit_varsel_200_2 <- LCR_Gibbs(X = sim1_X_200_2, Y = sim1_Y_200_2, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_200_3 <- LCR_Gibbs(X = sim1_X_200_3, Y = sim1_Y_200_3, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_200_4 <- LCR_Gibbs(X = sim1_X_200_4, Y = sim1_Y_200_4, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_200_5 <- LCR_Gibbs(X = sim1_X_200_5, Y = sim1_Y_200_5, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_200_6 <- LCR_Gibbs(X = sim1_X_200_6, Y = sim1_Y_200_6, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_200_7 <- LCR_Gibbs(X = sim1_X_200_7, Y = sim1_Y_200_7, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_200_8 <- LCR_Gibbs(X = sim1_X_200_8, Y = sim1_Y_200_8, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_200_9 <- LCR_Gibbs(X = sim1_X_200_9, Y = sim1_Y_200_9, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_200_10 <- LCR_Gibbs(X = sim1_X_200_10, Y = sim1_Y_200_10, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                        theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                        relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim1_LCR_pred_inclusion_mean_200 <- apply(cbind(sim1_LCR_fit_varsel_200_1$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_200_2$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_200_3$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_200_4$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_200_5$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_200_6$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_200_7$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_200_8$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_200_9$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_200_10$cov_inclusion_prob), 1, mean)

set.seed(128)
sim1_LCR_fit_varsel_300_1 <- LCR_Gibbs(X = sim1_X_300_1, Y = sim1_Y_300_1, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


set.seed(128)
sim1_LCR_fit_varsel_300_2 <- LCR_Gibbs(X = sim1_X_300_2, Y = sim1_Y_300_2, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_300_3 <- LCR_Gibbs(X = sim1_X_300_3, Y = sim1_Y_300_3, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_300_4 <- LCR_Gibbs(X = sim1_X_300_4, Y = sim1_Y_300_4, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_300_5 <- LCR_Gibbs(X = sim1_X_300_5, Y = sim1_Y_300_5, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_300_6 <- LCR_Gibbs(X = sim1_X_300_6, Y = sim1_Y_300_6, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_300_7 <- LCR_Gibbs(X = sim1_X_300_7, Y = sim1_Y_300_7, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_300_8 <- LCR_Gibbs(X = sim1_X_300_8, Y = sim1_Y_300_8, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_300_9 <- LCR_Gibbs(X = sim1_X_300_9, Y = sim1_Y_300_9, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_300_10 <- LCR_Gibbs(X = sim1_X_300_10, Y = sim1_Y_300_10, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                        theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                        relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim1_LCR_pred_inclusion_mean_300 <- apply(cbind(sim1_LCR_fit_varsel_300_1$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_300_2$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_300_3$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_300_4$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_300_5$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_300_6$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_300_7$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_300_8$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_300_9$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_300_10$cov_inclusion_prob), 1, mean)

set.seed(128)
sim1_LCR_fit_varsel_400_1 <- LCR_Gibbs(X = sim1_X_400_1, Y = sim1_Y_400_1, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


set.seed(128)
sim1_LCR_fit_varsel_400_2 <- LCR_Gibbs(X = sim1_X_400_2, Y = sim1_Y_400_2, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_400_3 <- LCR_Gibbs(X = sim1_X_400_3, Y = sim1_Y_400_3, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_400_4 <- LCR_Gibbs(X = sim1_X_400_4, Y = sim1_Y_400_4, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_400_5 <- LCR_Gibbs(X = sim1_X_400_5, Y = sim1_Y_400_5, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_400_6 <- LCR_Gibbs(X = sim1_X_400_6, Y = sim1_Y_400_6, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_400_7 <- LCR_Gibbs(X = sim1_X_400_7, Y = sim1_Y_400_7, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_400_8 <- LCR_Gibbs(X = sim1_X_400_8, Y = sim1_Y_400_8, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_400_9 <- LCR_Gibbs(X = sim1_X_400_9, Y = sim1_Y_400_9, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim1_LCR_fit_varsel_400_10 <- LCR_Gibbs(X = sim1_X_400_10, Y = sim1_Y_400_10, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                        theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                        relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim1_LCR_pred_inclusion_mean_400 <- apply(cbind(sim1_LCR_fit_varsel_400_1$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_400_2$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_400_3$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_400_4$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_400_5$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_400_6$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_400_7$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_400_8$cov_inclusion_prob,
                                                sim1_LCR_fit_varsel_400_9$cov_inclusion_prob, 
                                                sim1_LCR_fit_varsel_400_10$cov_inclusion_prob), 1, mean)


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

#We want to try create a ridgeline plot for the coefficients for only the iterations that the predictor was being included in the model 

sim1_varsel_beta_1_inclusion <- sim1_LCR_fit_varsel_full$samples$beta_samples[2,1,which(sim1_LCR_fit_varsel_full$samples$gamma_samples[2,] == 1)]
sim1_varsel_beta_2_inclusion <- sim1_LCR_fit_varsel_full$samples$beta_samples[3,1,which(sim1_LCR_fit_varsel_full$samples$gamma_samples[3,] == 1)]
sim1_varsel_beta_3_inclusion <- sim1_LCR_fit_varsel_full$samples$beta_samples[4,1,which(sim1_LCR_fit_varsel_full$samples$gamma_samples[4,] == 1)]
sim1_varsel_beta_4_inclusion <- sim1_LCR_fit_varsel_full$samples$beta_samples[5,1,which(sim1_LCR_fit_varsel_full$samples$gamma_samples[5,] == 1)]
sim1_varsel_beta_5_inclusion <- sim1_LCR_fit_varsel_full$samples$beta_samples[6,1,which(sim1_LCR_fit_varsel_full$samples$gamma_samples[6,] == 1)]
sim1_varsel_beta_6_inclusion <- sim1_LCR_fit_varsel_full$samples$beta_samples[7,1,which(sim1_LCR_fit_varsel_full$samples$gamma_samples[7,] == 1)]




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
#   plot = final_ridgeline_plot_beta_1_sim2,
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

#Want to create a plot with the two combined

param_levels <- levels(long_beta1_samples_data_sim2$Parameter)  

long_beta1 <- long_beta1_samples_data_sim2 %>% mutate(Group = "Group 1")
long_beta2 <- long_beta2_samples_data_sim2 %>% mutate(Group = "Group 2")
long_both  <- bind_rows(long_beta1, long_beta2) %>%
  mutate(
    Parameter = factor(Parameter, levels = param_levels),
    Group = factor(Group, levels = c("Group 1","Group 2")),
    Strip = interaction(Group, Parameter, sep = " • ", lex.order = TRUE)
  )

true_beta1 <- true_values_beta1_df_sim2 %>% mutate(Group = "Group 1")
true_beta2 <- true_values_beta2_df_sim2 %>% mutate(Group = "Group 2")
true_both  <- bind_rows(true_beta1, true_beta2) %>%
  mutate(
    Parameter = factor(Parameter, levels = param_levels),
    Group = factor(Group, levels = c("Group 1","Group 2")),
    Strip = interaction(Group, Parameter, sep = " • ", lex.order = TRUE)
  ) 

strip_df <- long_both |> dplyr::distinct(Strip)

xbreaks <- seq(-3, 3, by = 1)
xmin <- -3
xmax <- 3

# Define custom labels for the ridgelines
custom_labels <- rev(c(
  expression(beta[20]), expression(beta[21]), expression(beta[22]), 
  expression(beta[23]), expression(beta[24]), expression(beta[25]), 
  expression(beta[26]), expression(beta[10]), expression(beta[11]), 
  expression(beta[12]), expression(beta[13]), expression(beta[14]), 
  expression(beta[15]), expression(beta[16])
))

gg_stacked <- ggplot(long_both, aes(x = Value, y = Strip, fill = Group)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    scale = 0.9, alpha = 0.75, color = "gray30", linewidth = 0.25
  ) +
  geom_segment(
    data = strip_df,
    aes(x = xmin, xend = xmax, y = Strip, yend = Strip),
    inherit.aes = FALSE, color = "gray70", linewidth = 0.6
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  geom_point(
    data = true_both,
    aes(x = TrueValue, y = Strip, color = Group),
    shape = 18, size = 2.8, show.legend = FALSE
  ) +
  geom_text(
    data = true_both,
    aes(x = TrueValue, y = Strip, label = round(TrueValue, 2), color = Group),
    vjust = 2, size = 3, show.legend = FALSE
  ) +
  scale_fill_manual(values = cb_cols) +
  scale_color_manual(values = cb_cols) +
  scale_y_discrete(
    labels = custom_labels,  # Use the custom LaTeX labels
    expand = c(0.01, 0)
  ) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, expand = expansion(mult = 0)) +
  labs(
    x = "Logit Coefficient", y = "",
    title = "Ridgeline Plot of Posterior Distributions",
    subtitle = "Posterior density with median and 95% HDI. Diamonds mark true values.",
    fill = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top")


#ggsave(
#  './sim_study_plots/sim_study2_plots/beta_plots/sim2_combined_ridgeline.png',
#  plot = gg_stacked,
#  width = 8,
#  height = 10,
#  dpi = 600
#)

# ggsave(
#  './sim_study_plots/sim_study2_plots/beta_plots/sim2_combined_ridgeline.pdf',
#  plot = gg_stacked,
#  width = 8,
#  height = 10
# )



#Diagnostic plots

png("./sim_study_plots/sim_study2_plots/sim2_itemsel_log_post.png", width = 800, height = 600)
plot(sim2_LCR_fit_itemsel$samples$log_post_samples, type = 'l', main = 'Simulated Data 2', xlab = 'Iteration Number', ylab = 'log posterior')
dev.off()

png("./sim_study_plots/sim_study2_plots/sim2_itemsel_log_post_acf.png", width = 800, height = 600)
acf(sim2_LCR_fit_itemsel$samples$log_post_samples, main = 'Simulated Data 2')
dev.off()



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






sim2_X_150_1 <- sim2_data$X[sample_indices_150_1,]
sim2_X_150_2 <- sim2_data$X[sample_indices_150_2,]
sim2_X_150_3 <- sim2_data$X[sample_indices_150_3,]
sim2_X_150_4 <- sim2_data$X[sample_indices_150_4,]
sim2_X_150_5 <- sim2_data$X[sample_indices_150_5,]
sim2_X_150_6 <- sim2_data$X[sample_indices_150_6,]
sim2_X_150_7 <- sim2_data$X[sample_indices_150_7,]
sim2_X_150_8 <- sim2_data$X[sample_indices_150_8,]
sim2_X_150_9 <- sim2_data$X[sample_indices_150_9,]
sim2_X_150_10 <- sim2_data$X[sample_indices_150_10,]


sim2_X_200_1 <- sim2_data$X[sample_indices_200_1,]
sim2_X_200_2 <- sim2_data$X[sample_indices_200_2,]
sim2_X_200_3 <- sim2_data$X[sample_indices_200_3,]
sim2_X_200_4 <- sim2_data$X[sample_indices_200_4,]
sim2_X_200_5 <- sim2_data$X[sample_indices_200_5,]
sim2_X_200_6 <- sim2_data$X[sample_indices_200_6,]
sim2_X_200_7 <- sim2_data$X[sample_indices_200_7,]
sim2_X_200_8 <- sim2_data$X[sample_indices_200_8,]
sim2_X_200_9 <- sim2_data$X[sample_indices_200_9,]
sim2_X_200_10 <- sim2_data$X[sample_indices_200_10,]


sim2_X_300_1 <- sim2_data$X[sample_indices_300_1,]
sim2_X_300_2 <- sim2_data$X[sample_indices_300_2,]
sim2_X_300_3 <- sim2_data$X[sample_indices_300_3,]
sim2_X_300_4 <- sim2_data$X[sample_indices_300_4,]
sim2_X_300_5 <- sim2_data$X[sample_indices_300_5,]
sim2_X_300_6 <- sim2_data$X[sample_indices_300_6,]
sim2_X_300_7 <- sim2_data$X[sample_indices_300_7,]
sim2_X_300_8 <- sim2_data$X[sample_indices_300_8,]
sim2_X_300_9 <- sim2_data$X[sample_indices_300_9,]
sim2_X_300_10 <- sim2_data$X[sample_indices_300_10,]


sim2_X_400_1 <- sim2_data$X[sample_indices_400_1,]
sim2_X_400_2 <- sim2_data$X[sample_indices_400_2,]
sim2_X_400_3 <- sim2_data$X[sample_indices_400_3,]
sim2_X_400_4 <- sim2_data$X[sample_indices_400_4,]
sim2_X_400_5 <- sim2_data$X[sample_indices_400_5,]
sim2_X_400_6 <- sim2_data$X[sample_indices_400_6,]
sim2_X_400_7 <- sim2_data$X[sample_indices_400_7,]
sim2_X_400_8 <- sim2_data$X[sample_indices_400_8,]
sim2_X_400_9 <- sim2_data$X[sample_indices_400_9,]
sim2_X_400_10 <- sim2_data$X[sample_indices_400_10,]












sim2_Y_150_1 <- sim2_data$Y[sample_indices_150_1,]
sim2_Y_150_2 <- sim2_data$Y[sample_indices_150_2,]
sim2_Y_150_3 <- sim2_data$Y[sample_indices_150_3,]
sim2_Y_150_4 <- sim2_data$Y[sample_indices_150_4,]
sim2_Y_150_5 <- sim2_data$Y[sample_indices_150_5,]
sim2_Y_150_6 <- sim2_data$Y[sample_indices_150_6,]
sim2_Y_150_7 <- sim2_data$Y[sample_indices_150_7,]
sim2_Y_150_8 <- sim2_data$Y[sample_indices_150_8,]
sim2_Y_150_9 <- sim2_data$Y[sample_indices_150_9,]
sim2_Y_150_10 <- sim2_data$Y[sample_indices_150_10,]


sim2_Y_200_1 <- sim2_data$Y[sample_indices_200_1,]
sim2_Y_200_2 <- sim2_data$Y[sample_indices_200_2,]
sim2_Y_200_3 <- sim2_data$Y[sample_indices_200_3,]
sim2_Y_200_4 <- sim2_data$Y[sample_indices_200_4,]
sim2_Y_200_5 <- sim2_data$Y[sample_indices_200_5,]
sim2_Y_200_6 <- sim2_data$Y[sample_indices_200_6,]
sim2_Y_200_7 <- sim2_data$Y[sample_indices_200_7,]
sim2_Y_200_8 <- sim2_data$Y[sample_indices_200_8,]
sim2_Y_200_9 <- sim2_data$Y[sample_indices_200_9,]
sim2_Y_200_10 <- sim2_data$Y[sample_indices_200_10,]


sim2_Y_300_1 <- sim2_data$Y[sample_indices_300_1,]
sim2_Y_300_2 <- sim2_data$Y[sample_indices_300_2,]
sim2_Y_300_3 <- sim2_data$Y[sample_indices_300_3,]
sim2_Y_300_4 <- sim2_data$Y[sample_indices_300_4,]
sim2_Y_300_5 <- sim2_data$Y[sample_indices_300_5,]
sim2_Y_300_6 <- sim2_data$Y[sample_indices_300_6,]
sim2_Y_300_7 <- sim2_data$Y[sample_indices_300_7,]
sim2_Y_300_8 <- sim2_data$Y[sample_indices_300_8,]
sim2_Y_300_9 <- sim2_data$Y[sample_indices_300_9,]
sim2_Y_300_10 <- sim2_data$Y[sample_indices_300_10,]


sim2_Y_400_1 <- sim2_data$Y[sample_indices_400_1,]
sim2_Y_400_2 <- sim2_data$Y[sample_indices_400_2,]
sim2_Y_400_3 <- sim2_data$Y[sample_indices_400_3,]
sim2_Y_400_4 <- sim2_data$Y[sample_indices_400_4,]
sim2_Y_400_5 <- sim2_data$Y[sample_indices_400_5,]
sim2_Y_400_6 <- sim2_data$Y[sample_indices_400_6,]
sim2_Y_400_7 <- sim2_data$Y[sample_indices_400_7,]
sim2_Y_400_8 <- sim2_data$Y[sample_indices_400_8,]
sim2_Y_400_9 <- sim2_data$Y[sample_indices_400_9,]
sim2_Y_400_10 <- sim2_data$Y[sample_indices_400_10,]


sim2_subsample_indices <- list(
  n150 = list(
    sample_indices_150_1, sample_indices_150_2, sample_indices_150_3, sample_indices_150_4, sample_indices_150_5,
    sample_indices_150_6, sample_indices_150_7, sample_indices_150_8, sample_indices_150_9, sample_indices_150_10
  ),
  n200 = list(
    sample_indices_200_1, sample_indices_200_2, sample_indices_200_3, sample_indices_200_4, sample_indices_200_5,
    sample_indices_200_6, sample_indices_200_7, sample_indices_200_8, sample_indices_200_9, sample_indices_200_10
  ),
  n300 = list(
    sample_indices_300_1, sample_indices_300_2, sample_indices_300_3, sample_indices_300_4, sample_indices_300_5,
    sample_indices_300_6, sample_indices_300_7, sample_indices_300_8, sample_indices_300_9, sample_indices_300_10
  ),
  n400 = list(
    sample_indices_400_1, sample_indices_400_2, sample_indices_400_3, sample_indices_400_4, sample_indices_400_5,
    sample_indices_400_6, sample_indices_400_7, sample_indices_400_8, sample_indices_400_9, sample_indices_400_10
  )
)

saveRDS(
  sim2_subsample_indices,
  file.path("data", "sim2_subsample_indices.rds")
)






set.seed(128)

sim2_LCR_fit_varsel_full <- LCR_Gibbs(X = sim2_data$X, Y = sim2_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                      theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                      relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_150_1 <- LCR_Gibbs(X = sim2_X_150_1, Y = sim2_Y_150_1, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


set.seed(128)
sim2_LCR_fit_varsel_150_2 <- LCR_Gibbs(X = sim2_X_150_2, Y = sim2_Y_150_2, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_150_3 <- LCR_Gibbs(X = sim2_X_150_3, Y = sim2_Y_150_3, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_150_4 <- LCR_Gibbs(X = sim2_X_150_4, Y = sim2_Y_150_4, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_150_5 <- LCR_Gibbs(X = sim2_X_150_5, Y = sim2_Y_150_5, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_150_6 <- LCR_Gibbs(X = sim2_X_150_6, Y = sim2_Y_150_6, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_150_7 <- LCR_Gibbs(X = sim2_X_150_7, Y = sim2_Y_150_7, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_150_8 <- LCR_Gibbs(X = sim2_X_150_8, Y = sim2_Y_150_8, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_150_9 <- LCR_Gibbs(X = sim2_X_150_9, Y = sim2_Y_150_9, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_150_10 <- LCR_Gibbs(X = sim2_X_150_10, Y = sim2_Y_150_10, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                        theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                        relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


sim2_LCR_pred_inclusion_mean_150 <- apply(cbind(sim2_LCR_fit_varsel_150_1$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_150_2$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_150_3$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_150_4$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_150_5$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_150_6$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_150_7$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_150_8$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_150_9$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_150_10$cov_inclusion_prob), 1, mean)





set.seed(128)
sim2_LCR_fit_varsel_200_1 <- LCR_Gibbs(X = sim2_X_200_1, Y = sim2_Y_200_1, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


set.seed(128)
sim2_LCR_fit_varsel_200_2 <- LCR_Gibbs(X = sim2_X_200_2, Y = sim2_Y_200_2, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_200_3 <- LCR_Gibbs(X = sim2_X_200_3, Y = sim2_Y_200_3, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_200_4 <- LCR_Gibbs(X = sim2_X_200_4, Y = sim2_Y_200_4, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_200_5 <- LCR_Gibbs(X = sim2_X_200_5, Y = sim2_Y_200_5, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_200_6 <- LCR_Gibbs(X = sim2_X_200_6, Y = sim2_Y_200_6, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_200_7 <- LCR_Gibbs(X = sim2_X_200_7, Y = sim2_Y_200_7, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_200_8 <- LCR_Gibbs(X = sim2_X_200_8, Y = sim2_Y_200_8, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_200_9 <- LCR_Gibbs(X = sim2_X_200_9, Y = sim2_Y_200_9, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_200_10 <- LCR_Gibbs(X = sim2_X_200_10, Y = sim2_Y_200_10, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                        theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                        relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_pred_inclusion_mean_200 <- apply(cbind(sim2_LCR_fit_varsel_200_1$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_200_2$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_200_3$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_200_4$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_200_5$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_200_6$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_200_7$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_200_8$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_200_9$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_200_10$cov_inclusion_prob), 1, mean)

set.seed(128)
sim2_LCR_fit_varsel_300_1 <- LCR_Gibbs(X = sim2_X_300_1, Y = sim2_Y_300_1, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


set.seed(128)
sim2_LCR_fit_varsel_300_2 <- LCR_Gibbs(X = sim2_X_300_2, Y = sim2_Y_300_2, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_300_3 <- LCR_Gibbs(X = sim2_X_300_3, Y = sim2_Y_300_3, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_300_4 <- LCR_Gibbs(X = sim2_X_300_4, Y = sim2_Y_300_4, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_300_5 <- LCR_Gibbs(X = sim2_X_300_5, Y = sim2_Y_300_5, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_300_6 <- LCR_Gibbs(X = sim2_X_300_6, Y = sim2_Y_300_6, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_300_7 <- LCR_Gibbs(X = sim2_X_300_7, Y = sim2_Y_300_7, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_300_8 <- LCR_Gibbs(X = sim2_X_300_8, Y = sim2_Y_300_8, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_300_9 <- LCR_Gibbs(X = sim2_X_300_9, Y = sim2_Y_300_9, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_300_10 <- LCR_Gibbs(X = sim2_X_300_10, Y = sim2_Y_300_10, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                        theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                        relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_pred_inclusion_mean_300 <- apply(cbind(sim2_LCR_fit_varsel_300_1$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_300_2$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_300_3$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_300_4$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_300_5$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_300_6$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_300_7$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_300_8$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_300_9$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_300_10$cov_inclusion_prob), 1, mean)

set.seed(128)
sim2_LCR_fit_varsel_400_1 <- LCR_Gibbs(X = sim2_X_400_1, Y = sim2_Y_400_1, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


set.seed(128)
sim2_LCR_fit_varsel_400_2 <- LCR_Gibbs(X = sim2_X_400_2, Y = sim2_Y_400_2, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_400_3 <- LCR_Gibbs(X = sim2_X_400_3, Y = sim2_Y_400_3, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_400_4 <- LCR_Gibbs(X = sim2_X_400_4, Y = sim2_Y_400_4, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_400_5 <- LCR_Gibbs(X = sim2_X_400_5, Y = sim2_Y_400_5, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_400_6 <- LCR_Gibbs(X = sim2_X_400_6, Y = sim2_Y_400_6, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_400_7 <- LCR_Gibbs(X = sim2_X_400_7, Y = sim2_Y_400_7, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_400_8 <- LCR_Gibbs(X = sim2_X_400_8, Y = sim2_Y_400_8, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_400_9 <- LCR_Gibbs(X = sim2_X_400_9, Y = sim2_Y_400_9, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                       theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                       relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

set.seed(128)
sim2_LCR_fit_varsel_400_10 <- LCR_Gibbs(X = sim2_X_400_10, Y = sim2_Y_400_10, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                        theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                        relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

sim2_LCR_pred_inclusion_mean_400 <- apply(cbind(sim2_LCR_fit_varsel_400_1$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_400_2$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_400_3$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_400_4$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_400_5$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_400_6$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_400_7$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_400_8$cov_inclusion_prob,
                                                sim2_LCR_fit_varsel_400_9$cov_inclusion_prob, 
                                                sim2_LCR_fit_varsel_400_10$cov_inclusion_prob), 1, mean)


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


#Reading in the data frame for the survey
CSHQ_df <- read.csv('./CSHQ_df.csv')

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

saveRDS(CSHQ_4group_minVI_cluster$cl, "CSHQ_LCA_4group_minVI_cluster")

#VI distance from the 4 group cluster solution to the collapsed model cluster solution

CSHQ_dist_4group_to_point_estimate <- vi.dist(CSHQ_minVI_cluster$cl, CSHQ_4group_minVI_cluster$cl)

#So the 4 group LCA model is contained within the 95% credible ball, hence is a plausible clustering solution
#Checking the cross-classification table and ARI to compare clusterings

table_collapsed_with_4_group <- table(CSHQ_minVI_cluster$cl, CSHQ_4group_minVI_cluster$cl)
ari_collapsed_with_4_group <- adj.rand.index(CSHQ_minVI_cluster$cl, CSHQ_4group_minVI_cluster$cl)

#ARI of 0.86 indicates strong clustering correspondence


#Looking at the group characteristics

CSHQ_groupA <- CSHQ_df[which(CSHQ_4group_minVI_cluster$cl ==1),]
CSHQ_groupB <- CSHQ_df[which(CSHQ_4group_minVI_cluster$cl ==2),]
CSHQ_groupC <- CSHQ_df[which(CSHQ_4group_minVI_cluster$cl ==3),]
CSHQ_groupD <- CSHQ_df[which(CSHQ_4group_minVI_cluster$cl ==4),]

#mean total CSHQ scores
groupA_mean_total_score <- mean(CSHQ_groupA$sleep_data_total_score)
groupB_mean_total_score <- mean(CSHQ_groupB$sleep_data_total_score)
groupC_mean_total_score <- mean(CSHQ_groupC$sleep_data_total_score)
groupD_mean_total_score <- mean(CSHQ_groupD$sleep_data_total_score)

#sd of total CSHQ scores
groupA_sd_total_score <- sd(CSHQ_groupA$sleep_data_total_score)
groupB_sd_total_score <- sd(CSHQ_groupB$sleep_data_total_score)
groupC_sd_total_score <- sd(CSHQ_groupC$sleep_data_total_score)
groupD_sd_total_score <- sd(CSHQ_groupD$sleep_data_total_score)
  
#Proportion above the threshold of 41 indicating a clinical sleep disorder
groupA_sleep_disorder_proportion <- mean(CSHQ_groupA$sleep_disorders)
groupB_sleep_disorder_proportion <- mean(CSHQ_groupB$sleep_disorders)
groupC_sleep_disorder_proportion <- mean(CSHQ_groupC$sleep_disorders)
groupD_sleep_disorder_proportion <- mean(CSHQ_groupD$sleep_disorders)

#mean and sd age by group

groupA_mean_age <- mean(CSHQ_groupA$ChildAge)
groupB_mean_age <- mean(CSHQ_groupB$ChildAge)
groupC_mean_age <- mean(CSHQ_groupC$ChildAge)
groupD_mean_age <- mean(CSHQ_groupD$ChildAge)

groupA_sd_age <- sd(CSHQ_groupA$ChildAge)
groupB_sd_age <- sd(CSHQ_groupB$ChildAge)
groupC_sd_age <- sd(CSHQ_groupC$ChildAge)
groupD_sd_age <- sd(CSHQ_groupD$ChildAge)



#Creating profile plot for 4 group LCA groups

Y_CSHQ_reduced <-Y_CSHQ[,which(CSHQ_preliminary_4group_LCA$model.indicator == 1)] 

#setting the minimum possible subscale totals as the baseline
minBR <- 3
minSOD <- 1
minSD <- 3
minSA <- 4
minNW <- 3
minP <- 3
minDS <- 4

CSHQ_Y_reduced_BR <- Y_CSHQ_reduced[,1:3]
CSHQ_Y_reduced_SOD <- Y_CSHQ_reduced[,4]
CSHQ_Y_reduced_SD <- Y_CSHQ_reduced[,5:7]
CSHQ_Y_reduced_SA <- Y_CSHQ_reduced[,8:11]
CSHQ_Y_reduced_NW <- Y_CSHQ_reduced[,12:14]
CSHQ_Y_reduced_P <- Y_CSHQ_reduced[,15:17]
CSHQ_Y_reduced_DS <- Y_CSHQ_reduced[,18:21]

total_BR_relevelled <- rowSums(CSHQ_Y_reduced_BR) - minBR
total_SOD_relevelled <- CSHQ_Y_reduced_SOD - minSOD
total_SD_relevelled <- rowSums(CSHQ_Y_reduced_SD) - minSD
total_SA_relevelled <- rowSums(CSHQ_Y_reduced_SA) - minSA
total_NW_relevelled <- rowSums(CSHQ_Y_reduced_NW) - minNW
total_P_relevelled <- rowSums(CSHQ_Y_reduced_P) - minP
total_DS_relevelled <- rowSums(CSHQ_Y_reduced_DS) - minDS

total_BR_standardised <- scale(total_BR_relevelled)
total_SOD_standardised <- scale(total_SOD_relevelled)
total_SD_standardised <- scale(total_SD_relevelled)
total_SA_standardised <- scale(total_SA_relevelled)
total_NW_standardised <- scale(total_NW_relevelled)
total_P_standardised <- scale(total_P_relevelled)
total_DS_standardised <- scale(total_DS_relevelled)

subscale_total_mat_CSHQ_itemsel <- cbind(total_BR_relevelled, 
                                         total_SOD_relevelled, 
                                         total_SD_relevelled,
                                         total_SA_relevelled,
                                         total_NW_relevelled,
                                         total_P_relevelled,
                                         total_DS_relevelled) 

subscale_total_mat_CSHQ_itemsel_standardised <- cbind(total_BR_standardised, 
                                                      total_SOD_standardised, 
                                                      total_SD_standardised,
                                                      total_SA_standardised,
                                                      total_NW_standardised,
                                                      total_P_standardised,
                                                      total_DS_standardised) 

colnames(subscale_total_mat_CSHQ_itemsel) <- c('BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'DS')
colnames(subscale_total_mat_CSHQ_itemsel_standardised) <- c('BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'DS')

#subsetting into assigned clusters
subscale_total_mat_CSHQ_LCA_groupA <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_4group_minVI_cluster$cl == 1),]
subscale_total_mat_CSHQ_LCA_groupB <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_4group_minVI_cluster$cl == 2),]
subscale_total_mat_CSHQ_LCA_groupC <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_4group_minVI_cluster$cl == 3),]
subscale_total_mat_CSHQ_LCA_groupD <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_4group_minVI_cluster$cl == 4),]

subscale_total_mat_CSHQ_LCA_groupA_standardised <- subscale_total_mat_CSHQ_itemsel_standardised[which(CSHQ_4group_minVI_cluster$cl == 1),]
subscale_total_mat_CSHQ_LCA_groupB_standardised <- subscale_total_mat_CSHQ_itemsel_standardised[which(CSHQ_4group_minVI_cluster$cl == 2),]
subscale_total_mat_CSHQ_LCA_groupC_standardised <- subscale_total_mat_CSHQ_itemsel_standardised[which(CSHQ_4group_minVI_cluster$cl == 3),]
subscale_total_mat_CSHQ_LCA_groupD_standardised <- subscale_total_mat_CSHQ_itemsel_standardised[which(CSHQ_4group_minVI_cluster$cl == 4),]

#creating a matrix with groups as rows, and corresponding mean subscale totals as columns
CSHQ_LCA_profile_mat <- rbind(apply(subscale_total_mat_CSHQ_LCA_groupA, 2, mean), 
                                  apply(subscale_total_mat_CSHQ_LCA_groupB, 2, mean), 
                                  apply(subscale_total_mat_CSHQ_LCA_groupC, 2, mean),
                                  apply(subscale_total_mat_CSHQ_LCA_groupD, 2, mean))
CSHQ_LCA_profile_mat_standardised <- rbind(apply(subscale_total_mat_CSHQ_LCA_groupA_standardised, 2, mean), 
                              apply(subscale_total_mat_CSHQ_LCA_groupB_standardised, 2, mean), 
                              apply(subscale_total_mat_CSHQ_LCA_groupC_standardised, 2, mean),
                              apply(subscale_total_mat_CSHQ_LCA_groupD_standardised, 2, mean))

subscale_names <- c("Bedtime Resistance", "Sleep Onset Delay", "Sleep Duration", "Sleep Anxiety", "Night Waking", "Parasomnias", "Daytime Sleepiness")
CSHQ_LCA_profile_df <- tibble(
  Subscale = rep(subscale_names, times = 4),
  Score = c(t(CSHQ_LCA_profile_mat)),
  Group = rep(c('A', 'B', 'C', 'D'), each = length(subscale_names))
)

CSHQ_LCA_profile_df_standardised <- tibble(
  Subscale = rep(subscale_names, times = 4),
  Score = c(t(CSHQ_LCA_profile_mat_standardised)),
  Group = rep(c('A', 'B', 'C', 'D'), each = length(subscale_names))
)

CSHQ_LCA_profile_plot <- ggplot(CSHQ_LCA_profile_df,
                                aes(x = Subscale, y = Score, group = Group, color = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.4) +
  theme_minimal(base_size = 18) +
  labs(
    title = "Sleep Profile by Latent Class",
    subtitle = "Mean Relevelled Subscale Scores by Group",
    x = "CSHQ Subscale",
    y = "Mean Score"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18),
    legend.text = element_text(size = 14)
  )

CSHQ_LCA_profile_plot_standardised <- ggplot(CSHQ_LCA_profile_df_standardised,
                                aes(x = Subscale, y = Score, group = Group, color = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.4) +
  theme_minimal(base_size = 18) +
  labs(
    title = "Sleep Profile by Latent Class",
    subtitle = "Mean Relevelled Subscale Scores (standardised) by Group",
    x = "CSHQ Subscale",
    y = "Mean Score"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18),
    legend.text = element_text(size = 14)
  )


#ggsave(
 # filename = "./CSHQ_plots/CSHQ_sleep_profile_plot.png",
 # plot = CSHQ_LCA_profile_plot,
 # width = 16, height = 9, units = "in",
 # dpi = 600
#)

# ggsave(
# filename = "./CSHQ_plots/CSHQ_sleep_profile_plot.pdf",
# plot = CSHQ_LCA_profile_plot,
# width = 16, height = 8.5
# )
# 
# ggsave(
# filename = "./CSHQ_plots/CSHQ_sleep_profile_plot_standardised.pdf",
# plot = CSHQ_LCA_profile_plot_standardised,
# width = 16, height = 8.5
# )














#Creating Mosaic Plots 
#Using variables BR5, SOD1, SD3, SA1, SA3, NW1

CSHQ_mosaic_variables <- CSHQ_preliminary_4group_LCA$itemprob[which(CSHQ_preliminary_4group_LCA$model.indicator == 1)][c(3,4,7,8,10,12)]

mosaic_plot_6grid_selection <- plot_mosaic_gg_LCA(itemprob = CSHQ_mosaic_variables, classprob = CSHQ_preliminary_4group_LCA$classprob, show_y_axis_numbers = TRUE, show_title = FALSE)

mosaic_plot_6grid_titles <- c('Bedtime Resistance 5', 'Sleep Onset Delay', 'Sleep Duration 3', 'Sleep Anxiety 1', 'Sleep Anxiety 3', 'Night Waking 1')

cb_palette_mosaic <- c(
  "Response 1" = "#0072B2",
  "Response 2" = "#E69F00",
  "Response 3" = "#D55E00"
)

mosaic_plot_6grid_selection <- Map(function(p, ttl) {
  p +
    labs(title = ttl, fill = "Response") +
    scale_fill_manual(
      name   = "Response",
      values = cb_palette_mosaic,
      limits = c("Response 1","Response 2","Response 3"),
      labels = c("1","2","3"),
      drop   = FALSE
    ) +
    theme_minimal(base_size = 15) + 
    theme(plot.title = element_text(hjust = 0.5))
}, mosaic_plot_6grid_selection, mosaic_plot_6grid_titles)

blank_y <- theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)
ncol <- 3
mosaic_plots2_6grid_selection <- lapply(seq_along(mosaic_plot_6grid_selection), function(i) {
  col_idx <- ((i - 1) %% ncol) + 1
  if (col_idx != 1) mosaic_plot_6grid_selection[[i]] + blank_y else mosaic_plot_6grid_selection[[i]]
})

CSHQ_LCA_6grid_mosaic_plot <- wrap_plots(mosaic_plots2_6grid_selection, ncol = 3) +
  plot_layout(guides = "collect") +
  theme(
    text = element_text(size = 16),          
    plot.title = element_text(size = 18),    
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 12),     
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank()
  )



# ggsave(
#   filename = "./CSHQ_plots/CSHQ_6grid_LCA_mosaic.png",
#   plot = CSHQ_LCA_6grid_mosaic_plot ,
#   width = 16, height = 9, units = "in",
#   dpi = 600
# )

# ggsave(
#   filename = "./CSHQ_plots/CSHQ_6grid_LCA_mosaic.pdf",
#   plot = CSHQ_LCA_6grid_mosaic_plot ,
#   width = 16, height = 9,
# )


#Now creating a grid to display mosaic plots of all of the item variables. This will consist of a grid of 3 rows and a grid of 4 rows 
#Will keep the mosaic plot legend at the very top for clarity

CSHQ_LCA_full_mosaic_grid1_vars <- CSHQ_preliminary_4group_LCA$itemprob[which(CSHQ_preliminary_4group_LCA$model.indicator == 1)][1:12]
CSHQ_LCA_full_mosaic_grid2_vars <- CSHQ_preliminary_4group_LCA$itemprob[which(CSHQ_preliminary_4group_LCA$model.indicator == 1)][13:21]

#Now want to create each of the grid plots

#Starting with the first one

mosaic_plot_grid_full1 <- plot_mosaic_gg_LCA(itemprob = CSHQ_LCA_full_mosaic_grid1_vars, classprob = CSHQ_preliminary_4group_LCA$classprob, show_y_axis_numbers = TRUE, show_title = FALSE)

mosaic_plot_grid_full1_names <- c('Bedtime Resistance 2', 'Bedtime Resistance 3', 'Bedtime Resistance 5', 'Sleep Onset Delay', 'Sleep Duration 1', 'Sleep Duration 2', 'Sleep Duration 3', 'Sleep Anxiety 1', 'Sleep Anxiety 2', 'Sleep Anxiety 3', 'Sleep Anxiety 4', 'Night Waking 1')

mosaic_plot_grid_full1 <- Map(function(p, ttl) {
  p +
    labs(title = ttl, fill = "Response") +
    scale_fill_manual(
      name   = "Response",
      values = cb_palette_mosaic,
      limits = c("Response 1","Response 2","Response 3"),
      labels = c("1","2","3"),
      drop   = FALSE
    ) +
    theme_minimal(base_size = 15) + 
    theme(plot.title = element_text(hjust = 0.5))
}, mosaic_plot_grid_full1, mosaic_plot_grid_full1_names)

blank_y <- theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)
ncol <- 3
mosaic_plots2_grid_full1 <- lapply(seq_along(mosaic_plot_grid_full1), function(i) {
  col_idx <- ((i - 1) %% ncol) + 1
  if (col_idx != 1) mosaic_plot_grid_full1[[i]] + blank_y else mosaic_plot_grid_full1[[i]]
})



CSHQ_LCA_mosaic_full1_plot <- wrap_plots(mosaic_plots2_grid_full1, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(
    text = element_text(size = 16),          
    plot.title = element_text(size = 18),    
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 12),     
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 10),
    legend.position = "right",
    legend.key.width = unit(0.5, 'cm'),    
    legend.key.height = unit(1.2, 'cm'),   
    panel.grid = element_blank()
  )

#ggsave(
#  filename = './CSHQ_plots/CSHQ_full_grid_LCA_mosaic1.png',
#  plot = CSHQ_LCA_mosaic_full1_plot,
#  width = 16, height = 18, units = 'in',
#  dpi = 600
#)



mosaic_plot_grid_full2 <- plot_mosaic_gg_LCA(itemprob = CSHQ_LCA_full_mosaic_grid2_vars, classprob = CSHQ_preliminary_4group_LCA$classprob, show_y_axis_numbers = TRUE, show_title = FALSE)

mosaic_plot_grid_full2_names <- c('Night Waking 2', 'Night Waking 3', 'Parasomnias 2', 'Parasomnias 3', 'Parasomnias 5', 'Daytime Sleepiness 2', 'Daytime Sleepiness 4', 'Daytime Sleepiness 5', 'Daytime Sleepiness 6')

mosaic_plot_grid_full2 <- Map(function(p, ttl) {
  p +
    labs(title = ttl, fill = "Response") +
    scale_fill_manual(
      name   = "Response",
      values = cb_palette_mosaic,
      limits = c("Response 1","Response 2","Response 3"),
      labels = c("1","2","3"),
      drop   = FALSE
    ) +
    theme_minimal(base_size = 15) + 
    theme(plot.title = element_text(hjust = 0.5))
}, mosaic_plot_grid_full2, mosaic_plot_grid_full2_names)

blank_y <- theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)
ncol <- 3
mosaic_plots2_grid_full2 <- lapply(seq_along(mosaic_plot_grid_full2), function(i) {
  col_idx <- ((i - 1) %% ncol) + 1
  if (col_idx != 1) mosaic_plot_grid_full2[[i]] + blank_y else mosaic_plot_grid_full2[[i]]
})



CSHQ_LCA_mosaic_full2_plot <- wrap_plots(mosaic_plots2_grid_full2, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(
    text = element_text(size = 16),          
    plot.title = element_text(size = 18),    
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 12),     
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 10),
    legend.position = "right",
    legend.key.width = unit(0.5, 'cm'),    
    legend.key.height = unit(1.2, 'cm'),   
    panel.grid = element_blank()
  )

#ggsave(
#  filename = './CSHQ_plots/CSHQ_full_grid_LCA_mosaic2.png',
#  plot = CSHQ_LCA_mosaic_full2_plot,
#  width = 16, height = 13.5, units = 'in',
#  dpi = 600
#)

# ggsave(
#   filename = './CSHQ_plots/CSHQ_full_grid_LCA_mosaic1.pdf',
#   plot = CSHQ_LCA_mosaic_full1_plot,
#   width = 16, height = 18
# )
# 
# ggsave(
#  filename = './CSHQ_plots/CSHQ_full_grid_LCA_mosaic2.pdf',
#  plot = CSHQ_LCA_mosaic_full2_plot,
#  width = 16, height = 13.5
# )


















set.seed(135)

#Running the LCR model without variable selection first

CSHQ_LCR_no_varsel <- LCR_Gibbs(X = X_CSHQ, Y = Y_CSHQ, G = 4,
                                           beta_prior_cov = diag(c(rep(10^2,4),5^2,5^2)),
                                           beta_prior_mean = rep(0,6), theta_hyperparam = 1,
                                           clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE,
                                           verbose = TRUE, relabel = TRUE, n_samples = 5000,
                                           burnin = 1000, thinby = 10)

#We can see some signs of non-identifiability here, with the log posterior switching between modes

#png("./CSHQ_plots/CSHQ_LCR_log_post_trace_no_varsel.png", width = 800, height = 600)
plot(CSHQ_LCR_no_varsel$samples$log_post_samples, type = 'l', main = 'CSHQ Data LCR log posterior (no variable selection)', xlab = 'Iteration Number', ylab = 'log posterior')
#dev.off()

#png("./CSHQ_plots/CSHQ_LCR_log_post_acf_no_varsel.png", width = 800, height = 600)
acf(CSHQ_LCR_no_varsel$samples$log_post_samples, main = 'CSHQ Data LCR log posterior (no variable selection)')
#dev.off()

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

#Creating Mosaic plots for each of the item variables that was included more than 0.5 of the time
names(CSHQ_LCR_itemsel$itemprob) <- colnames(Y_CSHQ)

mosaic_plots_CSHQ_itemsel <- plot_mosaic_gg_LCA(itemprob = CSHQ_LCR_itemsel$itemprob, classprob = CSHQ_LCR_itemsel$pi, show_y_axis_numbers = TRUE)




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

final_ridgeline_plot_betaB_CSHQ
final_ridgeline_plot_betaC_CSHQ
final_ridgeline_plot_betaD_CSHQ


#Creating stacked ridgeline plots

#Want to use the colour palette such that the colours on the profile plot match those on the ridgeline plot
#Getting the Dark2 colour palette details from the documentation

group_colors <- c("Group B" = "#D95F02FF", "Group C" = "#7570B3FF", "Group D" = "#E7298AFF")

long_betaB_samples_data_CSHQ$Group <- "Group B"
long_betaC_samples_data_CSHQ$Group <- "Group C" 
long_betaD_samples_data_CSHQ$Group <- "Group D"

combined_beta_data <- rbind(
  long_betaB_samples_data_CSHQ,
  long_betaC_samples_data_CSHQ,
  long_betaD_samples_data_CSHQ
)

combined_beta_data$Parameter_Group <- interaction(combined_beta_data$Parameter, combined_beta_data$Group)


param_group_levels <- c()
groups_order <- c("Group B", "Group C", "Group D")
for(group in groups_order) {
  for(param in CSHQ_beta_names) { 
    param_group_levels <- c(param_group_levels, paste(param, group, sep = "."))
  }
}


combined_beta_data$Parameter_Group <- factor(combined_beta_data$Parameter_Group, levels = rev(param_group_levels))

custom_param_labels <- c("Intercept", "Age", "Sex", "ASD diagnosis", "ID Diagnosis", "Other diagnosis")


y_labels <- c()
groups_order <- c("Group B", "Group C", "Group D")
for(group in groups_order) {
  for(label in custom_param_labels) {
    y_labels <- c(y_labels, label)
  }
}


y_labels <- rev(y_labels)

CSHQ_itemsel <- ggplot(combined_beta_data, aes(x = Value, y = Parameter_Group, fill = Group)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_fill_manual(values = group_colors) +
  scale_y_discrete(labels = y_labels) +
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions by Group",
    subtitle = "Posterior density with median & 95% HDI"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top"
  ) + 
  coord_cartesian(xlim = c(-5, 5))

print(final_combined_ridgeline_plot_CSHQ_itemsel)

#ggsave(
#  './CSHQ_plots/CSHQ_ridgeline_plot_combined.png',
#  plot = final_combined_ridgeline_plot_CSHQ_itemsel,
#  width = 8,
#  height = 10,
#  dpi = 600
#)


#Diagnostic plots

png("./CSHQ_plots/CSHQ_LCR_itemsel_log_post.png", width = 800, height = 600)
plot(CSHQ_LCR_itemsel$samples$log_post_samples, type = 'l', main = 'CSHQ Data LCR log posterior (item selection)', xlab = 'Iteration Number', ylab = 'log posterior')
dev.off()

png("./CSHQ_plots/CSHQ_LCR_itemsel_log_post_acf.png", width = 800, height = 600)
acf(CSHQ_LCR_itemsel$samples$log_post_samples, main = 'CSHQ Data LCR log posterior (item selection)')
dev.off()


#formatting the clustering iterations into correct format for computing point estimate
CSHQ_LCR_itemsel_cluster_mat <- matrix(0, nrow = dim(CSHQ_LCR_itemsel$samples$z_samples)[3], ncol = nrow(Y_CSHQ))
for (slice in 1:dim(CSHQ_LCR_itemsel_cluster_mat)[1]) {
  CSHQ_LCR_itemsel_cluster_mat[slice, ] <- max.col(CSHQ_LCR_itemsel$samples$z_samples[, , slice])
}

#Computing the point estimate for the clustering solution (for item variable selection model)
CSHQ_itemsel_psm_mat <- comp.psm(CSHQ_LCR_itemsel_cluster_mat)
CSHQ_itemsel_minVI_cluster <- minVI(psm = CSHQ_itemsel_psm_mat, method = 'greedy', start.cl = max.col(CSHQ_LCR_itemsel$Z))

#Creating profile plot across the survey subscales
#The idea is to basically relevel the subscale totals in terms of the minimum possible scores
#since the minimum score for each item is 1, the minimum total for a given subscale will be the number of items for that subscale
#so to get the relevelled subscale total we subtract this minimum from the total for each observation
#We only consider the variables with a posterior inclusion probability of >=0.5

Y_CSHQ_reduced <-Y_CSHQ[,which(CSHQ_LCR_itemsel$item_inclusion_prob >= 0.5)] 

#setting the minimum possible subscale totals as the baseline
minBR <- 3
minSOD <- 1
minSD <- 3
minSA <- 4
minNW <- 3
minP <- 3
minDS <- 4

CSHQ_Y_reduced_BR <- Y_CSHQ_reduced[,1:3]
CSHQ_Y_reduced_SOD <- Y_CSHQ_reduced[,4]
CSHQ_Y_reduced_SD <- Y_CSHQ_reduced[,5:7]
CSHQ_Y_reduced_SA <- Y_CSHQ_reduced[,8:11]
CSHQ_Y_reduced_NW <- Y_CSHQ_reduced[,12:14]
CSHQ_Y_reduced_P <- Y_CSHQ_reduced[,15:17]
CSHQ_Y_reduced_DS <- Y_CSHQ_reduced[,18:21]

total_BR_relevelled <- rowSums(CSHQ_Y_reduced_BR) - minBR
total_SOD_relevelled <- CSHQ_Y_reduced_SOD - minSOD
total_SD_relevelled <- rowSums(CSHQ_Y_reduced_SD) - minSD
total_SA_relevelled <- rowSums(CSHQ_Y_reduced_SA) - minSA
total_NW_relevelled <- rowSums(CSHQ_Y_reduced_NW) - minNW
total_P_relevelled <- rowSums(CSHQ_Y_reduced_P) - minP
total_DS_relevelled <- rowSums(CSHQ_Y_reduced_DS) - minDS

subscale_total_mat_CSHQ_itemsel <- cbind(total_BR_relevelled, 
                                         total_SOD_relevelled, 
                                         total_SD_relevelled,
                                         total_SA_relevelled,
                                         total_NW_relevelled,
                                         total_P_relevelled,
                                         total_DS_relevelled) 

colnames(subscale_total_mat_CSHQ_itemsel) <- c('BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'DS')

#subsetting into assigned clusters
subscale_total_mat_CSHQ_itemsel_groupA <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_itemsel_minVI_cluster$cl == 1),]
subscale_total_mat_CSHQ_itemsel_groupB <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_itemsel_minVI_cluster$cl == 2),]
subscale_total_mat_CSHQ_itemsel_groupC <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_itemsel_minVI_cluster$cl == 3),]
subscale_total_mat_CSHQ_itemsel_groupD <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_itemsel_minVI_cluster$cl == 4),]

#creating a matrix with groups as rows, and corresponding mean subscale totals as columns
CSHQ_itemsel_profile_mat <- rbind(apply(subscale_total_mat_CSHQ_itemsel_groupA, 2, mean), 
                     apply(subscale_total_mat_CSHQ_itemsel_groupB, 2, mean), 
                     apply(subscale_total_mat_CSHQ_itemsel_groupC, 2, mean),
                     apply(subscale_total_mat_CSHQ_itemsel_groupD, 2, mean))

subscale_names <- c("Bedtime Resistance", "Sleep Onset Delay", "Sleep Duration", "Sleep Anxiety", "Night Waking", "Parasomnias", "Daytime Sleepiness")
CSHQ_itemsel_profile_df <- tibble(
  Subscale = rep(subscale_names, times = 4),
  Score = c(t(CSHQ_itemsel_profile_mat)),
  Group = rep(c('A', 'B', 'C', 'D'), each = length(subscale_names))
)

CSHQ_LCR_itemsel_profile_plot <- ggplot(CSHQ_itemsel_profile_df, aes(x = Subscale, y = Score, group = Group, color = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Sleep Profile by Latent Class",
    subtitle = "Mean Relevelled Subscale Scores by Group",
    x = "CSHQ Subscale",
    y = "Mean Score"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


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

#formatting the clustering iterations into correct format for computing point estimate
CSHQ_LCR_varsel_cluster_mat <- matrix(0, nrow = dim(CSHQ_LCR_varsel$samples$z_samples)[3], ncol = nrow(Y_CSHQ))
for (slice in 1:dim(CSHQ_LCR_varsel_cluster_mat)[1]) {
  CSHQ_LCR_varsel_cluster_mat[slice, ] <- max.col(CSHQ_LCR_varsel$samples$z_samples[, , slice])
}

#Computing the point estimate for the clustering solution (for simultaneous variable selection model)
CSHQ_varsel_psm_mat <- comp.psm(CSHQ_LCR_varsel_cluster_mat)
CSHQ_varsel_minVI_cluster <- minVI(psm = CSHQ_varsel_psm_mat, method = 'greedy', start.cl = max.col(CSHQ_LCR_varsel$Z))


#Classifying observations using minimum posterior expected variation of information for the variable selection model for comparison with the item selection model

itemsel_vs_varsel_cross_classification_table <- table(CSHQ_itemsel_minVI_cluster$cl, CSHQ_varsel_minVI_cluster$cl)
itemsel_vs_varsel_ari <- adj.rand.index(CSHQ_itemsel_minVI_cluster$cl, CSHQ_varsel_minVI_cluster$cl)


#Creating ridgeline plots for the beta coefficients for the variable selection model

CSHQ_betaB_samples_varsel <- CSHQ_LCR_varsel$samples$beta_samples[,2,]
CSHQ_betaC_samples_varsel <- CSHQ_LCR_varsel$samples$beta_samples[,3,]
CSHQ_betaD_samples_varsel <- CSHQ_LCR_varsel$samples$beta_samples[,4,]


CSHQ_beta_names <- paste0("beta", 0:5)
rownames(CSHQ_betaB_samples_varsel) <- CSHQ_beta_names
rownames(CSHQ_betaC_samples_varsel) <- CSHQ_beta_names
rownames(CSHQ_betaD_samples_varsel) <- CSHQ_beta_names

long_betaB_samples_data_CSHQ_varsel <- as.data.frame(t(CSHQ_betaB_samples_varsel)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )
long_betaC_samples_data_CSHQ_varsel <- as.data.frame(t(CSHQ_betaC_samples_varsel)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )
long_betaD_samples_data_CSHQ_varsel <- as.data.frame(t(CSHQ_betaD_samples_varsel)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )

long_betaB_samples_data_CSHQ_varsel$Parameter <- factor(long_betaB_samples_data_CSHQ_varsel$Parameter, levels = rev(CSHQ_beta_names))
long_betaC_samples_data_CSHQ_varsel$Parameter <- factor(long_betaC_samples_data_CSHQ_varsel$Parameter, levels = rev(CSHQ_beta_names))
long_betaD_samples_data_CSHQ_varsel$Parameter <- factor(long_betaD_samples_data_CSHQ_varsel$Parameter, levels = rev(CSHQ_beta_names))


math_labels_CSHQ <- expression(
  beta[5], beta[4], beta[3], beta[2], beta[1], beta[0]
)


final_ridgeline_plot_betaB_CSHQ_varsel <- ggplot(long_betaB_samples_data_CSHQ_varsel, aes(x = Value, y = Parameter)) +
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


final_ridgeline_plot_betaC_CSHQ_varsel <- ggplot(long_betaC_samples_data_CSHQ_varsel, aes(x = Value, y = Parameter)) +
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


final_ridgeline_plot_betaD_CSHQ_varsel <- ggplot(long_betaD_samples_data_CSHQ_varsel, aes(x = Value, y = Parameter)) +
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






#Creating Mosaic Plots for the variable selection LCR (using pi estimate for the widths)
#Using variables BR5, SOD1, SD3, SA1, SA3, NW1

CSHQ_mosaic_variables <- CSHQ_LCR_varsel$itemprob[which(CSHQ_LCR_varsel$item.ind == 1)][c(3,4,7,8,10,12)]

mosaic_plot_6grid_selection <- plot_mosaic_gg_LCA(itemprob = CSHQ_mosaic_variables, classprob = CSHQ_LCR_varsel$pi, show_y_axis_numbers = TRUE, show_title = FALSE)

mosaic_plot_6grid_titles <- c('Bedtime Resistance 5', 'Sleep Onset Delay', 'Sleep Duration 3', 'Sleep Anxiety 1', 'Sleep Anxiety 3', 'Night Waking 1')

cb_palette_mosaic <- c(
  "Response 1" = "#0072B2",
  "Response 2" = "#E69F00",
  "Response 3" = "#D55E00"
)

mosaic_plot_6grid_selection <- Map(function(p, ttl) {
  p +
    labs(title = ttl, fill = "Response") +
    scale_fill_manual(
      name   = "Response",
      values = cb_palette_mosaic,
      limits = c("Response 1","Response 2","Response 3"),
      labels = c("1","2","3"),
      drop   = FALSE
    ) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),  
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"       
    )
}, mosaic_plot_6grid_selection, mosaic_plot_6grid_titles)



blank_y <- theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)
ncol <- 3
mosaic_plots2_6grid_selection <- lapply(seq_along(mosaic_plot_6grid_selection), function(i) {
  col_idx <- ((i - 1) %% ncol) + 1
  if (col_idx != 1) mosaic_plot_6grid_selection[[i]] + blank_y else mosaic_plot_6grid_selection[[i]]
})

CSHQ_LCR_6grid_mosaic_plot <- wrap_plots(mosaic_plots2_6grid_selection, ncol = 3) +
  plot_layout(guides = "collect") +
  theme(
    text = element_text(size = 16),          
    plot.title = element_text(size = 18),    
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 12),     
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank()
  )



# ggsave(
#   filename = "./CSHQ_plots/CSHQ_6grid_LCR_mosaic.png",
#   plot = CSHQ_LCR_6grid_mosaic_plot ,
#   width = 16, height = 9, units = "in",
#   dpi = 600
# )
# 
# ggsave(
#   filename = "./CSHQ_plots/CSHQ_6grid_LCR_mosaic.pdf",
#   plot = CSHQ_LCR_6grid_mosaic_plot ,
#   width = 16, height = 9,
# )


#Now creating a grid to display mosaic plots of all of the item variables. This will consist of a grid of 3 rows and a grid of 4 rows 
#Will keep the mosaic plot legend at the very top for clarity

CSHQ_LCR_full_mosaic_grid1_vars <- CSHQ_LCR_varsel$itemprob[which(CSHQ_LCR_varsel$item.ind == 1)][1:12]
CSHQ_LCR_full_mosaic_grid2_vars <- CSHQ_LCR_varsel$itemprob[which(CSHQ_LCR_varsel$item.ind == 1)][13:21]

#Now want to create each of the grid plots

#Starting with the first one

mosaic_plot_grid_full1 <- plot_mosaic_gg_LCA(itemprob = CSHQ_LCR_full_mosaic_grid1_vars, classprob = CSHQ_LCR_varsel$pi, show_y_axis_numbers = TRUE, show_title = FALSE)

mosaic_plot_grid_full1_names <- c('Bedtime Resistance 2', 'Bedtime Resistance 3', 'Bedtime Resistance 5', 'Sleep Onset Delay', 'Sleep Duration 1', 'Sleep Duration 2', 'Sleep Duration 3', 'Sleep Anxiety 1', 'Sleep Anxiety 2', 'Sleep Anxiety 3', 'Sleep Anxiety 4', 'Night Waking 1')

mosaic_plot_grid_full1 <- Map(function(p, ttl) {
  p +
    labs(title = ttl, fill = "Response") +
    scale_fill_manual(
      name   = "Response",
      values = cb_palette_mosaic,
      limits = c("Response 1","Response 2","Response 3"),
      labels = c("1","2","3"),
      drop   = FALSE
    ) +
    theme_minimal(base_size = 15) + 
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",    
      panel.grid = element_blank()  
    )
}, mosaic_plot_grid_full1, mosaic_plot_grid_full1_names)


blank_y <- theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)
ncol <- 3
mosaic_plots2_grid_full1 <- lapply(seq_along(mosaic_plot_grid_full1), function(i) {
  col_idx <- ((i - 1) %% ncol) + 1
  if (col_idx != 1) mosaic_plot_grid_full1[[i]] + blank_y else mosaic_plot_grid_full1[[i]]
})



CSHQ_LCR_mosaic_full1_plot <- wrap_plots(mosaic_plots2_grid_full1, ncol = 3) +
  plot_layout(guides = "collect") +
  theme(
    text = element_text(size = 16),          
    plot.title = element_text(size = 18),    
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 12),     
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 10),
    legend.position = "right",
    legend.key.width = unit(0.5, 'cm'),    
    legend.key.height = unit(1.2, 'cm'),   
    panel.grid = element_blank()
  )


ggsave(
 filename = './CSHQ_plots/CSHQ_full_grid_LCR_mosaic1.png',
 plot = CSHQ_LCR_mosaic_full1_plot,
 width = 16, height = 18, units = 'in',
 dpi = 600
)
ggsave(
 filename = './CSHQ_plots/CSHQ_full_grid_LCR_mosaic1.pdf',
 plot = CSHQ_LCR_mosaic_full1_plot,
 width = 16, height = 18
)

#I accidentally overwrote the LCA mosaic plots here so need to run that again if I want that back.

mosaic_plot_grid_full2 <- plot_mosaic_gg_LCA(itemprob = CSHQ_LCR_full_mosaic_grid2_vars, classprob = CSHQ_LCR_varsel$pi, show_y_axis_numbers = TRUE, show_title = FALSE)

mosaic_plot_grid_full2_names <- c('Night Waking 2', 'Night Waking 3', 'Parasomnias 2', 'Parasomnias 3', 'Parasomnias 5', 'Daytime Sleepiness 2', 'Daytime Sleepiness 4', 'Daytime Sleepiness 5', 'Daytime Sleepiness 6')

mosaic_plot_grid_full2 <- Map(function(p, ttl) {
  p +
    labs(title = ttl, fill = "Response") +
    scale_fill_manual(
      name   = "Response",
      values = cb_palette_mosaic,
      limits = c("Response 1","Response 2","Response 3"),
      labels = c("1","2","3"),
      drop   = FALSE
    ) +
    theme_minimal(base_size = 15) + 
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",      
      panel.grid = element_blank()   
    )
}, mosaic_plot_grid_full2, mosaic_plot_grid_full2_names)


blank_y <- theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)
ncol <- 3
mosaic_plots2_grid_full2 <- lapply(seq_along(mosaic_plot_grid_full2), function(i) {
  col_idx <- ((i - 1) %% ncol) + 1
  if (col_idx != 1) mosaic_plot_grid_full2[[i]] + blank_y else mosaic_plot_grid_full2[[i]]
})



CSHQ_LCR_mosaic_full2_plot <- wrap_plots(mosaic_plots2_grid_full2, ncol = 3) +
  plot_layout(guides = "collect") +
  theme(
    text = element_text(size = 16),          
    plot.title = element_text(size = 18),    
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 12),     
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 10),
    legend.position = "right",
    legend.key.width = unit(0.5, 'cm'),    
    legend.key.height = unit(1.2, 'cm'),   
    panel.grid = element_blank()
  )

# ggsave(
#  filename = './CSHQ_plots/CSHQ_full_grid_LCR_mosaic2.png',
#  plot = CSHQ_LCR_mosaic_full2_plot,
#  width = 16, height = 13.5, units = 'in',
#  dpi = 600
# )
# ggsave(
#  filename = './CSHQ_plots/CSHQ_full_grid_LCR_mosaic2.pdf',
#  plot = CSHQ_LCR_mosaic_full2_plot,
#  width = 16, height = 13.5
# )





#Creating profile plot for 4 group LCR varsel groups

Y_CSHQ_reduced <- Y_CSHQ[,which(CSHQ_LCR_varsel$item.ind == 1)] 

#setting the minimum possible subscale totals as the baseline
minBR <- 3
minSOD <- 1
minSD <- 3
minSA <- 4
minNW <- 3
minP <- 3
minDS <- 4

CSHQ_Y_reduced_BR <- Y_CSHQ_reduced[,1:3]
CSHQ_Y_reduced_SOD <- Y_CSHQ_reduced[,4]
CSHQ_Y_reduced_SD <- Y_CSHQ_reduced[,5:7]
CSHQ_Y_reduced_SA <- Y_CSHQ_reduced[,8:11]
CSHQ_Y_reduced_NW <- Y_CSHQ_reduced[,12:14]
CSHQ_Y_reduced_P <- Y_CSHQ_reduced[,15:17]
CSHQ_Y_reduced_DS <- Y_CSHQ_reduced[,18:21]

total_BR_relevelled <- rowSums(CSHQ_Y_reduced_BR) - minBR
total_SOD_relevelled <- CSHQ_Y_reduced_SOD - minSOD
total_SD_relevelled <- rowSums(CSHQ_Y_reduced_SD) - minSD
total_SA_relevelled <- rowSums(CSHQ_Y_reduced_SA) - minSA
total_NW_relevelled <- rowSums(CSHQ_Y_reduced_NW) - minNW
total_P_relevelled <- rowSums(CSHQ_Y_reduced_P) - minP
total_DS_relevelled <- rowSums(CSHQ_Y_reduced_DS) - minDS

total_BR_standardised <- scale(total_BR_relevelled)
total_SOD_standardised <- scale(total_SOD_relevelled)
total_SD_standardised <- scale(total_SD_relevelled)
total_SA_standardised <- scale(total_SA_relevelled)
total_NW_standardised <- scale(total_NW_relevelled)
total_P_standardised <- scale(total_P_relevelled)
total_DS_standardised <- scale(total_DS_relevelled)

subscale_total_mat_CSHQ_itemsel <- cbind(total_BR_relevelled, 
                                         total_SOD_relevelled, 
                                         total_SD_relevelled,
                                         total_SA_relevelled,
                                         total_NW_relevelled,
                                         total_P_relevelled,
                                         total_DS_relevelled) 

subscale_total_mat_CSHQ_itemsel_standardised <- cbind(total_BR_standardised, 
                                                      total_SOD_standardised, 
                                                      total_SD_standardised,
                                                      total_SA_standardised,
                                                      total_NW_standardised,
                                                      total_P_standardised,
                                                      total_DS_standardised) 

colnames(subscale_total_mat_CSHQ_itemsel) <- c('BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'DS')
colnames(subscale_total_mat_CSHQ_itemsel_standardised) <- c('BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'DS')

#subsetting into assigned clusters
subscale_total_mat_CSHQ_LCR_groupA <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_varsel_minVI_cluster$cl == 1),]
subscale_total_mat_CSHQ_LCR_groupB <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_varsel_minVI_cluster$cl == 2),]
subscale_total_mat_CSHQ_LCR_groupC <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_varsel_minVI_cluster$cl == 3),]
subscale_total_mat_CSHQ_LCR_groupD <- subscale_total_mat_CSHQ_itemsel[which(CSHQ_varsel_minVI_cluster$cl == 4),]

subscale_total_mat_CSHQ_LCR_groupA_standardised <- subscale_total_mat_CSHQ_itemsel_standardised[which(CSHQ_varsel_minVI_cluster$cl == 1),]
subscale_total_mat_CSHQ_LCR_groupB_standardised <- subscale_total_mat_CSHQ_itemsel_standardised[which(CSHQ_varsel_minVI_cluster$cl == 2),]
subscale_total_mat_CSHQ_LCR_groupC_standardised <- subscale_total_mat_CSHQ_itemsel_standardised[which(CSHQ_varsel_minVI_cluster$cl == 3),]
subscale_total_mat_CSHQ_LCR_groupD_standardised <- subscale_total_mat_CSHQ_itemsel_standardised[which(CSHQ_varsel_minVI_cluster$cl == 4),]

#creating a matrix with groups as rows, and corresponding mean subscale totals as columns
CSHQ_LCR_profile_mat <- rbind(apply(subscale_total_mat_CSHQ_LCR_groupA, 2, mean), 
                              apply(subscale_total_mat_CSHQ_LCR_groupB, 2, mean), 
                              apply(subscale_total_mat_CSHQ_LCR_groupC, 2, mean),
                              apply(subscale_total_mat_CSHQ_LCR_groupD, 2, mean))
CSHQ_LCR_profile_mat_standardised <- rbind(apply(subscale_total_mat_CSHQ_LCR_groupA_standardised, 2, mean), 
                                           apply(subscale_total_mat_CSHQ_LCR_groupB_standardised, 2, mean), 
                                           apply(subscale_total_mat_CSHQ_LCR_groupC_standardised, 2, mean),
                                           apply(subscale_total_mat_CSHQ_LCR_groupD_standardised, 2, mean))

subscale_names <- c("Bedtime Resistance", "Sleep Onset Delay", "Sleep Duration", "Sleep Anxiety", "Night Waking", "Parasomnias", "Daytime Sleepiness")
CSHQ_LCR_profile_df <- tibble(
  Subscale = rep(subscale_names, times = 4),
  Score = c(t(CSHQ_LCR_profile_mat)),
  Group = rep(c('A', 'B', 'C', 'D'), each = length(subscale_names))
)

CSHQ_LCR_profile_df_standardised <- tibble(
  Subscale = rep(subscale_names, times = 4),
  Score = c(t(CSHQ_LCR_profile_mat_standardised)),
  Group = rep(c('A', 'B', 'C', 'D'), each = length(subscale_names))
)

CSHQ_LCR_profile_plot <- ggplot(CSHQ_LCR_profile_df,
                                aes(x = Subscale, y = Score, group = Group, color = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.4) +
  theme_minimal(base_size = 18) +
  labs(
    title = "Sleep Profile by Latent Class",
    subtitle = "Mean Relevelled Subscale Scores by Group",
    x = "CSHQ Subscale",
    y = "Mean Score"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18),
    legend.text = element_text(size = 14)
  )

CSHQ_LCR_profile_plot_standardised <- ggplot(CSHQ_LCR_profile_df_standardised,
                                             aes(x = Subscale, y = Score, group = Group, color = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.4) +
  theme_minimal(base_size = 18) +
  labs(
    title = "Sleep Profile by Latent Class",
    subtitle = "Mean Relevelled Subscale Scores (standardised) by Group",
    x = "CSHQ Subscale",
    y = "Mean Score"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18),
    legend.text = element_text(size = 14)
  )



# ggsave(
# filename = "./CSHQ_plots/CSHQ_LCR_sleep_profile_plot.png",
# plot = CSHQ_LCR_profile_plot,
# width = 16, height = 9, units = "in",
# dpi = 600
# )
# ggsave(
# filename = "./CSHQ_plots/CSHQ_LCR_sleep_profile_plot.pdf",
# plot = CSHQ_LCR_profile_plot,
# width = 16, height = 8.5
# )
# ggsave(
# filename = "./CSHQ_plots/CSHQ_LCR_sleep_profile_plot.png",
# plot = CSHQ_LCR_profile_plot_standardised,
# width = 16, height = 9, units = "in",
# dpi = 600
# )
# ggsave(
# filename = "./CSHQ_plots/CSHQ_LCR_sleep_profile_plot_standardised.pdf",
# plot = CSHQ_LCR_profile_plot_standardised,
# width = 16, height = 8.5
# )







#Computing the probabilities of giving a response of either sometimes or usually

sometimes_usually_combined_itemprob_LCR <- lapply(CSHQ_LCR_varsel$itemprob, function(mat) {
  cbind(mat[, 1], mat[, 2] + mat[, 3])
})

names(sometimes_usually_combined_itemprob_LCR) <- colnames(Y_CSHQ)
















#computing estimate for ASD beta coefficient for the variable selection model, by dropping all iterations which contain all zeros for the values of beta

ASD_beta_coeff_samples_varsel <- CSHQ_LCR_varsel$samples$beta_samples[4,,, drop = FALSE] 
ASD_beta_coeff_samples_varsel_with_intercept <- CSHQ_LCR_varsel$samples$beta_samples[c(1,4),,, drop = FALSE] 

zero_slices <- numeric(5000)
for (slice in 1:5000){
  zero_slices[slice] <- 1*ifelse(CSHQ_LCR_varsel$samples$gamma_samples[4,slice] == 0, 1, 0)
}

ASD_beta_coeff_nonzero_samples <- ASD_beta_coeff_samples_varsel[,,which(zero_slices == 0), drop = FALSE]
ASD_beta_coeff_nonzero_samples_with_intercept <- ASD_beta_coeff_samples_varsel_with_intercept[,,which(zero_slices == 0), drop = FALSE]

CSHQ_ASD_beta_estimate_varsel <- apply(ASD_beta_coeff_nonzero_samples, c(1,2), mean)
CSHQ_ASD_beta_estimate_varsel_with_intercept <- apply(ASD_beta_coeff_nonzero_samples_with_intercept, c(1,2), mean)
CSHQ_ASD_beta_sd_varsel_with_intercept <- apply(ASD_beta_coeff_nonzero_samples_with_intercept, c(1,2), sd)



#Creating a ridgeline plot with only the intercept and ASD variable

ridgeline_array_ASD_only <- ASD_beta_coeff_nonzero_samples_with_intercept[,2:4,, drop = FALSE]



CSHQ_betaB_samples_ASD_only <- ridgeline_array_ASD_only[,1,]
CSHQ_betaC_samples_ASD_only <- ridgeline_array_ASD_only[,2,]
CSHQ_betaD_samples_ASD_only <- ridgeline_array_ASD_only[,3,]


CSHQ_beta_names_ASD_only <- paste0("beta", 0:1)
rownames(CSHQ_betaB_samples_ASD_only) <- CSHQ_beta_names_ASD_only
rownames(CSHQ_betaC_samples_ASD_only) <- CSHQ_beta_names_ASD_only
rownames(CSHQ_betaD_samples_ASD_only) <- CSHQ_beta_names_ASD_only

long_betaB_samples_data_CSHQ_ASD_only <- as.data.frame(t(CSHQ_betaB_samples_ASD_only)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )
long_betaC_samples_data_CSHQ_ASD_only <- as.data.frame(t(CSHQ_betaC_samples_ASD_only)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )
long_betaD_samples_data_CSHQ_ASD_only <- as.data.frame(t(CSHQ_betaD_samples_ASD_only)) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Parameter", 
    values_to = "Value"
  )

long_betaB_samples_data_CSHQ_ASD_only$Parameter <- factor(long_betaB_samples_data_CSHQ_ASD_only$Parameter, levels = rev(CSHQ_beta_names_ASD_only))
long_betaC_samples_data_CSHQ_ASD_only$Parameter <- factor(long_betaC_samples_data_CSHQ_ASD_only$Parameter, levels = rev(CSHQ_beta_names_ASD_only))
long_betaD_samples_data_CSHQ_ASD_only$Parameter <- factor(long_betaD_samples_data_CSHQ_ASD_only$Parameter, levels = rev(CSHQ_beta_names_ASD_only))


math_labels_CSHQ_ASD_only <- expression(
  beta[1], beta[0]
)

final_ridgeline_plot_betaB_CSHQ_ASD_only <- ggplot(long_betaB_samples_data_CSHQ_ASD_only, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_y_discrete(labels = math_labels_CSHQ_ASD_only) + 
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


final_ridgeline_plot_betaC_CSHQ_ASD_only <- ggplot(long_betaC_samples_data_CSHQ_ASD_only, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_y_discrete(labels = math_labels_CSHQ_ASD_only) + 
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

final_ridgeline_plot_betaD_CSHQ_ASD_only <- ggplot(long_betaD_samples_data_CSHQ_ASD_only, aes(x = Value, y = Parameter)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    fill = cb_blue,
    alpha = 0.7,
    scale = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_y_discrete(labels = math_labels_CSHQ_ASD_only) + 
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

final_ridgeline_plot_betaB_CSHQ
final_ridgeline_plot_betaC_CSHQ
final_ridgeline_plot_betaD_CSHQ


#Creating stacked ridgeline plots

#Want to use the colour palette such that the colours on the profile plot match those on the ridgeline plot
#Getting the Dark2 colour palette details from the documentation

group_colors <- c("Group B" = "#D95F02FF", "Group C" = "#7570B3FF", "Group D" = "#E7298AFF")

long_betaB_samples_data_CSHQ_ASD_only$Group <- "Group B"
long_betaC_samples_data_CSHQ_ASD_only$Group <- "Group C" 
long_betaD_samples_data_CSHQ_ASD_only$Group <- "Group D"

combined_beta_data_ASD_only <- rbind(
  long_betaB_samples_data_CSHQ_ASD_only,
  long_betaC_samples_data_CSHQ_ASD_only,
  long_betaD_samples_data_CSHQ_ASD_only
)

combined_beta_data_ASD_only$Parameter_Group <- interaction(combined_beta_data_ASD_only$Parameter, combined_beta_data_ASD_only$Group)


param_group_levels <- c()
groups_order <- c("Group B", "Group C", "Group D")
for(group in groups_order) {
  for(param in CSHQ_beta_names_ASD_only) { 
    param_group_levels <- c(param_group_levels, paste(param, group, sep = "."))
  }
}


combined_beta_data_ASD_only$Parameter_Group <- factor(combined_beta_data_ASD_only$Parameter_Group, levels = rev(param_group_levels))

custom_param_labels_ASD_only <- c("Intercept", "ASD diagnosis")


y_labels <- c()
groups_order <- c("Group B", "Group C", "Group D")
for(group in groups_order) {
  for(label in custom_param_labels_ASD_only) {
    y_labels <- c(y_labels, label)
  }
}


y_labels <- rev(y_labels)

final_combined_ridgeline_plot_CSHQ_ASD_only <- ggplot(combined_beta_data_ASD_only, aes(x = Value, y = Parameter_Group, fill = Group)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    alpha = 0.7,
    scale = 0.9,
    trim = FALSE,           
    rel_min_height = 0,     
    from = -5, to = 5       
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_fill_manual(values = group_colors) +
  scale_y_discrete(labels = y_labels) +
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Ridgeline Plot of Posterior Distributions by Group",
    subtitle = "Posterior density with median & 95% HDI"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top"
  ) + 
  coord_cartesian(xlim = c(-5, 5))



# ggsave(
#  './CSHQ_plots/CSHQ_ridgeline_plot_combined_ASD_only.png',
#  plot = final_combined_ridgeline_plot_CSHQ_ASD_only,
#  width = 8,
#  height = 5,
#  dpi = 600
# )

# ggsave(
#   './CSHQ_plots/CSHQ_ridgeline_plot_combined_ASD_only.pdf',
#   plot = final_combined_ridgeline_plot_CSHQ_ASD_only,
#   width = 8,
#   height = 3.75
# )


















#Creating visualisations tailored to the LCR model with CSHQ data specifically



#Plotting the pmf of total CSHQ score for ASD vs non-ASD, with stacked plots for group membership 
asd_vs_nonasd_CSHQ_total_density <- plot_predictive_T(theta = CSHQ_LCR_varsel$itemprob, beta = CSHQ_ASD_beta_estimate_varsel_with_intercept)


# In order to get the expected subscale totals, we need to first extract point estimates for the theta values 
#on each of the sampler iterations.

theta_samples_collapsed_estimate <- compute_theta_samples_from_counts(N_gjk_samples = CSHQ_LCR_varsel$samples$N_gjk_samples, N_g_samples = CSHQ_LCR_varsel$samples$N_g_samples)

# Need to account for permutation of groups as well

theta_samples_collapsed_estimate_perm <- theta_samples_collapsed_estimate[c(4,1,3,2),,,]









#We'll plot with all iterations included as well as the ones where only ASD was included. Also should compare with standardising per subscale in some way

#We have the iterations for which ASD is not included as zero_slices, so we can easily subset using this




###CURRENTLY IN THE PROCESS OF CONVERTING THIS INTO A FUNCTION!!!


CSHQ_LCR_varsel_compute_expected_subscale_totals_asd_nonzero <- compute_expected_subscale_total_across_iterations(theta_samples = theta_samples_collapsed_estimate_perm[,,,which(zero_slices == 0)],
                                                                                                          asd = 1, 
                                                                                                          beta_samples = CSHQ_LCR_varsel$samples$beta_samples[,,which(zero_slices == 0)])


CSHQ_LCR_varsel_compute_expected_subscale_totals_non_asd_nonzero <- compute_expected_subscale_total_across_iterations(theta_samples = theta_samples_collapsed_estimate_perm[,,,which(zero_slices == 0)],
                                                                                                              asd = 0, 
                                                                                                              beta_samples = CSHQ_LCR_varsel$samples$beta_samples[,,which(zero_slices == 0)])





diff_CSHQ_LCR_varsel_compute_expected_subscale_totals_nonzero <- CSHQ_LCR_varsel_compute_expected_subscale_totals_asd_nonzero$expected_subscale_total - CSHQ_LCR_varsel_compute_expected_subscale_totals_non_asd_nonzero$expected_subscale_total

diff_nonzero_df <- as.data.frame(diff_CSHQ_LCR_varsel_compute_expected_subscale_totals_nonzero)
colnames(diff_nonzero_df) <- c('BR','SOD','SD','SA','NW','P','SDB','DS')

diff_df_long_nonzero <- diff_nonzero_df %>%
  pivot_longer(everything(), names_to = "Column", values_to = "Value")


ggplot(diff_df_long_nonzero, aes(x = Column, y = Value, fill = Column)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplots of Difference between Subscale Totals", x = "Columns", y = "Values")











label_map <- c(
  'BR' = 'Bedtime Resistance',
  'SOD' = 'Sleep Onset Delay', 
  'SD' = 'Sleep Duration',
  'SA' = 'Sleep Anxiety',
  'NW' = 'Night Waking',
  'P' = 'Parasomnias',
  'SDB' = 'Sleep Disordered Breathing',
  'DS' = 'Daytime Sleepiness'
)

diff_df_long_nonzero <- diff_df_long_nonzero %>%
  mutate(
    Column_long = label_map[Column],
    median_val = tapply(Value, Column, median)[Column]
  ) %>%
  arrange(median_val) %>%
  mutate(Column_long = factor(Column_long, levels = unique(Column_long)))

boxplot_expected_CSHQ_subscale_total_difference <- ggplot(diff_df_long_nonzero, aes(x = Column_long, y = Value, fill = Column_long)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.5) +
  stat_summary(fun.data = "median_hilow", geom = "errorbar", 
               fun.args = list(conf.int = 0.95), width = 0.2, size = 1.5) +
  stat_summary(fun = "median", geom = "point", size = 5, color = "white") +
  labs(
    title = "Posterior Differences in Expected CSHQ Subscale Totals (ASD - Non-ASD)",
    x = "CSHQ Subscales",
    y = "Differences in Expected Total"
  ) +
  theme_classic(base_size = 16) +  # Global base font size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),     
    axis.text.y = element_text(size = 14),                                
    axis.title = element_text(size = 18),                               
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  
    plot.margin = margin(20, 20, 20, 20)                                
  ) +
  scale_fill_brewer(palette = "Dark2")


ggsave(
  filename = "./CSHQ_plots/boxplot_expected_CSHQ_subscale_total_difference.png",
  plot = boxplot_expected_CSHQ_subscale_total_difference,
  width = 16, height = 9, units = "in",
  dpi = 600
)

ggsave(
  filename = "./CSHQ_plots/boxplot_expected_CSHQ_subscale_total_difference.pdf",
  plot = boxplot_expected_CSHQ_subscale_total_difference,
  width = 16, height = 9
)






#Plotting the predictive pmfs derived from item probabilities and regression coefficients
#along with stacked plots for group membership, and specified quantile based CI for each score.

combined_plot_density_predictive_T_autozoom_80_CI <- plot_combined_predictive_T(
  theta = CSHQ_LCR_varsel$itemprob,
  beta = CSHQ_ASD_beta_estimate_varsel_with_intercept,
  theta_samples = theta_samples_collapsed_estimate_perm,
  beta_samples = CSHQ_LCR_varsel$samples$beta_samples[c(1,4),,],
  auto_ylim = TRUE,
  ci_level = 0.80
)

combined_plot_density_predictive_T_autozoom_90_CI <- plot_combined_predictive_T(
  theta = CSHQ_LCR_varsel$itemprob,
  beta = CSHQ_ASD_beta_estimate_varsel_with_intercept,
  theta_samples = theta_samples_collapsed_estimate_perm,
  beta_samples = CSHQ_LCR_varsel$samples$beta_samples[c(1,4),,],
  auto_ylim = TRUE,
  ci_level = 0.90
)

combined_plot_density_predictive_T_autozoom_95_CI <- plot_combined_predictive_T(
  theta = CSHQ_LCR_varsel$itemprob,
  beta = CSHQ_ASD_beta_estimate_varsel_with_intercept,
  theta_samples = theta_samples_collapsed_estimate_perm,
  beta_samples = CSHQ_LCR_varsel$samples$beta_samples[c(1,4),,],
  auto_ylim = TRUE,
  ci_level = 0.95
)


#Creating a plot with the labels changed to neurodiverse and neurotypical

combined_plot_density_predictive_T_autozoom_95_CI_label_change <- plot_combined_predictive_T_label_change(
  theta = CSHQ_LCR_varsel$itemprob,
  beta = CSHQ_ASD_beta_estimate_varsel_with_intercept,
  theta_samples = theta_samples_collapsed_estimate_perm[,,,],
  beta_samples = CSHQ_LCR_varsel$samples$beta_samples[c(1,4),,],
  auto_ylim = TRUE,
  ci_level = 0.95,
  main_title = 'Posterior Predictive Distributions of Total CSHQ Score by Neurotype',
  plot_headings = c('Neurotypical', 'Neurodivergent'),
  legend_title = 'Neurotype'
)






ggsave(
  filename = "./CSHQ_plots/combined_plot_density_predictive_T_autozoom_80_CI.png",
  plot = combined_plot_density_predictive_T_autozoom_80_CI,
  width = 16, height = 9, units = "in",
  dpi = 600
)

ggsave(
  filename = "./CSHQ_plots/combined_plot_density_predictive_T_autozoom_80_CI.pdf",
  plot = combined_plot_density_predictive_T_autozoom_80_CI,
  width = 16, height = 9
)


ggsave(
  filename = "./CSHQ_plots/combined_plot_density_predictive_T_autozoom_90_CI.png",
  plot = combined_plot_density_predictive_T_autozoom_90_CI,
  width = 16, height = 9, units = "in",
  dpi = 600
)

ggsave(
  filename = "./CSHQ_plots/combined_plot_density_predictive_T_autozoom_90_CI.pdf",
  plot = combined_plot_density_predictive_T_autozoom_90_CI,
  width = 16, height = 9
)


ggsave(
  filename = "./CSHQ_plots/combined_plot_density_predictive_T_autozoom_95_CI.png",
  plot = combined_plot_density_predictive_T_autozoom_95_CI,
  width = 16, height = 9, units = "in",
  dpi = 600
)

ggsave(
  filename = "./CSHQ_plots/combined_plot_density_predictive_T_autozoom_95_CI.pdf",
  plot = combined_plot_density_predictive_T_autozoom_95_CI,
  width = 16, height = 9
)



ggsave(
  filename = "./CSHQ_plots/combined_plot_density_predictive_T_autozoom_95_CI_label_change.pdf",
  plot = combined_plot_density_predictive_T_autozoom_95_CI_label_change,
  width = 16, height = 9
)



