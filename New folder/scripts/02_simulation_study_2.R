# ---- Simulation Study 2 ----

cat("\n=============================\n")
cat("Running Simulation Study 2\n")
cat("=============================\n\n")

# ---- Setup ----
G <- 3
p <- 6
M <- 13
K <- c(2, 2, 2, 3, 3, 3, 4, 4, 3, 3, 5, 5, 5)

sim2_beta <- matrix(
  c(0, 0, 1, -1, -1, 1, 0.5, -0.5, -0.4, 0.4, 0, 0, 0, 0),
  nrow = p + 1,
  byrow = TRUE
)

sim2_theta1  <- matrix(c(0.15, 0.6, 0.8, 0.85, 0.4, 0.2), nrow = G, ncol = K[1])
sim2_theta2  <- matrix(c(0.25, 0.45, 0.7, 0.75, 0.55, 0.3), nrow = G, ncol = K[2])
sim2_theta3  <- matrix(c(0.7, 0.2, 0.65, 0.3, 0.8, 0.35), nrow = G, ncol = K[3])
sim2_theta4  <- matrix(c(0.1, 0.35, 0.7, 0.25, 0.4, 0.2, 0.65, 0.25, 0.1), nrow = G, ncol = K[4])
sim2_theta5  <- matrix(c(0.2, 0.25, 0.65, 0.15, 0.6, 0.25, 0.65, 0.15, 0.1), nrow = G, ncol = K[5])
sim2_theta6  <- matrix(c(0.15, 0.5, 0.75, 0.2, 0.35, 0.15, 0.65, 0.15, 0.1), nrow = G, ncol = K[6])
sim2_theta7  <- matrix(c(0.1, 0.25, 0.6, 0.15, 0.35, 0.25, 0.25, 0.25, 0.1, 0.5, 0.15, 0.05), nrow = G, ncol = K[7])
sim2_theta8  <- matrix(c(0.15, 0.2, 0.55, 0.2, 0.45, 0.2, 0.2, 0.25, 0.15, 0.45, 0.1, 0.1), nrow = G, ncol = K[8])
sim2_theta9  <- matrix(c(0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1), nrow = G, ncol = K[9])
sim2_theta10 <- matrix(c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2), nrow = G, ncol = K[10])
sim2_theta11 <- matrix(1 / K[11], nrow = G, ncol = K[11])
sim2_theta12 <- matrix(c(0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.2, 0.2, 0.2, 0.25, 0.25, 0.25, 0.3, 0.3, 0.3), nrow = G, ncol = K[12])
sim2_theta13 <- matrix(c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = G, ncol = K[13])

sim2_theta <- list(
  sim2_theta1, sim2_theta2, sim2_theta3, sim2_theta4, sim2_theta5,
  sim2_theta6, sim2_theta7, sim2_theta8, sim2_theta9, sim2_theta10,
  sim2_theta11, sim2_theta12, sim2_theta13
)

# ---- Simulate data ----
set.seed(129)
sim2_data <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 500)
sim2_truelabels <- sim2_data$class

# ---- True-parameter benchmark ----
true_posterior_cluster_probs_sim2 <- LCR_posterior_membership_prob(
  theta = sim2_theta,
  beta = cbind(sim2_beta, 0),
  X = sim2_data$X,
  Y = sim2_data$Y
)

sim2_true_model_ari <- adj.rand.index(
  sim2_truelabels,
  max.col(true_posterior_cluster_probs_sim2)
)

cat("Adjusted Rand Index for true model parameters:", sim2_true_model_ari, "\n")

write.csv(
  data.frame(metric = "ARI_true_model", value = sim2_true_model_ari),
  file.path("output", "tables", "sim2_true_model_benchmark.csv"),
  row.names = FALSE
)

# ---- LCR fit without variable selection ----
set.seed(131)

sim2_LCR_fit_no_varsel <- LCR_Gibbs(
  X = sim2_data$X,
  Y = sim2_data$Y,
  G = 3,
  beta_prior_cov = diag(10^2, 7),
  beta_prior_mean = rep(0, 7),
  theta_hyperparam = 1,
  clust_var_prior = 0.5,
  item.sel = FALSE,
  cov.sel = FALSE,
  verbose = TRUE,
  relabel = TRUE,
  n_samples = 5000,
  burnin = 1000,
  thinby = 10
)

sim2_beta_summary_table_no_varsel <- multinomial_logit_coefficient_summary_table_MCMC(
  sim2_LCR_fit_no_varsel$samples$beta_samples
)
print(sim2_beta_summary_table_no_varsel)

write.csv(
  sim2_beta_summary_table_no_varsel,
  file.path("output", "tables", "sim2_beta_summary_no_varsel.csv"),
  row.names = FALSE
)

sim2_cluster_mat_no_varsel <- apply(
  sim2_LCR_fit_no_varsel$samples$z_samples,
  3,
  function(matrix_slice) apply(matrix_slice, 1, which.max)
)

sim2_psm_mat_no_varsel <- comp.psm(t(sim2_cluster_mat_no_varsel))
sim2_cluster_minVI_no_varsel <- minVI(
  psm = sim2_psm_mat_no_varsel,
  method = "greedy",
  start.cl = max.col(sim2_LCR_fit_no_varsel$Z)
)

sim2_cross_no_varsel <- table(sim2_truelabels, sim2_cluster_minVI_no_varsel$cl)
cat(
  "Cross-classification table for true labels vs. LCR model (no variable selection):\n",
  paste(capture.output(print(sim2_cross_no_varsel)), collapse = "\n"),
  "\n"
)

write.csv(
  as.data.frame.matrix(sim2_cross_no_varsel),
  file.path("output", "tables", "sim2_cross_classification_no_varsel.csv")
)

# ---- LCR fit with item selection ----
set.seed(132)

sim2_LCR_fit_itemsel <- LCR_Gibbs(
  X = sim2_data$X,
  Y = sim2_data$Y,
  G = 3,
  beta_prior_cov = diag(10^2, 7),
  beta_prior_mean = rep(0, 7),
  theta_hyperparam = 1,
  clust_var_prior = 0.5,
  item.sel = TRUE,
  cov.sel = FALSE,
  verbose = TRUE,
  relabel = TRUE,
  n_samples = 5000,
  burnin = 1000,
  thinby = 10
)

sim2_item_inclusion_probs_itemsel <- sim2_LCR_fit_itemsel$item_inclusion_prob
names(sim2_item_inclusion_probs_itemsel) <- paste0("Y", seq_along(sim2_item_inclusion_probs_itemsel))

cat(
  "Posterior inclusion probabilities for item variables:\n",
  paste(names(sim2_item_inclusion_probs_itemsel), collapse = "    "), "\n",
  paste(sprintf("%.4f", sim2_item_inclusion_probs_itemsel), collapse = "  "), "\n"
)

write.csv(
  data.frame(
    Item = names(sim2_item_inclusion_probs_itemsel),
    InclusionProbability = sim2_item_inclusion_probs_itemsel
  ),
  file.path("output", "tables", "sim2_item_inclusion_probabilities.csv"),
  row.names = FALSE
)

# Relevel beta samples for consistency with simulated values
relevelled_beta_samples <- sim2_LCR_fit_itemsel$samples$beta_samples
for (sample in seq_len(dim(relevelled_beta_samples)[3])) {
  relevelled_beta_samples[, , sample] <-
    relevelled_beta_samples[, , sample] - relevelled_beta_samples[, 3, sample]
}

sim2_beta_summary_table_itemsel <- multinomial_logit_coefficient_summary_table_MCMC(
  relevelled_beta_samples
)
print(sim2_beta_summary_table_itemsel)

write.csv(
  sim2_beta_summary_table_itemsel,
  file.path("output", "tables", "sim2_beta_summary_itemsel.csv"),
  row.names = FALSE
)

# ---- Simulation 2 combined beta ridgeline plot ----
sim2_beta1_samples <- relevelled_beta_samples[, 2, ]
sim2_beta2_samples <- relevelled_beta_samples[, 1, ]
sim2_beta_names <- paste0("beta", 0:6)

long_beta1_samples_data_sim2 <- prepare_beta_long_df(
  beta_samples = sim2_beta1_samples,
  beta_names = sim2_beta_names,
  group = "Group 1"
)

long_beta2_samples_data_sim2 <- prepare_beta_long_df(
  beta_samples = sim2_beta2_samples,
  beta_names = sim2_beta_names,
  group = "Group 2"
)

true_values_beta1_df_sim2 <- prepare_true_beta_df(
  true_beta = sim2_beta[, 1],
  beta_names = sim2_beta_names,
  group = "Group 1"
)

true_values_beta2_df_sim2 <- prepare_true_beta_df(
  true_beta = sim2_beta[, 2],
  beta_names = sim2_beta_names,
  group = "Group 2"
)

custom_labels_sim2 <- rev(c(
  expression(beta[20]), expression(beta[21]), expression(beta[22]),
  expression(beta[23]), expression(beta[24]), expression(beta[25]),
  expression(beta[26]), expression(beta[10]), expression(beta[11]),
  expression(beta[12]), expression(beta[13]), expression(beta[14]),
  expression(beta[15]), expression(beta[16])
))

sim2_combined_beta_ridgeline_plot <- plot_combined_beta_ridgeline(
  long_dfs = list(long_beta1_samples_data_sim2, long_beta2_samples_data_sim2),
  true_dfs = list(true_values_beta1_df_sim2, true_values_beta2_df_sim2),
  group_order = c("Group 1", "Group 2"),
  custom_labels = custom_labels_sim2,
  fill_cols = c("Group 1" = cb_blue, "Group 2" = cb_orange),
  xlim = c(-3, 3),
  title = "Simulation 2: Posterior distributions of regression coefficients by group",
  subtitle = "Posterior density with median and 95% HDI. Diamonds mark true values."
)

ggsave(
  file.path("output", "figures", "sim2_beta_posterior_ridgeline_combined.pdf"),
  plot = sim2_combined_beta_ridgeline_plot,
  width = 8,
  height = 10
)

# ---- Item-selection clustering summary ----
sim2_cluster_mat_itemsel <- apply(
  sim2_LCR_fit_itemsel$samples$z_samples,
  3,
  function(matrix_slice) apply(matrix_slice, 1, which.max)
)

sim2_psm_mat_itemsel <- comp.psm(t(sim2_cluster_mat_itemsel))
sim2_cluster_minVI_itemsel <- minVI(
  psm = sim2_psm_mat_itemsel,
  method = "greedy",
  start.cl = max.col(sim2_LCR_fit_itemsel$Z)
)

sim2_cross_itemsel <- table(sim2_truelabels, sim2_cluster_minVI_itemsel$cl)
sim2_ari_itemsel <- adj.rand.index(sim2_truelabels, sim2_cluster_minVI_itemsel$cl)

cat(
  "Cross-classification table for true labels vs. LCR model with item selection:\n",
  paste(capture.output(print(sim2_cross_itemsel)), collapse = "\n"),
  "\n"
)
cat("Adjusted Rand Index for LCR model with item selection:", sim2_ari_itemsel, "\n")

write.csv(
  as.data.frame.matrix(sim2_cross_itemsel),
  file.path("output", "tables", "sim2_cross_classification_itemsel.csv")
)

write.csv(
  data.frame(metric = "ARI_item_selection", value = sim2_ari_itemsel),
  file.path("output", "tables", "sim2_item_selection_ari.csv"),
  row.names = FALSE
)

# ---- LCR fit with predictor selection ----
set.seed(133)

sim2_LCR_fit_covsel <- LCR_Gibbs(
  X = sim2_data$X,
  Y = sim2_data$Y,
  G = 3,
  beta_prior_cov = diag(10^2, 7),
  beta_prior_mean = rep(0, 7),
  theta_hyperparam = 1,
  clust_var_prior = 0.5,
  item.sel = FALSE,
  cov.sel = TRUE,
  verbose = TRUE,
  relabel = TRUE,
  n_samples = 5000,
  burnin = 1000,
  thinby = 10
)

sim2_predictor_inclusion_probs_covsel <- sim2_LCR_fit_covsel$cov_inclusion_prob
names(sim2_predictor_inclusion_probs_covsel) <- paste0("X", seq_along(sim2_predictor_inclusion_probs_covsel))

cat(
  "Posterior inclusion probabilities for predictor variables:\n",
  paste(names(sim2_predictor_inclusion_probs_covsel), collapse = "    "), "\n",
  paste(sprintf("%.4f", sim2_predictor_inclusion_probs_covsel), collapse = "  "), "\n"
)

write.csv(
  data.frame(
    Predictor = names(sim2_predictor_inclusion_probs_covsel),
    InclusionProbability = sim2_predictor_inclusion_probs_covsel
  ),
  file.path("output", "tables", "sim2_predictor_inclusion_probabilities.csv"),
  row.names = FALSE
)

# ---- LCR fit with simultaneous variable selection ----

set.seed(128)

sim2_LCR_fit_varsel <- LCR_Gibbs(
  X = sim2_data$X,
  Y = sim2_data$Y,
  G = 3,
  beta_prior_cov = diag(10^2, 7),
  beta_prior_mean = rep(0, 7),
  theta_hyperparam = 1,
  clust_var_prior = 0.5,
  item.sel = TRUE,
  cov.sel = TRUE,
  verbose = TRUE,
  relabel = TRUE,
  n_samples = 5000,
  burnin = 1000,
  thinby = 10
)

sim2_item_inclusion_probs_varsel <- sim2_LCR_fit_varsel$item_inclusion_prob
names(sim2_item_inclusion_probs_varsel) <- paste0("Y", seq_along(sim2_item_inclusion_probs_varsel))

sim2_predictor_inclusion_probs_varsel <- sim2_LCR_fit_varsel$cov_inclusion_prob
names(sim2_predictor_inclusion_probs_varsel) <- paste0("X", seq_along(sim2_predictor_inclusion_probs_varsel))

cat(
  "Posterior inclusion probabilities for item variables (simultaneous selection):\n",
  paste(names(sim2_item_inclusion_probs_varsel), collapse = "    "), "\n",
  paste(sprintf("%.4f", sim2_item_inclusion_probs_varsel), collapse = "  "), "\n"
)

cat(
  "Posterior inclusion probabilities for predictor variables (simultaneous selection):\n",
  paste(names(sim2_predictor_inclusion_probs_varsel), collapse = "    "), "\n",
  paste(sprintf("%.4f", sim2_predictor_inclusion_probs_varsel), collapse = "  "), "\n"
)

write.csv(
  data.frame(
    Item = names(sim2_item_inclusion_probs_varsel),
    InclusionProbability = sim2_item_inclusion_probs_varsel
  ),
  file.path("output", "tables", "sim2_item_inclusion_probabilities_varsel.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Predictor = names(sim2_predictor_inclusion_probs_varsel),
    InclusionProbability = sim2_predictor_inclusion_probs_varsel
  ),
  file.path("output", "tables", "sim2_predictor_inclusion_probabilities_varsel.csv"),
  row.names = FALSE
)

sim2_cluster_mat_varsel <- apply(
  sim2_LCR_fit_varsel$samples$z_samples,
  3,
  function(matrix_slice) apply(matrix_slice, 1, which.max)
)

sim2_psm_mat_varsel <- comp.psm(t(sim2_cluster_mat_varsel))
sim2_cluster_minVI_varsel <- minVI(
  psm = sim2_psm_mat_varsel,
  method = "greedy",
  start.cl = max.col(sim2_LCR_fit_varsel$Z)
)

sim2_cross_varsel_true <- table(sim2_truelabels, sim2_cluster_minVI_varsel$cl)
sim2_ari_varsel_true <- adj.rand.index(sim2_truelabels, sim2_cluster_minVI_varsel$cl)

cat(
  "Cross-classification table for true labels vs. LCR model with simultaneous variable selection:\n",
  paste(capture.output(print(sim2_cross_varsel_true)), collapse = "\n"),
  "\n"
)
cat(
  "Adjusted Rand Index for true labels vs. LCR model with simultaneous variable selection:",
  sim2_ari_varsel_true, "\n"
)

write.csv(
  as.data.frame.matrix(sim2_cross_varsel_true),
  file.path("output", "tables", "sim2_cross_classification_varsel_vs_true.csv")
)

write.csv(
  data.frame(metric = "ARI_varsel_vs_true", value = sim2_ari_varsel_true),
  file.path("output", "tables", "sim2_varsel_ari_vs_true.csv"),
  row.names = FALSE
)

sim2_cross_itemsel_vs_varsel <- table(
  sim2_cluster_minVI_itemsel$cl,
  sim2_cluster_minVI_varsel$cl
)
sim2_ari_itemsel_vs_varsel <- adj.rand.index(
  sim2_cluster_minVI_itemsel$cl,
  sim2_cluster_minVI_varsel$cl
)

cat(
  "Cross-classification table for item-selection vs. simultaneous-selection cluster labels:\n",
  paste(capture.output(print(sim2_cross_itemsel_vs_varsel)), collapse = "\n"),
  "\n"
)
cat(
  "Adjusted Rand Index for item-selection vs. simultaneous-selection:",
  sim2_ari_itemsel_vs_varsel, "\n"
)

write.csv(
  as.data.frame.matrix(sim2_cross_itemsel_vs_varsel),
  file.path("output", "tables", "sim2_cross_classification_itemsel_vs_varsel.csv")
)

write.csv(
  data.frame(metric = "ARI_itemsel_vs_varsel", value = sim2_ari_itemsel_vs_varsel),
  file.path("output", "tables", "sim2_itemsel_vs_varsel_ari.csv"),
  row.names = FALSE
)

cat("\nSimulation Study 2 completed.\n")