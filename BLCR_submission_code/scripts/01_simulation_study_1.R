# ---- Simulation Study 1 ----

cat("\n=============================\n")
cat("Running Simulation Study 1\n")
cat("=============================\n\n")


G <- 2
p <- 6
M <- 8
K <- rep(3, 8)

sim1_beta <- matrix(c(-0.5, 0.8, 1.2, 1, 0.4, 0, 0), nrow = p + 1)

sim1_theta1 <- matrix(c(0.15, 0.25, 0.6, 0.7, 0.2, 0.1), nrow = G, ncol = K[1], byrow = TRUE)
sim1_theta2 <- matrix(c(0.2, 0.35, 0.45, 0.55, 0.3, 0.15), nrow = G, ncol = K[2], byrow = TRUE)
sim1_theta3 <- matrix(c(0.1, 0.15, 0.75, 0.8, 0.15, 0.05), nrow = G, ncol = K[3], byrow = TRUE)
sim1_theta4 <- matrix(c(0.25, 0.4, 0.35, 0.45, 0.35, 0.2), nrow = G, ncol = K[4], byrow = TRUE)
sim1_theta5 <- matrix(c(0.4, 0.5, 0.1, 0.4, 0.5, 0.1), nrow = G, ncol = K[5], byrow = TRUE)
sim1_theta6 <- matrix(c(0.7, 0.1, 0.2, 0.7, 0.1, 0.2), nrow = G, ncol = K[6], byrow = TRUE)
sim1_theta7 <- matrix(c(0.1, 0.5, 0.4, 0.1, 0.5, 0.4), nrow = G, ncol = K[7], byrow = TRUE)
sim1_theta8 <- matrix(1 / K[8], nrow = G, ncol = K[8])

sim1_theta <- list(
  sim1_theta1, sim1_theta2, sim1_theta3, sim1_theta4,
  sim1_theta5, sim1_theta6, sim1_theta7, sim1_theta8
)

# ---- Simulate data ----
set.seed(123)
sim1_data <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 500)
sim1_truelabels <- sim1_data$class

# ---- True-parameter benchmark ----
true_posterior_cluster_probs_sim1 <- LCR_posterior_membership_prob(
  theta = sim1_theta,
  beta = cbind(sim1_beta, 0),
  X = sim1_data$X,
  Y = sim1_data$Y
)

sim1_true_model_ari <- adj.rand.index(
  sim1_truelabels,
  max.col(true_posterior_cluster_probs_sim1)
)

cat("Adjusted Rand Index for true model parameters:", sim1_true_model_ari, "\n")

write.csv(
  data.frame(metric = "ARI_true_model", value = sim1_true_model_ari),
  file.path("output", "tables", "sim1_true_model_benchmark.csv"),
  row.names = FALSE
)


# ---- LCR fit without variable selection ----
set.seed(125)

sim1_LCR_fit_no_varsel <- LCR_Gibbs(
  X = sim1_data$X,
  Y = sim1_data$Y,
  G = 2,
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

sim1_beta_summary_table_no_varsel <- multinomial_logit_coefficient_summary_table_MCMC(
  sim1_LCR_fit_no_varsel$samples$beta_samples
)
print(sim1_beta_summary_table_no_varsel)

write.csv(
  sim1_beta_summary_table_no_varsel,
  file.path("output", "tables", "sim1_beta_summary_no_varsel.csv"),
  row.names = FALSE
)

# ---- Simulation 1 beta ridgeline plot ----
sim1_beta_samples <- sim1_LCR_fit_no_varsel$samples$beta_samples[, 2, ]
sim1_beta_names <- paste0("beta", 0:6)

long_beta_samples_data_sim1 <- prepare_beta_long_df(
  beta_samples = sim1_beta_samples,
  beta_names = sim1_beta_names
)

true_values_df_sim1 <- prepare_true_beta_df(
  true_beta = sim1_beta[, 1],
  beta_names = sim1_beta_names
)

math_labels_sim1 <- expression(
  beta[6], beta[5], beta[4], beta[3], beta[2], beta[1], beta[0]
)

final_ridgeline_plot_sim1 <- plot_beta_ridgeline(
  long_df = long_beta_samples_data_sim1,
  true_df = true_values_df_sim1,
  math_labels = math_labels_sim1,
  fill_col = cb_blue,
  true_col = cb_orange,
  xlim = c(-3, 3),
  title = "Simulation 1: Posterior distributions of regression coefficients",
  subtitle = "Posterior density with median & 95% HDI. Diamonds mark true values."
)

ggsave(
  file.path("output", "figures", "sim1_beta_posterior_ridgeline.pdf"),
  plot = final_ridgeline_plot_sim1,
  width = 8,
  height = 3.5
)

# ---- No-selection clustering summary ----
sim1_cluster_mat_no_varsel <- apply(
  sim1_LCR_fit_no_varsel$samples$z_samples,
  3,
  function(matrix_slice) apply(matrix_slice, 1, which.max)
)

sim1_psm_mat_no_varsel <- comp.psm(t(sim1_cluster_mat_no_varsel))
sim1_cluster_minVI_no_varsel <- minVI(
  psm = sim1_psm_mat_no_varsel,
  method = "greedy",
  start.cl = max.col(sim1_LCR_fit_no_varsel$Z)
)

sim1_cross_no_varsel <- table(sim1_truelabels, sim1_cluster_minVI_no_varsel$cl)
sim1_ari_no_varsel <- adj.rand.index(sim1_truelabels, sim1_cluster_minVI_no_varsel$cl)

cat(
  "Cross-classification table for true labels vs. LCR model (no variable selection):\n",
  paste(capture.output(print(sim1_cross_no_varsel)), collapse = "\n"),
  "\n"
)
cat("Adjusted Rand Index for LCR model:", sim1_ari_no_varsel, "\n")

write.csv(
  as.data.frame.matrix(sim1_cross_no_varsel),
  file.path("output", "tables", "sim1_cross_classification_no_varsel.csv")
)

write.csv(
  data.frame(metric = "ARI_no_varsel", value = sim1_ari_no_varsel),
  file.path("output", "tables", "sim1_no_varsel_ari.csv"),
  row.names = FALSE
)

# ---- LCR fit with item selection ----
set.seed(126)

sim1_LCR_fit_itemsel <- LCR_Gibbs(
  X = sim1_data$X,
  Y = sim1_data$Y,
  G = 2,
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

sim1_item_inclusion_probs_itemsel <- sim1_LCR_fit_itemsel$item_inclusion_prob
names(sim1_item_inclusion_probs_itemsel) <- paste0("Y", seq_along(sim1_item_inclusion_probs_itemsel))

cat(
  "Posterior inclusion probabilities for item variables:\n",
  paste(names(sim1_item_inclusion_probs_itemsel), collapse = "    "), "\n",
  paste(sprintf("%.4f", sim1_item_inclusion_probs_itemsel), collapse = "  "), "\n"
)

write.csv(
  data.frame(
    Item = names(sim1_item_inclusion_probs_itemsel),
    InclusionProbability = sim1_item_inclusion_probs_itemsel
  ),
  file.path("output", "tables", "sim1_item_inclusion_probabilities.csv"),
  row.names = FALSE
)

sim1_LCR_no_varsel_vs_itemsel_mean_abs_diff <- mean(
  abs(sim1_LCR_fit_no_varsel$beta_estimate - sim1_LCR_fit_itemsel$beta_estimate)
)

write.csv(
  data.frame(
    metric = "mean_abs_diff_beta_no_varsel_vs_itemsel",
    value = sim1_LCR_no_varsel_vs_itemsel_mean_abs_diff
  ),
  file.path("output", "tables", "sim1_no_varsel_vs_itemsel_beta_difference.csv"),
  row.names = FALSE
)

sim1_cluster_mat_itemsel <- apply(
  sim1_LCR_fit_itemsel$samples$z_samples,
  3,
  function(matrix_slice) apply(matrix_slice, 1, which.max)
)

sim1_psm_mat_itemsel <- comp.psm(t(sim1_cluster_mat_itemsel))
sim1_cluster_minVI_itemsel <- minVI(
  psm = sim1_psm_mat_itemsel,
  method = "greedy",
  start.cl = max.col(sim1_LCR_fit_itemsel$Z)
)

sim1_cross_itemsel <- table(sim1_truelabels, sim1_cluster_minVI_itemsel$cl)
sim1_ari_itemsel <- adj.rand.index(sim1_truelabels, sim1_cluster_minVI_itemsel$cl)

cat(
  "Cross-classification table for true labels vs. LCR model with item selection:\n",
  paste(capture.output(print(sim1_cross_itemsel)), collapse = "\n"),
  "\n"
)
cat("Adjusted Rand Index for LCR model with item selection:", sim1_ari_itemsel, "\n")

write.csv(
  as.data.frame.matrix(sim1_cross_itemsel),
  file.path("output", "tables", "sim1_cross_classification_itemsel.csv")
)

write.csv(
  data.frame(metric = "ARI_itemsel", value = sim1_ari_itemsel),
  file.path("output", "tables", "sim1_itemsel_ari.csv"),
  row.names = FALSE
)

# ---- LCR fit with predictor selection ----
set.seed(127)

sim1_LCR_fit_covsel <- LCR_Gibbs(
  X = sim1_data$X,
  Y = sim1_data$Y,
  G = 2,
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

sim1_predictor_inclusion_probs_covsel <- sim1_LCR_fit_covsel$cov_inclusion_prob
names(sim1_predictor_inclusion_probs_covsel) <- paste0("X", seq_along(sim1_predictor_inclusion_probs_covsel))

cat(
  "Posterior inclusion probabilities for predictor variables:\n",
  paste(names(sim1_predictor_inclusion_probs_covsel), collapse = "    "), "\n",
  paste(sprintf("%.4f", sim1_predictor_inclusion_probs_covsel), collapse = "  "), "\n"
)

write.csv(
  data.frame(
    Predictor = names(sim1_predictor_inclusion_probs_covsel),
    InclusionProbability = sim1_predictor_inclusion_probs_covsel
  ),
  file.path("output", "tables", "sim1_predictor_inclusion_probabilities.csv"),
  row.names = FALSE
)

# ---- LCR fit with simultaneous variable selection ----

set.seed(128)

sim1_LCR_fit_varsel <- LCR_Gibbs(
  X = sim1_data$X,
  Y = sim1_data$Y,
  G = 2,
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

sim1_item_inclusion_probs_varsel <- sim1_LCR_fit_varsel$item_inclusion_prob
names(sim1_item_inclusion_probs_varsel) <- paste0("Y", seq_along(sim1_item_inclusion_probs_varsel))

sim1_predictor_inclusion_probs_varsel <- sim1_LCR_fit_varsel$cov_inclusion_prob
names(sim1_predictor_inclusion_probs_varsel) <- paste0("X", seq_along(sim1_predictor_inclusion_probs_varsel))

cat(
  "Posterior inclusion probabilities for item variables (simultaneous selection):\n",
  paste(names(sim1_item_inclusion_probs_varsel), collapse = "    "), "\n",
  paste(sprintf("%.4f", sim1_item_inclusion_probs_varsel), collapse = "  "), "\n"
)

cat(
  "Posterior inclusion probabilities for predictor variables (simultaneous selection):\n",
  paste(names(sim1_predictor_inclusion_probs_varsel), collapse = "    "), "\n",
  paste(sprintf("%.4f", sim1_predictor_inclusion_probs_varsel), collapse = "  "), "\n"
)

write.csv(
  data.frame(
    Item = names(sim1_item_inclusion_probs_varsel),
    InclusionProbability = sim1_item_inclusion_probs_varsel
  ),
  file.path("output", "tables", "sim1_item_inclusion_probabilities_varsel.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Predictor = names(sim1_predictor_inclusion_probs_varsel),
    InclusionProbability = sim1_predictor_inclusion_probs_varsel
  ),
  file.path("output", "tables", "sim1_predictor_inclusion_probabilities_varsel.csv"),
  row.names = FALSE
)

sim1_cluster_mat_varsel <- apply(
  sim1_LCR_fit_varsel$samples$z_samples,
  3,
  function(matrix_slice) apply(matrix_slice, 1, which.max)
)

sim1_psm_mat_varsel <- comp.psm(t(sim1_cluster_mat_varsel))
sim1_cluster_minVI_varsel <- minVI(
  psm = sim1_psm_mat_varsel,
  method = "greedy",
  start.cl = max.col(sim1_LCR_fit_varsel$Z)
)

sim1_cross_varsel_true <- table(sim1_truelabels, sim1_cluster_minVI_varsel$cl)
sim1_ari_varsel_true <- adj.rand.index(sim1_truelabels, sim1_cluster_minVI_varsel$cl)

cat(
  "Cross-classification table for true labels vs. simultaneous-selection LCR:\n",
  paste(capture.output(print(sim1_cross_varsel_true)), collapse = "\n"),
  "\n"
)
cat(
  "Adjusted Rand Index for true labels vs. simultaneous-selection LCR:",
  sim1_ari_varsel_true, "\n"
)

write.csv(
  as.data.frame.matrix(sim1_cross_varsel_true),
  file.path("output", "tables", "sim1_cross_classification_varsel_vs_true.csv")
)

write.csv(
  data.frame(metric = "ARI_varsel_vs_true", value = sim1_ari_varsel_true),
  file.path("output", "tables", "sim1_varsel_ari_vs_true.csv"),
  row.names = FALSE
)

sim1_cross_no_varsel_vs_varsel <- table(
  sim1_cluster_minVI_no_varsel$cl,
  sim1_cluster_minVI_varsel$cl
)
sim1_ari_no_varsel_vs_varsel <- adj.rand.index(
  sim1_cluster_minVI_no_varsel$cl,
  sim1_cluster_minVI_varsel$cl
)

cat(
  "Cross-classification table for no-selection vs. simultaneous-selection cluster labels:\n",
  paste(capture.output(print(sim1_cross_no_varsel_vs_varsel)), collapse = "\n"),
  "\n"
)
cat(
  "Adjusted Rand Index for no-selection vs. simultaneous-selection:",
  sim1_ari_no_varsel_vs_varsel, "\n"
)

write.csv(
  as.data.frame.matrix(sim1_cross_no_varsel_vs_varsel),
  file.path("output", "tables", "sim1_cross_classification_no_varsel_vs_varsel.csv")
)

write.csv(
  data.frame(metric = "ARI_no_varsel_vs_varsel", value = sim1_ari_no_varsel_vs_varsel),
  file.path("output", "tables", "sim1_no_varsel_vs_varsel_ari.csv"),
  row.names = FALSE
)

cat("\nSimulation Study 1 completed.\n")