# CSHQ Sleep Patterns Data 

cat("\n=============================\n")
cat("Running CSHQ analysis\n")
cat("=============================\n\n")

# Read data 
CSHQ_df <- read.csv(file.path("data", "CSHQ_df.csv"))
Y_CSHQ <- as.matrix(read.csv(file.path("data", "CSHQ_response_Y.csv")))
X_CSHQ <- as.matrix(read.csv(file.path("data", "CSHQ_predictor_X.csv")))[,-1]

# Collapsed model-selection run 
set.seed(134)

CSHQ_preliminary_collapsed_run <- blca.collapsed(
  X = Y_CSHQ,
  G = 1,
  iter = 50000,
  burn.in = 10000,
  thin = 1 / 10,
  G.sel = TRUE,
  var.sel = TRUE,
  post.hoc.run = TRUE,
  verbose = TRUE,
  relabel = TRUE
)

CSHQ_collapsed_group_number_post_prob <-
  table(CSHQ_preliminary_collapsed_run$samples$G) /
  sum(table(CSHQ_preliminary_collapsed_run$samples$G))

write.csv(
  data.frame(
    Groups = names(CSHQ_collapsed_group_number_post_prob),
    PosteriorProbability = as.numeric(CSHQ_collapsed_group_number_post_prob)
  ),
  file.path("output", "tables", "cshq_collapsed_group_number_posterior_probabilities.csv"),
  row.names = FALSE
)

LCA_item_inclusion_ind_CSHQ <- do.call(rbind, CSHQ_preliminary_collapsed_run$samples$var.ind)
LCA_item_inclusion_prob_CSHQ <- apply(LCA_item_inclusion_ind_CSHQ, 2, mean)

write.csv(
  data.frame(
    Item = colnames(Y_CSHQ),
    InclusionProbability = LCA_item_inclusion_prob_CSHQ
  ),
  file.path("output", "tables", "cshq_lca_item_inclusion_probabilities.csv"),
  row.names = FALSE
)

CSHQ_collapsed_cluster_mat <- do.call(rbind, CSHQ_preliminary_collapsed_run$samples$labels)
CSHQ_collapsed_psm_mat <- comp.psm(CSHQ_collapsed_cluster_mat)

CSHQ_minVI_cluster <- minVI(
  psm = CSHQ_collapsed_psm_mat,
  method = "greedy",
  start.cl = max.col(CSHQ_preliminary_collapsed_run$Z),
  suppress.comment = FALSE
)

CSHQ_cluster_95_cred_ball <- credibleball(
  c.star = CSHQ_minVI_cluster$cl,
  cls.draw = CSHQ_collapsed_cluster_mat,
  c.dist = "VI",
  alpha = 0.05
)

#Fixed 4-group LCA comparison
CSHQ_preliminary_4group_LCA <- blca.collapsed(
  X = Y_CSHQ,
  G = 4,
  iter = 50000,
  burn.in = 10000,
  thin = 1 / 10,
  G.sel = FALSE,
  var.sel = TRUE,
  post.hoc.run = TRUE,
  verbose = TRUE
)

CSHQ_4group_cluster_mat <- CSHQ_preliminary_4group_LCA$samples$labels
CSHQ_4group_psm_mat <- comp.psm(CSHQ_4group_cluster_mat$`G = 4`)

CSHQ_4group_minVI_cluster <- minVI(
  psm = CSHQ_4group_psm_mat,
  method = "greedy",
  suppress.comment = FALSE,
  start.cl = max.col(CSHQ_preliminary_4group_LCA$Z)
)

CSHQ_dist_4group_to_point_estimate <- vi.dist(
  CSHQ_minVI_cluster$cl,
  CSHQ_4group_minVI_cluster$cl
)

table_collapsed_with_4_group <- table(
  CSHQ_minVI_cluster$cl,
  CSHQ_4group_minVI_cluster$cl
)
ari_collapsed_with_4_group <- adj.rand.index(
  CSHQ_minVI_cluster$cl,
  CSHQ_4group_minVI_cluster$cl
)

write.csv(
  as.data.frame.matrix(table_collapsed_with_4_group),
  file.path("output", "tables", "cshq_collapsed_vs_4group_lca_cross_classification.csv")
)

write.csv(
  data.frame(
    metric = c("vi_distance_4group_to_point_estimate", "ari_collapsed_vs_4group"),
    value = c(CSHQ_dist_4group_to_point_estimate, ari_collapsed_with_4_group)
  ),
  file.path("output", "tables", "cshq_lca_clustering_comparison_metrics.csv"),
  row.names = FALSE
)

# LCR with item selection
set.seed(136)

CSHQ_LCR_itemsel <- LCR_Gibbs(
  X = X_CSHQ,
  Y = Y_CSHQ,
  G = 4,
  beta_prior_cov = diag(c(rep(10^2, 4), 5^2, 5^2)),
  beta_prior_mean = rep(0, 6),
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

itemsel_item_inclusion_prob_CSHQ <- CSHQ_LCR_itemsel$item_inclusion_prob
write.csv(
  data.frame(
    Item = colnames(Y_CSHQ),
    InclusionProbability = itemsel_item_inclusion_prob_CSHQ
  ),
  file.path("output", "tables", "cshq_itemsel_item_inclusion_probabilities.csv"),
  row.names = FALSE
)

CSHQ_beta_summary_table_itemsel <- multinomial_logit_coefficient_summary_table_MCMC(
  CSHQ_LCR_itemsel$samples$beta_samples
)
write.csv(
  CSHQ_beta_summary_table_itemsel,
  file.path("output", "tables", "cshq_itemsel_beta_summary.csv"),
  row.names = FALSE
)

CSHQ_LCR_itemsel_cluster_mat <- matrix(
  0,
  nrow = dim(CSHQ_LCR_itemsel$samples$z_samples)[3],
  ncol = nrow(Y_CSHQ)
)

for (slice in seq_len(nrow(CSHQ_LCR_itemsel_cluster_mat))) {
  CSHQ_LCR_itemsel_cluster_mat[slice, ] <- max.col(CSHQ_LCR_itemsel$samples$z_samples[, , slice])
}

CSHQ_itemsel_psm_mat <- comp.psm(CSHQ_LCR_itemsel_cluster_mat)
CSHQ_itemsel_minVI_cluster <- minVI(
  psm = CSHQ_itemsel_psm_mat,
  method = "greedy",
  start.cl = max.col(CSHQ_LCR_itemsel$Z)
)

# LCR with predictor selection
set.seed(137)

CSHQ_LCR_predsel <- LCR_Gibbs(
  X = X_CSHQ,
  Y = Y_CSHQ,
  G = 4,
  beta_prior_cov = diag(c(rep(10^2, 4), 5^2, 5^2)),
  beta_prior_mean = rep(0, 6),
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

predsel_predictor_inclusion_prob_CSHQ <- CSHQ_LCR_predsel$cov_inclusion_prob

predictor_names_CSHQ <- c("Intercept", colnames(X_CSHQ))

write.csv(
  data.frame(
    Predictor = predictor_names_CSHQ,
    InclusionProbability = predsel_predictor_inclusion_prob_CSHQ
  ),
  file.path("output", "tables", "cshq_predsel_predictor_inclusion_probabilities.csv"),
  row.names = FALSE
)

#LCR with simultaneous variable selection
set.seed(138)

CSHQ_LCR_varsel <- LCR_Gibbs(
  X = X_CSHQ,
  Y = Y_CSHQ,
  G = 4,
  beta_prior_cov = diag(c(rep(10^2, 4), 5^2, 5^2)),
  beta_prior_mean = rep(0, 6),
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

varsel_predictor_inclusion_prob_CSHQ <- CSHQ_LCR_varsel$cov_inclusion_prob
varsel_item_inclusion_prob_CSHQ <- CSHQ_LCR_varsel$item_inclusion_prob

write.csv(
  data.frame(
    Item = colnames(Y_CSHQ),
    InclusionProbability = varsel_item_inclusion_prob_CSHQ
  ),
  file.path("output", "tables", "cshq_varsel_item_inclusion_probabilities.csv"),
  row.names = FALSE
)

predictor_names_CSHQ <- c("Intercept", colnames(X_CSHQ))

write.csv(
  data.frame(
    Predictor = predictor_names_CSHQ,
    InclusionProbability = varsel_predictor_inclusion_prob_CSHQ
  ),
  file.path("output", "tables", "cshq_varsel_predictor_inclusion_probabilities.csv"),
  row.names = FALSE
)

CSHQ_LCR_varsel_cluster_mat <- matrix(
  0,
  nrow = dim(CSHQ_LCR_varsel$samples$z_samples)[3],
  ncol = nrow(Y_CSHQ)
)

for (slice in seq_len(nrow(CSHQ_LCR_varsel_cluster_mat))) {
  CSHQ_LCR_varsel_cluster_mat[slice, ] <- max.col(CSHQ_LCR_varsel$samples$z_samples[, , slice])
}

CSHQ_varsel_psm_mat <- comp.psm(CSHQ_LCR_varsel_cluster_mat)
CSHQ_varsel_minVI_cluster <- minVI(
  psm = CSHQ_varsel_psm_mat,
  method = "greedy",
  start.cl = max.col(CSHQ_LCR_varsel$Z)
)

itemsel_vs_varsel_cross_classification_table <- table(
  CSHQ_itemsel_minVI_cluster$cl,
  CSHQ_varsel_minVI_cluster$cl
)
itemsel_vs_varsel_ari <- adj.rand.index(
  CSHQ_itemsel_minVI_cluster$cl,
  CSHQ_varsel_minVI_cluster$cl
)

write.csv(
  as.data.frame.matrix(itemsel_vs_varsel_cross_classification_table),
  file.path("output", "tables", "cshq_itemsel_vs_varsel_cross_classification.csv")
)

write.csv(
  data.frame(metric = "ari_itemsel_vs_varsel", value = itemsel_vs_varsel_ari),
  file.path("output", "tables", "cshq_itemsel_vs_varsel_ari.csv"),
  row.names = FALSE
)

# Mosaic plots for BLCR retained items
cb_palette_mosaic <- c(
  "Response 1" = "#0072B2",
  "Response 2" = "#E69F00",
  "Response 3" = "#D55E00"
)

mosaic_global_theme <- theme(
  text = element_text(size = 16),
  plot.title = element_text(size = 18),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.position = "right",
  legend.key.width = grid::unit(0.5, "cm"),
  legend.key.height = grid::unit(1.2, "cm"),
  panel.grid = element_blank()
)

retained_itemprob <- CSHQ_LCR_varsel$itemprob[
  which(CSHQ_LCR_varsel$item.ind == 1)
]

mosaic_plot_6grid_titles <- c(
  "Bedtime Resistance 5",
  "Sleep Onset Delay",
  "Sleep Duration 3",
  "Sleep Anxiety 1",
  "Sleep Anxiety 3",
  "Night Waking 1"
)

CSHQ_LCR_6grid_mosaic_plot <- build_mosaic_grid(
  itemprob = retained_itemprob[c(3, 4, 7, 8, 10, 12)],
  classprob = CSHQ_LCR_varsel$pi,
  item_titles = mosaic_plot_6grid_titles,
  ncol = 3,
  show_y_axis_numbers = TRUE,
  palette = cb_palette_mosaic,
  base_size = 15,
  global_theme = mosaic_global_theme
)

ggsave(
  file.path("output", "figures", "cshq_lcr_mosaic_6grid.pdf"),
  plot = CSHQ_LCR_6grid_mosaic_plot,
  width = 16,
  height = 9
)

mosaic_plot_grid_full1_names <- c(
  "Bedtime Resistance 2",
  "Bedtime Resistance 3",
  "Bedtime Resistance 5",
  "Sleep Onset Delay",
  "Sleep Duration 1",
  "Sleep Duration 2",
  "Sleep Duration 3",
  "Sleep Anxiety 1",
  "Sleep Anxiety 2",
  "Sleep Anxiety 3",
  "Sleep Anxiety 4",
  "Night Waking 1"
)

CSHQ_LCR_mosaic_full1_plot <- build_mosaic_grid(
  itemprob = retained_itemprob[1:12],
  classprob = CSHQ_LCR_varsel$pi,
  item_titles = mosaic_plot_grid_full1_names,
  ncol = 3,
  show_y_axis_numbers = TRUE,
  palette = cb_palette_mosaic,
  base_size = 15,
  global_theme = mosaic_global_theme
)

ggsave(
  file.path("output", "figures", "cshq_lcr_mosaic_full_grid1.pdf"),
  plot = CSHQ_LCR_mosaic_full1_plot,
  width = 16,
  height = 18
)

mosaic_plot_grid_full2_names <- c(
  "Night Waking 2",
  "Night Waking 3",
  "Parasomnias 2",
  "Parasomnias 3",
  "Parasomnias 5",
  "Daytime Sleepiness 2",
  "Daytime Sleepiness 4",
  "Daytime Sleepiness 5",
  "Daytime Sleepiness 6"
)

CSHQ_LCR_mosaic_full2_plot <- build_mosaic_grid(
  itemprob = retained_itemprob[13:21],
  classprob = CSHQ_LCR_varsel$pi,
  item_titles = mosaic_plot_grid_full2_names,
  ncol = 3,
  show_y_axis_numbers = TRUE,
  palette = cb_palette_mosaic,
  base_size = 15,
  global_theme = mosaic_global_theme
)

ggsave(
  file.path("output", "figures", "cshq_lcr_mosaic_full_grid2.pdf"),
  plot = CSHQ_LCR_mosaic_full2_plot,
  width = 16,
  height = 13.5
)

#Standardised BLCR profile plot
Y_CSHQ_reduced <- Y_CSHQ[, which(CSHQ_LCR_varsel$item.ind == 1)]

CSHQ_LCR_profiles <- compute_cshq_profiles(
  Y_reduced = Y_CSHQ_reduced,
  cluster_labels = CSHQ_varsel_minVI_cluster$cl,
  standardise = TRUE
)

cshq_profile_cols <- c(
  "Group 1" = "#1B9E77",
  "Group 2" = "#D95F02",
  "Group 3" = "#7570B3",
  "Group 4" = "#E7298A"
)

CSHQ_LCR_profile_plot_standardised <- plot_cshq_profiles(
  profile_mat = CSHQ_LCR_profiles$profile_standardised,
  title = "BLCR CSHQ subscale profiles (standardised)",
  subtitle = "Standardised reduced subscale means by group",
  ylab = "Standardised mean score",
  palette = cshq_profile_cols
)

ggsave(
  file.path("output", "figures", "cshq_blcr_profile_plot_standardised.pdf"),
  plot = CSHQ_LCR_profile_plot_standardised,
  width = 8,
  height = 5
)

# ASD-only coefficient summary and stacked ridgelin
ASD_beta_coeff_samples_varsel <- CSHQ_LCR_varsel$samples$beta_samples[4, , , drop = FALSE]
ASD_beta_coeff_samples_varsel_with_intercept <- CSHQ_LCR_varsel$samples$beta_samples[c(1, 4), , , drop = FALSE]

zero_slices <- numeric(dim(CSHQ_LCR_varsel$samples$gamma_samples)[2])
for (slice in seq_along(zero_slices)) {
  zero_slices[slice] <- 1 * ifelse(CSHQ_LCR_varsel$samples$gamma_samples[4, slice] == 0, 1, 0)
}

ASD_beta_coeff_nonzero_samples <- ASD_beta_coeff_samples_varsel[, , which(zero_slices == 0), drop = FALSE]
ASD_beta_coeff_nonzero_samples_with_intercept <- ASD_beta_coeff_samples_varsel_with_intercept[, , which(zero_slices == 0), drop = FALSE]

CSHQ_ASD_beta_estimate_varsel <- apply(ASD_beta_coeff_nonzero_samples, c(1, 2), mean)
CSHQ_ASD_beta_estimate_varsel_with_intercept <- apply(ASD_beta_coeff_nonzero_samples_with_intercept, c(1, 2), mean)
CSHQ_ASD_beta_sd_varsel_with_intercept <- apply(ASD_beta_coeff_nonzero_samples_with_intercept, c(1, 2), sd)

write.csv(
  data.frame(
    Parameter = c("Intercept", "ASD"),
    GroupB_Mean = CSHQ_ASD_beta_estimate_varsel_with_intercept[, 2],
    GroupC_Mean = CSHQ_ASD_beta_estimate_varsel_with_intercept[, 3],
    GroupD_Mean = CSHQ_ASD_beta_estimate_varsel_with_intercept[, 4],
    GroupB_SD = CSHQ_ASD_beta_sd_varsel_with_intercept[, 2],
    GroupC_SD = CSHQ_ASD_beta_sd_varsel_with_intercept[, 3],
    GroupD_SD = CSHQ_ASD_beta_sd_varsel_with_intercept[, 4]
  ),
  file.path("output", "tables", "cshq_asd_only_beta_summary.csv"),
  row.names = FALSE
)

ridgeline_array_ASD_only <- ASD_beta_coeff_nonzero_samples_with_intercept[, 2:4, , drop = FALSE]

CSHQ_betaB_samples_ASD_only <- ridgeline_array_ASD_only[, 1, ]
CSHQ_betaC_samples_ASD_only <- ridgeline_array_ASD_only[, 2, ]
CSHQ_betaD_samples_ASD_only <- ridgeline_array_ASD_only[, 3, ]

CSHQ_beta_names_ASD_only <- paste0("beta", 0:1)

long_betaB_samples_data_CSHQ_ASD_only <- prepare_beta_long_df(
  beta_samples = CSHQ_betaB_samples_ASD_only,
  beta_names = CSHQ_beta_names_ASD_only,
  group = "Group B"
)

long_betaC_samples_data_CSHQ_ASD_only <- prepare_beta_long_df(
  beta_samples = CSHQ_betaC_samples_ASD_only,
  beta_names = CSHQ_beta_names_ASD_only,
  group = "Group C"
)

long_betaD_samples_data_CSHQ_ASD_only <- prepare_beta_long_df(
  beta_samples = CSHQ_betaD_samples_ASD_only,
  beta_names = CSHQ_beta_names_ASD_only,
  group = "Group D"
)

group_colors <- c(
  "Group B" = "#D95F02FF",
  "Group C" = "#7570B3FF",
  "Group D" = "#E7298AFF"
)

groups_order <- c("Group B", "Group C", "Group D")
custom_param_labels_ASD_only <- c("Intercept", "ASD diagnosis")

combined_beta_data_ASD_only <- bind_rows(
  long_betaB_samples_data_CSHQ_ASD_only,
  long_betaC_samples_data_CSHQ_ASD_only,
  long_betaD_samples_data_CSHQ_ASD_only
) %>%
  mutate(
    Group = factor(Group, levels = groups_order),
    Parameter = factor(Parameter, levels = CSHQ_beta_names_ASD_only),
    Parameter_Group = interaction(Parameter, Group, sep = ".", lex.order = FALSE)
  )

param_group_levels <- rev(as.vector(outer(
  CSHQ_beta_names_ASD_only,
  groups_order,
  paste,
  sep = "."
)))
y_labels <- rev(rep(custom_param_labels_ASD_only, times = length(groups_order)))

combined_beta_data_ASD_only$Parameter_Group <- factor(
  combined_beta_data_ASD_only$Parameter_Group,
  levels = param_group_levels
)

final_combined_ridgeline_plot_CSHQ_ASD_only <- ggplot(
  combined_beta_data_ASD_only,
  aes(x = Value, y = Parameter_Group, fill = Group)
) +
  geom_density_ridges(
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    alpha = 0.7,
    scale = 0.9,
    rel_min_height = 0,
    from = -5,
    to = 5
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_fill_manual(values = group_colors) +
  scale_y_discrete(labels = y_labels) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(
    x = "Logit Coefficient",
    y = "",
    title = "Posterior distributions of regression coefficients by group",
    subtitle = "Posterior density with median and 95% HDI.",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top"
  )

ggsave(
  file.path("output", "figures", "cshq_asd_only_beta_posterior_ridgeline_combined.pdf"),
  plot = final_combined_ridgeline_plot_CSHQ_ASD_only,
  width = 8,
  height = 3.75
)

# ASD-effect visualisation
theta_samples_collapsed_estimate <- compute_theta_samples_from_counts(
  N_gjk_samples = CSHQ_LCR_varsel$samples$N_gjk_samples,
  N_g_samples = CSHQ_LCR_varsel$samples$N_g_samples
)

theta_samples_collapsed_estimate_perm <- theta_samples_collapsed_estimate[c(4, 1, 3, 2), , , ]

asd_included_idx <- which(zero_slices == 0)

theta_samples_asd_nonzero <- theta_samples_collapsed_estimate_perm[, , , asd_included_idx, drop = FALSE]
beta_samples_asd_nonzero <- CSHQ_LCR_varsel$samples$beta_samples[, , asd_included_idx, drop = FALSE]

boxplot_expected_CSHQ_subscale_total_difference <- plot_subscale_diff_boxplot(
  theta_samples = theta_samples_asd_nonzero,
  beta_samples = beta_samples_asd_nonzero
)

ggsave(
  file.path("output", "figures", "cshq_subscale_total_difference_asd_vs_non_asd.pdf"),
  plot = boxplot_expected_CSHQ_subscale_total_difference,
  width = 16,
  height = 9
)

predictive_total_cshq_by_asd_plot <- plot_combined_predictive_T(
  theta = CSHQ_LCR_varsel$itemprob,
  beta = CSHQ_ASD_beta_estimate_varsel_with_intercept,
  theta_samples = theta_samples_collapsed_estimate_perm,
  beta_samples = CSHQ_LCR_varsel$samples$beta_samples[c(1, 4), , ],
  auto_ylim = TRUE,
  ci_level = 0.95
)

ggsave(
  file.path("output", "figures", "cshq_total_score_predictive_distribution_by_asd.pdf"),
  plot = predictive_total_cshq_by_asd_plot,
  width = 12,
  height = 8
)

cat("\nCSHQ analysis completed.\n")