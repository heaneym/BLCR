# ---- Simulation 2 subsample sensitivity analysis ----
# Uses pre-saved exact subsample indices for exact reproducibility.

cat("\n===============================================\n")
cat("Running Simulation 2 subsample study\n")
cat("===============================================\n\n")

if (!exists("sim2_data")) {
  stop("sim2_data not found. Run the main Simulation 2 script first.")
}

indices_path <- file.path("data", "subsample_indices.rds")
if (!file.exists(indices_path)) {
  stop("Missing exact subsample index file: ", indices_path)
}


sim2_subsample_indices <- readRDS(indices_path)

subsample_sizes <- c(150, 200, 300, 400)
n_reps <- 10
fit_seed <- 128

gibbs_args <- list(
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

fit_sim2_subsample_varsel <- function(X_sub, Y_sub, seed = 128) {
  set.seed(seed)
  do.call(
    LCR_Gibbs,
    c(list(X = X_sub, Y = Y_sub), gibbs_args)
  )
}

subsample_fit_results <- list()
predictor_inclusion_summary <- list()

for (n in subsample_sizes) {
  cat("\n--- Subsample size:", n, "---\n")
  
  idx_list <- sim2_subsample_indices[[paste0("n", n)]]
  if (length(idx_list) != n_reps) {
    stop("Expected ", n_reps, " replicates for subsample size ", n)
  }
  
  fit_list_n <- vector("list", length = n_reps)
  
  for (r in seq_len(n_reps)) {
    cat("Fitting replicate", r, "of", n_reps, "for n =", n, "\n")
    
    idx <- idx_list[[r]]
    X_sub <- sim2_data$X[idx, , drop = FALSE]
    Y_sub <- sim2_data$Y[idx, , drop = FALSE]
    
    fit_list_n[[r]] <- fit_sim2_subsample_varsel(
      X_sub = X_sub,
      Y_sub = Y_sub,
      seed = fit_seed
    )
  }
  
  names(fit_list_n) <- paste0("rep", seq_len(n_reps))
  subsample_fit_results[[paste0("n", n)]] <- fit_list_n
  
  cov_prob_mat <- do.call(
    cbind,
    lapply(fit_list_n, function(fit) fit$cov_inclusion_prob)
  )
  
  rownames(cov_prob_mat) <- paste0("X", seq_len(nrow(cov_prob_mat)))
  colnames(cov_prob_mat) <- paste0("rep", seq_len(ncol(cov_prob_mat)))
  
  mean_probs <- rowMeans(cov_prob_mat)
  
  predictor_inclusion_summary[[paste0("n", n)]] <- list(
    cov_prob_mat = cov_prob_mat,
    mean_probs = mean_probs
  )
  
  write.csv(
    data.frame(
      Predictor = names(mean_probs),
      MeanInclusionProbability = unname(mean_probs)
    ),
    file.path("output", "tables", paste0("sim2_subsample_mean_predictor_inclusion_n", n, ".csv")),
    row.names = FALSE
  )
}

combined_predictor_summary <- do.call(
  rbind,
  lapply(subsample_sizes, function(n) {
    mean_probs <- predictor_inclusion_summary[[paste0("n", n)]]$mean_probs
    data.frame(
      SubsampleSize = n,
      Predictor = names(mean_probs),
      MeanInclusionProbability = unname(mean_probs)
    )
  })
)

write.csv(
  combined_predictor_summary,
  file.path("output", "tables", "sim2_subsample_mean_predictor_inclusion_all_sizes.csv"),
  row.names = FALSE
)

saveRDS(
  subsample_fit_results,
  file.path("output", "rds", "sim2_subsample_varsel_fits.rds")
)

saveRDS(
  predictor_inclusion_summary,
  file.path("output", "rds", "sim2_subsample_predictor_inclusion_summary.rds")
)

cat("\nSimulation 2 subsample study completed.\n")