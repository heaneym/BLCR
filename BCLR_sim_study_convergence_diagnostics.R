library(posterior)
library(coda)

# Here we want to compute all of the relevant diagnostics for assessing convergence

# In particular, for each of the simulations and for the CSHQ data, we want to report

# multiple chains 
# effective sample sizes
# trace plots
# autocorrelation diagnostics


library(posterior)
library(coda)
library(bayesplot)


relabel_fit_draws <- function(fit, perm, G, K, M) {
  
  n_iter <- dim(fit$samples$beta_samples)[3]
  beta_relab <- array(NA, dim = dim(fit$samples$beta_samples))
  for (i in 1:n_iter) {
    beta_relab[, , i] <- relabel_beta_mat(fit$samples$beta_samples[, , i], perm)
  }
  
  n_iter_theta <- dim(fit$samples$theta_samples)[4]
  theta_relab <- vector("list", M)
  for (j in 1:M) {
    Kj <- K[j]
    arr_j <- fit$samples$theta_samples[, j, 1:Kj, , drop = TRUE]
    if (Kj == 1) dim(arr_j) <- c(G, 1, n_iter_theta)
    arr_relab <- array(NA, dim = dim(arr_j))
    for (i in 1:n_iter_theta) arr_relab[, , i] <- relabel_theta_mat(arr_j[, , i], perm)
    theta_relab[[j]] <- arr_relab
  }
  
  list(beta = beta_relab, theta = theta_relab,
       log_post = fit$log_post_trace)  
}

# ---------------------------------------------------------------------------
# Main diagnostic function
# ---------------------------------------------------------------------------
# fits        : list of n_chains fit objects, each an independent LCR_Gibbs run
#               on the SAME data with different seeds/starting values
# true_z      : reference partition to align every chain's labeling to
#               (use the true simulated z for a simulation study)
# beta_idx    : which (predictor, class-contrast) pairs of beta to diagnose,
#               as a data.frame with columns j, g -- default: all of them
# theta_idx   : which (item, class, category) triples of theta to diagnose,
#               as a data.frame with columns j, g, k -- default: all of them
# make_plots  : if FALSE, skip trace/acf plot generation (fast path for
#               looping diagnostics over many replicates; turn on for spot checks)

LCR_convergence_diagnose <- function(fits, true_z, G, K, p, M,
                                     beta_idx = NULL, theta_idx = NULL,
                                     make_plots = TRUE, max_lag = 40) {
  
  n_chains <- length(fits)
  
  #relabelling runs in terms of common baseline
  relabeled <- vector("list", n_chains)
  for (c in 1:n_chains) {
    perm <- get_relabel_perm(true_labels = true_z,
                             est_labels = max.col(fits[[c]]$assignment_prob),
                             G = G)
    relabeled[[c]] <- relabel_fit_draws(fits[[c]], perm, G = G, K = K, M = M)
  }
  
  n_iter <- length(relabeled[[1]]$log_post)

  if (is.null(beta_idx)) {
    beta_idx <- expand.grid(j = 1:(p + 1), g = 1:(G - 1))
  }
  if (is.null(theta_idx)) {
    theta_idx <- do.call(rbind, lapply(1:M, function(j) {
      expand.grid(j = j, g = 1:G, k = 1:K[j])
    }))
  }
  
  var_names <- c(
    "log_post",
    paste0("beta[", beta_idx$j, ",", beta_idx$g, "]"),
    paste0("theta[", theta_idx$j, ",", theta_idx$g, ",", theta_idx$k, "]")
  )
  
  draws_arr <- array(NA, dim = c(n_iter, n_chains, length(var_names)),
                     dimnames = list(NULL, NULL, var_names))
  
  for (c in 1:n_chains) {
    draws_arr[, c, "log_post"] <- relabeled[[c]]$log_post
    
    for (r in 1:nrow(beta_idx)) {
      j <- beta_idx$j[r]; g <- beta_idx$g[r]
      vname <- paste0("beta[", j, ",", g, "]")
      draws_arr[, c, vname] <- relabeled[[c]]$beta[j, g, ] - relabeled[[c]]$beta[j, G, ]
    }
    
    for (r in 1:nrow(theta_idx)) {
      j <- theta_idx$j[r]; g <- theta_idx$g[r]; k <- theta_idx$k[r]
      vname <- paste0("theta[", j, ",", g, ",", k, "]")
      draws_arr[, c, vname] <- relabeled[[c]]$theta[[j]][g, k, ]
    }
  }
  
  draws_obj <- posterior::as_draws_array(draws_arr)
  
  diag_summary <- posterior::summarise_draws(
    draws_obj, "mean", "sd", "rhat", "ess_bulk", "ess_tail"
  )

  plots <- NULL
  if (make_plots) {
    plots <- list(
      trace_log_post = bayesplot::mcmc_trace(draws_obj, pars = "log_post"),
      trace_beta     = bayesplot::mcmc_trace(draws_obj, regex_pars = "^beta\\["),
      acf_log_post   = bayesplot::mcmc_acf(draws_obj, pars = "log_post", lags = max_lag),
      acf_beta       = bayesplot::mcmc_acf(draws_obj, regex_pars = "^beta\\[", lags = max_lag)
    )
  }
  
  list(summary = diag_summary, draws = draws_obj, plots = plots)
}









#creating list of initial values


# n_chains distinct, overdispersed beta inits, each (p+1) x (G-1)
init_beta_list <- lapply(1:n_chains, function(c) {
  matrix(rnorm((p+1)*(G-1), mean = 0, sd = 3), nrow = p+1, ncol = G-1)
})











