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
      trace_beta = bayesplot::mcmc_trace(draws_obj, regex_pars = "^beta\\["),
      acf_log_post = bayesplot::mcmc_acf(draws_obj, pars = "log_post", lags = max_lag),
      acf_beta = bayesplot::mcmc_acf(draws_obj, regex_pars = "^beta\\[", lags = max_lag)
    )
  }
  list(summary = diag_summary, draws = draws_obj, plots = plots)
}









#creating list of initial values

n_chains <- 4

init_params_diagnostics <- function(n_chains = 4, p, G, M, K, n, 
                                    theta_dispersion = c(0.1, 0.5, 5, 50)){
  
  init_beta_list <- lapply(1:n_chains, function(c) {
    matrix(rnorm((p+1)*(G-1), mean = 0, sd = 3), nrow = p+1, ncol = G-1)
  })
  theta_dispersion <- c(0.1, 0.5, 5, 50)
  
  init_theta_list <- lapply(1:n_chains, function(c) {
    theta_init <- array(0, dim = c(G, M, max(K)))
    for (g in 1:G) {
      for (j in 1:M) {
        theta_init[g, j, 1:K[j]] <- rdirichlet(1, rep(theta_dispersion[c], K[j]))
      }
    }
    theta_init
  })
  
  init_nu_list <- lapply(1:n_chains, function(c) {
    if (c == 1){
      rep(0,M)
    } else if (c == n_chains){
      rep(1,M)
    } else {
      sample(c(0,1), M, replace = TRUE)
    }
  })
  
  init_gamma_list <- lapply(1:n_chains, function(c) {
    if (c == 1){
      c(1, rep(0,p))
    } else if (c == n_chains){
      c(1, rep(1,p))
    } else {
      c(1,sample(c(0,1), p, replace = TRUE))
    }
  })
  
  init_z_list <- lapply(1:n_chains, function(c) {
    if (c == 1) {
      Z <- matrix(0, nrow = n, ncol = G)
      Z[, 1] <- 1
    } else {
      labels <- sample(1:G, n, replace = TRUE)
      Z <- matrix(0, nrow = n, ncol = G)
      Z[cbind(1:n, labels)] <- 1
    }
    Z
  })
  
  return(
    list(
    init_beta_list = init_beta_list,
    init_theta_list = init_theta_list,
    init_nu_list = init_nu_list,
    init_gamma_list = init_gamma_list,
    init_z_list = init_z_list
    )
  )
}


# We're gonna run this for simulation 2 and CSHQ data for the moment
# Then possibly expand to simulation 1

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


# Running the full LCR first (no selection)
set.seed(123)
sim2_init_params <- init_params_diagnostics(p = p, G = G, M = M, K = K, n = 500)

sim2_fits_no_sel <- list()
set.seed(124)
for (t in 1:n_chains){
  sim2_fits_no_sel[[t]] <- LCR_Gibbs(X = sim2_data$X, Y = sim2_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                                     theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = FALSE, verbose = TRUE, 
                                     relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10, init_seed = TRUE, 
                                     beta_init = cbind(sim2_init_params$init_beta_list[[t]],0),
                                     z_init = sim2_init_params$init_z_list[[t]], 
                                     theta_init = sim2_init_params$init_theta_list[[t]],
                                     nu_init = sim2_init_params$init_nu_list[[t]], 
                                     gamma_init = sim2_init_params$init_gamma_list[[t]])
} 










