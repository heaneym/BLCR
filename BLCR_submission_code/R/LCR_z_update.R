
Rcpp::sourceCpp(file.path("R", "z_update_collapsed.cpp"))
Rcpp::sourceCpp(file.path("R", "z_update_uncollapsed.cpp"))


get_z_update_function <- function(item.sel) {
  if (item.sel) {
    return(function(params) {
      with(params, z_update_collapsed(
        z, nu, mu, K, theta_hyperparam, 
        as.vector(N_gjk),          
        N_g,
        as.vector(Y_indicator),     
        Y_missing_indicator,        
        n, G, M, omega, C, p,
        beta, beta_cov_inv, beta_cov_inv_chol,
        beta_prior_mean, beta_prior_cov_inv,
        gamma, X_current))
    })
  } else {
    return(function(params) {
      with(params, z_update_uncollapsed(log_logit_probs, log_theta, Y, 
                                        Y_missing_indicator, G, n, M, K))
    })
  }
}








