LCR_log_post_compute <- function(mu, Y_indicator, log_theta, z, theta_hyperparam, beta_prior_mean, beta_prior_cov_inv, beta, G){
  term1_mat <- log_softmax(mu)
  term2_array <- einsum("ijk,gjk -> igj", Y_indicator, log_theta)
  term2_mat <- apply(term2_array, c(1, 2), sum)
  combined_terms <- z*(term1_mat + term2_mat)
  sum_combined <- sum(combined_terms)
  theta_prior_term1 <- theta_hyperparam - 1 
  #log_theta_scaled <- log_theta * array(theta_prior_term1, dim = c(1, 1, length(theta_prior_term1)))
  log_theta_scaled <- log_theta * array(theta_prior_term1, dim = dim(log_theta))
  prior_term1 <- sum(log_theta_scaled)
  beta_prior_term1 <- (-1/2)*sum(colSums(beta * (beta_prior_cov_inv %*% beta)))
  beta_prior_term2 <- sum(crossprod(beta_prior_mean, beta_prior_cov_inv %*% beta))
  log_like <- sum_combined
  log_post <- sum_combined + prior_term1 + beta_prior_term1 + beta_prior_term2
  return(list(log_like = log_like, log_post = log_post))
}


LCR_collapsed_log_post_compute <- function(beta, beta_prior_cov_inv, beta_prior_mean, nu, mu, M, clust_var_prior, z, K, N_jk, N_gjk, N_g, n, theta_hyperparam, G){
  beta_prior_term1 <- (-1/2)*sum(colSums(beta * (beta_prior_cov_inv %*% beta)))
  beta_prior_term2 <- sum(crossprod(beta_prior_mean, beta_prior_cov_inv %*% beta))
  beta_prior_term <- beta_prior_term1 + beta_prior_term2
  current_indices <- which(nu == 1)
  excl_indices <- which(nu == 0)
  sum_current_indices <- sum(nu)
  nu_prior_term <- sum_current_indices*log(clust_var_prior) + (M-sum_current_indices)*log(1-clust_var_prior)
  logit_term <- sum(z*log_softmax(mu))
  K_current <- K[current_indices]
  K_excl <- K[excl_indices]
  N_jk_excl <- N_jk[excl_indices,]
  N_jk_current <- N_jk[current_indices,]
  N_gjk_current <- N_gjk[,current_indices,]
  term1 <- sum(lgamma(K_excl*theta_hyperparam))
  term2 <- sum(K_excl*lgamma(theta_hyperparam))
  term3 <- sum(lgamma(N_jk_excl + theta_hyperparam))
  term4 <- sum(lgamma(n + K_excl*theta_hyperparam))
  log_like_term1 <- term1 - term2 + term3 - term4
  term5 <- G*sum(lgamma(K_current*theta_hyperparam))
  term6 <- G*sum(K_current*lgamma(theta_hyperparam))
  term7 <- sum(lgamma(N_gjk_current + theta_hyperparam))
  term8 <- sum_current_indices*sum(lgamma(N_g+theta_hyperparam))
  log_like_term2 <- term5 - term6 + term7 - term8
  log_like <- log_like_term1 + log_like_term2
  log_post <- log_like + beta_prior_term + nu_prior_term + logit_term
  return(list(log_like = log_like, log_post = log_post))
}


get_log_post_function <- function(item.sel) {
  if (item.sel) {
    return(function(params) {
      with(params, LCR_collapsed_log_post_compute(beta, beta_prior_cov_inv, beta_prior_mean, nu, mu, M, clust_var_prior, z, K, N_jk, N_gjk, N_g, n, theta_hyperparam, G))
    })
  } else {
    return(function(params) {
      with(params, LCR_log_post_compute(mu, Y_indicator, log_theta, z, theta_hyperparam, beta_prior_mean, beta_prior_cov_inv, beta, G))
    })
  }
}