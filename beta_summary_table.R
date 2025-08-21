multinomial_logit_coefficient_summary_table_MCMC <- function(beta_samples){
  #beta_samples is an array of dimension (p+1)xGxn_samples (note the beta samples are assumed to already be relevelled to the desired baseline)
  p <- dim(beta_samples)[1] - 1
  G <- dim(beta_samples)[2]
  n_samples <- dim(beta_samples)[3]
  #We want to create a table/dataframe/matrix with 4-5 columns (parameter name, posterior mean, posterior sd, 95% HDI lower bound, 95% HDI upper bound) 
  #the table should have (p+1)*(G-1) rows for beta_10, beta_11,...,beta_1p, beta_20,...,beta_2p,...beta_(G-1)0,...beta_(G-1)p
  beta_post_mean <- apply(beta_samples, c(1,2), mean)
  beta_post_sd <- apply(beta_samples, c(1,2), sd)
  ci_lower_mat <- matrix(0, nrow = p+1, ncol = G)
  ci_upper_mat <- matrix(0, nrow = p+1, ncol = G)
  for (l in 1:(p+1)){
    for (g in 1:G){
      hdi_temp <- hdi(beta_samples[l,g,], ci = 0.95)
      ci_upper_mat[l,g] <- hdi_temp$CI_high
      ci_lower_mat[l,g] <- hdi_temp$CI_low
    } 
  }
  beta_names <- c()
  for(g in 1:(G)) {
    for(j in 0:(p)) {
      beta_names <- c(beta_names, paste0("beta_", g, j))
    }
  }
  summary_table <- data.frame(beta_names = beta_names, beta_posterior_mean = c(beta_post_mean), beta_posterior_sd = c(beta_post_sd), beta_ci_lower = c(ci_lower_mat), beta_ci_upper = c(ci_upper_mat))
  return(summary_table)
}
