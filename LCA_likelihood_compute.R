LCA_row_likelihood_compute <- function(y_i, pi, theta){
  G <- length(pi)
  M <- length(y_i)
  K <- sapply(theta, ncol)
  log_likelihoods <- numeric(G)
  for (g in 1:G) {
    log_prob_sum <- 0
    for (j in 1:M) {
      log_prob_sum <- log_prob_sum + log(theta[[j]][g, y_i[j]])
    }
    log_likelihoods[g] <- log(pi[g]) + log_prob_sum
  }
  max_log_lik <- max(log_likelihoods)
  log_likelihood <- max_log_lik + log(sum(exp(log_likelihoods - max_log_lik)))
  return(exp(log_likelihood))
}


LCA_row_log_likelihood_compute <- function(y_i, pi, theta) {
  G <- length(pi)
  M <- length(y_i)
  K <- sapply(theta, ncol)
  
  log_likelihoods <- numeric(G)
  for (g in 1:G) {
    log_prob_sum <- sum(log(sapply(1:M, function(j) theta[[j]][g, y_i[j]])))
    log_likelihoods[g] <- log(pi[g]) + log_prob_sum
  }
  
  max_log_lik <- max(log_likelihoods)
  return(max_log_lik + log(sum(exp(log_likelihoods - max_log_lik))))
}
