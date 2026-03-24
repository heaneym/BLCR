#Function for computing the value of the mixture components across all iterations, given ASD status (ASD/non-ASD)

compute_pi_asd_across_iterations <- function(beta_samples, asd){
  n_iter <- dim(beta_samples)[3]
  pi_mat <- matrix(0, nrow = n_iter, ncol = dim(beta_samples)[2])
  for (t in 1:n_iter){
    beta <- beta_samples[c(1,4),,t]
    pi_mat[t,] <- compute_pi_asd(asd = asd, beta = beta)
  }
  return(pi_mat)
}