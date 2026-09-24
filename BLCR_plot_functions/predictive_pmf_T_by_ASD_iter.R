#Function for computing the predictive pmf on a given iteration given ASD/non-ASD, returning for each iteration:

# Iteration number
# ASD status
# T_int - possible score values
# pmf_T - corresponding probability for each score

predictive_pmf_T_by_ASD_iter <- function(t, asd, theta_samples, beta_samples, scores = c(1,2,3)) {
  G      <- dim(theta_samples)[1]
  M      <- dim(theta_samples)[2]
  Kmax   <- dim(theta_samples)[3]
  theta_iter <- vector("list", M)
  for (j in 1:M) {
    theta_iter[[j]] <- theta_samples[, j, , t]
  }
  beta_iter <- beta_samples[,,t]
  res <- predictive_pmf_T_by_ASD(asd = asd, theta = theta_iter, beta = beta_iter, scores = scores)
  data.frame(
    iter  = t,
    ASD   = if (asd == 1) "ASD" else "No ASD",
    T_int = res$totals,
    pmf_T = res$pmf_T
  )
}
