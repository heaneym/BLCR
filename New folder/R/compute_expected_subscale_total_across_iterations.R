# Function for computing the expected CSHQ subscale totals for given ASD status (ASD/non-ASD)
# across iterations. The funciton returns a list containing:

# expected_subscale_total - matrix of dimension n_iter x n_subscale with the expected subscale totals on each iteration
# expected_subscale_per_group - same as above but an array stratifying this into latent groups.
# pi_across_iterations - the value of the mixture proportions given ASD status across iterations (n_iter x G matrix)



compute_expected_subscale_total_across_iterations <- function(theta_samples, asd, beta_samples, 
                                                              subscale_item_numbers = c(4,1,3,4,3,7,3,8),
                                                              levels = c(1,2,3)) {
  pi_across_iterations <- compute_pi_asd_across_iterations(beta_samples = beta_samples, asd = asd)
  
  G      <- dim(theta_samples)[1]
  M      <- dim(theta_samples)[2]
  Kmax   <- dim(theta_samples)[3]
  n_iter <- dim(theta_samples)[4]
  
  if (n_iter != nrow(pi_across_iterations) || G != ncol(pi_across_iterations)) {
    stop("Dimensions of theta_samples and pi_across_iterations do not match.")
  }
  
  n_subscale <- length(subscale_item_numbers)
  
  expected_subscale_per_group_across_iterations <- array(NA_real_, dim = c(n_iter, n_subscale, G))
  expected_subscale_total_across_iterations     <- matrix(NA_real_, nrow = n_iter, ncol = n_subscale)
  
  subscale_item_indices <- vector("list", n_subscale)
  start_item <- 1
  for (l in 1:n_subscale) {
    end_item <- start_item + subscale_item_numbers[l] - 1
    if (end_item > M) stop("subscale_item_numbers implies more items than available in theta_samples.")
    subscale_item_indices[[l]] <- start_item:end_item
    start_item <- end_item + 1
  }
  
  for (t in 1:n_iter) {
    for (g in 1:G) {
      for (l in 1:n_subscale) {
        items_l <- subscale_item_indices[[l]]
        
        expected_subscale_l_g <- 0
        
        for (m_idx in items_l) {
          for (k in levels) {
            expected_subscale_l_g <- expected_subscale_l_g + k * theta_samples[g, m_idx, k, t]
          }
        }
        
        expected_subscale_per_group_across_iterations[t, l, g] <- expected_subscale_l_g
      }
    }
    
    for (l in 1:n_subscale) {
      expected_subscale_total_across_iterations[t, l] <-
        sum(pi_across_iterations[t, ] * expected_subscale_per_group_across_iterations[t, l, ])
    }
  }
  
  list(
    expected_subscale_total     = expected_subscale_total_across_iterations,
    expected_subscale_per_group = expected_subscale_per_group_across_iterations,
    pi_across_iterations        = pi_across_iterations
  )
}