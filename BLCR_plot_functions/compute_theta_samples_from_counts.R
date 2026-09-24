#This is a function to compute the values of the theta parameters from the 'empirical' (imputed/estimated from z) quantities at each iteration 

#This is from the collapsed sampler, where our ourput is N_gjk and N_g rather than Z

#beta is the hyperparameter for the Dirichlet prior on theta

compute_theta_samples_from_counts <- function(N_gjk_samples, N_g_samples, 
                                              beta = 1) {
  
  dims <- dim(N_gjk_samples)
  G      <- dims[1]
  M      <- dims[2]
  Kmax   <- dims[3]
  n_iter <- dims[4]
  
  if (!all(dim(N_g_samples) == c(G, n_iter))) {
    stop("Dimensions of N_g_samples must be G x n_iter and match N_gjk_samples.")
  }
  
  
  if (length(beta) == 1L) {
    beta_vec <- rep(beta, Kmax)
  } else if (length(beta) == Kmax) {
    beta_vec <- beta
  } else {
    stop("beta must be scalar or length Kmax.")
  }
  
  theta_samples <- array(NA_real_, dim = c(G, M, Kmax, n_iter))
  
  for (t in 1:n_iter) {
    for (g in 1:G) {
      Ng_t <- N_g_samples[g, t]  
      
      for (j in 1:M) {
        counts_gj <- N_gjk_samples[g, j, , t]
        
        alpha_post <- counts_gj + beta_vec
        
        denom <- Ng_t + sum(beta_vec)  
        
        
        theta_samples[g, j, , t] <- alpha_post / denom
      }
    }
  }
  
  theta_samples
}
