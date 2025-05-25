relabel_outputs <- function(beta_samples, z_samples, w_samples, theta_samples,
                            N_g_samples, N_gjk_samples, n_samples, n, G) {
  z_mat <- matrix(0, nrow = n_samples, ncol = n)
  for (iter in 1:n_samples) {
    temp     <- z_samples[, , iter]
    temp_row <- apply(temp, 1, which.max)
    z_mat[iter, ] <- temp_row
  }
  w_samples_reshape <- aperm(w_samples, c(3,1,2))
  ls <- suppressMessages(
    suppressWarnings({
      .ls_tmp <- NULL
      invisible(
        capture.output(
          .ls_tmp <- label.switching(
            method = "STEPHENS",
            z      = z_mat,
            K      = G,
            p      = w_samples_reshape
          ),
          file = NULL
        )
      )
      .ls_tmp 
    })
  )
  ls_perm <- ls$permutations$STEPHENS
  reordered_beta_samples   <- array(0, dim = dim(beta_samples))
  reordered_z_samples      <- array(0, dim = dim(z_samples))
  reordered_w_samples      <- array(0, dim = dim(w_samples))
  reordered_theta_samples  <- array(0, dim = dim(theta_samples))
  reordered_N_gjk_samples  <- array(0, dim = dim(N_gjk_samples))
  reordered_N_g_samples    <- array(0, dim = dim(N_g_samples))
  for (i in 1:n_samples) {
    reordered_beta_samples[ , , i] <- beta_samples[ , ls_perm[i, ], i]
    baseline <- reordered_beta_samples[ , G, i]
    reordered_beta_samples[ , , i] <- reordered_beta_samples[ , , i] - baseline
    reordered_z_samples[ , , i]     <- z_samples[ , ls_perm[i, ], i]
    reordered_w_samples[ , ls_perm[i, ], i] <- w_samples[ , ls_perm[i, ], i]
    reordered_theta_samples[ls_perm[i, ], , , i]  <- theta_samples[ls_perm[i, ], , , i]
    reordered_N_gjk_samples[ls_perm[i, ], , , i]  <- N_gjk_samples[ls_perm[i, ], , , i]
    reordered_N_g_samples[ls_perm[i, ], i]        <- N_g_samples[ls_perm[i, ], i]
  }
  reordered_samples <- list(
    beta_samples   = reordered_beta_samples,
    z_samples      = reordered_z_samples,
    w_samples      = reordered_w_samples,
    theta_samples  = reordered_theta_samples,
    N_gjk_samples  = reordered_N_gjk_samples,
    N_g_samples    = reordered_N_g_samples
  )
  return(reordered_samples)
}