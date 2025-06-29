# z_update_uncollapsed <- function(log_logit_probs,log_theta,Y,G,n,M,K){
#   log_theta_list <- lapply(1:M, function(j) {
#     log_theta[, j, Y[, j]]
#   })
#   sum_log_theta <- Reduce(`+`, log_theta_list)
#   log_w <- log_logit_probs + t(sum_log_theta)
#   w <- exp(log_w)
#   z <- t(apply(w, 1, function(row) rmultinom(1, 1, row)))
#   return(list(w = w, z = z))
# }

setwd('C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel')

sourceCpp("z_update_collapsed.cpp")
sourceCpp("z_update_uncollapsed.cpp")




# z_update_collapsed <- function(z, nu, mu, K, theta_hyperparam, N_gjk, N_g, Y_indicator, n, G,
#                      omega, C, p, beta, beta_cov_inv, beta_cov_inv_chol,
#                      beta_prior_mean, beta_prior_cov_inv, gamma, X_current) {
#   w <- matrix(0, nrow = n, ncol = G)
#   which_item_var <- which(nu == 1)
#   #K_current <- K[nu == 1]
#   K_current <- K[which_item_var]
#   for (i in 1:n) {
#     curr_z_i <- z[i,]  # Store current assignment
# 
#     # Remove contribution of observation i from counts
#     N_g_minus_i <- N_g - curr_z_i
#     
#     mask_vec <- as.logical(Y_indicator[i, , ] == 1)
# 
#     for (g in 1:G) {
#       # Calculate new counts if observation i is assigned to group g
#       new_z_i <- rep(0, G)
#       new_z_i[g] <- 1
# 
#       # Update counts
#       N_g_temp <- N_g_minus_i + new_z_i
# 
#       N_gjk_temp <- N_gjk
#       
#       
#       N_gjk_mat      <- matrix(N_gjk, nrow = G)
#       N_gjk_mat_temp <- matrix(N_gjk_temp, nrow = G)
#       delta <- -curr_z_i + as.numeric(seq_len(G) == g)
#       N_gjk_mat_temp[ , mask_vec] <- N_gjk_mat[ , mask_vec] + delta
#       dim(N_gjk_mat_temp) <- c(G,M,max(K))
#       N_gjk_temp <- N_gjk_mat_temp
# 
# 
#       log_gamma_N_gjk_temp <- lgamma(N_gjk_temp[,which_item_var,, drop=FALSE] + theta_hyperparam)
#       term2 <- sum(log_gamma_N_gjk_temp)
# 
#       temp_mat <- outer(N_g_temp, K_current * theta_hyperparam, "+")
#       term3 <- sum(lgamma(temp_mat))
#       
#       term0 <- mu[i,g]
# 
# 
#       w[i,g] <- term0 + term2 - term3
#     }
# 
# 
#     max_w_i <- max(w[i,])
#     exp_w_i <- exp(w[i,] - max_w_i)
#     w[i,] <- exp_w_i/sum(exp_w_i)
# 
# 
#     z[i,] <- rmultinom(1, 1, w[i,])
# 
#     # Update overall counts after this reassignment
#     N_g <- N_g_minus_i + z[i,]
# 
#     # Update N_gjk
#     delta <- z[i, ] - curr_z_i
#     N_gjk_mat[ , mask_vec] <- N_gjk_mat[ , mask_vec] + delta
#     dim(N_gjk_mat) <- c(G,M,max(K))
#     N_gjk <- N_gjk_mat
#   }
# 
#   # Calculate final kappa
#   kappa <- z - 0.5
# 
#   return(list(z = z, w = w, kappa = kappa, N_g = N_g, N_gjk = N_gjk))
# }



get_z_update_function <- function(item.sel) {
  if (item.sel) {
    return(function(params) {
      with(params, z_update_collapsed(z, nu, mu, K, theta_hyperparam, N_gjk, N_g,
                                      Y_indicator, n, G, M,  omega, C, p,
                                      beta, beta_cov_inv, beta_cov_inv_chol,
                                      beta_prior_mean, beta_prior_cov_inv,
                                      gamma, X_current))
    })
  } else {
    return(function(params) {
      with(params, z_update_uncollapsed(log_logit_probs, log_theta, Y, G, n, M, K))
    })
  }
}

