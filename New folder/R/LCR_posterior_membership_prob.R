
LCR_posterior_membership_prob <- function(theta, beta, X, Y){
  n <- nrow(X)
  M <- ncol(Y)
  G <- ncol(beta)
  
  S <- matrix(0, nrow = n, ncol = G)
  
  for (j in 1:M){
    log_theta_j <- log(theta[[j]])
    log_theta_j[is.infinite(log_theta_j)] <- 0
    
    y_j <- Y[,j]
    obs_idx <- which(!is.na(y_j))
    
    for (i in obs_idx){
      S[i,] <- S[i,] + log_theta_j[, y_j[i]]
    }
  }
  

  X_1 <- cbind(1, X)
  mu <- X_1 %*% beta
  

  log_post_prop <- mu + S
  

  row_maxes <- apply(log_post_prop, 1, max)
  post_prop <- exp(log_post_prop - row_maxes)
  post <- post_prop / rowSums(post_prop)
  
  return(post)
}








