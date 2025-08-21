LCR_posterior_membership_prob <- function(theta, beta, X, Y){
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X)
  G <- ncol(beta)
  K <- vapply(1:M, function(j) length(unique(na.omit(Y[,j]))), numeric(1))
  max_K <- max(K)
  Y_indicator <- array(0, dim = c(n, M, max_K))
  for (i in 1:n){
    for (j in 1:M){
      for (k in 1:K[j]){
        Y_indicator[i,j,k] <- 1*(Y[i,j] == k)
      }
    }
  }
  padded_theta <- lapply(theta, function(mat) {
    if (ncol(mat) < max_K) {
      zeros <- matrix(0, nrow = G, ncol = max_K - ncol(mat))
      mat <- cbind(mat, zeros)
    }
    return(mat)
  })
  theta_array <- abind::abind(padded_theta, along = 3)
  theta_array <- aperm(theta_array, c(1, 3, 2))
  log_theta <- log(theta_array)
  dim(Y_indicator) <- c(n, M * max_K)
  dim(log_theta) <- c(G, M * max_K)
  log_theta[log_theta == -Inf] <- 0
  S <- Y_indicator %*% t(log_theta)
  X_1 <- cbind(1, X)
  mu <- X_1 %*% beta
  log_post_prop <- mu + S
  row_maxes <- apply(log_post_prop, 1, max)
  post_prop <- exp(log_post_prop - row_maxes)
  post <- post_prop / rowSums(post_prop)
  return(post)
}











