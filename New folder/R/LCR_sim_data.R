LCR_sim_data <- function(theta,beta,n_samples){
  #theta is a list of matrices (same as in poLCA) where #(list items) = #items, #rows = #classes, #cols = #(responses to that variable)
  #beta is a matrix of dimension (p+1)x(G-1)
  #n_samples is the number of samples we want to produce
  #n_noise_var is the number of useless/noise covariates we want to produce
  G <- ncol(beta)+1
  M <- length(theta)
  K <- sapply(theta, ncol)
  p <- nrow(beta)-1 
  
  #First we produce the X matrix from N(0,1)
  X <- cbind(1,matrix(rnorm(n_samples*p),nrow = n_samples))
  Y <- matrix(0,nrow = n_samples,ncol = M)
  beta_0 <- cbind(beta, 0)
  
  #Calculate the group membership probabilities and then decide class membership
  mu <- X%*%beta_0
  exp_mu <- exp(mu)
  gamma <- exp_mu/rowSums(exp_mu)
  class_vec <- apply(gamma, 1, function(x) sample(1:G, 1, prob = x))
  #generate responses based on theta
  for(i in 1:n_samples){
    class <- class_vec[i]
    for (j in 1:M){
      theta_mat_row <- theta[[j]][class,]
      Y[i,j] <- sample(1:K[j], 1, prob = theta_mat_row)
    }
  }
  pi <- numeric(G)
  for (g in 1:G){
    pi[g] <- length(class_vec[class_vec == g])/length(class_vec)
  }
  data <- list(X = X[,-1], Y = Y, class = class_vec, G = G, M = M, K = K, p = p, class_prob = gamma, pi = pi, real_beta = beta, real_theta = theta)
  return(data)
}




