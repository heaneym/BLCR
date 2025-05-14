LCR_sim_data_mixed_predictor <- function(theta, beta, n_samples, n_disc_unif_var = 0, disc_unif_range = 0, n_bernoulli_var = 0, bernoulli_params = 0, n_cat_var = 0,cat_outcomes = 0, cat_params = 0, n_ind_norm_var = 0){
  
  G <- ncol(beta) + 1
  M <- length(theta)
  K <- sapply(theta, ncol)
  p <- nrow(beta)-1 
  if(n_disc_unif_var + n_bernoulli_var + n_cat_var + n_ind_norm_var != p){
    stop("length of coefficient vector doesnt match the number of chosen variables")
  }
  
  if (n_disc_unif_var > 0){
    X_disc_unif <- matrix(, nrow = n_samples, ncol = n_disc_unif_var)
    for (i in 1:n_disc_unif_var){
      X_disc_unif[,i] <- sample(disc_unif_range[i,1]:disc_unif_range[i,2], size = n_samples, replace = TRUE)
    }
  } else {
    X_disc_unif <- NULL
  }
  
  #I think this Bernoulli is redundant since the categorical one covers it as a special case - I think this is the case for the discrete uniform as well
  if (n_bernoulli_var > 0){
    X_bernoulli <- matrix(, nrow = n_samples, ncol = n_bernoulli_var)
    for (i in 1:n_bernoulli_var){
      X_bernoulli[,i] <- rbinom(n_samples, 1, bernoulli_params[i])
    }
  } else {
    X_bernoulli <- NULL
  }
  
    #FINISH THIS
    
    if (n_cat_var > 0) {
      X_cat_vec <- matrix(NA, nrow = n_samples, ncol = n_cat_var)
      for (i in 1:n_cat_var) {
        n_levels <- sum(cat_params[i, ] != 0)
        X_cat_vec[, i] <- sample(1:n_levels, size = n_samples, replace = TRUE, prob = cat_params[i, 1:n_levels])
      }
      X_cat_df <- as.data.frame(X_cat_vec)
      names(X_cat_df) <- paste0("Var", 1:n_cat_var)
      X_cat_df[] <- lapply(X_cat_df, as.factor)
      X_cat <- model.matrix(~ . - 1, data = X_cat_df)
    } else {
      X_cat <- NULL
    }
    
    if(n_ind_norm_var > 0){
      X_norm <- matrix(rnorm(n_samples*n_ind_norm_var), nrow = n_samples, ncol = n_ind_norm_var)
    } else {
      X_norm <- NULL
    }
    
    X <- cbind(1,X_disc_unif, X_bernoulli, X_cat, X_norm)
    
    Y <- matrix(0,nrow = n_samples,ncol = M)
    
    print(dim(X))
    print(dim(beta))
    
    mu <- cbind(X%*%beta, 0)
    
    exp_mu <- exp(mu)
    
    logit_prob <- exp_mu/rowSums(exp_mu)
    
    class_vec <- apply(logit_prob, 1, function(x) sample(1:G, 1, prob = x))
    
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
    data <- list(X = X[,-1], Y = Y, class = class_vec, G = G, M = M, K = K, p = p, class_prob = logit_prob, pi = pi, real_beta = beta, real_theta = theta)
    return(data)
    
}
