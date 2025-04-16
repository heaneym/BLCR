#This function just simulates some data using the poLCA function and then prepares it for the BLCR sampler

prepare_data_2group_poLCA_covariates <- function(seed, n, M, p, beta,theta = NULL, G = 2){
  #seed is the seed we want to simulate the data with, n is the number of observations
  #M is the number of clustering variables, p is the number of covariates, and beta is the vector of these covariates
  set.seed(seed)
  if(is.null(beta)){
    b <- NULL
  } else {
    b <- matrix(beta,nrow = p+1, ncol = G-1)
  }
  dat <- poLCA.simdata(N = n, probs = theta, nclass = G, ndv = M, nresp = NULL, x = NULL, niv = p, b = b, P = NULL, missval = FALSE, pctmiss = NULL)
  Y_columns <- lapply(1:M, function(i) dat[["dat"]][[paste0("Y", i)]])
  Y <- do.call(cbind, Y_columns)
  X_columns <- lapply(1:p, function(i) dat[["dat"]][[paste0("X", i)]])
  X <- do.call(cbind, X_columns)
  tc <- dat[["trueclass"]]
  tc[tc == 1] <- 0
  tc[tc == 2] <- 1
  return(list(tc = tc, Y = Y, X = X, data = dat))
}


#This function simulates data without poLCA and then prepares it for the BLCR sampler
#It is not fully finished!!!


prepare_data_2group_covariates_not_poLCA <- function(seed, n, p, beta = NULL, theta=NULL, G = 2){
  set.seed(seed)
  if(is.null(beta)){
    b <- matrix(rnorm(p+1,0,5),nrow = p+1, ncol = G-1)
  } else {
    b <- matrix(beta,nrow = p+1, ncol = G-1)
  }
  dat <- generate_BLCR_probit_binary(beta = b, theta = theta, n = n, G)
  M <- ncol(dat$Y)
  Y_columns <- lapply(1:M, function(i) dat[["dat"]][[paste0("Y", i)]])
  Y <- do.call(cbind, Y_columns)
  X_columns <- lapply(1:p, function(i) dat[["dat"]][[paste0("X", i)]])
  X <- do.call(cbind, X_columns)
  
  tc <- dat[["membership"]]
  return(list(tc = tc, Y = Y, X = X, data = dat, theta = dat$itemprobs, beta = dat$beta ))
}










#This function takes response data Y and covariates X and initialises variables for the Gibbs Sampler

initialize_variables_2group_BLCR <- function(X,Y,G = 2, a = NULL){
  # This function initialises the parameters - I wanna see about making a function for each of these maybe
  p <- ncol(X)
  if (is.vector(X)){
    ones <- rep(1,length(X))
    X <- cbind(ones,X)
  } else{
    ones <- rep(1,nrow(X))
    X <- cbind(ones,X)
  }
  M <- ncol(Y)
  n <- nrow(Y)
  K <- apply(Y, 2, function(x) length(unique(x)))
  
  #-----------------------------------------------------------------------------
  #-------------------------theta initialisation--------------------------------
  #-----------------------------------------------------------------------------
  
  theta <- array(0, dim = c(G, M, max(K)))  #Initialising theta
  if (is.null(a)){
    alpha <- rep(1, max(K))
  } else {
    if(length(a) != max(K)){
      print("The parameter a has the wrong length")
    } else{
      alpha <- a
    }
  }
  for (g in 1:G) { #generating a sample from theta_{gj vectors for each g and j}
    for (j in 1:M) {
      theta[g, j, 1:K[j]] <- rdirichlet(1, alpha[1:K[j]])  
    }
  }
  
  #-----------------------------------------------------------------------------
  #-------------------------beta initialisation---------------------------------
  #-----------------------------------------------------------------------------
  
  beta <- array(rnorm((G-1)*(p+1),2,0.5),dim = c(G-1,p+1))
  
  #-----------------------------------------------------------------------------
  #-------------------------z initialisation------------------------------------
  #-----------------------------------------------------------------------------
  
  w <- array(0, dim = c(n,G)) 
  z <- array(0, dim = c(n,G)) 
  for(i in 1:n){
    w[i,] <- rep(1/G,G)
    z[i,] <- rmultinom(1,1,w[i,])
  }
  
  #-----------------------------------------------------------------------------
  #-------------------------C initialisation------------------------------------
  #-----------------------------------------------------------------------------
  
  C <- array(0, dim = c(n,1))
  for (i in 1:n){
    mu <-  beta %*% X[i,]
    C[i] <- (all(z[i,] == c(1,0)))*rtruncnorm(1,a=0,b=Inf,mean = mu,sd=1) + (all(z[i,] == c(0,1)))*rtruncnorm(1,a=-Inf,b=0,mean = mu,sd=1)
  }
  
  
  #Precalculating various quantities needed in the sampler
  
  Y_flat <- as.vector(t(Y))
  I <- rep(1:n, each = M * G)
  J <- rep(rep(1:M, each = G), n)
  G_0 <- rep(1:G, times = n * M)
  K_0 <- rep(Y_flat, each = G)
  
  #Precalculating the covariance for beta update
  beta_cov <- chol2inv(chol(crossprod(X)))
  
  return(list(X = X, Y = Y, theta = theta, alpha = alpha, beta = beta, w = w, z = z, C = C, n = n, M = M, K = K, p = p, Y_flat = Y_flat, I = I, J = J, G_0 = G_0, K_0 = K_0, beta_cov = beta_cov))
}




#This function runs the Gibbs sampler

run_gibbs_BLCR_2class <- function(Y,X,G = 2, n_iter, burnin, thin = 1, verbose = FALSE){
  #-----------------------------------------------------------------------------
  #-------------Initialising variables first and storing them-------------------
  #-----------------------------------------------------------------------------
  init <- initialize_variables_2group_BLCR(X,Y,G)
  X <- init$X
  Y <- init$Y
  theta <- init$theta
  alpha <- init$alpha
  beta <- init$beta
  w <- init$beta
  z <- init$z
  C <- init$C
  n <- init$n
  M <- init$M
  K <- init$K
  p <- init$p
  Y_flat <- init$Y_flat
  I <- init$I
  J <- init$J
  G_0 <- init$G_0
  K_0 <- init$K_0
  beta_cov <- init$beta_cov
  
  
  beta_samples <- array(0,dim = c(G-1,p+1,n_iter))
  z_samples <- array(0, dim = c(n,G,n_iter))
  C_samples <- array(0,dim = c(n,1,n_iter))
  theta_samples <- array(0, dim = c(G,M,max(K),n_iter))
  
  #-----------------------------------------------------------------------------
  #-------------------------------Update Loop-----------------------------------
  #-----------------------------------------------------------------------------
  
  
  for(count in 1:n_iter){
    
    #Update S and then theta
    S <- S_update(Y_flat,I,J,G_0,K_0,M,z)
    theta <- theta_update(alpha,G,M,K,S,theta)
    theta_samples[,,,count] <- theta
    
    #Update mu, Phi, then w, then z
    mu_vec <- X%*%t(matrix(beta, nrow=1))
    Phi_vec <- pnorm(mu_vec,0,1)
    not_Phi_vec <- 1 - Phi_vec
    log_Phi_vec <- log(Phi_vec)
    log_not_Phi_vec <- log(not_Phi_vec)
    log_theta <- log(theta)
    
    w <- w_update(Y,n,M,K,log_theta = log_theta, log_Phi_vec = log_Phi_vec, log_not_Phi_vec = log_not_Phi_vec)
    z <- z_update(w)
    z_samples[,,count] <- z
    
    #Update C
    C <- C_update(z,mu_vec)
    C_samples[,,count] <- C
    
    #Update beta
    beta <- beta_update(beta_cov,X,C)
    beta_samples[,,count] <- beta
    
    #print iteration number after every 500 iterations
    if(verbose == TRUE && count %% 500 == 0){
      cat("\rIteration:", count)
      flush.console()
    }
  }
  
  #-----------------------------------------------------------------------------
  #--------------------------Post Processing------------------------------------
  #-----------------------------------------------------------------------------
  
  
  beta_samples_burned <- beta_samples[,,burnin:(n_iter-1), drop = FALSE]
  theta_samples_burned <- theta_samples[,,,burnin:(n_iter-1), drop = FALSE]
  C_samples_burned <- C_samples[,,burnin:(n_iter-1), drop = FALSE]
  z_samples_burned <- z_samples[,,burnin:(n_iter-1), drop = FALSE]
  
  beta_samples_thin <- beta_samples_burned[,,seq(1, dim(beta_samples_burned)[3], by = thin), drop = FALSE]
  theta_samples_thin <- theta_samples_burned[,,,seq(1, dim(theta_samples_burned)[4], by = thin), drop = FALSE]
  C_samples_thin <- C_samples_burned[,,seq(1, dim(C_samples_burned)[2], by = thin), drop = FALSE]
  z_samples_thin <- z_samples[,,seq(1, dim(z_samples_burned)[3], by = thin), drop = FALSE]
  
  beta_estimate <- apply(beta_samples_thin, c(1,2), mean)
  beta_sd <- apply(beta_samples_thin, c(1,2), sd)
  theta_estimate <- apply(theta_samples_thin, c(1,2,3), mean)
  theta_sd <- apply(theta_samples_thin, c(1,2), sd)
  
  return(list(beta_estimate = beta_estimate, theta_estimate = theta_estimate, beta_sd = beta_sd, theta_sd = theta_sd, beta_samples = beta_samples_thin, theta_samples = theta_samples_thin, C_samples = C_samples_thin, z_samples = z_samples_thin, init = init))
}





















#Update function for S

S_update <- function(Y_flat,I,J,G_0,K_0,M,z){
  z_binded_M <- do.call(cbind, replicate(M, z, simplify = FALSE))
  z_flat <- as.vector(t(z_binded_M))
  S_df <- data.frame(
    i = I,         
    j = J,      
    g = G_0,         
    k = K_0,       
    z_value = z_flat 
  )
  S_df <- S_df[S_df$z_value == 1, ]
  S <- with(S_df, table(g, j, k))
  return(S)
}

#Update function for theta

theta_update <- function(alpha,G,M,K,S,theta){
  sapply(1:G, function(g) {
    sapply(1:M, function(j) {
      # Generate theta values for each g and j
      theta[g, j, 1:K[j]] <<- rdirichlet(1, alpha[1:K[j]] + S[g, j, 1:K[j]])
    })
  })
  return(theta)
}

#Update function for w

w_update <- function(Y,n,M,K,log_theta,log_Phi_vec,log_not_Phi_vec){
  log_sum1 <- numeric(n)
  log_sum2 <- numeric(n)
  for (i in 1:n) {
    for (j in 1:M) {
      y_ij <- Y[i,j]
      log_sum1[i] <- log_sum1[i] + (y_ij<=K[j])*log_theta[1, j, y_ij]
      log_sum2[i] <- log_sum2[i] + (y_ij<=K[j])*log_theta[2, j, y_ij]
    }
  }
  log_w1 <- log_Phi_vec + log_sum1
  log_w0 <- log_not_Phi_vec + log_sum2
  w <- exp(cbind(log_w1,log_w0))
  return(w)
}

#Update function for C

C_update <- function(z,mu_vec){
  a <- ifelse(z[,1] == 1, 0, -Inf)
  b <- ifelse(z[,1] == 1, Inf, 0)
  C <- rtruncnorm(1,a = a,b = b, mean = mu_vec)
  return(C)
}

beta_update <- function(beta_cov,X,C){
  beta <- mvrnorm(1,mu = beta_cov%*%(t(X)%*%C), Sigma = beta_cov)
  return(beta)
}

z_update <- function(w){
  z <- t(apply(w, 1, function(row) rmultinom(1, 1, row)))
}


#I want to try make a function here that plots the values of theta, and one for plotting the traces and densities of the beta parameters
#plot_theta <- function(){
  
#}


#compute_w <- function(beta,X,Y,n,M,K){
  #mu_vec <- X%*%beta
  #Phi_vec <- pnorm(mu_vec,0,1)
  #not_Phi_vec <- rep(1,n) - Phi_vec
  #log_Phi_vec <- log(Phi_vec)
  #log_not_Phi_vec <- log(not_Phi_vec)
  #log_theta <- log(theta)
  
  #This code works (?) but i could get it a bit more efficient without the loops maybe
  #log_sum1 <- numeric(n)
  #log_sum2 <- numeric(n)
  #for (i in 1:n) {
    #for (j in 1:M) {
      #y_ij <- Y[i,j]
      #log_sum1[i] <- log_sum1[i] + (y_ij<=K[j])*log_theta[1, j, y_ij]
      #log_sum2[i] <- log_sum2[i] + (y_ij<=K[j])*log_theta[2, j, y_ij]
    #}
  #}
  #log_w1 <- log_Phi_vec + log_sum1
  #log_w0 <- log_not_Phi_vec + log_sum2
  #w <- exp(cbind(log_w1,log_w0))
  #return(w)
#}















#initialize_variables <- function(Y, X, G = 2) {
  #n <- nrow(Y)
  #M <- ncol(Y)
  #K <- apply(Y, 2, function(x) length(unique(x)))
  #p <- ncol(as.matrix(X)) 
  
  #if (is.vector(X)) {
    #X <- cbind(1, X)
  #} else {
    #X <- cbind(1, X)
  #}
  
  #theta <- array(0, dim = c(G, M, max(K)))
  #alpha <- rep(1, max(K))
  
  #for (g in 1:G) {
    #for (j in 1:M) {
      #theta[g, j, 1:K[j]] <- rdirichlet(1, alpha[1:K[j]])
    #}
  #}
  
  #beta <- array(rnorm((G - 1) * (p + 1), 2, 0.5), dim = c(G - 1, p + 1))
  #w <- array(0, dim = c(n, G))
  #z <- array(0, dim = c(n, G))
  
  #for (i in 1:n) {
    #z[i, ] <- rmultinom(1, 1, c(0.5, 0.5))
  #}
  
  #C <- array(0, dim = c(n, 1))
  #for (i in 1:n) {
    #mu <- beta %*% X[i, ]
    #C[i] <- (all(z[i, ] == c(1, 0)))*rtruncnorm(1, a = 0, b = Inf, mean = mu, sd = 1) + 
      #(all(z[i, ] == c(0, 1)))*rtruncnorm(1, a = -Inf, b = 0, mean = mu, sd = 1)
  #}
  
  #beta <- mvrnorm(1, mu = solve(t(X) %*% X) %*% (t(X) %*% C), Sigma = solve(t(X) %*% X))
  
  #return(list(theta = theta, alpha = alpha, beta = beta, w = w, z = z, C = C, n = n, M = M, K = K, p = p))
#}




#run_gibbs_sampling <- function(tc, Y, X, n_iter = 100000, burnin = 10000, thin = 10) {
  #init <- initialize_variables(Y, X)
  #theta <- init$theta
  #alpha <- init$alpha
  #beta <- init$beta
  #w <- init$w
  #z <- init$z
  #C <- init$C
  #n <- init$n
  #M <- init$M
  #K <- init$K
  #p <- init$p
  
  #beta_samples <- array(0, dim = c(G - 1, p + 1, n_iter))
  #z_samples <- array(0, dim = c(n, G, n_iter))
  #C_samples <- array(0, dim = c(n, 1, n_iter))
  #theta_samples <- array(0, dim = c(G, M, max(K), n_iter))
  
  #for (count in 1:n_iter) {
    # Theta update
    #S <- S_create(Y_flat, I, J, G_0, K_0, M, z)
    #for (g in 1:G) {
      #for (j in 1:M) {
        #theta[g, j, 1:K[j]] <- rdirichlet(1, alpha[1:K[j]] + S[g, j, 1:K[j]])
      #}
    #}
    
    # Store samples
    #theta_samples[,,,count] <- theta
    
    # Compute posterior for z and w
    #w <- compute_w(beta,X,Y,n,M,K)
    #for (i in 1:n) {
      #z[i, ] <- rmultinom(1, 1, w[i, ])
    #}
    #z_samples[,,count] <- z
    
    # Sample C and beta
    #C <- sample_truncated_normal(z, beta, X, n, mu_vec)
    #C_samples[,,count] <- C
    #beta <- update_beta(X, C, p)
    #beta_samples[,,count] <- beta
    #if(count %% 500 == 0){
      #print(paste("Iteration ", count, " completed"))
    #}
  #}
  
  #return(list(beta_samples = beta_samples, theta_samples = theta_samples, C_samples = C_samples, z_samples = z_samples))
#}





