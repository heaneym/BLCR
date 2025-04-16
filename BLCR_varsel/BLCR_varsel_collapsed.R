
library(BayesLogit)

n_cov <- 1
n_noise <- 0
G <- 2
n <- 500

#noise_item_var <- sample(1:2,250, replace = TRUE)


#theta1 <- matrix(c(0.8,0.2,0.3,0.7,0.2,0.8), nrow = G, byrow = TRUE)
#theta3 <- matrix(c(0.7,0.3,0.4,0.6,0.2,0.8), nrow = G, byrow = TRUE)
#theta2 <- matrix(c(0.9,0.1,0.5,0.5,0.1,0.9), nrow = G, byrow = TRUE)
#theta4 <- matrix(c(0.6,0.4,0.3,0.7,0.8,0.2), nrow = G, byrow = TRUE)


#Dean and Raftery Binary data - I'll simulate some data with with the inclusion of covariates, which will vary in significance

theta1 <- matrix(c(0.6,0.4,0.2,0.8), nrow = G, byrow = TRUE)
theta2 <- matrix(c(0.8,0.2,0.5,0.5), nrow = G, byrow = TRUE)
theta3 <- matrix(c(0.7,0.3,0.4,0.6), nrow = G, byrow = TRUE)
theta4 <- matrix(c(0.6,0.4,0.9,0.1), nrow = G, byrow = TRUE)
theta5 <- matrix(c(0.5,0.5,0.5,0.5), nrow = G, byrow = TRUE)
theta6 <- matrix(c(0.4,0.6,0.4,0.6), nrow = G, byrow = TRUE)
theta7 <- matrix(c(0.3,0.7,0.3,0.7), nrow = G, byrow = TRUE)
theta8 <- matrix(c(0.2,0.8,0.2,0.8), nrow = G, byrow = TRUE)
theta9 <- matrix(c(0.9,0.1,0.9,0.1), nrow = G, byrow = TRUE)
theta10 <- matrix(c(0.6,0.4,0.6,0.4), nrow = G, byrow = TRUE)
theta11 <- matrix(c(0.7,0.3,0.7,0.3), nrow = G, byrow = TRUE)
theta12 <- matrix(c(0.8,0.2,0.8,0.2), nrow = G, byrow = TRUE)
theta13 <- matrix(c(0.1,0.9,0.1,0.9), nrow = G, byrow = TRUE)

realtheta <- list(theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10, theta11, theta12, theta13)
M <- length(realtheta)
K <- sapply(realtheta, ncol)

#Strongly informative predictor
X1 <- rnorm(n)


#X2 correlated to X1 
X2 <- 0.7*X1 + rnorm(n,0,0.2) 


#Pure noise variables
X3 <- rnorm(n)
X4 <- rnorm(n)

#weakly informative
X5 <- rnorm(n)

#for the coefficient of X3 here, i seems like if beta<0.7, then it excludes all the time, and when its bigger than 0.7, it can sometimes exclude it all the time and sometimes always include it.
realbeta <- c(0.4, 1, 0.5, 0, 0, 0.5)
p <- length(realbeta)-1

X <- cbind(1, X1, X2, X3, X4, X5)

mu_temp <- X%*%realbeta

logit_prob_temp <- exp(mu_temp)/(1+exp(mu_temp))

logit_prob_temp_mat <- cbind(logit_prob_temp, 1-logit_prob_temp)

class_vec <- apply(logit_prob_temp_mat, 1, function(x) sample(1:G, 1, prob = x))

Y <- matrix(0,nrow = n,ncol = M)

for(i in 1:n){
  class <- class_vec[i]
  for (j in 1:M){
    theta_mat_row <- realtheta[[j]][class,]
    Y[i,j] <- sample(1:K[j], 1, prob = theta_mat_row)
  }
}


X <- X[,2:6]









#realtheta <- list(theta1,theta2,theta3,theta4)


#realbeta <- matrix(c(-0.5,1,0.5,-1), nrow = n_cov+1, ncol = G-1, byrow = FALSE)
#data <- LCR_sim_data(theta = realtheta, beta = realbeta, n_samples = n, n_noise_var = n_noise)
#list2env(data, envir = .GlobalEnv)


#Y <- cbind(Y,noise_item_var)

init <- initialise_variables_BLCR_collapsed(G = G, X = X, Y = Y, beta_prior_mean = c(0,0), beta_prior_var_factor = 10)

a = 1
b = 1

list2env(init, envir = .GlobalEnv)


n_iter <- 50000
beta_samples <- array(0, dim = c(dim(beta), n_iter))
nu_samples <- array(0, dim = c(length(nu), n_iter))
z_samples <- array(0, dim = c(dim(z), n_iter))
w_samples <- array(0, dim = c(dim(w), n_iter))
gamma_samples <- array(0, dim = c(length(gamma),n_iter))
N_g_samples <- array(0, dim = c(length(N_g), n_iter))
N_gjk_samples <- array(0, dim = c(dim(N_gjk), n_iter))


for (count in 1:n_iter){
  
  gamma_up <- gamma_update(gamma = gamma, p = p, beta_prior_cov_inv = beta_prior_cov_inv, beta_prior_mean = beta_prior_mean, kappa = kappa , G = G, X_current = X_current, omega = omega, tau =  tau, X = X)
  
  list2env(gamma_up, envir = .GlobalEnv)
  
  #print(gamma)
  
  omega <- omega_update(eta = eta, mu = mu, G = G, n = n, omega = omega)
  
  exp_mu <- exp(mu)
  
  logit_probs <- logit_probs_update(eta = eta)
  
  log_logit_probs <- log(logit_probs)
  
  C <- C_update(mu = mu, exp_mu = exp_mu, G = G, C = C)
  
  kappa <- kappa_update(z)
  
  eta <- eta_update(mu = mu, C = C)
  
  A <- A_update(kappa = kappa, omega = omega, C = C)
  
  beta_up <- beta_update(gamma = gamma, X_current = X_current, G = G, p = p, beta_prior_cov_inv = beta_prior_cov_inv, A = A, beta_prior_mean = beta_prior_mean, omega = omega)
  
  list2env(beta_up, envir = .GlobalEnv)
  
  beta_samples[,,count] <- beta
  
  mu <- mu_update(X = X, beta = beta)
  
  tau <- tau_update(gamma = gamma, a = a, b = b)
  
  N_up <- N_updates(z = z, Y_indicator = Y_indicator, N_jk, nu)
  
  N_g <- N_up$N_g
  
  N_gjk <- N_up$N_gjk
  
  nu <- nu_update(nu = nu, M = M, K = K, G = G, alpha = alpha, N_g, N_gjk = N_gjk, inclusion_sum = inclusion_sum, exclusion_sum = exclusion_sum, count = count)
  
  #z_up <- z_update(z = z, nu = nu, mu = mu, K = K, alpha = alpha, N_gjk = N_gjk, N_g = N_g, n = n, G = G, Y_indicator = Y_indicator, log_logit_probs = log_logit_probs, gamma = gamma, X_current = X_current, A = A, beta = beta, kappa, beta_mean = beta_mean, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, omega = omega, C = C)
  z_up <- z_update(z = z, nu = nu, mu = mu, K = K, alpha = alpha, N_gjk = N_gjk, N_g = N_g, n = n, G = G, Y_indicator = Y_indicator, gamma = gamma, X_current = X_current,  beta = beta, beta_cov_inv = beta_cov_inv, beta_cov_inv_chol = beta_cov_inv_chol, omega = omega, C = C)
  
  
  list2env(z_up, envir = .GlobalEnv)
  
  nu_samples[,count] <- nu
  
  gamma_samples[,count] <- gamma
  
  z_samples[,,count] <- z
  
  w_samples[,,count] <- w
  
  N_g_samples[,count] <- N_g
  
  N_gjk_samples[,,,count] <- N_gjk
  if (count %% 500 == 0){
    print(count) 
  }
}



z_mat <- matrix(0, nrow = dim(z_samples)[3], ncol = n)
for (iter in 1:dim(z_samples)[3]){
  temp <- z_samples[,,iter]
  temp_row <- apply(temp,1,which.max)
  z_mat[iter,] <- temp_row 
}
w_samples_reshape <- aperm(w_samples,c(3,1,2))
ls <- label.switching(method = "STEPHENS", z = z_mat, K = G, p = w_samples_reshape)

ls_perm <- ls$permutations$STEPHENS








