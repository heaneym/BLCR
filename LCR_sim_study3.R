# Creating a simulation 3 with many more items to more accurately emulate the CSHQ data

# We'll do a simulation with 50 items, 5 predictors of varying type (binary, continuous)

# We can have 20-25 of the items be informative, and 1-2 of the predictors informative

# All items having 3 response levels, and vary the proporions of each group (try to have a small enough group as well) 

# G = 4


G <- 4
n <- 150

#We consider 5 predictors, one representing age (N(0,1)), and 4 representing clinical dianoses 
# (binary, proportions 0.2, 0,15, 0.1, 0.5 - replicating ASD, other, ID diagnoses, gender respectively)
p <- 5
sim3_X <- matrix(NA, ncol = p, nrow = n)

# Continuous, representing age
sim3_X[,1] <- rnorm(n)
# Binary (0.5 proportion) representing gender
sim3_X[,2] <- sample(c(0,1), n, replace = TRUE, prob = c(0.5,0.5))
# Binary (0.2 proportion) representing ASD
sim3_X[,3] <- sample(c(0,1), n, replace = TRUE, prob = c(0.8,0.2))
# Binary (0.15 proportion) representing other diagnoses
sim3_X[,4] <- sample(c(0,1), n, replace = TRUE, prob = c(0.85,0.15))
# Binary (0.06 proportion) representing ID diagnosis
sim3_X[,5] <- sample(c(0,1), n, replace = TRUE, prob = c(0.94,0.06))

M <- 40
K <- rep(3, 40)

source('./sim_study_parameters.R')

G <- 4
n <- 500
M <- 40
K <- rep(3,M)

p <- 5

# covariate effects matrix
# sim3_beta <- matrix(c(
#   -0.8,  -0.6, -1.5, 
#    0.2,     0, 0.15, 
#      0,     0,    0, 
#      1,     0,  1.5, 
#      0,   0.1,  0.2, 
#      0,     0,  0.4 
# ), nrow = p+1, byrow = TRUE)




#We let items 1-20 be informative with distinct parameters, and 21-40 non-informative

# sim3_theta1 <- matrix(c(
#   0.4, 0.3, 0.3,
#   0.2, 0.3, 0.5,
#   0.1, 0.2, 0.7,
#   0.6, 0.2, 0.2
# ), 
# nrow = G, ncol = 3, byrow = TRUE)
# sim3_theta2 <- matrix(c(
#   0.2, 0.4, 0.4,
#   0.5, 0.3, 0.2,
#   0.1, 0.3, 0.6,
#   0.7, 0.2, 0.1
# ), 
# nrow = G, ncol = 3, byrow = TRUE)
# sim3_theta3 <- matrix(c(
#   0.5, 0.2, 0.3,
#   0.1, 0.4, 0.5,
#   0.2, 0.2, 0.6,
#   0.4, 0.3, 0.3
# ), 
# nrow = G, ncol = 3, byrow = TRUE)
# sim3_theta4 <- matrix(c(
#   0.15, 0.3, 0.55,
#   0.6, 0.25, 0.15,
#   0.1, 0.3, 0.6,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta5 <- matrix(c(
#   0.2, 0.3, 0.5,
#   0.55, 0.3, 0.15,
#   0.15, 0.25, 0.6,
#   0.6, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta6 <- matrix(c(
#   0.15, 0.35, 0.5,
#   0.65, 0.2, 0.15,
#   0.1, 0.3, 0.6,
#   0.45, 0.35, 0.2
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta7 <- matrix(c(
#   0.2, 0.25, 0.55,
#   0.6, 0.25, 0.15,
#   0.4, 0.35, 0.25,
#   0.7, 0.2, 0.1
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta8 <- matrix(c(
#   0.15, 0.30, 0.55,
#   0.55, 0.25, 0.20,
#   0.15, 0.30, 0.55,
#   0.60, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta9 <- matrix(c(
#   0.2, 0.3, 0.5,
#   0.65, 0.2, 0.15,
#   0.1, 0.25, 0.65,
#   0.55, 0.3, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta10 <- matrix(c(
#   0.15, 0.35, 0.5,
#   0.6, 0.2, 0.2,
#   0.15, 0.3, 0.55,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta11 <- matrix(c(
#   0.2, 0.3, 0.5,
#   0.55, 0.3, 0.15,
#   0.45, 0.3, 0.25,
#   0.6, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta12 <- matrix(c(
#   0.6, 0.25, 0.15,
#   0.15, 0.3, 0.55,
#   0.1, 0.3, 0.6,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta13 <- matrix(c(
#   0.55, 0.3, 0.15,
#   0.2, 0.3, 0.5,
#   0.15, 0.25, 0.6,
#   0.45, 0.35, 0.2
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta14 <- matrix(c(
#   0.65, 0.2, 0.15,
#   0.15, 0.35, 0.5,
#   0.1, 0.3, 0.6,
#   0.6, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta15 <- matrix(c(
#   0.6, 0.25, 0.15,
#   0.2, 0.3, 0.5,
#   0.4, 0.35, 0.25,
#   0.7, 0.2, 0.1
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta16 <- matrix(c(
#   0.55, 0.25, 0.2,
#   0.15, 0.3, 0.55,
#   0.15, 0.3, 0.55,
#   0.55, 0.3, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta17 <- matrix(c(
#   0.6, 0.2, 0.2,
#   0.2, 0.3, 0.5,
#   0.1, 0.25, 0.65,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta18 <- matrix(c(
#   0.6, 0.2, 0.2,
#   0.15, 0.35, 0.5,
#   0.15, 0.3, 0.55,
#   0.6, 0.25, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta19 <- matrix(c(
#   0.55, 0.3, 0.15,
#   0.2, 0.3, 0.5,
#   0.45, 0.3, 0.25,
#   0.65, 0.2, 0.15
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# sim3_theta20 <- matrix(c(
#   0.5, 0.3, 0.2,
#   0.5, 0.3, 0.2,
#   0.15, 0.30, 0.55,
#   0.45, 0.35, 0.2
# ), nrow = G, ncol = 3, byrow = TRUE)
# 
# 

set.seed(129)

#Simulatiing dataset from using the above parameters
sim3_data <- LCR_sim_data_given_X(theta = sim3_theta, beta = sim3_beta, n_samples = n, X = sim3_X)




sim3_LCR_fit_varsel <- LCR_Gibbs(X = sim3_data$X, Y = sim3_data$Y, G = 4, beta_prior_cov = diag(10^2,6), beta_prior_mean = rep(0,6), 
                          theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                          relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)
