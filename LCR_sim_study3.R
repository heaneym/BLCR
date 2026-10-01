# Creating a simulation 3 with many more items to more accurately emulate the CSHQ data

# We'll do a simulation with 50 items, 5 predictors of varying type (binary, continuous)

# We can have 20-25 of the items be informative, and 1-2 of the predictors informative

# All items having 3 response levels, and vary the proporions of each group (try to have a small enough group as well) 

# G = 4


G <- 4
n <- 500

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

M <- 50
K <- rep(3, 50)

# covariate effects matrix
sim3_beta <- matrix(c(
  -0.8,  -0.6, -1.5, 
   0.2,     0, 0.15, 
     0,     0,    0, 
     1,     0,  1.5, 
     0,   0.1,  0.2, 
     0,     0,  0.4, 
), nrow = p+1, byrow = TRUE)

#We let items 1-25 be informative with distinct parameters, and 25-50 non-informative

sim3_theta1 <- matrix(c(
  0.4, 0.3, 0.3,
  0.2, 0.3, 0.5,
  0.1, 0.2, 0.7,
  0.6, 0.2, 0.2
), 
nrow = G, ncol = 3, byrow = TRUE)
sim3_theta2 <- matrix(c(
  0.2, 0.4, 0.4,
  0.5, 0.3, 0.2,
  0.1, 0.3, 0.6,
  0.7, 0.2, 0.1
), 
nrow = G, ncol = 3, byrow = TRUE)
sim3_theta3 <- matrix(c(
  0.5, 0.2, 0.3,
  0.1, 0.4, 0.5,
  0.2, 0.2, 0.6,
  0.4, 0.3, 0.3
), 
nrow = G, ncol = 3, byrow = TRUE)
sim3_theta4 <- matrix(c(
  0.15, 0.3, 0.55,
  0.6, 0.25, 0.15,
  0.1, 0.3, 0.6,
  0.65, 0.2, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta5 <- matrix(c(
  0.2, 0.3, 0.5,
  0.55, 0.3, 0.15,
  0.15, 0.25, 0.6,
  0.6, 0.25, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta6 <- matrix(c(
  0.15, 0.35, 0.5,
  0.65, 0.2, 0.15,
  0.1, 0.3, 0.6,
  0.45, 0.35, 0.2
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta7 <- matrix(c(
  0.2, 0.25, 0.55,
  0.6, 0.25, 0.15,
  0.4, 0.35, 0.25,
  0.7, 0.2, 0.1
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta8 <- matrix(c(
  0.15, 0.30, 0.55,
  0.55, 0.25, 0.20,
  0.15, 0.30, 0.55,
  0.60, 0.25, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta9 <- matrix(c(
  0.2, 0.3, 0.5,
  0.65, 0.2, 0.15,
  0.1, 0.25, 0.65,
  0.55, 0.3, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta10 <- matrix(c(
  0.15, 0.35, 0.5,
  0.6, 0.2, 0.2,
  0.15, 0.3, 0.55,
  0.65, 0.2, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta11 <- matrix(c(
  0.2, 0.3, 0.5,
  0.55, 0.3, 0.15,
  0.45, 0.3, 0.25,
  0.6, 0.25, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta12 <- matrix(c(
  0.6, 0.25, 0.15,
  0.15, 0.3, 0.55,
  0.1, 0.3, 0.6,
  0.65, 0.2, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta13 <- matrix(c(
  0.55, 0.3, 0.15,
  0.2, 0.3, 0.5,
  0.15, 0.25, 0.6,
  0.45, 0.35, 0.2
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta14 <- matrix(c(
  0.65, 0.2, 0.15,
  0.15, 0.35, 0.5,
  0.1, 0.3, 0.6,
  0.6, 0.25, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta15 <- matrix(c(
  0.6, 0.25, 0.15,
  0.2, 0.3, 0.5,
  0.4, 0.35, 0.25,
  0.7, 0.2, 0.1
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta16 <- matrix(c(
  0.55, 0.25, 0.2,
  0.15, 0.3, 0.55,
  0.15, 0.3, 0.55,
  0.55, 0.3, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta17 <- matrix(c(
  0.6, 0.2, 0.2,
  0.2, 0.3, 0.5,
  0.1, 0.25, 0.65,
  0.65, 0.2, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta18 <- matrix(c(
  0.6, 0.2, 0.2,
  0.15, 0.35, 0.5,
  0.15, 0.3, 0.55,
  0.6, 0.25, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta19 <- matrix(c(
  0.55, 0.3, 0.15,
  0.2, 0.3, 0.5,
  0.45, 0.3, 0.25,
  0.65, 0.2, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta20 <- matrix(c(
  0.5, 0.3, 0.2,
  0.5, 0.3, 0.2,
  0.15, 0.30, 0.55,
  0.45, 0.35, 0.2
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta21 <- matrix(c(
  0.55, 0.3, 0.15,
  0.55, 0.25, 0.2,
  0.1, 0.3, 0.6,
  0.65, 0.2, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta22 <- matrix(c(
  0.5, 0.35, 0.15,
  0.55, 0.3, 0.15,
  0.4, 0.35, 0.25,
  0.6, 0.25, 0.15
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta23 <- matrix(c(
  0.5, 0.3, 0.2,
  0.5, 0.3, 0.2,
  0.15, 0.25, 0.6,
  0.7, 0.2, 0.1
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta24 <- matrix(c(
  0.55, 0.25, 0.2,
  0.5, 0.3, 0.2,
  0.45, 0.3, 0.25,
  0.45, 0.35, 0.2
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta25 <- matrix(c(
  0.5, 0.3, 0.2,
  0.55, 0.25, 0.2,
  0.4, 0.35, 0.25,
  0.65, 0.2, 0.15
), nrow = G, ncol = 3, byrow = TRUE)



#Non-informative variables 
sim3_theta26 <- matrix(rep(c(0.5,0.4,0.1),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta27 <- matrix(rep(c(0.3,0.1,0.6),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta28 <- matrix(rep(c(0.5, 0.2, 0.3),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta29 <- matrix(rep(c(0.2, 0.05, 0.75),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta30 <- matrix(rep(c(0.7, 0.05, 0.25),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta31 <- matrix(rep(c(0.25, 0.5, 0.25),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta32 <- matrix(rep(c(0.2, 0.3, 0.5),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta33 <- matrix(rep(c(0.6, 0.2, 0.2),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta34 <- matrix(rep(c(0.7, 0.1, 0.2),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta35 <- matrix(rep(c(0.65, 0.15, 0.2),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta36 <- matrix(rep(c(0.15, 0.2, 0.65),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta37 <- matrix(rep(c(0.45, 0.25, 0.3),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta38 <- matrix(rep(c(0.05, 0.35, 0.6),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta39 <- matrix(rep(c(0.85, 0.1, 0.05),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta40 <- matrix(rep(c(0.1, 0.35, 0.55),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta41 <- matrix(rep(c(0.3, 0.67, 0.03),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta42 <- matrix(rep(c(0.8, 0.1, 0.1),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta43 <- matrix(rep(c(1/G, 1/G, 1/G),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta44 <- matrix(rep(c(0.1, 0.5, 0.4),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta45 <- matrix(rep(c(0.4, 0.3, 0.3),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta46 <- matrix(rep(c(0.5, 0.2, 0.3),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta47 <- matrix(rep(c(0.1, 0.8, 0.1),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta48 <- matrix(rep(c(0.9, 0.05, 0.05),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta49 <- matrix(rep(c(0.1, 0.05, 0.85),G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta50 <- matrix(rep(c(0.99,0.005,0.005),G), nrow = G, ncol = 3, byrow = TRUE)


sim3_theta <- list(
  sim3_theta1,
  sim3_theta2,
  sim3_theta3,
  sim3_theta4,
  sim3_theta5,
  sim3_theta6,
  sim3_theta7,
  sim3_theta8,
  sim3_theta9,
  sim3_theta10,
  sim3_theta11,
  sim3_theta12,
  sim3_theta13,
  sim3_theta14,
  sim3_theta15,
  sim3_theta16,
  sim3_theta17,
  sim3_theta18,
  sim3_theta19,
  sim3_theta20,
  sim3_theta21,
  sim3_theta22,
  sim3_theta23,
  sim3_theta24,
  sim3_theta25,
  sim3_theta26,
  sim3_theta27,
  sim3_theta28,
  sim3_theta29,
  sim3_theta30,
  sim3_theta31,
  sim3_theta32,
  sim3_theta33,
  sim3_theta34,
  sim3_theta35,
  sim3_theta36,
  sim3_theta37,
  sim3_theta38,
  sim3_theta39,
  sim3_theta40,
  sim3_theta41,
  sim3_theta42,
  sim3_theta43,
  sim3_theta44,
  sim3_theta45,
  sim3_theta46,
  sim3_theta47,
  sim3_theta48,
  sim3_theta49,
  sim3_theta50
)


set.seed(129)

#Simulatiing dataset from using the above parameters
sim3_data <- LCR_sim_data_given_X(theta = sim3_theta, beta = sim3_beta, n_samples = 500, X = sim3_X)