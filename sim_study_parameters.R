
#Simulation 1 parameters

G <- 2
p <- 6          
M <- 8
K <- rep(3, 8)  

sim1_beta <- matrix(c(-0.5, 0.8, 1.2, 1, 0.4, 0, 0), nrow = p + 1)

sim1_theta1 <- matrix(c(0.15,0.25,0.6,0.7,0.2,0.1), nrow = G, ncol = K[1], byrow = TRUE)
sim1_theta2 <- matrix(c(0.2,0.35,0.45,0.55,0.3,0.15), nrow = G, ncol = K[2], byrow = TRUE)
sim1_theta3 <- matrix(c(0.1,0.15,0.75,0.8,0.15,0.05), nrow = G, ncol = K[3], byrow = TRUE)
sim1_theta4 <- matrix(c(0.25,0.4,0.35,0.45,0.35,0.2), nrow = G, ncol = K[4], byrow = TRUE)
sim1_theta5 <- matrix(c(0.4, 0.5, 0.1, 0.4, 0.5, 0.1), nrow = G, ncol = K[5], byrow = TRUE)
sim1_theta6 <- matrix(c(0.7, 0.1, 0.2, 0.7, 0.1, 0.2), nrow = G, ncol = K[6], byrow = TRUE)
sim1_theta7 <- matrix(c(0.1, 0.5, 0.4, 0.1, 0.5, 0.4), nrow = G, ncol = K[7], byrow = TRUE)
sim1_theta8 <- matrix(1/K[8], nrow = G, ncol = K[8])

sim1_theta <- list(sim1_theta1, sim1_theta2, sim1_theta3, sim1_theta4,
                   sim1_theta5, sim1_theta6, sim1_theta7, sim1_theta8)





#Simulation 2 parameters


G <- 3
p <- 6
M <- 13
K <- c(2,2,2,3,3,3,4,4,3,3,5,5,5)

# covariate effects matrix
sim2_beta <- matrix(c(0, 0, 1, -1, -1, 1, 0.5, -0.5, -0.4, 0.4, 0, 0, 0, 0), nrow = p+1, byrow = TRUE)



sim2_theta1 <- matrix(c(0.15, 0.6, 0.8, 0.85, 0.4, 0.2), nrow = G, ncol = K[1])
sim2_theta2 <- matrix(c(0.25, 0.45, 0.7, 0.75, 0.55, 0.3), nrow = G, ncol = K[2])
sim2_theta3 <- matrix(c(0.7, 0.2, 0.65, 0.3, 0.8, 0.35), nrow = G, ncol = K[3])
sim2_theta4 <- matrix(c(0.1, 0.35, 0.7, 0.25, 0.4, 0.2, 0.65, 0.25, 0.1), nrow = G, ncol = K[4])
sim2_theta5 <- matrix(c(0.2, 0.25, 0.65, 0.15, 0.6, 0.25, 0.65, 0.15, 0.1), nrow = G, ncol = K[5])
sim2_theta6 <- matrix(c(0.15, 0.5, 0.75, 0.2, 0.35, 0.15, 0.65, 0.15, 0.1), nrow = G, ncol = K[6])
sim2_theta7 <- matrix(c(0.1, 0.25, 0.6, 0.15, 0.35, 0.25, 0.25, 0.25, 0.1, 0.5, 0.15, 0.05), nrow = G, ncol = K[7])
sim2_theta8 <- matrix(c(0.15, 0.2, 0.55, 0.2, 0.45, 0.2, 0.2, 0.25, 0.15, 0.45, 0.1, 0.1), nrow = G, ncol = K[8])
sim2_theta9 <- matrix(c(0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1), nrow = G, ncol = K[9])
sim2_theta10 <- matrix(c(0.7, 0.7, 0.7, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2), nrow = G, ncol = K[10])
sim2_theta11 <- matrix(1/K[11], nrow = G, ncol = K[11])
sim2_theta12 <- matrix(c(0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.2, 0.2, 0.2, 0.25, 0.25, 0.25, 0.3, 0.3, 0.3), nrow = G, ncol = K[12])
sim2_theta13 <- matrix(c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = G, ncol = K[13])


sim2_theta <- list(sim2_theta1, sim2_theta2, sim2_theta3, sim2_theta4, sim2_theta5, sim2_theta6, sim2_theta7, sim2_theta8, sim2_theta9, sim2_theta10, sim2_theta11, sim2_theta12, sim2_theta13)







#Simulation 3 parameters



G <- 4
n <- 500

p <- 5

M <- 40
K <- rep(3, M)


sim3_beta <- matrix(c(
  -1.6, -1.2, -1.6,   # intercept
  1.0, -0.8,  0.6,   # standard normal
  1.2,  0.0, -1.0,   # binary 
  -1.5,  1.2,  0.8,   # binary 
  1.0, -1.2,  1.5,   # binary 
  0.0,  0.0,  0.0    # binary - non-informative
), nrow = p + 1, byrow = TRUE)



sim3_theta1 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta2 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.05, 0.25, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta3 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.20, 0.10, 0.70,
  0.55, 0.30, 0.15     # G4 slightly elevated
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta4 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,    # G3 only moderate
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta5 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.10, 0.20, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta6 <- matrix(c(
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05,
  0.10, 0.15, 0.75,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

# Items 7-12: G1 low, G2 high, G3 high, G4 low
sim3_theta7 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta8 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta9 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta10 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.15, 0.35, 0.50,
  0.55, 0.30, 0.15     # G4 slightly elevated
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta11 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.05, 0.05, 0.90,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta12 <- matrix(c(
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.15, 0.15, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

# Items 13-16: G1 moderate, G2 low, G3 high, G4 low
sim3_theta13 <- matrix(c(
  0.35, 0.40, 0.25,
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta14 <- matrix(c(
  0.35, 0.40, 0.25,
  0.80, 0.15, 0.05,
  0.10, 0.10, 0.80,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta15 <- matrix(c(
  0.35, 0.40, 0.25,
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,    # G3 only moderate
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta16 <- matrix(c(
  0.35, 0.40, 0.25,
  0.80, 0.15, 0.05,
  0.10, 0.30, 0.60,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

# Items 17-20: G1 low, G2 moderate, G3 high, G4 low
sim3_theta17 <- matrix(c(
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,
  0.10, 0.20, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta18 <- matrix(c(
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,
  0.10, 0.30, 0.60,
  0.55, 0.30, 0.15     # G4 slightly elevated
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta19 <- matrix(c(
  0.55, 0.30, 0.15,    # G1 slightly elevated
  0.35, 0.40, 0.25,
  0.10, 0.10, 0.80,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta20 <- matrix(c(
  0.80, 0.15, 0.05,
  0.35, 0.40, 0.25,
  0.10, 0.20, 0.70,
  0.80, 0.15, 0.05
), nrow = G, ncol = 3, byrow = TRUE)

# ---------------- Non-informative items 21-40 ----------------
sim3_theta21 <- matrix(rep(c(0.5, 0.4, 0.1), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta22 <- matrix(rep(c(0.5, 0.4, 0.1), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta23 <- matrix(rep(c(0.3, 0.1, 0.6), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta24 <- matrix(rep(c(0.5, 0.2, 0.3), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta25 <- matrix(rep(c(0.2, 0.05, 0.75), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta26 <- matrix(rep(c(0.7, 0.05, 0.25), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta27 <- matrix(rep(c(0.25, 0.5, 0.25), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta28 <- matrix(rep(c(0.2, 0.3, 0.5), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta29 <- matrix(rep(c(0.6, 0.2, 0.2), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta30 <- matrix(rep(c(0.7, 0.1, 0.2), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta31 <- matrix(rep(c(0.65, 0.15, 0.2), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta32 <- matrix(rep(c(0.15, 0.2, 0.65), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta33 <- matrix(rep(c(0.45, 0.25, 0.3), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta34 <- matrix(rep(c(0.05, 0.35, 0.6), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta35 <- matrix(rep(c(0.85, 0.1, 0.05), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta36 <- matrix(rep(c(0.1, 0.35, 0.55), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta37 <- matrix(rep(c(0.3, 0.67, 0.03), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta38 <- matrix(rep(c(0.8, 0.1, 0.1), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta39 <- matrix(rep(c(1/3, 1/3, 1/3), G), nrow = G, ncol = 3, byrow = TRUE)
sim3_theta40 <- matrix(rep(c(0.1, 0.5, 0.4), G), nrow = G, ncol = 3, byrow = TRUE)

sim3_theta <- list(
  sim3_theta1,  sim3_theta2,  sim3_theta3,  sim3_theta4,  sim3_theta5,
  sim3_theta6,  sim3_theta7,  sim3_theta8,  sim3_theta9,  sim3_theta10,
  sim3_theta11, sim3_theta12, sim3_theta13, sim3_theta14, sim3_theta15,
  sim3_theta16, sim3_theta17, sim3_theta18, sim3_theta19, sim3_theta20,
  sim3_theta21, sim3_theta22, sim3_theta23, sim3_theta24, sim3_theta25,
  sim3_theta26, sim3_theta27, sim3_theta28, sim3_theta29, sim3_theta30,
  sim3_theta31, sim3_theta32, sim3_theta33, sim3_theta34, sim3_theta35,
  sim3_theta36, sim3_theta37, sim3_theta38, sim3_theta39, sim3_theta40
)
