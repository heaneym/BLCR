working_dir <- getwd()
source('LCR_Gibbs.R')
source('LCR_sim_data.R')

#Simulation Study 1

# 2 groups
# 8 item variables, 4 of which are useful for clustering (all item variables have 3 possible outcomes)
# 6 covariates, of varying effect size
# 300 observations

# structure parameters
G <- 2
p <- 6
M <- 8
K <- rep(3,8)

# covariate effects vector
sim1_beta <- matrix(c(-0.5, 0.8, 1.2, 1, 0.4, 0, 0), nrow = p+1)

#item probability parameters
sim1_theta1 <- matrix(c(0.15,0.25,0.6,0.7,0.2,0.1), nrow = G, ncol = K[1], byrow = TRUE)
sim1_theta2 <- matrix(c(0.2,0.35,0.45,0.55,0.3,0.15), nrow = G, ncol = K[2], byrow = TRUE)
sim1_theta3 <- matrix(c(0.1,0.15,0.75,0.8,0.15,0.05), nrow = G, ncol = K[3], byrow = TRUE)
sim1_theta4 <- matrix(c(0.25,0.4,0.35,0.45,0.35,0.2), nrow = G, ncol = K[4], byrow = TRUE)
sim1_theta5 <- matrix(c(0.4, 0.5, 0.1, 0.4, 0.5, 0.1), nrow = G, ncol = K[5], byrow = TRUE)
sim1_theta6 <- matrix(c(0.7, 0.1, 0.2, 0.7, 0.1, 0.2), nrow = G, ncol = K[6], byrow = TRUE)
sim1_theta7 <- matrix(c(0.1, 0.5, 0.4, 0.1, 0.5, 0.4), nrow = G, ncol = K[7], byrow = TRUE)
sim1_theta8 <- matrix(1/K[8], nrow = G, ncol = K[8])

sim1_theta <- list(sim1_theta1, sim1_theta2, sim1_theta3, sim1_theta4, sim1_theta5, sim1_theta6, sim1_theta7, sim1_theta8)

sim1_data <- LCR_sim_data(theta = sim1_theta, beta = sim1_beta, n_samples = 300, n_noise_var = 0)

sim1_LCR_fit <- LCR_Gibbs(X = sim1_data$X, Y = sim1_data$Y, G = 2, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                          theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                          relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


#Simulation Study 2

# 3 groups
# 13 item variables, 8 of which are useful for clustering 
# variables 1,2,3 have 2 possible responses, 4,5,6 have 3 possible responses, and 7 and 8 have 4 possible responses
# variables 9,10 have 3 possible responses, 11,12,13 have 5 possible responses
# 6 covariates, of varying effect size for different groups
# 300 observations

# structure parameters
G <- 3
p <- 6
M <- 13
K <- c(2,2,2,3,3,3,4,4,3,3,5,5,5)

# covariate effects matrix
sim2_beta <- matrix(c(-1, 0.8, 1.1, 1.2, 0.3, 0, 0, -0.3, 0.5, 0, 0.4, 0, 0, 0), nrow = p+1)


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


sim2_data <- LCR_sim_data(theta = sim2_theta, beta = sim2_beta, n_samples = 300, n_noise_var = 0)


sim2_LCR_fit <- LCR_Gibbs(X = sim2_data$X, Y = sim2_data$Y, G = 3, beta_prior_cov = diag(10^2,7), beta_prior_mean = rep(0,7), 
                          theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                          relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)





#Running with the CSHQ sleep survey data, we consider as covariates the variables 
setwd(working_dir)
setwd("..") 
CSHQ_data <- read.spss('./CSHQ_analysis/LynnWalsh thesis data files/Thesis Cleaned File.sav')
setwd(working_dir) 

data <- as.data.frame(data)
sleep_data_responses <- data[,65:99]

mapping <- c(
  "Rarely (0-1 times per week)" = 1,
  "Sometimes (2-4 times per week)" = 2,
  "Usually (5-7 times per week)" = 3
)

sleep_data_responses_numeric <- sleep_data_responses %>% mutate_all(~ mapping[as.character(.)])
#removing bedtime resistance 4 and 6 since they are double counted from sleep anxiety
sleep_data_responses_numeric <- sleep_data_responses_numeric[,-c(4,6)]


df <- cbind(data[,c(35:48,56)],sleep_data_responses_numeric)
df$Diagnosis_ID <- ifelse(is.na(df$Diagnosis_ID),0,1)
df$Diagnosis_CerebralPalsy <- ifelse(is.na(df$Diagnosis_CerebralPalsy),0,1)
df$Diagnosis_ASD <- ifelse(is.na(df$Diagnosis_ASD),0,1)
df$Diagnosis_ADHD <- ifelse(is.na(df$Diagnosis_ADHD),0,1)
df$Diagnosis_DS <- ifelse(is.na(df$Diagnosis_DS),0,1)
df$Diagnosis_Other <- ifelse(is.na(df$Diagnosis_Other),0,1)
df$Diagnosis_None <- ifelse(is.na(df$Diagnosis_None),0,1)


#I'm gonna remove rows that are almost all NA (all indicator variables are missing)

df <- df[complete.cases(df[,16]),]
colnames(df)[16:48] <-c('BR1','BR2','BR3','BR5','SOD1','SD1','SD2','SD3','SA1','SA2','SA3','SA4','NW1','NW2','NW3','P1','P2','P3','P4','P5','P6','P7','SDB1', 'SDB2', 'SDB3', 'DS1','DS2','DS3','DS4','DS5','DS6','DS7','DS8')

#reverse scoring the negative items

df$BR1 <- 4 - df$BR1
df$BR2 <- 4 - df$BR2
df$SOD1 <- 4 - df$SOD1
df$SD2 <- 4 - df$SD2
df$SD3 <- 4 - df$SD3
df$DS1 <- 4 - df$DS1

sleep_data_total_score <- rowSums(df[,16:48])
sleep_disorders <- ifelse(sleep_data_total_score >= 41,1,0)
df <- cbind(df,sleep_data_total_score,sleep_disorders)
df <- df[,-c(8,11)]

# We want to look at a variable for ASD (diagnosis_ASD), one for Intellectual Disability (diagnosis_ID), and then other, which includes anything else.
df$other_diagnoses <- df$Diagnosis_Other + df$Diagnosis_ADHD
df$other_diagnoses <- ifelse(df$other_diagnoses>0, 1, 0)

#Need to check but I'm not sure how relevant the asthma one is here
df$other_diagnoses[107] <- 0

df$interaction_ASD_ID <- (df$Diagnosis_ASD)*(df$Diagnosis_ID)
df$interaction_ASD_other <- (df$Diagnosis_ASD)*(df$other_diagnoses)
df$interaction_ID_other <- (df$Diagnosis_ID)*(df$other_diagnoses)

nd_vec <- rowSums(df[,8:10])
nd_vec[c(83,107,143)] <- 0
nd_vec[which(nd_vec==2 | nd_vec == 3) ] <- 1
nd_vec <- as.numeric(nd_vec)
df <- cbind(df, nd_vec)

CSHQ_Y <- as.matrix(df[,14:46])
#should I include interaction terms here?
df$scaled_age <- scale(df$ChildAge)
CSHQ_X <- as.matrix(cbind(df$Diagnosis_ID, df$Diagnosis_ASD, df$other_diagnoses, df$scaled_age))



CSHQ_LCR_fit_6group <- LCR_Gibbs(X = CSHQ_X, Y = CSHQ_Y, G = 6, beta_prior_cov = diag(10^2,5), beta_prior_mean = rep(0,5), 
                          theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                          relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

CSHQ_LCR_fit_5group <- LCR_Gibbs(X = CSHQ_X, Y = CSHQ_Y, G = 5, beta_prior_cov = diag(10^2,5), beta_prior_mean = rep(0,5), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)

CSHQ_LCR_fit_4group <- LCR_Gibbs(X = CSHQ_X, Y = CSHQ_Y, G = 4, beta_prior_cov = diag(10^2,5), beta_prior_mean = rep(0,5), 
                                 theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = TRUE, cov.sel = TRUE, verbose = TRUE, 
                                 relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)



