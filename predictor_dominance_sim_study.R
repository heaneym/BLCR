source('LCR_Gibbs.R')
source('LCR_sim_data.R')
library(ggplot2)
library(RColorBrewer)
library(dominanceanalysis)
library(reshape2) 
library(BayesLCA)

kl_divergence_rowwise <- function(p, q) {
  if (!identical(dim(p), dim(q))) {
    stop("Matrices must have the same dimensions")
  }
  epsilon <- 1e-10
  p <- pmax(p, epsilon)
  q <- pmax(q, epsilon)
  kl_values <- numeric(nrow(p))
  for (i in 1:nrow(p)) {
    kl_values[i] <- sum(p[i, ] * log(p[i, ] / q[i, ]))
  }
  kl_values[kl_values < 0 & abs(kl_values) < 1e-10] <- 0
  return(kl_values)
}


#This is a simulation study to assess how the Bayesian LCR variable selection handles cases where a 
#predictor has 0 coefficient for some groups but nonzero for others


n1 <- 75
n2 <- 150
n3 <- 500
n4 <- 1000

n <- 1000

G <- 4
M <- 21
sim_theta3 <- readRDS(file = 'C:/Users/matth/Documents/AIM CP Project/LCA/CSHQ_analysis/LCA_reduced_4group_itemprob')

#Creating a beta coefficient matrix with 7 predictors such that the variables have no effect for 
#some groups but varying nonzero effects for others.

# #Predictor X1: a binary variable such that all groups have nonzero effect (like the ASD predictor).
# #We set the proportion of observations with this variable to be. 
# 
X1 <- rbinom(n, 1, 0.5)

# #Predictor X2: continuous variable such that all groups have nonzero effect. 
# #Simulated from a standard normal
# 
X2 <- rnorm(n)

# #Predictor X3: binary variable which is informative for groups 2 and 4 but not 3
# 
X3 <- rbinom(n, 1, 0.5)

# #Predictor X4: continuous variable informative for groups 2 and 3 but not 4
# 
X4 <- rnorm(n)

# 
# #Predictor X5: continuous variable informative for group 2 only
# 
X5 <- rnorm(n)

# 
# #Predictor X6: binary variable informative for group 3 only
# 
X6 <- rbinom(n, 1, 0.5)

# 
# #predictor X7: continuous variable which is not informative at all.
# 
X7 <- rnorm(n)

# 
sim_X <- cbind(X1, X2, X3, X4, X5, X6, X7)


sim_beta <- matrix(c(0, -0.5, -0.5,   -1, 
                     0,    0.8,   -0.8,   0.8, 
                     0,    0.8,    0.8,   -0.8, 
                     0, 1,    0,    -1,
                     0, -1,   -1,    0,
                     0,  1,    0,    0,
                     0,    0,  0.8,    0,
                     0, 0, 0, 0), 
                   byrow = TRUE, nrow = 8, ncol = 4)

sim_true_logit_class_probs <- exp(cbind(1,sim_X)%*%sim_beta)/rowSums(exp(cbind(1,sim_X)%*%sim_beta))

sim_true_class_vec <- apply(sim_true_logit_class_probs, 1, function(x) sample(1:4, 1, prob = x))

sim_Y <- matrix(0, nrow = n, ncol = length(sim_theta3))

for(i in 1:n){
  class <- sim_true_class_vec[i]
  for (j in 1:M){
    theta_mat_row <- sim_theta3[[j]][class,]
    sim_Y[i,j] <- sample(1:K[j], 1, prob = theta_mat_row)
  }
}

small_subset_simulation_indices <- sample(1:n, 75, replace = FALSE) 
medium_subset_simulation_indices <- sample(1:n, 150, replace = FALSE) 
large_subset_simulation_indices <- sample(1:n, 500, replace = FALSE) 
verylarge_subset_simulation_indices <- sample(1:n, 1000, replace = FALSE) 

sim_X_small <- sim_X[small_subset_simulation_indices,]
sim_Y_small <- sim_Y[small_subset_simulation_indices,]
  
sim_X_medium <- sim_X[medium_subset_simulation_indices,]
sim_Y_medium <- sim_Y[medium_subset_simulation_indices,]
  
sim_X_large <- sim_X[large_subset_simulation_indices,]
sim_Y_large <- sim_Y[large_subset_simulation_indices,]
  
sim_X_verylarge <- sim_X[verylarge_subset_simulation_indices,]
sim_Y_verylarge <- sim_Y[verylarge_subset_simulation_indices,]

sim_coefficient_test_cov_sel_LCR_small <- LCR_Gibbs(X = sim_X_small, 
                                              Y = sim_Y_small, 
                                              G = 4,
                                              beta_prior_cov = diag(10^2, 8), 
                                              beta_prior_mean = rep(0,8), 
                                              theta_hyperparam = 1,
                                              clust_var_prior = 0.5,
                                              item.sel = FALSE,
                                              cov.sel = TRUE,
                                              verbose = TRUE,
                                              relabel = TRUE,
                                              n_samples = 5000,
                                              burnin = 1000,
                                              thinby = 10)

sim_coefficient_test_cov_sel_LCR_medium <- LCR_Gibbs(X = sim_X_medium, 
                                                    Y = sim_Y_medium, 
                                                    G = 4,
                                                    beta_prior_cov = diag(10^2, 8), 
                                                    beta_prior_mean = rep(0,8), 
                                                    theta_hyperparam = 1,
                                                    clust_var_prior = 0.5,
                                                    item.sel = FALSE,
                                                    cov.sel = TRUE,
                                                    verbose = TRUE,
                                                    relabel = TRUE,
                                                    n_samples = 5000,
                                                    burnin = 1000,
                                                    thinby = 10)


sim_coefficient_test_cov_sel_LCR_large <- LCR_Gibbs(X = sim_X_large, 
                                                    Y = sim_Y_large, 
                                                    G = 4,
                                                    beta_prior_cov = diag(10^2, 8), 
                                                    beta_prior_mean = rep(0,8), 
                                                    theta_hyperparam = 1,
                                                    clust_var_prior = 0.5,
                                                    item.sel = FALSE,
                                                    cov.sel = TRUE,
                                                    verbose = TRUE,
                                                    relabel = TRUE,
                                                    n_samples = 5000,
                                                    burnin = 1000,
                                                    thinby = 10)



sim_coefficient_test_cov_sel_LCR_verylarge <- LCR_Gibbs(X = sim_X_verylarge, 
                                                    Y = sim_Y_verylarge, 
                                                    G = 4,
                                                    beta_prior_cov = diag(10^2, 8), 
                                                    beta_prior_mean = rep(0,8), 
                                                    theta_hyperparam = 1,
                                                    clust_var_prior = 0.5,
                                                    item.sel = FALSE,
                                                    cov.sel = TRUE,
                                                    verbose = TRUE,
                                                    relabel = TRUE,
                                                    n_samples = 5000,
                                                    burnin = 1000,
                                                    thinby = 10)










#Carrying out dominance analysis for the logit component for comparison

data_df <- data.frame(
  outcome = as.factor(sim_true_class_vec),
  X1 = X1,
  X2 = X2,
  X3 = X3,
  X4 = X4,
  X5 = X5,
  X6 = X6,
  X7 = X7
)

full_model_ml_est <- multinom(outcome ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, data = data_df, trace = FALSE)
null_model_ml_est <- multinom(outcome ~ 1, data = data_df, trace = FALSE)
baseline_pseudo_r2 <- 1 - (logLik(full_model_ml_est) / logLik(null_model_ml_est))

predictor_names <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7")
n_predictors <- length(predictor_names)


calc_pseudo_r2 <- function(predictors_to_include) {
  if (length(predictors_to_include) == 0) {
    return(0)
  }
  formula_str <- paste("outcome ~", paste(predictors_to_include, collapse = " + "))
  model_formula <- as.formula(formula_str)
  model <- multinom(model_formula, data = data_df, trace = FALSE)
  pseudo_r2 <- 1 - (logLik(model) / logLik(null_model_ml_est))
  return(as.numeric(pseudo_r2))
}

all_subsets <- list()
subset_r2 <- c()


for (i in 1:n_predictors) {
  combinations <- combn(predictor_names, i, simplify = FALSE)
  all_subsets <- c(all_subsets, combinations)
}


for (i in 1:length(all_subsets)) {
  subset_r2[i] <- calc_pseudo_r2(all_subsets[[i]])
  if (i %% 10 == 0) cat("Completed", i, "of", length(all_subsets), "models\n")
}


dominance_scores <- rep(0, n_predictors)
names(dominance_scores) <- predictor_names


for (p in 1:n_predictors) {
  predictor <- predictor_names[p]
  total_contribution <- 0
  count <- 0
  for (i in 1:length(all_subsets)) {
    if (predictor %in% all_subsets[[i]]) {
      subset_without <- setdiff(all_subsets[[i]], predictor)
      if (length(subset_without) == 0) {
        contribution <- subset_r2[i] - 0
      } else {
        without_index <- which(sapply(all_subsets, function(x) setequal(x, subset_without)))
        if (length(without_index) > 0) {
          contribution <- subset_r2[i] - subset_r2[without_index[1]]
        } else {
          contribution <- 0
        }
      }
      
      total_contribution <- total_contribution + contribution
      count <- count + 1
    }
  }
  dominance_scores[p] <- total_contribution / count
}



#Testing and comparing to the KL divergence


assignment_probs_full_mu_sim_predsel_test <- cbind(1,sim_X)%*%sim_beta

assignment_probs_X1_mu_sim_predsel_test <- cbind(1,sim_X[,-1])%*%sim_beta[-2,]
assignment_probs_X2_mu_sim_predsel_test <- cbind(1,sim_X[,-2])%*%sim_beta[-3,]
assignment_probs_X3_mu_sim_predsel_test <- cbind(1,sim_X[,-3])%*%sim_beta[-4,]
assignment_probs_X4_mu_sim_predsel_test <- cbind(1,sim_X[,-4])%*%sim_beta[-5,]
assignment_probs_X5_mu_sim_predsel_test <- cbind(1,sim_X[,-5])%*%sim_beta[-6,]
assignment_probs_X6_mu_sim_predsel_test <- cbind(1,sim_X[,-6])%*%sim_beta[-7,]
assignment_probs_X7_mu_sim_predsel_test <- cbind(1,sim_X[,-7])%*%sim_beta[-8,]



assignment_matrix_full_predsel_test <- exp(assignment_probs_full_mu_sim_predsel_test)/rowSums(exp(assignment_probs_full_mu_sim_predsel_test))
assignment_matrix_X1_predsel_test <- exp(assignment_probs_X1_mu_sim_predsel_test)/rowSums(exp(assignment_probs_X1_mu_sim_predsel_test))
assignment_matrix_X2_predsel_test <- exp(assignment_probs_X2_mu_sim_predsel_test)/rowSums(exp(assignment_probs_X2_mu_sim_predsel_test))
assignment_matrix_X3_predsel_test <- exp(assignment_probs_X3_mu_sim_predsel_test)/rowSums(exp(assignment_probs_X3_mu_sim_predsel_test))
assignment_matrix_X4_predsel_test <- exp(assignment_probs_X4_mu_sim_predsel_test)/rowSums(exp(assignment_probs_X4_mu_sim_predsel_test))
assignment_matrix_X5_predsel_test <- exp(assignment_probs_X5_mu_sim_predsel_test)/rowSums(exp(assignment_probs_X5_mu_sim_predsel_test))
assignment_matrix_X6_predsel_test <- exp(assignment_probs_X6_mu_sim_predsel_test)/rowSums(exp(assignment_probs_X6_mu_sim_predsel_test))
assignment_matrix_X7_predsel_test <- exp(assignment_probs_X7_mu_sim_predsel_test)/rowSums(exp(assignment_probs_X7_mu_sim_predsel_test))






X1_kl_divergence_predsel_test <- kl_divergence_rowwise(assignment_matrix_full_predsel_test, assignment_matrix_X1_predsel_test)
X2_kl_divergence_predsel_test <- kl_divergence_rowwise(assignment_matrix_full_predsel_test, assignment_matrix_X2_predsel_test)
X3_kl_divergence_predsel_test <- kl_divergence_rowwise(assignment_matrix_full_predsel_test, assignment_matrix_X3_predsel_test)
X4_kl_divergence_predsel_test <- kl_divergence_rowwise(assignment_matrix_full_predsel_test, assignment_matrix_X4_predsel_test)
X5_kl_divergence_predsel_test <- kl_divergence_rowwise(assignment_matrix_full_predsel_test, assignment_matrix_X5_predsel_test)
X6_kl_divergence_predsel_test <- kl_divergence_rowwise(assignment_matrix_full_predsel_test, assignment_matrix_X6_predsel_test)
X7_kl_divergence_predsel_test <- kl_divergence_rowwise(assignment_matrix_full_predsel_test, assignment_matrix_X7_predsel_test)








