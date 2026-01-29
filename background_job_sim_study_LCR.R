library(MCMCpack)
library(BayesLogit)
library(label.switching)
source('LCR_Gibbs.R')

# sim_coefficient_test_cov_sel_LCR <- LCR_Gibbs(X = sim_X, 
#                                               Y = sim_Y, 
#                                               G = 4,
#                                               beta_prior_cov = diag(10^2, 8), 
#                                               beta_prior_mean = rep(0,8), 
#                                               theta_hyperparam = 1,
#                                               clust_var_prior = 0.5,
#                                               item.sel = FALSE,
#                                               cov.sel = TRUE,
#                                               verbose = TRUE,
#                                               relabel = TRUE,
#                                               n_samples = 5000,
#                                               burnin = 1000,
#                                               thinby = 10)




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