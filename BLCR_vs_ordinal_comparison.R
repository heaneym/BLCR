#Script for comparing clustering performance of LCA with the ordinal model

load("C:/Users/matth/Documents/AIM CP Project/LCA/ordinal_clustering_model/ordinal_sim1_data.RData")
load("C:/Users/matth/Documents/AIM CP Project/LCA/ordinal_clustering_model/ordinal_sim2_data.RData")


BLCR_ordinal_comparison_sim1_fit <- LCR_Gibbs(X = ordinal_sim1$X, Y = ordinal_sim1$Y, G = 2, beta_prior_cov = diag(10^2,5), beta_prior_mean = rep(0,5), 
                                              theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                              relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)
  
  
BLCR_ordinal_comparison_sim2_fit <- LCR_Gibbs(X = ordinal_sim2$X, Y = ordinal_sim2$Y, G = 4, beta_prior_cov = diag(10^2,5), beta_prior_mean = rep(0,5), 
                                              theta_hyperparam = 1, clust_var_prior = 0.5, item.sel = FALSE, cov.sel = TRUE, verbose = TRUE, 
                                              relabel = TRUE, n_samples = 5000, burnin = 1000, thinby = 10)


#Computing point estimates using minVI

ordinal_sim1_LCR_cluster_mat <- apply(BLCR_ordinal_comparison_sim1_fit$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
}) 

ordinal_sim1_LCR_psm_mat <- comp.psm(t(ordinal_sim1_LCR_cluster_mat))

ordinal_sim1_LCR_minVI_clust <- minVI(psm = ordinal_sim1_LCR_psm_mat, method = 'greedy', start.cl = max.col(BLCR_ordinal_comparison_sim1_fit$assignment_prob))






ordinal_sim2_LCR_cluster_mat <- apply(BLCR_ordinal_comparison_sim2_fit$samples$z_samples, 3, function(matrix_slice) {
  apply(matrix_slice, 1, which.max)
}) 

ordinal_sim2_LCR_psm_mat <- comp.psm(t(ordinal_sim2_LCR_cluster_mat))

ordinal_sim2_LCR_minVI_clust <- minVI(psm = ordinal_sim2_LCR_psm_mat, method = 'greedy', start.cl = max.col(BLCR_ordinal_comparison_sim2_fit$assignment_prob))



#Computing classification entropy

mean_entropy_LCR_ordinal_sim1 <- mean(-rowSums(BLCR_ordinal_comparison_sim1_fit$Z * log(BLCR_ordinal_comparison_sim1_fit$Z + 1e-12)))
mean_entropy_LCR_ordinal_sim2 <- mean(-rowSums(BLCR_ordinal_comparison_sim2_fit$Z * log(BLCR_ordinal_comparison_sim2_fit$Z + 1e-12)))