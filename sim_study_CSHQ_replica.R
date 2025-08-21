library(BayesLCA)
library(poLCA)

replica_sim_theta_4group_reduced <- readRDS(file = 'C:/Users/matth/Documents/AIM CP Project/LCA/CSHQ_analysis/sim_theta_4group_reduced')
replica_sim_theta_6group_reduced <- readRDS(file = 'C:/Users/matth/Documents/AIM CP Project/LCA/CSHQ_analysis/sim_theta_6group_reduced')

replica_sim_pi_4group_reduced <- readRDS(file = 'C:/Users/matth/Documents/AIM CP Project/LCA/CSHQ_analysis/sim_pi_4group_reduced') 
replica_sim_pi_6group_reduced <- readRDS(file = 'C:/Users/matth/Documents/AIM CP Project/LCA/CSHQ_analysis/sim_pi_6group_reduced') 


# replica_CSHQ_sim_data_6group1 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# replica_CSHQ_sim_data_6group2 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# replica_CSHQ_sim_data_6group3 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# replica_CSHQ_sim_data_6group4 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# replica_CSHQ_sim_data_6group5 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# replica_CSHQ_sim_data_6group6 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# replica_CSHQ_sim_data_6group7 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# replica_CSHQ_sim_data_6group8 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# replica_CSHQ_sim_data_6group9 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# replica_CSHQ_sim_data_6group10 <- poLCA.simdata(N = 150, probs = replica_sim_theta_6group_reduced, nclass = 6, P = replica_sim_pi_6group_reduced)
# 
# 
# replica_CSHQ_sim_data_4group1 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# replica_CSHQ_sim_data_4group2 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# replica_CSHQ_sim_data_4group3 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# replica_CSHQ_sim_data_4group4 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# replica_CSHQ_sim_data_4group5 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# replica_CSHQ_sim_data_4group6 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# replica_CSHQ_sim_data_4group7 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# replica_CSHQ_sim_data_4group8 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# replica_CSHQ_sim_data_4group9 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# replica_CSHQ_sim_data_4group10 <- poLCA.simdata(N = 150, probs = replica_sim_theta_4group_reduced, nclass = 4, P = replica_sim_pi_4group_reduced)
# 
# CSHQ_Y_replica_4group1 <- as.matrix(replica_CSHQ_sim_data_4group1$dat)
# CSHQ_Y_replica_4group2 <- as.matrix(replica_CSHQ_sim_data_4group2$dat)
# CSHQ_Y_replica_4group3 <- as.matrix(replica_CSHQ_sim_data_4group3$dat)
# CSHQ_Y_replica_4group4 <- as.matrix(replica_CSHQ_sim_data_4group4$dat)
# CSHQ_Y_replica_4group5 <- as.matrix(replica_CSHQ_sim_data_4group5$dat)
# CSHQ_Y_replica_4group6 <- as.matrix(replica_CSHQ_sim_data_4group6$dat)
# CSHQ_Y_replica_4group7 <- as.matrix(replica_CSHQ_sim_data_4group7$dat)
# CSHQ_Y_replica_4group8 <- as.matrix(replica_CSHQ_sim_data_4group8$dat)
# CSHQ_Y_replica_4group9 <- as.matrix(replica_CSHQ_sim_data_4group9$dat)
# CSHQ_Y_replica_4group10 <- as.matrix(replica_CSHQ_sim_data_4group10$dat)
# 
# CSHQ_Y_replica_6group1 <- as.matrix(replica_CSHQ_sim_data_6group1$dat)
# CSHQ_Y_replica_6group2 <- as.matrix(replica_CSHQ_sim_data_6group2$dat)
# CSHQ_Y_replica_6group3 <- as.matrix(replica_CSHQ_sim_data_6group3$dat)
# CSHQ_Y_replica_6group4 <- as.matrix(replica_CSHQ_sim_data_6group4$dat)
# CSHQ_Y_replica_6group5 <- as.matrix(replica_CSHQ_sim_data_6group5$dat)
# CSHQ_Y_replica_6group6 <- as.matrix(replica_CSHQ_sim_data_6group6$dat)
# CSHQ_Y_replica_6group7 <- as.matrix(replica_CSHQ_sim_data_6group7$dat)
# CSHQ_Y_replica_6group8 <- as.matrix(replica_CSHQ_sim_data_6group8$dat)
# CSHQ_Y_replica_6group9 <- as.matrix(replica_CSHQ_sim_data_6group9$dat)
# CSHQ_Y_replica_6group10 <- as.matrix(replica_CSHQ_sim_data_6group10$dat)
# 
# 
# 
# 
# LCA_collapsed_CSHQ_replica_4group1 <- blca.collapsed(X = CSHQ_Y_replica_4group1, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# LCA_collapsed_CSHQ_replica_4group2 <- blca.collapsed(X = CSHQ_Y_replica_4group2, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# LCA_collapsed_CSHQ_replica_4group3 <- blca.collapsed(X = CSHQ_Y_replica_4group3, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# LCA_collapsed_CSHQ_replica_4group4 <- blca.collapsed(X = CSHQ_Y_replica_4group4, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# LCA_collapsed_CSHQ_replica_4group5 <- blca.collapsed(X = CSHQ_Y_replica_4group5, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# LCA_collapsed_CSHQ_replica_4group6 <- blca.collapsed(X = CSHQ_Y_replica_4group6, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# LCA_collapsed_CSHQ_replica_4group7 <- blca.collapsed(X = CSHQ_Y_replica_4group7, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# LCA_collapsed_CSHQ_replica_4group8 <- blca.collapsed(X = CSHQ_Y_replica_4group8, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# LCA_collapsed_CSHQ_replica_4group9 <- blca.collapsed(X = CSHQ_Y_replica_4group9, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# LCA_collapsed_CSHQ_replica_4group10 <- blca.collapsed(X = CSHQ_Y_replica_4group10, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
# 
# 
# 
# 
# LCA_collapsed_CSHQ_replica_6group1 <- blca.collapsed(X = CSHQ_Y_replica_6group1, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# LCA_collapsed_CSHQ_replica_6group2 <- blca.collapsed(X = CSHQ_Y_replica_6group2, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# LCA_collapsed_CSHQ_replica_6group3 <- blca.collapsed(X = CSHQ_Y_replica_6group3, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# LCA_collapsed_CSHQ_replica_6group4 <- blca.collapsed(X = CSHQ_Y_replica_6group4, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# LCA_collapsed_CSHQ_replica_6group5 <- blca.collapsed(X = CSHQ_Y_replica_6group5, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# LCA_collapsed_CSHQ_replica_6group6 <- blca.collapsed(X = CSHQ_Y_replica_6group6, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# LCA_collapsed_CSHQ_replica_6group7 <- blca.collapsed(X = CSHQ_Y_replica_6group7, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# LCA_collapsed_CSHQ_replica_6group8 <- blca.collapsed(X = CSHQ_Y_replica_6group8, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# LCA_collapsed_CSHQ_replica_6group9 <- blca.collapsed(X = CSHQ_Y_replica_6group9, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# LCA_collapsed_CSHQ_replica_6group10 <- blca.collapsed(X = CSHQ_Y_replica_6group10, G = 1, G.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
# 
# 
# cluster_mat_collapsed_CSHQ_replica_4group1 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group1$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_4group2 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group2$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_4group3 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group3$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_4group4 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group4$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_4group5 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group5$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_4group6 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group6$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_4group7 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group7$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_4group8 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group8$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_4group9 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group9$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_4group10 <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group10$samples$labels) 
# 
# 
# 
# cluster_mat_collapsed_CSHQ_replica_6group1 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group1$samples$labels)
# cluster_mat_collapsed_CSHQ_replica_6group2 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group2$samples$labels)
# cluster_mat_collapsed_CSHQ_replica_6group3 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group3$samples$labels)
# cluster_mat_collapsed_CSHQ_replica_6group4 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group4$samples$labels)
# cluster_mat_collapsed_CSHQ_replica_6group5 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group5$samples$labels)
# cluster_mat_collapsed_CSHQ_replica_6group6 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group6$samples$labels)
# cluster_mat_collapsed_CSHQ_replica_6group7 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group7$samples$labels)
# cluster_mat_collapsed_CSHQ_replica_6group8 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group8$samples$labels)
# cluster_mat_collapsed_CSHQ_replica_6group9 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group9$samples$labels)
# cluster_mat_collapsed_CSHQ_replica_6group10 <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group10$samples$labels)
# 
# 
# psm_collapsed_CSHQ_replica_4group1 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group1)
# psm_collapsed_CSHQ_replica_4group2 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group2)
# psm_collapsed_CSHQ_replica_4group3 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group3)
# psm_collapsed_CSHQ_replica_4group4 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group4)
# psm_collapsed_CSHQ_replica_4group5 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group5)
# psm_collapsed_CSHQ_replica_4group6 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group6)
# psm_collapsed_CSHQ_replica_4group7 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group7)
# psm_collapsed_CSHQ_replica_4group8 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group8)
# psm_collapsed_CSHQ_replica_4group9 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group9)
# psm_collapsed_CSHQ_replica_4group10 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group10)
# 
# 
# 
# psm_collapsed_CSHQ_replica_6group1 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group1)
# psm_collapsed_CSHQ_replica_6group2 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group2)
# psm_collapsed_CSHQ_replica_6group3 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group3)
# psm_collapsed_CSHQ_replica_6group4 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group4)
# psm_collapsed_CSHQ_replica_6group5 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group5)
# psm_collapsed_CSHQ_replica_6group6 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group6)
# psm_collapsed_CSHQ_replica_6group7 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group7)
# psm_collapsed_CSHQ_replica_6group8 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group8)
# psm_collapsed_CSHQ_replica_6group9 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group9)
# psm_collapsed_CSHQ_replica_6group10 <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group10)
# 
# minVI_collapsed_CSHQ_replica_4group1 <- minVI(psm_collapsed_CSHQ_replica_4group1, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group1$Z))
# minVI_collapsed_CSHQ_replica_4group2 <- minVI(psm_collapsed_CSHQ_replica_4group2, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group2$Z))
# minVI_collapsed_CSHQ_replica_4group3 <- minVI(psm_collapsed_CSHQ_replica_4group3, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group3$Z))
# minVI_collapsed_CSHQ_replica_4group4 <- minVI(psm_collapsed_CSHQ_replica_4group4, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group4$Z))
# minVI_collapsed_CSHQ_replica_4group5 <- minVI(psm_collapsed_CSHQ_replica_4group5, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group5$Z))
# minVI_collapsed_CSHQ_replica_4group6 <- minVI(psm_collapsed_CSHQ_replica_4group6, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group6$Z))
# minVI_collapsed_CSHQ_replica_4group7 <- minVI(psm_collapsed_CSHQ_replica_4group7, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group7$Z))
# minVI_collapsed_CSHQ_replica_4group8 <- minVI(psm_collapsed_CSHQ_replica_4group8, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group8$Z))
# minVI_collapsed_CSHQ_replica_4group9 <- minVI(psm_collapsed_CSHQ_replica_4group9, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group9$Z))
# minVI_collapsed_CSHQ_replica_4group10 <- minVI(psm_collapsed_CSHQ_replica_4group10, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group10$Z))
#  
# minVI_collapsed_CSHQ_replica_6group1 <- minVI(psm_collapsed_CSHQ_replica_6group1, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group1$Z))
# minVI_collapsed_CSHQ_replica_6group2 <- minVI(psm_collapsed_CSHQ_replica_6group2, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group2$Z))
# minVI_collapsed_CSHQ_replica_6group3 <- minVI(psm_collapsed_CSHQ_replica_6group3, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group3$Z))
# minVI_collapsed_CSHQ_replica_6group4 <- minVI(psm_collapsed_CSHQ_replica_6group4, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group4$Z))
# minVI_collapsed_CSHQ_replica_6group5 <- minVI(psm_collapsed_CSHQ_replica_6group5, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group5$Z))
# minVI_collapsed_CSHQ_replica_6group6 <- minVI(psm_collapsed_CSHQ_replica_6group6, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group6$Z))
# minVI_collapsed_CSHQ_replica_6group7 <- minVI(psm_collapsed_CSHQ_replica_6group7, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group7$Z))
# minVI_collapsed_CSHQ_replica_6group8 <- minVI(psm_collapsed_CSHQ_replica_6group8, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group8$Z))
# minVI_collapsed_CSHQ_replica_6group9 <- minVI(psm_collapsed_CSHQ_replica_6group9, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group9$Z))
# minVI_collapsed_CSHQ_replica_6group10 <- minVI(psm_collapsed_CSHQ_replica_6group10, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group10$Z))
#  
#  

#creating a function that simulates from a specified number of LCA models with theta and pi parameters, with a specified number of observations for each dataset
#It then fits a collapsed model to each of these datasets, then clusters using the minVI method of Wade and Ghahramani, then also fits the credible ball
#We want to return as much info as possible, so I'll return a list containing the collapsed outputs, the minVI outputs, as well as the credible ball outputs. We also compute summary things like average number of groups etc. 

sim_study_replica_data_minVI_clusters <- function(obs_per_dataset, theta, pi){
  n_datasets <- length(obs_per_dataset)
  total_datapoints <- sum(obs_per_dataset)
  n_group <- length(pi)
  total_data <- poLCA.simdata(N = total_datapoints, probs = theta, nclass = n_group, P = pi)
  
  # Initialize lists
  Y_list <- list()
  collapsed_output_list <- list()
  collapsed_G_proportions_list <- list()
  cluster_mat_list <- list()
  psm_mat_list <- list()
  minVI_output_list <- list()
  minVI_clustering_list <- list()
  minVI_cluster_table_list <- list()
  minVI_cluster_proportions_list <- list()
  fitted_G_val_list <- list()
  cred_ball_list <- list()
  
  subset_indices <- cumsum(c(0,obs_per_dataset))
  failed_datasets <- c()  # Track which datasets failed
  
  for(i in 1:n_datasets){
    print(paste('Dataset ', i))
    
    tryCatch({
      Y_list[[i]] <- as.matrix(total_data$dat[(subset_indices[i] + 1):subset_indices[i+1],])
      collapsed_output_list[[i]] <- blca.collapsed(X = Y_list[[i]], G = 1, G.sel = TRUE, var.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE) 
      collapsed_G_proportions_list[[i]] <- table(collapsed_output_list[[i]]$samples$G)/sum(table(collapsed_output_list[[i]]$samples$G))
      cluster_mat_list[[i]] <- do.call(rbind, collapsed_output_list[[i]]$samples$labels)
      psm_mat_list[[i]] <- comp.psm(cluster_mat_list[[i]])
      minVI_output_list[[i]] <- minVI(psm_mat_list[[i]], method = 'greedy', suppress.comment = FALSE, start.cl = max.col(collapsed_output_list[[i]]$Z))
      minVI_clustering_list[[i]] <- minVI_output_list[[i]]$cl
      minVI_cluster_table_list[[i]] <- table(minVI_clustering_list[[i]])
      minVI_cluster_proportions_list[[i]] <- minVI_cluster_table_list[[i]]/sum(minVI_cluster_table_list[[i]])
      fitted_G_val_list[[i]] <- length(minVI_cluster_table_list[[i]])
      cred_ball_list[[i]] <- credibleball(c.star = minVI_clustering_list[[i]], cls.draw = cluster_mat_list[[i]], c.dist = 'VI')
      
    }, error = function(e) {
      print(paste("ERROR in dataset", i, ":", e$message))
      print("Skipping this dataset and continuing...")
      failed_datasets <<- c(failed_datasets, i)  
      Y_list[[i]] <<- NULL
      collapsed_output_list[[i]] <<- NULL
      collapsed_G_proportions_list[[i]] <<- NULL
      cluster_mat_list[[i]] <<- NULL
      psm_mat_list[[i]] <<- NULL
      minVI_output_list[[i]] <<- NULL
      minVI_clustering_list[[i]] <<- NULL
      minVI_cluster_table_list[[i]] <<- NULL
      minVI_cluster_proportions_list[[i]] <<- NULL
      fitted_G_val_list[[i]] <<- NULL
      cred_ball_list[[i]] <<- NULL
    })
  }
  successful_G_vals <- fitted_G_val_list[!sapply(fitted_G_val_list, is.null)]
  mean_G_across_datasets <- mean(unlist(successful_G_vals))
  
  print(paste("Total datasets processed:", n_datasets))
  print(paste("Failed datasets:", length(failed_datasets)))
  if(length(failed_datasets) > 0) {
    print(paste("Failed dataset indices:", paste(failed_datasets, collapse = ", ")))
  }
  
  return(list(Y_list = Y_list, collapsed_output = collapsed_output_list, 
              collapsed_G_proportions = collapsed_G_proportions_list, cluster_matrices = cluster_mat_list, 
              psm_matrices = psm_mat_list, minVI_outputs = minVI_output_list, 
              minVI_clusterings = minVI_clustering_list, minVI_cluster_sizes = minVI_cluster_table_list,
              minVI_cluster_proportions = minVI_cluster_proportions_list,
              clustered_G_vals = fitted_G_val_list, cred_ball_outputs = cred_ball_list,
              mean_G_minVI = mean_G_across_datasets, failed_datasets = failed_datasets))
}

 

sim_study_replica_CSHQ_4group <- sim_study_replica_data_minVI_clusters(obs_per_dataset = rep(150, 150), theta = replica_sim_theta_4group_reduced, pi = replica_sim_pi_4group_reduced)
sim_study_replica_CSHQ_6group <- sim_study_replica_data_minVI_clusters(obs_per_dataset = rep(150, 150), theta = replica_sim_theta_6group_reduced, pi = replica_sim_pi_6group_reduced)



# LCA_4groupfit_CSHQ_replica_4group <- blca.gibbs(X = CSHQ_Y_replica_4group, G = 4, verbose = TRUE)
# LCA_5groupfit_CSHQ_replica_4group <- blca.gibbs(X = CSHQ_Y_replica_4group, G = 5, verbose = TRUE)
# LCA_6groupfit_CSHQ_replica_4group <- blca.gibbs(X = CSHQ_Y_replica_4group, G = 6, verbose = TRUE)
#   
# LCA_4groupfit_CSHQ_replica_6group <- blca.gibbs(X = CSHQ_Y_replica_6group, G = 4, verbose = TRUE)
# LCA_5groupfit_CSHQ_replica_6group <- blca.gibbs(X = CSHQ_Y_replica_6group, G = 5, verbose = TRUE)
# LCA_6groupfit_CSHQ_replica_6group <- blca.gibbs(X = CSHQ_Y_replica_6group, G = 6, verbose = TRUE)
# 
# cluster_mat_collapsed_CSHQ_replica_4group <- do.call(rbind, LCA_collapsed_CSHQ_replica_4group$samples$labels) 
# cluster_mat_collapsed_CSHQ_replica_6group <- do.call(rbind, LCA_collapsed_CSHQ_replica_6group$samples$labels) 
# 
# 
# cluster_mat_4groupfit_CSHQ_4group <- LCA_4groupfit_CSHQ_replica_4group$samples$labels
# cluster_mat_5groupfit_CSHQ_4group <- LCA_5groupfit_CSHQ_replica_4group$samples$labels
# cluster_mat_6groupfit_CSHQ_4group <- LCA_6groupfit_CSHQ_replica_4group$samples$labels
# 
# cluster_mat_4groupfit_CSHQ_6group <- LCA_4groupfit_CSHQ_replica_6group$samples$labels
# cluster_mat_5groupfit_CSHQ_6group <- LCA_5groupfit_CSHQ_replica_6group$samples$labels
# cluster_mat_6groupfit_CSHQ_6group <- LCA_6groupfit_CSHQ_replica_6group$samples$labels
# 
# 
# 
# psm_collapsed_CSHQ_replica_4group <- comp.psm(cluster_mat_collapsed_CSHQ_replica_4group)
# psm_collapsed_CSHQ_replica_6group <- comp.psm(cluster_mat_collapsed_CSHQ_replica_6group) 
# 
# psm_4groupfit_CSHQ_replica_4group <- comp.psm(cluster_mat_4groupfit_CSHQ_4group)
# psm_5groupfit_CSHQ_replica_4group <- comp.psm(cluster_mat_5groupfit_CSHQ_4group)
# psm_6groupfit_CSHQ_replica_4group <- comp.psm(cluster_mat_6groupfit_CSHQ_4group)
# 
# psm_4groupfit_CSHQ_replica_6group <- comp.psm(cluster_mat_4groupfit_CSHQ_6group)
# psm_5groupfit_CSHQ_replica_6group <- comp.psm(cluster_mat_5groupfit_CSHQ_6group)
# psm_6groupfit_CSHQ_replica_6group <- comp.psm(cluster_mat_6groupfit_CSHQ_6group)
# 
# 
# minVI_collapsed_CSHQ_replica_4group <- minVI(psm_collapsed_CSHQ_replica_4group, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_4group$Z))
# minVI_collapsed_CSHQ_replica_6group <- minVI(psm_collapsed_CSHQ_replica_6group, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_collapsed_CSHQ_replica_6group$Z))
# 
# 
# minVI_4groupfit_CSHQ_replica_4group <- minVI(psm_4groupfit_CSHQ_replica_4group, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_4groupfit_CSHQ_replica_4group$Z))
# minVI_5groupfit_CSHQ_replica_4group <- minVI(psm_5groupfit_CSHQ_replica_4group, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_5groupfit_CSHQ_replica_4group$Z))
# minVI_6groupfit_CSHQ_replica_4group <- minVI(psm_6groupfit_CSHQ_replica_4group, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_6groupfit_CSHQ_replica_4group$Z))
# 
# minVI_4groupfit_CSHQ_replica_6group <- minVI(psm_4groupfit_CSHQ_replica_6group, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_4groupfit_CSHQ_replica_6group$Z))
# minVI_5groupfit_CSHQ_replica_6group <- minVI(psm_5groupfit_CSHQ_replica_6group, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_5groupfit_CSHQ_replica_6group$Z))
# minVI_6groupfit_CSHQ_replica_6group <- minVI(psm_6groupfit_CSHQ_replica_6group, method = 'greedy', suppress.comment = FALSE, start.cl = max.col(LCA_6groupfit_CSHQ_replica_6group$Z))
# 
# 
# cred_ball_collapsed_CSHQ_replica_4group <- credibleball(c.star = minVI_collapsed_CSHQ_replica_4group$cl, cls.draw = cluster_mat_collapsed_CSHQ_replica_4group, c.dist = 'VI' )
# cred_ball_collapsed_CSHQ_replica_6group <- credibleball(c.star = minVI_collapsed_CSHQ_replica_6group$cl, cls.draw = cluster_mat_collapsed_CSHQ_replica_6group, c.dist = 'VI' )
# 
# 
# cred_ball_4groupfit_CSHQ_replica_4group <- credibleball(c.star = minVI_4groupfit_CSHQ_replica_4group$cl, cls.draw = cluster_mat_4groupfit_CSHQ_4group, c.dist = 'VI')
# cred_ball_5groupfit_CSHQ_replica_4group <- credibleball(c.star = minVI_5groupfit_CSHQ_replica_4group$cl, cls.draw = cluster_mat_5groupfit_CSHQ_4group, c.dist = 'VI')
# cred_ball_6groupfit_CSHQ_replica_4group <- credibleball(c.star = minVI_6groupfit_CSHQ_replica_4group$cl, cls.draw = cluster_mat_6groupfit_CSHQ_4group, c.dist = 'VI')
# 
# cred_ball_4groupfit_CSHQ_replica_6group <- credibleball(c.star = minVI_4groupfit_CSHQ_replica_6group$cl, cls.draw = cluster_mat_4groupfit_CSHQ_6group, c.dist = 'VI')
# cred_ball_5groupfit_CSHQ_replica_6group <- credibleball(c.star = minVI_5groupfit_CSHQ_replica_6group$cl, cls.draw = cluster_mat_5groupfit_CSHQ_6group, c.dist = 'VI')
# cred_ball_6groupfit_CSHQ_replica_6group <- credibleball(c.star = minVI_6groupfit_CSHQ_replica_6group$cl, cls.draw = cluster_mat_6groupfit_CSHQ_6group, c.dist = 'VI')












 