#Creating a simulation study which contains a number of very small clusters similar to the CSHQ dataset


#The structure is:
#Group 1: High probability for outcome 1 on item 1-7, uniform on other items
#Group 2: High probability for outcome 2 on items 8-14, uniform on other items
#Group 3: High probability for outcome 3 on items 15-21, uniform on other items

#Group 4: Alternating outcome dominant (say 0.8 for item 1 category 1, 0.8 for item 2 category 2,...)
#Group 5: reversed pattern from Group 4
#Group 6 strong probability for even outcomes, uniform for odd ones

#The groups 1-3 are simulated in the standard way
#The groups 4-6 are added to the dataset manually, so we can control the number of observations that we are adding


M <- 21
G <- 6
K <- c(rep(3, M))
noisyclass_pi <- c(0.45, 0.25, 0.21, 0.06, 0.02, 0.01)


noisyclass_theta1 <- matrix(c(0.8, 0.1, 0.1, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.8, 0.1, 0.1, 0.1, 0.1, 0.8, 0.33, 0.33, 0.33), byrow = TRUE, nrow = 6, ncol = K[1])
noisyclass_theta2 <- matrix(c(0.8, 0.1, 0.1, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1), nrow = 6, ncol = K[2])
noisyclass_theta3 <- matrix(c(0.8, 0.1, 0.1, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.33, 0.33, 0.33), nrow = 6, ncol = K[3])
noisyclass_theta4 <- matrix(c(0.8, 0.1, 0.1, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.8, 0.1, 0.1, 0.1, 0.1, 0.8, 0.1, 0.8, 0.1), nrow = 6, ncol = K[4])
noisyclass_theta5 <- matrix(c(0.8, 0.1, 0.1, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33), nrow = 6, ncol = K[5])
noisyclass_theta6 <- matrix(c(0.8, 0.1, 0.1, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.1, 0.8, 0.1), nrow = 6, ncol = K[6])
noisyclass_theta7 <- matrix(c(0.8, 0.1, 0.1, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.8, 0.1, 0.1, 0.1, 0.1, 0.8, 0.33, 0.33, 0.33), nrow = 6, ncol = K[7])
noisyclass_theta8 <- matrix(c(0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1), nrow = 6, ncol = K[8])
noisyclass_theta9 <- matrix(c(0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.33, 0.33, 0.33), nrow = 6, ncol = K[9])
noisyclass_theta10 <- matrix(c(0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33, 0.8, 0.1, 0.1, 0.1, 0.1, 0.8, 0.1, 0.8, 0.1), nrow = 6, ncol = K[10])
noisyclass_theta11 <- matrix(c(0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33), nrow = 6, ncol = K[11])
noisyclass_theta12 <- matrix(c(0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.1, 0.8, 0.1), nrow = 6, ncol = K[12])
noisyclass_theta13 <- matrix(c(0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33, 0.8, 0.1, 0.1, 0.1, 0.1, 0.8, 0.33, 0.33, 0.33), nrow = 6, ncol = K[13])
noisyclass_theta14 <- matrix(c(0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1), nrow = 6, ncol = K[14])
noisyclass_theta15 <- matrix(c(0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.33, 0.33, 0.33), nrow = 6, ncol = K[15])
noisyclass_theta16 <- matrix(c(0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.1, 0.1, 0.8, 0.1, 0.8, 0.1), nrow = 6, ncol = K[16])
noisyclass_theta17 <- matrix(c(0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1, 0.33, 0.33, 0.33), nrow = 6, ncol = K[17])
noisyclass_theta18 <- matrix(c(0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.1, 0.8, 0.1), nrow = 6, ncol = K[18])
noisyclass_theta19 <- matrix(c(0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.1, 0.1, 0.8, 0.33, 0.33, 0.33), nrow = 6, ncol = K[19])
noisyclass_theta20 <- matrix(c(0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1, 0.1, 0.8, 0.1), nrow = 6, ncol = K[20])
noisyclass_theta21 <- matrix(c(0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.1, 0.1, 0.8, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.33, 0.33, 0.33), nrow = 6, ncol = K[21])

noisyclass_theta <- list(noisyclass_theta1, noisyclass_theta2, noisyclass_theta3, noisyclass_theta4, noisyclass_theta5, noisyclass_theta6, noisyclass_theta7, noisyclass_theta8, noisyclass_theta9, noisyclass_theta10, noisyclass_theta11, noisyclass_theta12, noisyclass_theta13, noisyclass_theta14, noisyclass_theta15, noisyclass_theta16, noisyclass_theta17, noisyclass_theta18, noisyclass_theta19, noisyclass_theta20, noisyclass_theta21)




sim_data_noisy_classes_verysmall <- poLCA.simdata(N = 75, probs = noisyclass_theta, nclass = G, ndv = 21, P = noisyclass_pi)
sim_data_noisy_classes_small <- poLCA.simdata(N = 150, probs = noisyclass_theta, nclass = G, ndv = 21, P = noisyclass_pi)
sim_data_noisy_classes_medium <- poLCA.simdata(N = 300, probs = noisyclass_theta, nclass = G, ndv = 21, P = noisyclass_pi)
sim_data_noisy_classes_large <- poLCA.simdata(N = 1000, probs = noisyclass_theta, nclass = G, ndv = 21, P = noisyclass_pi)
sim_data_noisy_classes_verylarge <- poLCA.simdata(N = 1500, probs = noisyclass_theta, nclass = G, ndv = 21, P = noisyclass_pi)


noisyclass_Y_verysmall <- as.matrix(sim_data_noisy_classes_verysmall$dat)
noisyclass_Y_small <- as.matrix(sim_data_noisy_classes_small$dat)
noisyclass_Y_medium <- as.matrix(sim_data_noisy_classes_medium$dat)
noisyclass_Y_large <- as.matrix(sim_data_noisy_classes_large$dat)
noisyclass_Y_verylarge <- as.matrix(sim_data_noisy_classes_verylarge$dat)




noisyclass_collapsed_fit_verysmall <- blca.collapsed(X = noisyclass_Y_verysmall, G = 1, iter = 50000, burn.in = 2000, thin = 1/10, G.sel = TRUE, var.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
noisyclass_collapsed_fit_small <- blca.collapsed(X = noisyclass_Y_small, G = 1, iter = 50000, burn.in = 2000, thin = 1/10, G.sel = TRUE, var.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
noisyclass_collapsed_fit_medium <- blca.collapsed(X = noisyclass_Y_medium, G = 1, iter = 50000, burn.in = 2000, thin = 1/10, G.sel = TRUE, var.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
noisyclass_collapsed_fit_large <- blca.collapsed(X = noisyclass_Y_large, G = 1, iter = 50000, burn.in = 2000, thin = 1/10, G.sel = TRUE, var.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)
noisyclass_collapsed_fit_verylarge <- blca.collapsed(X = noisyclass_Y_verylarge, G = 1, iter = 50000, burn.in = 2000, thin = 1/10, G.sel = TRUE, var.sel = TRUE, verbose = TRUE, post.hoc.run = TRUE)

collapsed_group_proportions_verysmall <- table(noisyclass_collapsed_fit_verysmall$samples$G)/sum(table(noisyclass_collapsed_fit_verysmall$samples$G)) 
collapsed_group_proportions_small <- table(noisyclass_collapsed_fit_small$samples$G)/sum(table(noisyclass_collapsed_fit_small$samples$G))
collapsed_group_proportions_medium <- table(noisyclass_collapsed_fit_medium$samples$G)/sum(table(noisyclass_collapsed_fit_medium$samples$G))
collapsed_group_proportions_large <- table(noisyclass_collapsed_fit_large$samples$G)/sum(table(noisyclass_collapsed_fit_large$samples$G))
collapsed_group_proportions_verylarge <- table(noisyclass_collapsed_fit_verylarge$samples$G)/sum(table(noisyclass_collapsed_fit_verylarge$samples$G))


noisyclass_4group_LCA_fit_verysmall <- blca.gibbs(X = noisyclass_Y_verysmall, G = 4, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_4group_LCA_fit_small <- blca.gibbs(X = noisyclass_Y_small, G = 4, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_4group_LCA_fit_medium <- blca.gibbs(X = noisyclass_Y_medium, G = 4, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_4group_LCA_fit_large <- blca.gibbs(X = noisyclass_Y_large, G = 4, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_4group_LCA_fit_verylarge <- blca.gibbs(X = noisyclass_Y_verylarge, G = 4, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)

noisyclass_5group_LCA_fit_verysmall <- blca.gibbs(X = noisyclass_Y_verysmall, G = 5, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_5group_LCA_fit_small <- blca.gibbs(X = noisyclass_Y_small, G = 5, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_5group_LCA_fit_medium <- blca.gibbs(X = noisyclass_Y_medium, G = 5, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_5group_LCA_fit_large <- blca.gibbs(X = noisyclass_Y_large, G = 5, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_5group_LCA_fit_verylarge <- blca.gibbs(X = noisyclass_Y_verylarge, G = 5, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)

noisyclass_6group_LCA_fit_verysmall <- blca.gibbs(X = noisyclass_Y_verysmall, G = 6, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_6group_LCA_fit_small <- blca.gibbs(X = noisyclass_Y_small, G = 6, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_6group_LCA_fit_medium <- blca.gibbs(X = noisyclass_Y_medium, G = 6, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_6group_LCA_fit_large <- blca.gibbs(X = noisyclass_Y_large, G = 6, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)
noisyclass_6group_LCA_fit_verylarge <- blca.gibbs(X = noisyclass_Y_verylarge, G = 6, iter = 50000, burn.in = 2000, thin = 1/10, verbose = TRUE)


noisyclass_cluster_matrix_verysmall <- do.call(rbind, noisyclass_collapsed_fit_verysmall$samples$labels) 
noisyclass_cluster_matrix_small <- do.call(rbind, noisyclass_collapsed_fit_small$samples$labels) 
noisyclass_cluster_matrix_medium <- do.call(rbind, noisyclass_collapsed_fit_medium$samples$labels) 
noisyclass_cluster_matrix_large <- do.call(rbind, noisyclass_collapsed_fit_large$samples$labels) 
noisyclass_cluster_matrix_verylarge <- do.call(rbind, noisyclass_collapsed_fit_verylarge$samples$labels) 



psm_noisyclass_verysmall <- comp.psm(noisyclass_cluster_matrix_verysmall)
psm_noisyclass_small <- comp.psm(noisyclass_cluster_matrix_small)
psm_noisyclass_medium <- comp.psm(noisyclass_cluster_matrix_medium)
psm_noisyclass_large <- comp.psm(noisyclass_cluster_matrix_large)
psm_noisyclass_verylarge <- comp.psm(noisyclass_cluster_matrix_verylarge)

psm_noisyclass_4group_verysmall <- comp.psm(noisyclass_4group_LCA_fit_verysmall$samples$labels)
psm_noisyclass_4group_small <- comp.psm(noisyclass_4group_LCA_fit_small$samples$labels)
psm_noisyclass_4group_medium <- comp.psm(noisyclass_4group_LCA_fit_medium$samples$labels)
psm_noisyclass_4group_large <- comp.psm(noisyclass_4group_LCA_fit_large$samples$labels)
psm_noisyclass_4group_verylarge <- comp.psm(noisyclass_4group_LCA_fit_verylarge$samples$labels)

psm_noisyclass_5group_verysmall <- comp.psm(noisyclass_5group_LCA_fit_verysmall$samples$labels)
psm_noisyclass_5group_small <- comp.psm(noisyclass_5group_LCA_fit_small$samples$labels)
psm_noisyclass_5group_medium <- comp.psm(noisyclass_5group_LCA_fit_medium$samples$labels)
psm_noisyclass_5group_large <- comp.psm(noisyclass_5group_LCA_fit_large$samples$labels)
psm_noisyclass_5group_verylarge <- comp.psm(noisyclass_5group_LCA_fit_verylarge$samples$labels)

psm_noisyclass_6group_verysmall <- comp.psm(noisyclass_6group_LCA_fit_verysmall$samples$labels)
psm_noisyclass_6group_small <- comp.psm(noisyclass_6group_LCA_fit_small$samples$labels)
psm_noisyclass_6group_medium <- comp.psm(noisyclass_6group_LCA_fit_medium$samples$labels)
psm_noisyclass_6group_large <- comp.psm(noisyclass_6group_LCA_fit_large$samples$labels)
psm_noisyclass_6group_verylarge <- comp.psm(noisyclass_6group_LCA_fit_verylarge$samples$labels)

minVI_noisyclass_verysmall <- minVI(psm = psm_noisyclass_verysmall, method = 'greedy', start.cl = max.col(noisyclass_collapsed_fit_verysmall$Z), suppress.comment = FALSE)
minVI_noisyclass_small <- minVI(psm = psm_noisyclass_small, method = 'greedy', start.cl = max.col(noisyclass_collapsed_fit_small$Z), suppress.comment = FALSE)
minVI_noisyclass_medium <- minVI(psm = psm_noisyclass_medium, method = 'greedy', start.cl = max.col(noisyclass_collapsed_fit_medium$Z), suppress.comment = FALSE)
minVI_noisyclass_large <- minVI(psm = psm_noisyclass_large, method = 'greedy', start.cl = max.col(noisyclass_collapsed_fit_large$Z), suppress.comment = FALSE)
minVI_noisyclass_verylarge <- minVI(psm = psm_noisyclass_verylarge, method = 'greedy', start.cl = max.col(noisyclass_collapsed_fit_verylarge$Z), suppress.comment = FALSE)

minVI_noisyclass_4group_verysmall <- minVI(psm = psm_noisyclass_4group_verysmall, method = 'greedy', start.cl = max.col(noisyclass_4group_LCA_fit_verysmall$Z), suppress.comment = FALSE)
minVI_noisyclass_4group_small <- minVI(psm = psm_noisyclass_4group_small, method = 'greedy', start.cl = max.col(noisyclass_4group_LCA_fit_small$Z), suppress.comment = FALSE)
minVI_noisyclass_4group_medium <- minVI(psm = psm_noisyclass_4group_medium, method = 'greedy', start.cl = max.col(noisyclass_4group_LCA_fit_medium$Z), suppress.comment = FALSE)
minVI_noisyclass_4group_large <- minVI(psm = psm_noisyclass_4group_large, method = 'greedy', start.cl = max.col(noisyclass_4group_LCA_fit_large$Z), suppress.comment = FALSE)
minVI_noisyclass_4group_verylarge <- minVI(psm = psm_noisyclass_4group_verylarge, method = 'greedy', start.cl = max.col(noisyclass_4group_LCA_fit_verylarge$Z), suppress.comment = FALSE)

minVI_noisyclass_5group_verysmall <- minVI(psm = psm_noisyclass_5group_verysmall, method = 'greedy', start.cl = max.col(noisyclass_5group_LCA_fit_verysmall$Z), suppress.comment = FALSE)
minVI_noisyclass_5group_small <- minVI(psm = psm_noisyclass_5group_small, method = 'greedy', start.cl = max.col(noisyclass_5group_LCA_fit_small$Z), suppress.comment = FALSE)
minVI_noisyclass_5group_medium <- minVI(psm = psm_noisyclass_5group_medium, method = 'greedy', start.cl = max.col(noisyclass_5group_LCA_fit_medium$Z), suppress.comment = FALSE)
minVI_noisyclass_5group_large <- minVI(psm = psm_noisyclass_5group_large, method = 'greedy', start.cl = max.col(noisyclass_5group_LCA_fit_large$Z), suppress.comment = FALSE)
minVI_noisyclass_5group_verylarge <- minVI(psm = psm_noisyclass_5group_verylarge, method = 'greedy', start.cl = max.col(noisyclass_5group_LCA_fit_verylarge$Z), suppress.comment = FALSE)

minVI_noisyclass_6group_verysmall <- minVI(psm = psm_noisyclass_6group_verysmall, method = 'greedy', start.cl = max.col(noisyclass_6group_LCA_fit_verysmall$Z), suppress.comment = FALSE)
minVI_noisyclass_6group_small <- minVI(psm = psm_noisyclass_6group_small, method = 'greedy', start.cl = max.col(noisyclass_6group_LCA_fit_small$Z), suppress.comment = FALSE)
minVI_noisyclass_6group_medium <- minVI(psm = psm_noisyclass_6group_medium, method = 'greedy', start.cl = max.col(noisyclass_6group_LCA_fit_medium$Z), suppress.comment = FALSE)
minVI_noisyclass_6group_large <- minVI(psm = psm_noisyclass_6group_large, method = 'greedy', start.cl = max.col(noisyclass_6group_LCA_fit_large$Z), suppress.comment = FALSE)
minVI_noisyclass_6group_verylarge <- minVI(psm = psm_noisyclass_6group_verylarge, method = 'greedy', start.cl = max.col(noisyclass_6group_LCA_fit_verylarge$Z), suppress.comment = FALSE)

noisyclass_cred_ball_verysmall <- credibleball(c.star = minVI_noisyclass_verysmall$cl, cls.draw = noisyclass_cluster_matrix_verysmall , c.dist = 'VI', alpha = 0.05)
noisyclass_cred_ball_small <- credibleball(c.star = minVI_noisyclass_small$cl, cls.draw = noisyclass_cluster_matrix_small , c.dist = 'VI', alpha = 0.05)
noisyclass_cred_ball_medium <- credibleball(c.star = minVI_noisyclass_medium$cl, cls.draw = noisyclass_cluster_matrix_medium , c.dist = 'VI', alpha = 0.05)
noisyclass_cred_ball_large <- credibleball(c.star = minVI_noisyclass_large$cl, cls.draw = noisyclass_cluster_matrix_large , c.dist = 'VI', alpha = 0.05)
noisyclass_cred_ball_verylarge <- credibleball(c.star = minVI_noisyclass_verylarge$cl, cls.draw = noisyclass_cluster_matrix_verylarge , c.dist = 'VI', alpha = 0.05)






