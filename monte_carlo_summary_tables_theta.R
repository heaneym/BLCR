source('./sim_study_monte_carlo_functions.R')
source('./sim_study_monte_carlo_summary_analysis.R')

library(xtable)

n_replicates <- 300
sample_sizes <- c(150, 300, 500)


sim1_results_N150 <- readRDS('./sim1_results_N150.rds')
sim1_results_N300 <- readRDS('./sim1_results_N300.rds')
sim1_results_N500 <- readRDS('./sim1_results_N500.rds')

sim2_results_N150 <- readRDS('./sim2_results_N150.rds')
sim2_results_N300 <- readRDS('./sim2_results_N300.rds')
sim2_results_N500 <- readRDS('./sim2_results_N500.rds')

sim3_results_N150 <- readRDS('./sim3_results_N150.rds')
sim3_results_N300 <- readRDS('./sim3_results_N300.rds')
sim3_results_N500 <- readRDS('./sim3_results_N500.rds')



sim1_theta_N150 <- theta_summarise(sim1_results_N150, true_theta = sim1_theta)
sim1_theta_N300 <- theta_summarise(sim1_results_N300, true_theta = sim1_theta)
sim1_theta_N500 <- theta_summarise(sim1_results_N500, true_theta = sim1_theta)

sim2_theta_N150 <- theta_summarise(sim2_results_N150, true_theta = sim2_theta)
sim2_theta_N300 <- theta_summarise(sim2_results_N300, true_theta = sim2_theta)
sim2_theta_N500 <- theta_summarise(sim2_results_N500, true_theta = sim2_theta)

sim3_theta_N150 <- theta_summarise(sim3_results_N150, true_theta = sim3_theta)
sim3_theta_N300 <- theta_summarise(sim3_results_N300, true_theta = sim3_theta)
sim3_theta_N500 <- theta_summarise(sim3_results_N500, true_theta = sim3_theta)

# 
# 
# 
# #Creating a summary table for simulation 1 first
# sim1_theta_bias_g1_N150 <- lapply(sim1_theta_N150$mean_bias, function(x) x[1,])
# sim1_theta_bias_g2_N150 <- lapply(sim1_theta_N150$mean_bias, function(x) x[2,])
# 
# sim1_theta_bias_g1_N300 <- lapply(sim1_theta_N300$mean_bias, function(x) x[1,])
# sim1_theta_bias_g2_N300 <- lapply(sim1_theta_N300$mean_bias, function(x) x[2,])
# 
# sim1_theta_bias_g1_N500 <- lapply(sim1_theta_N500$mean_bias, function(x) x[1,])
# sim1_theta_bias_g2_N500 <- lapply(sim1_theta_N500$mean_bias, function(x) x[2,])
# 
# sim1_theta_mse_g1_N150 <- lapply(sim1_theta_N150$mean_mse, function(x) x[1,])
# sim1_theta_mse_g2_N150 <- lapply(sim1_theta_N150$mean_mse, function(x) x[2,])
# 
# sim1_theta_mse_g1_N300 <- lapply(sim1_theta_N300$mean_mse, function(x) x[1,])
# sim1_theta_mse_g2_N300 <- lapply(sim1_theta_N300$mean_mse, function(x) x[2,])
# 
# sim1_theta_mse_g1_N500 <- lapply(sim1_theta_N500$mean_mse, function(x) x[1,])
# sim1_theta_mse_g2_N500 <- lapply(sim1_theta_N500$mean_mse, function(x) x[2,])
# 
# 
# sim1_theta_coverage_g1_N150 <- lapply(sim1_theta_N150$mean_coverage, function(x) x[1,])
# sim1_theta_coverage_g2_N150 <- lapply(sim1_theta_N150$mean_coverage, function(x) x[2,])
# 
# sim1_theta_coverage_g1_N300 <- lapply(sim1_theta_N300$mean_coverage, function(x) x[1,])
# sim1_theta_coverage_g2_N300 <- lapply(sim1_theta_N300$mean_coverage, function(x) x[2,])
# 
# sim1_theta_coverage_g1_N500 <- lapply(sim1_theta_N500$mean_coverage, function(x) x[1,])
# sim1_theta_coverage_g2_N500 <- lapply(sim1_theta_N500$mean_coverage, function(x) x[2,])
# 
# 
# 
# sim1_theta_mean_bias_g1 <- c(mean(unlist(sim1_theta_bias_g1_N150)), mean(unlist(sim1_theta_bias_g1_N300)), mean(unlist(sim1_theta_bias_g1_N500)))
# sim1_theta_mean_bias_g2 <- c(mean(unlist(sim1_theta_bias_g2_N150)), mean(unlist(sim1_theta_bias_g2_N300)), mean(unlist(sim1_theta_bias_g2_N500)))
# 
# 
# sim1_theta_mean_abs_bias_g1 <- c(mean(abs(unlist(sim1_theta_bias_g1_N150))), mean(abs(unlist(sim1_theta_bias_g1_N300))), mean(abs(unlist(sim1_theta_bias_g1_N500))))
# sim1_theta_mean_abs_bias_g2 <- c(mean(abs(unlist(sim1_theta_bias_g2_N150))), mean(abs(unlist(sim1_theta_bias_g2_N300))), mean(abs(unlist(sim1_theta_bias_g2_N500))))
# 
# sim1_theta_max_abs_bias_g1 <- c(max(abs(unlist(sim1_theta_bias_g1_N150))), max(abs(unlist(sim1_theta_bias_g1_N300))), max(abs(unlist(sim1_theta_bias_g1_N500))))
# sim1_theta_max_abs_bias_g2 <- c(max(abs(unlist(sim1_theta_bias_g2_N150))), max(abs(unlist(sim1_theta_bias_g2_N300))), max(abs(unlist(sim1_theta_bias_g2_N500))))
# 
# sim1_theta_mean_mse_g1 <- c(mean(unlist(sim1_theta_mse_g1_N150)), mean(unlist(sim1_theta_mse_g1_N300)), mean(unlist(sim1_theta_mse_g1_N500)))
# sim1_theta_mean_mse_g2 <- c(mean(unlist(sim1_theta_mse_g2_N150)), mean(unlist(sim1_theta_mse_g2_N300)), mean(unlist(sim1_theta_mse_g2_N500)))
# 
# sim1_theta_mean_coverage_g1 <- c(mean(unlist(sim1_theta_coverage_g1_N150)), mean(unlist(sim1_theta_coverage_g1_N300)), mean(unlist(sim1_theta_coverage_g1_N500)))
# sim1_theta_mean_coverage_g2 <- c(mean(unlist(sim1_theta_coverage_g2_N150)), mean(unlist(sim1_theta_coverage_g2_N300)), mean(unlist(sim1_theta_coverage_g2_N500)))
# 
# sim1_theta_min_coverage_g1 <- c(min(unlist(sim1_theta_coverage_g1_N150)), min(unlist(sim1_theta_coverage_g1_N300)), min(unlist(sim1_theta_coverage_g1_N500)))
# sim1_theta_min_coverage_g2 <- c(min(unlist(sim1_theta_coverage_g2_N150)), min(unlist(sim1_theta_coverage_g2_N300)), min(unlist(sim1_theta_coverage_g2_N500)))
# 
# sim1_theta_max_coverage_g1 <- c(max(unlist(sim1_theta_coverage_g1_N150)), max(unlist(sim1_theta_coverage_g1_N300)), max(unlist(sim1_theta_coverage_g1_N500)))
# sim1_theta_max_coverage_g2 <- c(max(unlist(sim1_theta_coverage_g2_N150)), max(unlist(sim1_theta_coverage_g2_N300)), max(unlist(sim1_theta_coverage_g2_N500)))
# 


library(xtable)

make_sim_table <- function(theta_list, G, Ns = c(150, 300, 500)) {
  
  # One block of rows (one row per N) for each group
  blocks <- lapply(seq_len(G), function(g) {
    t(sapply(seq_along(theta_list), function(i) {
      res  <- theta_list[[i]]
      bias <- unlist(lapply(res$mean_bias,     function(x) x[g, ]))
      mse  <- unlist(lapply(res$mean_mse,      function(x) x[g, ]))
      cov  <- unlist(lapply(res$mean_coverage, function(x) x[g, ]))
      
      c(N             = Ns[i],
        mean_bias     = mean(bias),
        mean_abs_bias = mean(abs(bias)),
        max_abs_bias  = max(abs(bias)),
        mean_mse      = mean(mse),
        mean_coverage = mean(cov),
        min_coverage  = min(cov),
        max_coverage  = max(cov))
    }))
  })
  
  tab <- as.data.frame(do.call(rbind, blocks))
  n_N <- length(Ns)
  
  
  g_col <- unlist(lapply(seq_len(G), function(g) {
    c(paste0("$g = ", g, "$"), rep("", n_N - 1))
  }))
  
  tab <- data.frame(g = g_col, tab, check.names = FALSE, stringsAsFactors = FALSE)
  
  colnames(tab) <- c(
    "",
    "$N$",
    "\\shortstack{Mean\\\\Bias}",
    "\\shortstack{Mean\\\\$|\\mathrm{Bias}|$}",
    "\\shortstack{Max\\\\$|\\mathrm{Bias}|$}",
    "\\shortstack{Mean\\\\MSE}",
    "\\shortstack{Mean\\\\Coverage}",
    "\\shortstack{Min\\\\Coverage}",
    "\\shortstack{Max\\\\Coverage}"
  )
  
  tab_xt <- xtable(
    tab,
    digits = c(0, 0, 0, rep(3, 7)),
    align  = c("l", "l", rep("c", 8))
  )
  
  # A midrule after every group except the last
  mid_pos <- as.list(n_N * seq_len(G - 1))
  
  print(
    tab_xt,
    include.rownames = FALSE,
    sanitize.colnames.function = identity,
    sanitize.text.function = identity,
    booktabs = TRUE,
    size = "\\small",
    add.to.row = list(pos = mid_pos, command = rep("\\midrule\n", G - 1)),
    hline.after = c(-1, 0, nrow(tab))
  )
}

# Simulation 1 (2 groups)
make_sim_table(list(sim1_theta_N150, sim1_theta_N300, sim1_theta_N500), G = 2)

# Simulation 2 (3 groups)
make_sim_table(list(sim2_theta_N150, sim2_theta_N300, sim2_theta_N500), G = 3)

# Simulation 3 (4 groups)
make_sim_table(list(sim3_theta_N150, sim3_theta_N300, sim3_theta_N500), G = 4)
