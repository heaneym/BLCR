#script for reading in monte carlo results and getting summaries

source('./sim_study_monte_carlo_functions.R')
source('./sim_study_monte_carlo_summary_analysis.R')




sim1_results_N150 <- readRDS('./sim1_results_N150.rds')
sim1_results_N300 <- readRDS('./sim1_results_N300.rds')
sim1_results_N500 <- readRDS('./sim1_results_N500.rds')
  
sim2_results_N150 <- readRDS('./sim2_results_N150.rds')
sim2_results_N300 <- readRDS('./sim2_results_N300.rds')
sim2_results_N500 <- readRDS('./sim2_results_N500.rds')


#Simulation 1

sim1_beta_N150 <- beta_summarise(sim1_results_N150)
sim1_beta_N300 <- beta_summarise(sim1_results_N300)
sim1_beta_N500 <- beta_summarise(sim1_results_N500)

sim1_theta_N150 <- theta_summarise(sim1_results_N150, true_theta = sim1_theta)
sim1_theta_N300 <- theta_summarise(sim1_results_N300, true_theta = sim1_theta)
sim1_theta_N500 <- theta_summarise(sim1_results_N500, true_theta = sim1_theta)

sim1_pred_sel_N150 <- pred_var_sel_summarise(fit = sim1_results_N150, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_pred_sel_N300 <- pred_var_sel_summarise(fit = sim1_results_N300, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_pred_sel_N500 <- pred_var_sel_summarise(fit = sim1_results_N500, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)


sim1_item_sel_N150 <- item_var_sel_summarise(fit = sim1_results_N150, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_item_sel_N300 <- item_var_sel_summarise(fit = sim1_results_N300, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_item_sel_N500 <- item_var_sel_summarise(fit = sim1_results_N500, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)










#Simulation 2

sim2_beta_N150 <- beta_summarise(sim2_results_N150)
sim2_beta_N300 <- beta_summarise(sim2_results_N300)
sim2_beta_N500 <- beta_summarise(sim2_results_N500)

sim2_theta_N150 <- theta_summarise(fit = sim2_results_N150, true_theta = sim2_theta)
sim2_theta_N300 <- theta_summarise(fit = sim2_results_N300, true_theta = sim2_theta)
sim2_theta_N500 <- theta_summarise(fit = sim2_results_N500, true_theta = sim2_theta)

sim2_pred_sel_N150 <- pred_var_sel_summarise(fit = sim2_results_N150, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)
sim2_pred_sel_N300 <- pred_var_sel_summarise(fit = sim2_results_N300, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)
sim2_pred_sel_N500 <- pred_var_sel_summarise(fit = sim2_results_N500, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)


sim2_item_sel_N150 <- item_var_sel_summarise(fit = sim2_results_N150, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim2_item_sel_N300 <- item_var_sel_summarise(fit = sim2_results_N300, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim2_item_sel_N500 <- item_var_sel_summarise(fit = sim2_results_N500, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)