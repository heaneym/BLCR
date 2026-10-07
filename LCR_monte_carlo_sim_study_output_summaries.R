#script for reading in monte carlo results and getting summaries
# We'll keep this for only the beta stuff for the moment and make
# separate files for theta and variable selection stuff


source('./sim_study_monte_carlo_functions.R')
source('./sim_study_monte_carlo_summary_analysis.R')

source('./sim_study_parameters.R')


get_row_sim1_beta <- function(metric, k) {
  c(sim1_beta_N150$conditional[[metric]][k, ],
    sim1_beta_N300$conditional[[metric]][k, ],
    sim1_beta_N500$conditional[[metric]][k, ])
}


make_block_sim1_beta <- function(title, metric) {
  header <- data.frame(
    sample_sizes = title,
    intercept = NA, X_1 = NA, X_2 = NA, X_3 = NA, X_4 = NA,
    stringsAsFactors = FALSE
  )
  body <- data.frame(
    sample_sizes = c('$N = 150$', '$N = 300$', '$N = 500$'),
    intercept = get_row_sim1_beta(metric, 1),
    X_1 = get_row_sim1_beta(metric, 2),
    X_2 = get_row_sim1_beta(metric, 3),
    X_3 = get_row_sim1_beta(metric, 4),
    X_4 = get_row_sim1_beta(metric, 5),
    stringsAsFactors = FALSE
  )
  rbind(header, body)
}

get_row_sim2_beta <- function(metric, k, j) {
  c(sim2_beta_N150$conditional[[metric]][k, j],
    sim2_beta_N300$conditional[[metric]][k, j],
    sim2_beta_N500$conditional[[metric]][k, j])
}

get_row_sim3_beta <- function(metric, k, j) {
  c(sim3_beta_N150$conditional[[metric]][k, j],
    sim3_beta_N300$conditional[[metric]][k, j],
    sim3_beta_N500$conditional[[metric]][k, j])
}

make_block_sim2_beta <- function(title, metric, j) {
  header <- data.frame(
    sample_sizes = title,
    intercept = NA, X_1 = NA, X_2 = NA, X_3 = NA, X_4 = NA,
    stringsAsFactors = FALSE
  )
  body <- data.frame(
    sample_sizes = c('$N = 150$', '$N = 300$', '$N = 500$'),
    intercept = get_row_sim2_beta(metric, 1, j),
    X_1 = get_row_sim2_beta(metric, 2, j),
    X_2 = get_row_sim2_beta(metric, 3, j),
    X_3 = get_row_sim2_beta(metric, 4, j),
    X_4 = get_row_sim2_beta(metric, 5, j),
    stringsAsFactors = FALSE
  )
  rbind(header, body)
}

make_block_sim3_beta <- function(title, metric, j) {
  header <- data.frame(
    sample_sizes = title,
    intercept = NA, X_1 = NA, X_2 = NA, X_3 = NA, X_4 = NA,
    stringsAsFactors = FALSE
  )
  body <- data.frame(
    sample_sizes = c('$N = 150$', '$N = 300$', '$N = 500$'),
    intercept = get_row_sim3_beta(metric, 1, j),
    X_1 = get_row_sim3_beta(metric, 2, j),
    X_2 = get_row_sim3_beta(metric, 3, j),
    X_3 = get_row_sim3_beta(metric, 4, j),
    X_4 = get_row_sim3_beta(metric, 5, j),
    stringsAsFactors = FALSE
  )
  rbind(header, body)
}

make_panel_sim2_beta <- function(j) {
  section <- data.frame(
    sample_sizes = sprintf('$g = %d$', j),
    intercept = NA, X_1 = NA, X_2 = NA, X_3 = NA, X_4 = NA,
    stringsAsFactors = FALSE
  )
  true_row <- data.frame(
    sample_sizes = 'True Coefficient',
    intercept = sim2_beta[1, j],
    X_1 = sim2_beta[2, j],
    X_2 = sim2_beta[3, j],
    X_3 = sim2_beta[4, j],
    X_4 = sim2_beta[5, j],
    stringsAsFactors = FALSE
  )
  rbind(
    section,
    make_block_sim2_beta('Coverage', 'mean_coverage', j),
    make_block_sim2_beta('Bias', 'mean_bias', j),
    make_block_sim2_beta('MSE', 'mean_mse', j),
    true_row
  )
}

make_panel_sim3_beta <- function(j) {
  section <- data.frame(
    sample_sizes = sprintf('$g = %d$', j),
    intercept = NA, X_1 = NA, X_2 = NA, X_3 = NA, X_4 = NA,
    stringsAsFactors = FALSE
  )
  true_row <- data.frame(
    sample_sizes = 'True Coefficient',
    intercept = sim3_beta[1, j],
    X_1 = sim3_beta[2, j],
    X_2 = sim3_beta[3, j],
    X_3 = sim3_beta[4, j],
    X_4 = sim3_beta[5, j],
    stringsAsFactors = FALSE
  )
  rbind(
    section,
    make_block_sim3_beta('Coverage', 'mean_coverage', j),
    make_block_sim3_beta('Bias', 'mean_bias', j),
    make_block_sim3_beta('MSE', 'mean_mse', j),
    true_row
  )
}

get_n <- function(fit) {
  rowSums(!is.na(fit$conditional$coverage[1:5, , , drop = FALSE]), dims = 1)
}

# per-replicate values for coefficient k from one fit
get_reps <- function(fit, metric, k) as.vector(fit$conditional[[metric]][k, , ])







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
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_item_sel_N300 <- item_var_sel_summarise(fit = sim1_results_N300, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_item_sel_N500 <- item_var_sel_summarise(fit = sim1_results_N500, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)


# Want to output our desired quantities using xtable so we can put the tables directly into latex

# We want coverage, bias, MSE for informative predictors and mean posterior inclusion for non-informative ones  

# Firstly for coverage we can include a table across the informative variables
# Then create some kind of plot with the binomial confidence intervals?




true_row <- data.frame(
  sample_sizes = 'True Coefficient',
  intercept = sim1_beta[1, ],
  X_1 = sim1_beta[2, ],
  X_2 = sim1_beta[3, ],
  X_3 = sim1_beta[4, ],
  X_4 = sim1_beta[5, ],
  stringsAsFactors = FALSE
)

sim1_beta_summary_frame <- rbind(
  make_block_sim1_beta('Coverage', 'mean_coverage'),
  make_block_sim1_beta('Bias', 'mean_bias'),
  make_block_sim1_beta('MSE', 'mean_mse'),
  true_row
)

colnames(sim1_beta_summary_frame) <- c(
  '',
  '\\multicolumn{1}{c}{\\phantom{$-$}$\\beta_0$}',
  sprintf('\\multicolumn{1}{c}{$\\beta_%d$}', 1:4)
)

print(
  xtable(
    sim1_beta_summary_frame,
    digits = c(0, 0, 3, 3, 3, 3, 3)
  ),
  sanitize.text.function = identity,
  sanitize.colnames.function = identity,
  math.style.negative = TRUE,
  include.rownames = FALSE,
  booktabs = TRUE,
  na.print = '',
  hline.after = c(-1, 0, 4, 8, 12, nrow(sim1_beta_summary_frame))
)




# Want to plot the coverage along with binomial confidence intervals

n_mat <- cbind(
  get_n(sim1_beta_N150),
  get_n(sim1_beta_N300),
  get_n(sim1_beta_N500)
)

n_replicates <- 300  

sample_sizes <- c(150, 300, 500)
coef_names <- c("beta_0", "beta_1", "beta_2", "beta_3", "beta_4")


cov_mat <- cbind(
  sim1_beta_N150$conditional$mean_coverage[1:5, ],
  sim1_beta_N300$conditional$mean_coverage[1:5, ],
  sim1_beta_N500$conditional$mean_coverage[1:5, ]
)

# Wald 95% interval for each coverage estimate
se_mat <- sqrt(cov_mat * (1 - cov_mat) / n_mat)
lower_mat <- pmax(cov_mat - 1.96 * se_mat, 0)
upper_mat <- pmin(cov_mat + 1.96 * se_mat, 1)

# One panel per coefficient
pdf("./sim_study_plots/sim_study1_plots/sim1_beta_coverage_plot.pdf", width = 7, height = 2.4)
par(mfrow = c(1, 5),
    oma = c(2, 3, 2, 0),     
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 1)

for (k in 1:5) {
  plot(1:3, cov_mat[k, ],
       ylim = c(min(lower_mat, na.rm = TRUE), 1),
       xaxt = "n",
       yaxt = "n",
       xlim = c(0.5, 3.5),
       pch = 16,
       xlab = "",
       ylab = "",
       main = bquote(beta[.(k - 1)]))
  axis(1, at = 1:3, labels = sample_sizes)
  if (k == 1) axis(2)
  arrows(1:3, lower_mat[k, ], 1:3, upper_mat[k, ],
         angle = 90, code = 3, length = 0.03)
  abline(h = 0.95, lty = 2, col = "red")
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("Coverage", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Empirical coverage for regression coefficients with 95% Wald intervals",
      side = 3, outer = TRUE, line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))




#FOR THE BIAS WE MAKE A BOX PLOT
fits <- list(sim1_beta_N150, sim1_beta_N300, sim1_beta_N500)
sample_sizes <- c(150, 300, 500)


metric <- "bias"   


y_lim <- range(unlist(lapply(1:5, function(k)
  lapply(fits, get_reps, metric = metric, k = k))), na.rm = TRUE)

pdf("./sim_study_plots/sim_study1_plots/sim1_beta_bias_plot.pdf", width = 7, height = 4)
par(mfrow = c(1, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 1)

for (k in 1:5) {
  dat <- lapply(fits, get_reps, metric = metric, k = k)
  names(dat) <- sample_sizes
  boxplot(dat,
          ylim = y_lim,
          yaxt = "n",
          pch = 16, cex = 0.4,          
          main = bquote(beta[.(k - 1)]))
  if (k == 1) axis(2)
  abline(h = 0, lty = 2, col = "red")
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("Bias", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Bias of coefficient estimates", side = 3, outer = TRUE,
      line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))




#Boxplots with outliers removed

metric <- "bias"   


y_lim <- c(-1.5,2)

pdf("./sim_study_plots/sim_study1_plots/sim1_beta_bias_plot_no_outliers.pdf", width = 7, height = 2.4)
par(mfrow = c(1, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 1)

for (k in 1:5) {
  dat <- lapply(fits, get_reps, metric = metric, k = k)
  names(dat) <- sample_sizes
  boxplot(dat,
          ylim = y_lim,
          yaxt = "n",
          pch = 16, cex = 0.4,          
          main = bquote(beta[.(k - 1)]),
          outline = FALSE)
  if (k == 1) axis(2)
  abline(h = 0, lty = 2, col = "red")
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("Bias", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Bias of coefficient estimates", side = 3, outer = TRUE,
      line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))





#NOW DOING MSE

mse_mat <- cbind(
  sim1_beta_N150$conditional$mean_mse[1:5, ],
  sim1_beta_N300$conditional$mean_mse[1:5, ],
  sim1_beta_N500$conditional$mean_mse[1:5, ]
)

y_lim <- range(mse_mat, na.rm = TRUE)

pdf("./sim_study_plots/sim_study1_plots/sim1_beta_mse_plot.pdf", width = 7, height = 2.7)
par(mfrow = c(1, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 1)

for (k in 1:5) {
  plot(1:3, mse_mat[k, ],
       #log = "y",
       ylim = y_lim,
       xlim = c(0.5, 3.5),
       xaxt = "n", yaxt = "n",
       type = "b", pch = 16,
       xlab = "", ylab = "",
       main = bquote(beta[.(k - 1)]))
  axis(1, at = 1:3, labels = sample_sizes)
  if (k == 1) axis(2)
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("MSE", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Mean squared error of coefficient estimates", side = 3, outer = TRUE,
      line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))















# 
# sim1_beta_coverage_frame_uncond <- data.frame(
#   sample_sizes = c('$N = 150$', '$N = 300$', '$N = 500$', 'True Coefficient'),
#   intercept = c(sim1_beta_N150$unconditional$mean_coverage[1,],
#                 sim1_beta_N300$unconditional$mean_coverage[1,],
#                 sim1_beta_N500$unconditional$mean_coverage[1,],
#                 sim1_beta[1,]),
#   X_1 = c(sim1_beta_N150$unconditional$mean_coverage[2,],
#           sim1_beta_N300$unconditional$mean_coverage[2,],
#           sim1_beta_N500$unconditional$mean_coverage[2,],
#           sim1_beta[2,]),
#   X_2 = c(sim1_beta_N150$unconditional$mean_coverage[3,],
#           sim1_beta_N300$unconditional$mean_coverage[3,],
#           sim1_beta_N500$unconditional$mean_coverage[3,],
#           sim1_beta[3,]),
#   X_3 = c(sim1_beta_N150$unconditional$mean_coverage[4,],
#           sim1_beta_N300$unconditional$mean_coverage[4,],
#           sim1_beta_N500$unconditional$mean_coverage[4,],
#           sim1_beta[4,]),
#   X_4 = c(sim1_beta_N150$unconditional$mean_coverage[5,],
#           sim1_beta_N300$unconditional$mean_coverage[5,],
#           sim1_beta_N500$unconditional$mean_coverage[5,],
#           sim1_beta[5,])
#)












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









sim2_beta_summary_frame <- rbind(
  make_panel_sim2_beta(1),
  make_panel_sim2_beta(2)
)

colnames(sim2_beta_summary_frame) <- c(
  '',
  '\\multicolumn{1}{c}{\\phantom{$-$}$\\beta_{0 g}$}',
  sprintf('\\multicolumn{1}{c}{$\\beta_{%d g}$}', 1:4)
)

panel_rows <- nrow(sim2_beta_summary_frame) / 2   

print(
  xtable(
    sim2_beta_summary_frame,
    digits = c(0, 0, 3, 3, 3, 3, 3)
  ),
  sanitize.text.function = identity,
  sanitize.colnames.function = identity,
  math.style.negative = TRUE,
  include.rownames = FALSE,
  booktabs = TRUE,
  na.print = '',
  hline.after = c(-1, 0, panel_rows, nrow(sim2_beta_summary_frame))
)




#Creating plots for the coverage


get_n_sim2 <- function(fit) {
  apply(!is.na(fit$conditional$coverage[1:5, , , drop = FALSE]), c(1, 2), sum)
}

n_mat <- cbind(
  get_n_sim2(sim2_beta_N150),
  get_n_sim2(sim2_beta_N300),
  get_n_sim2(sim2_beta_N500)
)

n_replicates <- 300

sample_sizes <- c(150, 300, 500)
coef_names <- c("beta_0", "beta_1", "beta_2", "beta_3", "beta_4")


cov_mat <- cbind(
  sim2_beta_N150$conditional$mean_coverage[1:5, ],
  sim2_beta_N300$conditional$mean_coverage[1:5, ],
  sim2_beta_N500$conditional$mean_coverage[1:5, ]
)

# Wald 95% interval for each coverage estimate
se_mat <- sqrt(cov_mat * (1 - cov_mat) / n_mat)
lower_mat <- pmax(cov_mat - 1.96 * se_mat, 0)
upper_mat <- pmin(cov_mat + 1.96 * se_mat, 1)

# One panel per coefficient (columns), one row per g
pdf("./sim_study_plots/sim_study2_plots/sim2_beta_coverage_plot.pdf", width = 7, height = 4.4)
par(mfrow = c(2, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 0.9)

for (g in 1:2) {
  cols <- seq(g, by = 2, length.out = 3)   
  for (k in 1:5) {
    plot(1:3, cov_mat[k, cols],
         ylim = c(min(lower_mat, na.rm = TRUE), 1),
         xaxt = "n",
         yaxt = "n",
         xlim = c(0.5, 3.5),
         pch = 16,
         xlab = "",
         ylab = "",
         main = bquote(beta[.(k-1) * .(g)]))
    axis(1, at = 1:3, labels = sample_sizes)
    if (k == 1) axis(2)
    arrows(1:3, lower_mat[k, cols], 1:3, upper_mat[k, cols],
           angle = 90, code = 3, length = 0.03)
    abline(h = 0.95, lty = 2, col = "red")
  }
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("Coverage", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Empirical coverage for regression coefficients with 95% Wald intervals",
      side = 3, outer = TRUE, line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))




get_reps_sim2 <- function(fit, metric, k, g) as.vector(fit$conditional[[metric]][k, g, ])


#FOR THE BIAS WE MAKE A BOX PLOT
fits <- list(sim2_beta_N150, sim2_beta_N300, sim2_beta_N500)
sample_sizes <- c(150, 300, 500)


metric <- "bias"


# shared y-axis range across both g
y_lim <- range(unlist(lapply(1:2, function(g)
  lapply(1:5, function(k)
    lapply(fits, get_reps_sim2, metric = metric, k = k, g = g)))), na.rm = TRUE)

pdf("./sim_study_plots/sim_study2_plots/sim2_beta_bias_plot.pdf", width = 7, height = 6)
par(mfrow = c(2, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 0.9)

for (g in 1:2) {
  for (k in 1:5) {
    dat <- lapply(fits, get_reps_sim2, metric = metric, k = k, g = g)
    names(dat) <- sample_sizes
    boxplot(dat,
            ylim = y_lim,
            yaxt = "n",
            pch = 16, cex = 0.4,
            main = bquote(beta[.(k-1) * .(g)]))
    if (k == 1) axis(2)
    abline(h = 0, lty = 2, col = "red")
  }
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("Bias", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Bias of coefficient estimates", side = 3, outer = TRUE,
      line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))




#Boxplots with outliers removed

metric <- "bias"


y_lim <- c(-4, 4.5)   

pdf("./sim_study_plots/sim_study2_plots/sim2_beta_bias_plot_no_outliers.pdf", width = 7, height = 4.4)
par(mfrow = c(2, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 0.9)

for (g in 1:2) {
  for (k in 1:5) {
    dat <- lapply(fits, get_reps_sim2, metric = metric, k = k, g = g)
    names(dat) <- sample_sizes
    boxplot(dat,
            ylim = y_lim,
            yaxt = "n",
            pch = 16, cex = 0.4,
            main = bquote(beta[.(k-1) * .(g)]),
            outline = FALSE)
    if (k == 1) axis(2)
    abline(h = 0, lty = 2, col = "red")
  }
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("Bias", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Bias of coefficient estimates", side = 3, outer = TRUE,
      line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))





#NOW DOING MSE

# columns ordered N150-g1, N150-g2, N300-g1, N300-g2, N500-g1, N500-g2
mse_mat <- cbind(
  sim2_beta_N150$conditional$mean_mse[1:5, ],
  sim2_beta_N300$conditional$mean_mse[1:5, ],
  sim2_beta_N500$conditional$mean_mse[1:5, ]
)

y_lim <- range(mse_mat, na.rm = TRUE)

pdf("./sim_study_plots/sim_study2_plots/sim2_beta_mse_plot.pdf", width = 7, height = 4.8)
par(mfrow = c(2, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 0.9)

for (g in 1:2) {
  cols <- seq(g, by = 2, length.out = 3)   # the three sample sizes for this g
  for (k in 1:5) {
    plot(1:3, mse_mat[k, cols],
         #log = "y",
         ylim = y_lim,
         xlim = c(0.5, 3.5),
         xaxt = "n", yaxt = "n",
         type = "b", pch = 16,
         xlab = "", ylab = "",
         main = bquote(beta[.(k-1) * .(g)]))
    axis(1, at = 1:3, labels = sample_sizes)
    if (k == 1) axis(2)
  }
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("MSE", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Mean squared error of coefficient estimates", side = 3, outer = TRUE,
      line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))




















#Looking now at simulation 3 (4 group scenario similar to CSHQ data )

sim3_results_N150 <- readRDS('./sim3_results_N150.rds')
sim3_results_N300 <- readRDS('./sim3_results_N300.rds')
sim3_results_N500 <- readRDS('./sim3_results_N500.rds')




#Simulation 3

sim3_beta_N150 <- beta_summarise(sim3_results_N150)
sim3_beta_N300 <- beta_summarise(sim3_results_N300)
sim3_beta_N500 <- beta_summarise(sim3_results_N500)

sim3_theta_N150 <- theta_summarise(sim3_results_N150, true_theta = sim3_theta)
sim3_theta_N300 <- theta_summarise(sim3_results_N300, true_theta = sim3_theta)
sim3_theta_N500 <- theta_summarise(sim3_results_N500, true_theta = sim3_theta)

sim3_pred_sel_N150 <- pred_var_sel_summarise(fit = sim3_results_N150, 
                                             true_pred_active = rep(TRUE, 5), 
                                             threshold = 0.5)
sim3_pred_sel_N300 <- pred_var_sel_summarise(fit = sim3_results_N300, 
                                             true_pred_active = rep(TRUE, 5), 
                                             threshold = 0.5)
sim3_pred_sel_N500 <- pred_var_sel_summarise(fit = sim3_results_N500, 
                                             true_pred_active = rep(TRUE, 5), 
                                             threshold = 0.5)


sim3_item_sel_N150 <- item_var_sel_summarise(fit = sim3_results_N150, 
                                             true_item_active = c(rep(TRUE, 20), rep(FALSE, 20)), 
                                             threshold = 0.5)
sim3_item_sel_N300 <- item_var_sel_summarise(fit = sim3_results_N300, 
                                             true_item_active = c(rep(TRUE, 20), rep(FALSE, 20)), 
                                             threshold = 0.5)
sim3_item_sel_N500 <- item_var_sel_summarise(fit = sim3_results_N500, 
                                             true_item_active = c(rep(TRUE, 20), rep(FALSE, 20)), 
                                             threshold = 0.5)






sim3_beta_summary_frame <- rbind(
  make_panel_sim3_beta(1),
  make_panel_sim3_beta(2),
  make_panel_sim3_beta(3)
)

colnames(sim3_beta_summary_frame) <- c(
  '',
  '\\multicolumn{1}{c}{\\phantom{$-$}$\\beta_{0g}$}',
  sprintf('\\multicolumn{1}{c}{$\\beta_{%d g}$}', 1:4)
)

panel_rows1 <- nrow(sim3_beta_summary_frame) / 3   
panel_rows2 <- 2*nrow(sim3_beta_summary_frame) / 3  

print(
  xtable(
    sim3_beta_summary_frame,
    digits = c(0, 0, 3, 3, 3, 3, 3)
  ),
  sanitize.text.function = identity,
  sanitize.colnames.function = identity,
  math.style.negative = TRUE,
  include.rownames = FALSE,
  booktabs = TRUE,
  na.print = '',
  hline.after = c(-1, 0, panel_rows1, panel_rows2, nrow(sim3_beta_summary_frame))
)




#Creating plots for the coverage


get_n_sim3 <- function(fit) {
  apply(!is.na(fit$conditional$coverage[1:5, , , drop = FALSE]), c(1, 2), sum)
}

n_mat <- cbind(
  get_n_sim3(sim3_beta_N150),
  get_n_sim3(sim3_beta_N300),
  get_n_sim3(sim3_beta_N500)
)

n_replicates <- 300

sample_sizes <- c(150, 300, 500)
coef_names <- c("beta_0", "beta_1", "beta_2", "beta_3", "beta_4")


cov_mat <- cbind(
  sim3_beta_N150$conditional$mean_coverage[1:5, ],
  sim3_beta_N300$conditional$mean_coverage[1:5, ],
  sim3_beta_N500$conditional$mean_coverage[1:5, ]
)

# Wald 95% interval for each coverage estimate
se_mat <- sqrt(cov_mat * (1 - cov_mat) / n_mat)
lower_mat <- pmax(cov_mat - 1.96 * se_mat, 0)
upper_mat <- pmin(cov_mat + 1.96 * se_mat, 1)

# One panel per coefficient (columns), one row per g
pdf("./sim_study_plots/sim_study3_plots/sim3_beta_coverage_plot.pdf", width = 7, height = 6.5)
par(mfrow = c(3, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 0.9)

for (g in 1:3) {
  cols <- seq(g, by = 3, length.out = 3)
  for (k in 1:5) {
    plot(1:3, cov_mat[k, cols],
         ylim = c(min(lower_mat, na.rm = TRUE), 1),
         xaxt = "n",
         yaxt = "n",
         xlim = c(0.5, 3.5),
         pch = 16,
         xlab = "",
         ylab = "",
         main = bquote(beta[.(k-1) * .(g)]))
    axis(1, at = 1:3, labels = sample_sizes)
    if (k == 1) axis(2)
    arrows(1:3, lower_mat[k, cols], 1:3, upper_mat[k, cols],
           angle = 90, code = 3, length = 0.03)
    abline(h = 0.95, lty = 2, col = "red")
  }
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("Coverage", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Empirical coverage for regression coefficients with 95% Wald intervals",
      side = 3, outer = TRUE, line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))




get_reps_sim3 <- function(fit, metric, k, g) as.vector(fit$conditional[[metric]][k, g, ])


#FOR THE BIAS WE MAKE A BOX PLOT
fits <- list(sim3_beta_N150, sim3_beta_N300, sim3_beta_N500)
sample_sizes <- c(150, 300, 500)


metric <- "bias"


# shared y-axis range across both g
y_lim <- range(unlist(lapply(1:3, function(g)
  lapply(1:5, function(k)
    lapply(fits, get_reps_sim3, metric = metric, k = k, g = g)))), na.rm = TRUE)

pdf("./sim_study_plots/sim_study3_plots/sim3_beta_bias_plot.pdf", width = 7, height = 6)
par(mfrow = c(3, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 0.9)

for (g in 1:3) {
  for (k in 1:5) {
    dat <- lapply(fits, get_reps_sim3, metric = metric, k = k, g = g)
    names(dat) <- sample_sizes
    boxplot(dat,
            ylim = y_lim,
            yaxt = "n",
            pch = 16, cex = 0.4,
            main = bquote(beta[.(k-1) * .(g)]))
    if (k == 1) axis(2)
    abline(h = 0, lty = 2, col = "red")
  }
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("Bias", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Bias of coefficient estimates", side = 3, outer = TRUE,
      line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))




#Boxplots with outliers removed

metric <- "bias"


y_lim <- range(unlist(lapply(1:3, function(g)
  lapply(1:5, function(k)
    lapply(fits, function(f)
      boxplot.stats(get_reps_sim3(f, metric, k, g))$stats)))), na.rm = TRUE)
y_lim <- y_lim + c(-1, 1) * 0.05 * diff(y_lim)  

pdf("./sim_study_plots/sim_study3_plots/sim3_beta_bias_plot_no_outliers.pdf", width = 7, height = 4.4)
par(mfrow = c(3, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 0.9)

for (g in 1:3) {
  for (k in 1:5) {
    dat <- lapply(fits, get_reps_sim3, metric = metric, k = k, g = g)
    names(dat) <- sample_sizes
    boxplot(dat,
            ylim = y_lim,
            yaxt = "n",
            pch = 16, cex = 0.4,
            main = bquote(beta[.(k-1) * .(g)]),
            outline = FALSE)
    if (k == 1) axis(2)
    abline(h = 0, lty = 2, col = "red")
  }
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("Bias", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Bias of coefficient estimates", side = 3, outer = TRUE,
      line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))





#NOW DOING MSE

# columns ordered N150-g1, N150-g2, N300-g1, N300-g2, N500-g1, N500-g2
mse_mat <- cbind(
  sim3_beta_N150$conditional$mean_mse[1:5, ],
  sim3_beta_N300$conditional$mean_mse[1:5, ],
  sim3_beta_N500$conditional$mean_mse[1:5, ]
)

# Log scale needs positive values
y_lim <- range(mse_mat[mse_mat > 0], na.rm = TRUE)

pdf("./sim_study_plots/sim_study3_plots/sim3_beta_mse_plot.pdf",
    width = 7, height = 4.8)
par(mfrow = c(3, 5),
    oma = c(2, 3, 2, 0),
    mar = c(2, 0.6, 2, 0.6),
    mgp = c(2, 0.7, 0),
    cex.axis = 0.8,
    cex.main = 0.9)

for (g in 1:3) {
  cols <- seq(g, by = 3, length.out = 3)   # was by = 2
  for (k in 1:5) {
    plot(1:3, mse_mat[k, cols],
         log = "y",
         ylim = y_lim,
         xlim = c(0.5, 3.5),
         xaxt = "n", yaxt = "n",
         type = "b", pch = 16,
         xlab = "", ylab = "",
         main = bquote(beta[.(k - 1) * .(g)]))
    axis(1, at = 1:3, labels = sample_sizes)
    if (k == 1) axis(2, las = 1)
  }
}

mtext("Sample size", side = 1, outer = TRUE, line = 0.8, cex = 0.9)
mtext("MSE", side = 2, outer = TRUE, line = 1.8, cex = 0.9)
mtext("Mean squared error of coefficient estimates", side = 3, outer = TRUE,
      line = 0.5, font = 2, cex = 1)

dev.off()
par(mfrow = c(1, 1))





