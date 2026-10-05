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


make_theta_metric_table <- function(fit, true_theta, metric,
                                    metric_label, digits_metric = 3,
                                    caption, label,
                                    file = "", digits_true = 2,
                                    longtable = FALSE) {
  n_items   <- length(true_theta)
  n_classes <- nrow(true_theta[[1]])
  
  rows <- list()
  for (j in seq_len(n_items)) {
    m_true <- true_theta[[j]]
    m_est  <- fit[[metric]][[j]]
    ncat   <- sum(colSums(!is.na(m_true)) > 0)   # drop padded NA categories
    for (cc in seq_len(ncat)) {
      r <- list(Var = if (cc == 1) as.character(j) else "",
                Cat = as.character(cc))
      for (k in seq_len(n_classes)) {
        r[[paste0("True", k)]]   <- m_true[k, cc]
        r[[paste0("Metric", k)]] <- m_est[k, cc]
      }
      rows[[length(rows) + 1]] <- as.data.frame(r, stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows)
  
  colnames(tab) <- c(
    "\\textbf{Var.}", "\\textbf{Cat.}",
    unlist(lapply(seq_len(n_classes), function(k)
      c(sprintf("\\textbf{True %d}", k),
        sprintf("\\textbf{%s %d}", metric_label, k))))
  )
  
  dig <- c(0, 0, 0, rep(c(digits_true, digits_metric), n_classes))
  
  xt <- xtable(tab,
               digits  = dig,
               align   = c("l", "l", "l", rep("c", 2 * n_classes)),
               caption = caption,
               label   = label)
  
  print(xt,
        file = file,
        include.rownames = FALSE,
        booktabs = TRUE,
        floating = !longtable,
        tabular.environment = if (longtable) "longtable" else "tabular",
        table.placement = if (longtable) NULL else "H",
        caption.placement = "bottom",
        sanitize.text.function = identity,
        sanitize.colnames.function = identity,
        math.style.negative = TRUE,
        na.print = "")
}


sim_specs <- list(
  sim1 = list(fits = list("150" = sim1_theta_N150,
                          "300" = sim1_theta_N300,
                          "500" = sim1_theta_N500),
              true = sim1_theta, longtable = FALSE, num = 1),
  sim2 = list(fits = list("150" = sim2_theta_N150,
                          "300" = sim2_theta_N300,
                          "500" = sim2_theta_N500),
              true = sim2_theta, longtable = FALSE, num = 2),
  sim3 = list(fits = list("150" = sim3_theta_N150,
                          "300" = sim3_theta_N300,
                          "500" = sim3_theta_N500),
              true = sim3_theta, longtable = TRUE, num = 3)
)

dir.create("./sim_study_tables", showWarnings = FALSE)

for (s in names(sim_specs)) {
  sp <- sim_specs[[s]]
  for (n in names(sp$fits)) {
    for (m in names(metrics)) {
      spec <- metrics[[m]]
      make_theta_metric_table(
        fit           = sp$fits[[n]],
        true_theta    = sp$true,
        metric        = spec$slot,
        metric_label  = spec$label,
        digits_metric = spec$digits,
        longtable     = sp$longtable,
        caption = sprintf("True item probability parameters and %s of their estimates for simulation %d ($N = %s$).",
                          spec$long, sp$num, n),
        label   = sprintf("table:%s_item_prob_%s_N%s", s, m, n),
        file    = sprintf("./sim_study_tables/%s_item_prob_%s_N%s.tex", s, m, n)
      )
    }
  }
}