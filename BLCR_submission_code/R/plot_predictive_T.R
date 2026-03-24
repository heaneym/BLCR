#This function plots the predictive pmf for ASD and non-ASD, without any uncertainty quantification (point estimate only).

plot_predictive_T <- function(theta, beta, scores = c(1, 2, 3), cutoff = 41) {
  G <- ncol(beta)
  cols <- brewer.pal(max(G, 3), "Set2")[1:G]
  
  res <- build_plot_df(theta, beta, scores)
  df_long <- res$df_long
  pmf_df  <- res$pmf_df
  
  ggplot() +
    geom_area(
      data = df_long,
      aes(x = T_int, y = p, fill = factor(class)),
      position = "stack",
      alpha = 0.5,
      colour = NA
    ) +
    geom_line(
      data = pmf_df,
      aes(x = T_int, y = pmf_T / max(pmf_T), group = ASD),  
      colour = "black",
      linewidth = 0.6
    ) +
    geom_vline(xintercept = cutoff, linetype = "dashed", colour = "red") +
    facet_wrap(~ ASD, ncol = 1) +
    scale_fill_manual(values = cols, name = "Latent class") +
    xlab("Total CSHQ score") +
    ylab("Probability (stacked by class)") +
    coord_cartesian(xlim = c(33, 80)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}