


# Tailored visual for CSHQ analysis - but likely generalisable

# Function for plotting the prdictive distribution of the total CSHQ score with specified ci
# and background stacked plot indicating group membership. Also stratified by ASD

plot_combined_predictive_T <- function(theta, beta, theta_samples, beta_samples,
                                       scores = c(1, 2, 3), cutoff = 41, 
                                       ci_level = 0.95,
                                       y_zoom = c(0, 0.25),
                                       auto_ylim = TRUE,
                                       show_progress = TRUE) {  
  
  res0 <- predictive_pmf_T_by_ASD(0, theta, beta, scores)
  res1 <- predictive_pmf_T_by_ASD(1, theta, beta, scores)
  mean_pmf_df <- rbind(
    data.frame(ASD = "No ASD", T_int = res0$totals, pmf_mean = res0$pmf_T),
    data.frame(ASD = "ASD",   T_int = res1$totals, pmf_mean = res1$pmf_T)
  )
  
  
  n_iter <- dim(theta_samples)[4]
  
  if (show_progress) {
    cat("Computing posterior predictive PMFs for", n_iter, "iterations...\n")
    pb <- progress::progress_bar$new(
      format = "  [:bar] :percent eta: :eta",
      total = n_iter * 2,
      clear = FALSE,
      width = 60
    )
  }
  
  pmf_all <- data.frame()
  for (t in 1:n_iter) {
    for (asd in c(0, 1)) {
      df_iter <- predictive_pmf_T_by_ASD_iter(t, asd, theta_samples, beta_samples, scores)
      pmf_all <- rbind(pmf_all, df_iter)
      
      if (show_progress) {
        pb$tick()
      }
    }
  }
  
  if (show_progress) {
    cat("Computing credible intervals...\n")
  }
  
  
  alpha <- 1 - ci_level
  ci_df <- pmf_all %>%
    group_by(ASD, T_int) %>%
    summarise(
      pmf_lower = quantile(pmf_T, probs = alpha/2, na.rm = TRUE),
      pmf_upper = quantile(pmf_T, probs = 1 - alpha/2, na.rm = TRUE),
      .groups = "drop"
    )
  
  mean_pmf_df <- mean_pmf_df %>%
    left_join(ci_df, by = c("ASD", "T_int"))
  
  
  if (auto_ylim) {
    max_pmf <- max(mean_pmf_df$pmf_upper, na.rm = TRUE)
    y_zoom <- c(0, max_pmf * 1.05)  
  }
  
  if (show_progress) {
    cat("Building plot...\n")
  }
  
  
  res <- build_plot_df(theta, beta, scores)
  df_long <- res$df_long
  max_p_original <- max(df_long$p, na.rm = TRUE)
  df_long_scaled <- df_long %>%
    mutate(p_scaled = p * (y_zoom[2] / max_p_original)) %>%
    mutate(group = interaction(ASD, class))
  
  G <- ncol(beta)
  cols <- RColorBrewer::brewer.pal(max(G, 3), "Set2")[1:G]
  
  pmf_color <- "#2c3e50"  
  prob_color <- "#7f8c8d"  
  
  p <- ggplot() +
    
    geom_area(data = df_long_scaled, 
              aes(x = T_int, y = p_scaled, fill = factor(class), group = group),
              position = "stack", alpha = 0.25, colour = NA) +
    
    geom_ribbon(data = mean_pmf_df,
                aes(x = T_int, ymin = pmf_lower, ymax = pmf_upper, group = ASD),
                fill = "grey70", alpha = 0.3, inherit.aes = FALSE) +
    
    geom_line(data = mean_pmf_df, 
              aes(x = T_int, y = pmf_mean, colour = ASD),
              linewidth = 2.2, alpha = 1) +
    
    geom_vline(xintercept = cutoff, linetype = "dashed", colour = "red", linewidth = 1.2) +
    facet_wrap(~ ASD, ncol = 1) +
    
    scale_fill_brewer(palette = "Set2", 
                      name = "Latent Class",
                      labels = c("A", "B", "C", "D")) +
    scale_colour_manual(values = c("No ASD" = "#1f78b4", "ASD" = "#e31a1c"), 
                        name = "ASD Status") +
    coord_cartesian(xlim = c(33, 80), ylim = y_zoom) +
    scale_y_continuous(
      name = "Predictive PMF (red/blue)",
      sec.axis = sec_axis(~ . / y_zoom[2], 
                          name = "P(group g | score T) (stacked areas)",
                          breaks = seq(0, 1, by = 0.25))
    ) +
    labs(
      x = "Total CSHQ score T",
      title = "Posterior Predictive Distributions of Total CSHQ Score by ASD Status",
      subtitle = sprintf("Point estimates with %d%% credible intervals and stacked latent class PMFs", 
                         round(ci_level * 100))
    ) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      legend.position = "bottom",
      axis.title.y.left = element_text(colour = pmf_color, face = "bold", size = 16),
      axis.text.y.left = element_text(colour = pmf_color, size = 14),
      axis.text.y.right = element_text(colour = prob_color, size = 14),
      axis.title.y.right = element_text(colour = prob_color, face = "bold", size = 16),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.ticks.y.left = element_line(colour = pmf_color),
      axis.ticks.y.right = element_line(colour = prob_color),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 13),
      strip.text = element_text(size = 14, face = "bold"),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  if (show_progress) {
    cat("Done!\n")
  }
  
  p
}