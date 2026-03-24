#Creating a function to plot the subscale total difference boxplots

plot_subscale_diff_boxplot <- function(theta_samples, beta_samples, subscale_names = c('BR','SOD','SD','SA','NW','P','SDB','DS')){
  compute_expected_subscale_totals_asd <- compute_expected_subscale_total_across_iterations(theta_samples = theta_samples_collapsed_estimate_perm,
                                                                                                            asd = 1, 
                                                                                                            beta_samples = CSHQ_LCR_varsel$samples$beta_samples)
  
  compute_expected_subscale_totals_non_asd <- compute_expected_subscale_total_across_iterations(theta_samples = theta_samples_collapsed_estimate_perm,
                                                                                                                asd = 0, 
                                                                                                                beta_samples = CSHQ_LCR_varsel$samples$beta_samples)
  
  diff_expected_subscale_totals <- compute_expected_subscale_totals_asd$expected_subscale_total - compute_expected_subscale_totals_non_asd$expected_subscale_total
  diff_df <- as.data.frame(diff_expected_subscale_totals)
  colnames(diff_df) <- subscale_names
  
  
  diff_df_long <- diff_df %>%
    pivot_longer(everything(), names_to = "Column", values_to = "Value")
  
  
  plot_result <- ggplot(diff_df_long, aes(x = Column, y = Value, fill = Column)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Boxplots of Difference between Subscale Totals", x = "Columns", y = "Values")
  
  return(plot_result)
}