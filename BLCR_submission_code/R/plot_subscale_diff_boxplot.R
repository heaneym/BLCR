#Creating a function to plot the subscale total difference boxplots

plot_subscale_diff_boxplot <- function(
    theta_samples,
    beta_samples,
    subscale_names = c("BR", "SOD", "SD", "SA", "NW", "P", "SDB", "DS")
) {
  expected_asd <- compute_expected_subscale_total_across_iterations(
    theta_samples = theta_samples,
    asd = 1,
    beta_samples = beta_samples
  )
  
  expected_non_asd <- compute_expected_subscale_total_across_iterations(
    theta_samples = theta_samples,
    asd = 0,
    beta_samples = beta_samples
  )
  
  diff_df_long <- as.data.frame(
    expected_asd$expected_subscale_total - expected_non_asd$expected_subscale_total
  ) %>%
    setNames(subscale_names) %>%
    pivot_longer(everything(), names_to = "Column", values_to = "Value")
  
  label_map <- c(
    "BR" = "Bedtime Resistance",
    "SOD" = "Sleep Onset Delay",
    "SD" = "Sleep Duration",
    "SA" = "Sleep Anxiety",
    "NW" = "Night Waking",
    "P" = "Parasomnias",
    "SDB" = "Sleep Disordered Breathing",
    "DS" = "Daytime Sleepiness"
  )
  
  diff_df_long <- diff_df_long %>%
    mutate(
      Column_long = label_map[Column],
      median_val = tapply(Value, Column, median)[Column]
    ) %>%
    arrange(median_val) %>%
    mutate(Column_long = factor(Column_long, levels = unique(Column_long)))
  
  ggplot(diff_df_long, aes(x = Column_long, y = Value, fill = Column_long)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "white") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1.2) +
    stat_summary(
      fun.data = "median_hilow",
      geom = "errorbar",
      fun.args = list(conf.int = 0.95),
      width = 0.2,
      linewidth = 1
    ) +
    stat_summary(fun = "median", geom = "point", size = 4, color = "white") +
    labs(
      title = "Posterior Differences in Expected CSHQ Subscale Totals (ASD - Non-ASD)",
      x = "CSHQ Subscales",
      y = "Difference in Expected Total"
    ) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_fill_brewer(palette = "Dark2")
}