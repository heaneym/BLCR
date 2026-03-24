prepare_beta_long_df <- function(beta_samples, beta_names, group = NULL) {
  rownames(beta_samples) <- beta_names
  
  out <- as.data.frame(t(beta_samples)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Parameter",
      values_to = "Value"
    ) %>%
    mutate(
      Parameter = factor(Parameter, levels = rev(beta_names))
    )
  
  if (!is.null(group)) {
    out <- out %>% mutate(Group = group)
  }
  
  out
}


prepare_true_beta_df <- function(true_beta, beta_names, group = NULL) {
  out <- data.frame(
    Parameter = beta_names,
    TrueValue = true_beta
  ) %>%
    mutate(
      Parameter = factor(Parameter, levels = rev(beta_names))
    )
  
  if (!is.null(group)) {
    out <- out %>% mutate(Group = group)
  }
  
  out
}


plot_beta_ridgeline <- function(
    long_df,
    true_df = NULL,
    math_labels,
    fill_col = cb_blue,
    true_col = cb_orange,
    xlim = NULL,
    title = "Ridgeline Plot of Posterior Distributions",
    subtitle = "Posterior density with median & 95% HDI"
) {
  p <- ggplot(long_df, aes(x = Value, y = Parameter)) +
    geom_density_ridges(
      quantile_lines = TRUE,
      quantiles = c(0.025, 0.5, 0.975),
      fill = fill_col,
      alpha = 0.7,
      scale = 0.9
    ) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
    scale_y_discrete(labels = math_labels) +
    labs(
      x = "Logit Coefficient",
      y = "",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  if (!is.null(true_df)) {
    p <- p +
      geom_point(
        data = true_df,
        aes(x = TrueValue, y = Parameter),
        color = true_col,
        shape = 18,
        size = 3,
        inherit.aes = FALSE
      ) +
      geom_text(
        data = true_df,
        aes(x = TrueValue, y = Parameter, label = round(TrueValue, 2)),
        color = true_col,
        vjust = 2,
        size = 3,
        inherit.aes = FALSE
      )
  }
  
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim, clip = "off")
  }
  
  p
}


plot_combined_beta_ridgeline <- function(
    long_dfs,
    true_dfs = NULL,
    group_order,
    custom_labels,
    fill_cols,
    xlim = c(-3, 3),
    title = "Ridgeline Plot of Posterior Distributions",
    subtitle = "Posterior density with median and 95% HDI. Diamonds mark true values."
) {
  param_levels <- levels(long_dfs[[1]]$Parameter)
  
  long_both <- bind_rows(long_dfs) %>%
    mutate(
      Parameter = factor(Parameter, levels = param_levels),
      Group = factor(Group, levels = group_order),
      Strip = interaction(Group, Parameter, sep = " • ", lex.order = TRUE)
    )
  
  strip_df <- distinct(long_both, Strip)
  
  xmin <- xlim[1]
  xmax <- xlim[2]
  xbreaks <- seq(xmin, xmax, by = 1)
  
  p <- ggplot(long_both, aes(x = Value, y = Strip, fill = Group)) +
    geom_density_ridges(
      quantile_lines = TRUE,
      quantiles = c(0.025, 0.5, 0.975),
      scale = 0.9,
      alpha = 0.75,
      color = "gray30",
      linewidth = 0.25
    ) +
    geom_segment(
      data = strip_df,
      aes(x = xmin, xend = xmax, y = Strip, yend = Strip),
      inherit.aes = FALSE,
      color = "gray70",
      linewidth = 0.6
    ) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
    scale_fill_manual(values = fill_cols) +
    scale_y_discrete(
      labels = custom_labels,
      expand = c(0.01, 0)
    ) +
    scale_x_continuous(
      limits = c(xmin, xmax),
      breaks = xbreaks,
      expand = expansion(mult = 0)
    ) +
    labs(
      x = "Logit Coefficient",
      y = "",
      title = title,
      subtitle = subtitle,
      fill = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "top"
    )
  
  if (!is.null(true_dfs)) {
    true_both <- bind_rows(true_dfs) %>%
      mutate(
        Parameter = factor(Parameter, levels = param_levels),
        Group = factor(Group, levels = group_order),
        Strip = interaction(Group, Parameter, sep = " • ", lex.order = TRUE)
      )
    
    p <- p +
      geom_point(
        data = true_both,
        aes(x = TrueValue, y = Strip, color = Group),
        shape = 18,
        size = 2.8,
        show.legend = FALSE
      ) +
      geom_text(
        data = true_both,
        aes(x = TrueValue, y = Strip, label = round(TrueValue, 2), color = Group),
        vjust = 2,
        size = 3,
        show.legend = FALSE
      ) +
      scale_color_manual(values = fill_cols)
  }
  
  p
}