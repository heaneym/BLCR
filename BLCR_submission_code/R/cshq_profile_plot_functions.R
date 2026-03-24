compute_cshq_subscale_totals <- function(Y_reduced, standardise = TRUE) {
  mins <- c(BR = 3, SOD = 1, SD = 3, SA = 4, NW = 3, P = 3, DS = 4)
  
  totals <- cbind(
    BR  = rowSums(Y_reduced[, 1:3, drop = FALSE])   - mins["BR"],
    SOD = Y_reduced[, 4]                            - mins["SOD"],
    SD  = rowSums(Y_reduced[, 5:7, drop = FALSE])   - mins["SD"],
    SA  = rowSums(Y_reduced[, 8:11, drop = FALSE])  - mins["SA"],
    NW  = rowSums(Y_reduced[, 12:14, drop = FALSE]) - mins["NW"],
    P   = rowSums(Y_reduced[, 15:17, drop = FALSE]) - mins["P"],
    DS  = rowSums(Y_reduced[, 18:21, drop = FALSE]) - mins["DS"]
  )
  
  out <- list(raw = totals)
  
  if (standardise) {
    out$standardised <- scale(totals)
  }
  
  out
}

compute_profile_matrix <- function(subscale_mat, cluster_labels) {
  groups <- sort(unique(cluster_labels))
  
  profile_mat <- do.call(
    rbind,
    lapply(groups, function(g) {
      colMeans(subscale_mat[cluster_labels == g, , drop = FALSE])
    })
  )
  
  rownames(profile_mat) <- paste0("Group ", groups)
  profile_mat
}

compute_cshq_profiles <- function(Y_reduced, cluster_labels, standardise = TRUE) {
  totals <- compute_cshq_subscale_totals(Y_reduced, standardise = standardise)
  
  out <- list(
    profile = compute_profile_matrix(totals$raw, cluster_labels)
  )
  
  if (standardise && !is.null(totals$standardised)) {
    out$profile_standardised <- compute_profile_matrix(totals$standardised, cluster_labels)
  }
  
  out
}

profile_mat_to_long <- function(profile_mat) {
  df <- as.data.frame(profile_mat)
  df$Group <- rownames(df)
  
  df_long <- df %>%
    pivot_longer(
      cols = -Group,
      names_to = "Subscale",
      values_to = "Value"
    )
  
  df_long$Group <- factor(df_long$Group, levels = rownames(profile_mat))
  df_long$Subscale <- factor(df_long$Subscale, levels = colnames(profile_mat))
  
  df_long
}

plot_cshq_profiles <- function(
    profile_mat,
    title = "CSHQ Profile Plot",
    subtitle = NULL,
    ylab = "Standardised mean score",
    palette = NULL
) {
  profile_long <- profile_mat_to_long(profile_mat)
  
  p <- ggplot(
    profile_long,
    aes(x = Subscale, y = Value, group = Group, colour = Group)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Subscale",
      y = ylab,
      colour = "Group"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  if (!is.null(palette)) {
    p <- p + scale_colour_manual(values = palette)
  }
  
  p
}