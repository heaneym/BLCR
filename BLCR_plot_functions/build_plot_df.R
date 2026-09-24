
#Function for constructing data frame for plotting the pmf for ASD vs. non-ASD.

build_plot_df <- function(theta, beta, scores = c(1, 2, 3)) {
  res0 <- predictive_pmf_T_by_ASD(0, theta, beta, scores)
  res1 <- predictive_pmf_T_by_ASD(1, theta, beta, scores)
  
  G <- ncol(beta)
  
  df0 <- as.data.frame(res0$contrib)
  df0$T_int <- res0$totals
  df0$ASD   <- "No ASD"
  
  df1 <- as.data.frame(res1$contrib)
  df1$T_int <- res1$totals
  df1$ASD   <- "ASD"
  
  df <- rbind(df0, df1)
  
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = dplyr::matches("^[0-9]+$"),  
    names_to = "class",
    values_to = "p_joint"
  )
  
  
  df_long <- df_long %>%
    dplyr::group_by(ASD, T_int) %>%
    dplyr::mutate(
      p = p_joint / sum(p_joint)
    ) %>%
    dplyr::ungroup()
  
  
  pmf_df <- rbind(
    data.frame(
      ASD   = "No ASD",
      T_int = res0$totals,
      pmf_T = res0$pmf_T
    ),
    data.frame(
      ASD   = "ASD",
      T_int = res1$totals,
      pmf_T = res1$pmf_T
    )
  )
  
  list(df_long = df_long, pmf_df = pmf_df)
}
