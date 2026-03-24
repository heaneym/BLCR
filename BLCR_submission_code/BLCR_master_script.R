### ------- BLCR Master Script ------- ###

#This file gets the relevant packages, then the files and runs through simulation 1, 2 and the CSHQ data analysis

# ---- Checking and loading packages ----
required_packages <- c(
  "ggplot2","RColorBrewer","reshape2","BayesLCA","bayestestR",
  "mcclust.ext","dplyr","tidyr","ggridges","MCMCpack",
  "BayesLogit","label.switching","einsum","Rcpp",
  "foreign","fossil","patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))


# ---- Sourcing functions ----
helper_files <- c(
  "LCR_Gibbs.R",
  "LCR_sim_data.R",
  "beta_summary_table.R",
  "LCR_posterior_membership_prob.R",
  "LCA_mosaic_plot.R",
  "plot_combined_predictive_T.R",
  "compute_expected_subscale_total_across_iterations.R",
  "predictive_pmf_T_by_ASD_iter.R",
  "compute_pi_asd.R",
  "pmf_T_given_class.R",
  "predictive_pmf_T_by_ASD.R",
  "build_plot_df.R",
  "plot_predictive_T.R",
  "compute_pi_asd_across_iterations.R",
  "compute_theta_samples_from_counts.R",
  "plot_subscale_diff_boxplot.R",
  "combined_beta_ridgeline_plot.R",
  "cshq_profile_plot_functions.R"
)
invisible(lapply(file.path("R", helper_files), source))


# ---- Getting Plot colours ----
cb_blue <- "#4477AA"
cb_orange <- "#EE7733"
cb_cols <- c(cb_blue, cb_orange)

# ---- Run scripts ----
source(file.path("scripts", "01_simulation_study_1.R"))
source(file.path("scripts", "01_1_simulation_study_1_sample_sensitivity.R"))
source(file.path("scripts", "02_simulation_study_2.R"))
source(file.path("scripts", "02_1_simulation_study_2_sample_sensitivity.R"))
source(file.path("scripts", "03_cshq_analysis.R"))











