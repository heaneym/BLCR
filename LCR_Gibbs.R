library(MCMCpack)
library(BayesLogit)
library(label.switching)
library(einsum)
#source z_update functions
source("./LCR_z_update.R")
#source init function
source("./LCR_init.R")
#source collapsed sampler functions
source("./LCR_collapsed_functions.R")
#source uncollapsed sampler functions
source("./LCR_uncollapsed_functions.R")
#source gamma functions
source("./LCR_gamma_update.R")
#source general functions
source("./LCR_general_functions.R")
#source post processing functions
source("./LCR_post_process.R")
#source relabelling functions
source("./LCR_relabel_outputs.R")
#source relabelling functions
source("./LCR_log_post_compute.R")
#source main function
source("./LCR_Gibbs_function.R")










