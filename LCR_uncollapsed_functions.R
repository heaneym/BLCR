S_update_uncollapsed <- function(Y_indicator,z,M,G,K,S){
  S <- array(0, dim = c(G, M, max(K)))
  for (k in 1:max(K)) {
    S[,,k] <- t(z)%*%Y_indicator[,,k]
  }
  return(S)
}

theta_update_uncollapsed <- function(S,theta_hyperparam,K,G,M,theta){
  alpha_S <- array(theta_hyperparam, dim = dim(S)) + S
  gam_samples <- array(rgamma(prod(dim(S)), shape = as.vector(alpha_S), scale = 1),dim = dim(S))
  gam_sums <- apply(gam_samples, c(1, 2), sum)
  theta <- sweep(gam_samples, MARGIN = c(1, 2), gam_sums, FUN = "/")
  return(theta)
}

nu_update_uncollapsed <- function(nu, M, K, G, theta_hyperparam, N_g, N_gjk, inclusion_sum, exclusion_sum, count){
  return(nu)
}

get_S_update <- function(item.sel){
  if (item.sel){
    S_update <- S_update_collapsed
  } else {
    S_update <- S_update_uncollapsed
  }
  return(S_update)
}

get_theta_update <- function(item.sel){
  if (item.sel){
    theta_update <- theta_update_collapsed
  } else {
    theta_update <- theta_update_uncollapsed
  }
  return(theta_update)
}

