N_updates <- function(z, Y_indicator, N_jk, nu){
  N_g <- colSums(z)
  N_gjk <- apply(Y_indicator, c(2, 3), function(S_jk) t(z) %*% S_jk)
  return(list(N_g = N_g, N_gjk = N_gjk))
}


nu_update_collapsed <- function(nu, M, K, G, theta_hyperparam, N_g, N_gjk, inclusion_sum, exclusion_sum, count){
  nu_prop <- nu
  j_prop <- sample(1:M,1)
  if(nu[j_prop] == 1){
    nu_prop[j_prop] <- 0
    log_gamma_N_alpha <- lgamma(N_gjk[,j_prop,] + theta_hyperparam)
    log_sum1 <- sum(log_gamma_N_alpha)
    diff_log <- sum(lgamma(N_g + K[j_prop]*theta_hyperparam)) - log_sum1
    log_accept_ratio <- exclusion_sum[j_prop] + diff_log
  } else {
    nu_prop[j_prop] <- 1
    log_gamma_N_alpha <- lgamma(N_gjk[,j_prop,] + theta_hyperparam)
    log_sum1 <- sum(log_gamma_N_alpha)
    diff_log <- log_sum1 - sum(lgamma(N_g + K[j_prop]*theta_hyperparam))
    log_accept_ratio <- inclusion_sum[j_prop] + diff_log
    
  }
  accept_ratio <- exp(log_accept_ratio)
  accept_ratio <- as.numeric(accept_ratio)
  if(runif(1) < min(1,accept_ratio)){
    nu <- nu_prop
  }
  return(nu)
}

get_nu_update <- function(item.sel){
  if (item.sel){
    nu_update <- nu_update_collapsed
  } else {
    nu_update <- nu_update_uncollapsed
  }
  return(nu_update)
}

S_update_collapsed <- function(Y_indicator,z,M,G,K,S){
  return(S)
}

theta_update_collapsed <- function(S,theta_hyperparam,K,G,M,theta){
  return(theta)
}