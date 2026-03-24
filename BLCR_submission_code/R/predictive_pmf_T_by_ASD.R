#This is a function for computing the predictive pmf of the CSHQ total, given ASD status (ASD/non-ASD) 

predictive_pmf_T_by_ASD <- function(asd, theta, beta, scores = c(1, 2, 3)) {
  M <- length(theta)
  min_score <- M * min(scores)
  max_score <- M * max(scores)
  totals <- min_score:max_score
  
  G <- ncol(beta)
  pi <- compute_pi_asd(asd, beta)
  
  
  pmf_classes <- lapply(1:G, function(g) pmf_T_given_class(g, theta, scores))
  
  pmf_mat <- matrix(0, nrow = length(totals), ncol = G,
                    dimnames = list(as.character(totals), 1:G))
  for (g in 1:G) {
    v <- pmf_classes[[g]]
    pmf_mat[names(v), as.character(g)] <- v
  }
  contrib <- sweep(pmf_mat, 2, pi, `*`)         
  pmf_T   <- rowSums(contrib)                  
  
  list(
    totals = totals,
    pmf_T  = pmf_T,
    contrib = contrib, 
    pi = pi
  )
}