#This is a function for computing the predictive pmf across CSHQ scores given membership in group g

#returns a vector of probabilities corresponding to the mass at each score

pmf_T_given_class <- function(g, theta, scores = c(1, 2, 3)) {
  M <- length(theta)
  prob_mat1 <- theta[[1]]       
  K1 <- ncol(prob_mat1)
  s1_vals <- sort(unique(scores[1:K1]))
  pmf <- setNames(numeric(length(s1_vals)), s1_vals)
  for (k in 1:K1) {
    s <- scores[k]
    pmf[as.character(s)] <- pmf[as.character(s)] + prob_mat1[g, k]
  }
  for (j in 2:M) {
    prob_mat <- theta[[j]]
    Kj <- ncol(prob_mat)
    old_totals <- as.integer(names(pmf))
    new_min <- min(old_totals) + min(scores[1:Kj])
    new_max <- max(old_totals) + max(scores[1:Kj])
    new_totals <- new_min:new_max
    pmf_new <- setNames(numeric(length(new_totals)), new_totals)
    
    for (t_str in names(pmf)) {
      p_t <- pmf[t_str]
      if (p_t == 0) next
      t <- as.integer(t_str)
      for (k in 1:Kj) {
        s <- scores[k]
        t_new <- t + s
        pmf_new[as.character(t_new)] <-
          pmf_new[as.character(t_new)] + p_t * prob_mat[g, k]
      }
    }
    pmf <- pmf_new
  }
  pmf
}