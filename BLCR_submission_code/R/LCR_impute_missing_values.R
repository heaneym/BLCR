impute_missing_values <- function(Y, Y_missing_indicator, z, theta, M, K) {
  n <- nrow(Y)
  Y_imputed <- Y
  
  for (i in 1:n) {
    class_i <- which(z[i,] == 1)
    
    for (j in 1:M) {
      if (Y_missing_indicator[i, j] == 1) {  
        probs <- theta[[j]][class_i, ]
        Y_imputed[i, j] <- sample(1:K[j], size = 1, prob = probs)
      }
    }
  }
  
  return(Y_imputed)
}
