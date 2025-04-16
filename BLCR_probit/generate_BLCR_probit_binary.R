library(MCMCpack)
library(extraDistr)
#Need to fix this!
generate_BLCR_probit_2group <- function(beta, theta, n, G = 2){
  p <- length(beta) - 1
  M <- length(theta)
  K <- unlist(lapply(theta, ncol))
  X_mat <- cbind(1,matrix(rnorm(p*n), nrow = n))
  column_names <- c("intercept", paste("X", 1:p, sep = ""))
  X <- as.data.frame(X_mat)
  colnames(X) <- column_names
  mu <- as.matrix(X)%*%beta
  probs <- pnorm(mu)
  class <- ifelse(probs>=0.5,1,2)
  class1size <- sum(class==1)
  class2size <- n-class1size
  groupsize <- c(class1size,class2size)
  pi <- groupsize/n
  Y <- data.frame(matrix(ncol = M, nrow = 0))
  colnames(Y) <- paste("Y", 1:M, sep = "")
  for (i in 1:n){
    cl_row <- class[i]
    response <- numeric(M)
    for (j in 1:M){
        levels <- 1:K[j]
        probabilities <- theta[[j]][cl_row,]
        response[j] <- sample(levels, size = 1, replace = TRUE, prob = probabilities)
    }
    response <- as.data.frame(t(response))
    colnames(response) <- colnames(Y) 
    Y <- rbind(Y,response)
  }
}

