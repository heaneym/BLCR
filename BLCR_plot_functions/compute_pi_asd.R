#Function for computing the mixture proportions for a given ASD status (ASD/non-ASD) and estimate for beta

#Note that the beta needs to be formatted as a 2 x G matrix with the 1st row being intercept and 2nd ASD coefficient.

#This function assumes that the baseline group corresponds to column 1 of beta. 

compute_pi_asd <- function(asd, beta) {
  G <- ncol(beta)
  
  eta <- numeric(G)
  eta[1] <- 0  
  
  if (G > 1) {
    for (g in 2:G) {
      eta[g] <- beta[1, g] + beta[2, g] * asd
    }
  }
  
  pi <- exp(eta) / sum(exp(eta))
  pi
}
