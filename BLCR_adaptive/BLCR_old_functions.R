#Just Putting the old code in here for BLCR Probit

#This is for getting S in the loop - obviously the triple nested for loop is not ideal
#for (j in 1:M){
#  for (g in 1:G){
#    for (k in 1:K[j]){
#      S[j,g,k] = sum((Y[,j] == k)*z[,g])
#    }
#  }
#}


#Old code for generating theta
#for (g in 1:G) { #generating a sample from theta_{gj} vectors for each g and j}
#  for (j in 1:M) {
#    theta[g, j, 1:K[j]] <- rdirichlet(1, alpha[1:K[j]] + S[g,j,1:K[j]])  
#  }
#}


#These are various ones for part of getting the double product for the vector w - I dont think they work

#-Inf in the theta will cause issues so for the purpose of the sum Im going to put zeroes in their place
#log_theta[log_theta == -Inf] <- 0 
#log_sum1 <- vapply(1:n, function(i) sum(log_theta[1,,Y[i,]]), numeric(1))
#log_sum2 <- vapply(1:n, function(i) sum(log_theta[2,,Y[i,]]), numeric(1))
#for (i in 1:n){
#  log_sum1[i] <- sum(log_theta[1, 1:M, Y[i, 1:M]])
#  log_sum2[i] <- sum(log_theta[2, 1:M, Y[i, 1:M]])
#}


#This one is a slow way to do the w vector, it does work im pretty sure
#for (i in 1:n){
#  log_sum1[i] <- 0
#  log_sum2[i] <- 0
#  for (j in 1:M){
#    for (k in 1:K[j]){
#      log_sum1[i] <- log_sum1[i] + (Y[i,j] == k)*log_theta[1,j,Y[i,j]]
#      log_sum2[i] <- log_sum2[i] + (Y[i,j] == k)*log_theta[2,j,Y[i,j]]
#    }
#  }
#}



#Old code for getting z
#for (i in 1:n){
#  I1 <- Phi_vec[i]
#  I0 <- not_Phi_vec[i]
#  p1 <- 1
#  p0 <- 1
#  for(j in 1:M){
#    for(k in 1:K[j]){
#      p1<- p1*theta[1,j,k]^(Y[i,j] == k)
#      p0<- p0*theta[2,j,k]^(Y[i,j] == k)
#    }
#  }
#  w[i, ] <- c(I1*p1,I0*p0)
#  z[i,] <- rmultinom(1,1,w[i,])
#}


