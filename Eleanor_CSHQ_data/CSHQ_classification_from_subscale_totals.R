#Creating a function to compute the probability of group membership in the CSHQ LCA/LCR model for a given observation, given that we only know the subscale totals.

generate_patterns_for_sum <- function(sum_score, n_items, n_categories = 3,
                                      min_response = 1, max_response = 3) {
  patterns <- list()
  
  generate_recursive <- function(remaining_sum, remaining_items, current_pattern) {
    if (remaining_items == 0) {
      if (remaining_sum == 0) {
        patterns <<- c(patterns, list(current_pattern))
      }
      return()
    }
    for (response in min_response:max_response) {
      if (response <= remaining_sum) {
        generate_recursive(
          remaining_sum - response,
          remaining_items - 1,
          c(current_pattern, response)
        )
      }
    }
  }
  
  generate_recursive(sum_score, n_items, c())
  
  if (length(patterns) > 0) {
    return(do.call(rbind, patterns))
  } else {
    return(matrix(nrow = 0, ncol = n_items))
  }
}



calc_prob_sum_given_class <- function(sum_score, class_num, item_indices,
                                      item_probs, min_response = 1,
                                      max_response = 3) {
  n_items <- length(item_indices)
  
  patterns <- generate_patterns_for_sum(sum_score, n_items, 
                                        n_categories = 3,
                                        min_response = min_response,
                                        max_response = max_response)
  
  if (nrow(patterns) == 0) {
    return(0)  
  }
  
  prob_sum <- 0
  
  for (i in 1:nrow(patterns)) {
    pattern <- patterns[i, ]
    prob_pattern <- 1
    valid_pattern <- TRUE
    
    for (j in 1:n_items) {
      item_idx <- item_indices[j]
      response <- pattern[j]
      
      item_matrix <- item_probs[[item_idx]]
      if (response > ncol(item_matrix) || response < 1) {
        valid_pattern <- FALSE
        break
      }
      
      prob_pattern <- prob_pattern * item_matrix[class_num, response]
    }
    
    if (valid_pattern) {
      prob_sum <- prob_sum + prob_pattern
    }
  }
  
  return(prob_sum)
}





classify_observation <- function(subscale_sums, subscale_items, pi, item_probs,
                                 min_response = 1, max_response = 3) {
  G <- length(pi)
  n_subscales <- length(subscale_sums)
  
  log_likelihoods <- rep(0, G)  
  
  # For each subscale, calculate log P(S_k | C) for each class
  for (l in 1:n_subscales) {
    sum_l <- subscale_sums[l]
    items_l <- subscale_items[[l]]
    
    for (g in 1:G) {
      prob_sum_given_class <- calc_prob_sum_given_class(
        sum_score = sum_l,
        class_num = g,
        item_indices = items_l,
        item_probs = item_probs,
        min_response = min_response,
        max_response = max_response
      )

      if (prob_sum_given_class > 0) {
        log_likelihoods[g] <- log_likelihoods[g] + log(prob_sum_given_class)
      } else {
        log_likelihoods[g] <- -Inf
      }
    }
  }
  
  # Add log prior: log P(C | S) = log P(S | C) + log P(C) - log P(S)
  log_posteriors <- log_likelihoods + log(pi)
  posteriors <- exp(log_posteriors - max(log_posteriors))
  posteriors <- posteriors / sum(posteriors)
  
  return(posteriors)
}


classify_observation_regression <- function(subscale_sums, subscale_items, beta, x, item_probs,
                                 min_response = 1, max_response = 3) {
  G <- ncol(beta)
  p <- length(x)
  n_subscales <- length(subscale_sums)
  
  log_likelihoods <- rep(0, G)  
  pi_unnormalised <- numeric(G)
  
  # For each subscale, calculate log P(S_k | C) for each class
  for (l in 1:n_subscales) {
    sum_l <- subscale_sums[l]
    items_l <- subscale_items[[l]]
    
    for (g in 1:G) {
      pi_unnormalised[g] <- exp(c(1,x)%*%beta[,g])
      prob_sum_given_class <- calc_prob_sum_given_class(
        sum_score = sum_l,
        class_num = g,
        item_indices = items_l,
        item_probs = item_probs,
        min_response = min_response,
        max_response = max_response
      )
      
      if (prob_sum_given_class > 0) {
        log_likelihoods[g] <- log_likelihoods[g] + log(prob_sum_given_class)
      } else {
        log_likelihoods[g] <- -Inf
      }
    }
  }
  pi <- pi_unnormalised/sum(pi_unnormalised)
  
  # Add log prior: log P(C | S) = log P(S | C) + log P(C) - log P(S)
  log_posteriors <- log_likelihoods + log(pi)
  
  # Convert back to probability scale using log-sum-exp trick
  posteriors <- exp(log_posteriors - max(log_posteriors))
  posteriors <- posteriors / sum(posteriors)
  
  return(posteriors)
}






#Attempting to apply to CSHQ fitted model (4 groups)

#need to append 0 to any of the item probability matrices of differing dimension

CSHQ_LCR_itemsel$itemprob$P4 <- cbind(CSHQ_LCR_itemsel$itemprob$P4, 0)
CSHQ_LCR_itemsel$itemprob$P6 <- cbind(CSHQ_LCR_itemsel$itemprob$P6, 0)
CSHQ_LCR_itemsel$itemprob$SDB2 <- cbind(CSHQ_LCR_itemsel$itemprob$SDB2, 0)
CSHQ_LCR_itemsel$itemprob$SDB3 <- cbind(CSHQ_LCR_itemsel$itemprob$SDB3, 0)

#Test case
classify_observation(subscale_sums = c(8,2,6,8,6,14,6,16), subscale_items = list(1:4, 5, 6:8, 9:12, 13:15, 16:22, 23:25, 26:33), pi = CSHQ_LCR_varsel$pi, item_probs = CSHQ_LCR_varsel$itemprob, min_response = 1, max_response = 3)


#Function for getting subscale totals from items

#Per observation
CSHQ_compute_subscale_totals_obs <- function(y){
  return(c(sum(y[1:4]), y[5], sum(y[6:8]), sum(y[9:12]), sum(y[13:15]), sum(y[16:22]), sum(y[23:25]), sum(y[26:33])))
}

#Per dataset
CSHQ_compute_subscale_totals_mat <- function(Y){
  n <- nrow(Y)
  result_mat <- matrix(0, nrow = n, ncol = 8)
  for (i in 1:n){
    result_mat[i,] <- CSHQ_compute_subscale_totals_obs(Y[i,])
  }
  return(result_mat)
}

CSHQ_Y_subscale_total_mat <- CSHQ_compute_subscale_totals_mat(Y_CSHQ)


classify_observation_CSHQ <- function(subscale_totals, pi, item_probs){
  return(classify_observation(subscale_sums = subscale_totals, subscale_items = list(1:4, 5, 6:8, 9:12, 13:15, 16:22, 23:25, 26:33), pi = pi, item_probs = item_probs, min_response = 1, max_response = 3))
}

classify_observation_regression_CSHQ <- function(subscale_totals, beta, x, item_probs){
  return(classify_observation_regression(subscale_sums = subscale_totals, subscale_items = list(1:4, 5, 6:8, 9:12, 13:15, 16:22, 23:25, 26:33), beta = beta, x = x, item_probs = item_probs, min_response = 1, max_response = 3))
}

classify_observation_CSHQ(subscale_totals = CSHQ_Y_subscale_total_mat[1,], pi = CSHQ_LCR_varsel$pi, item_probs = CSHQ_LCR_varsel$itemprob)
classify_observation_regression_CSHQ(subscale_totals = CSHQ_Y_subscale_total_mat[1,], beta = CSHQ_LCR_varsel$beta_estimate, x = X_CSHQ[1,], item_probs = CSHQ_LCR_varsel$itemprob)

classify_from_subscale_matrix_CSHQ <- function(subscale_totals_mat, pi, item_probs){
  n <- nrow(subscale_totals_mat)
  G <- length(pi)
  result_matrix <- matrix(0, nrow = n, ncol = G)
  for (i in 1:n){
    result_matrix[i,] <- classify_observation_CSHQ(subscale_totals = subscale_totals_mat[i,], pi = pi, item_probs = item_probs)
  }
  return(result_matrix)
}

classify_from_subscale_matrix_CSHQ_regression <- function(subscale_totals_mat, beta, X, item_probs){
  n <- nrow(subscale_totals_mat)
  G <- ncol(beta)
  result_matrix <- matrix(0, nrow = n, ncol = G)
  for (i in 1:n){
    result_matrix[i,] <- classify_observation_regression_CSHQ(subscale_totals = subscale_totals_mat[i,], beta = beta, x = X[i,], item_probs = item_probs)
  }
  return(result_matrix)
}


#Testing the correspondence between classification by this subscale method vs via modal assignment
subscale_classification_mat_CSHQ <- classify_from_subscale_matrix_CSHQ(subscale_totals_mat = CSHQ_Y_subscale_total_mat, pi = CSHQ_LCR_varsel$pi, item_probs = CSHQ_LCR_varsel$itemprob)

subscale_assignment_CSHQ <- max.col(subscale_classification_mat_CSHQ)

modal_assignment_CSHQ <- max.col(CSHQ_LCR_varsel$Z)

minVI_assignment_4group_CSHQ <- CSHQ_4group_minVI_cluster$cl

table(modal_assignment_CSHQ, subscale_assignment_CSHQ)
table(minVI_assignment_4group_CSHQ, subscale_assignment_CSHQ)







