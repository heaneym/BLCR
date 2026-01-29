LCR_inject_missingness <- function(Y, prop_missing = 0.05, max_missing_per_row = 3) {
  n <- nrow(Y)
  M <- ncol(Y)
  Y_missing <- Y
  
  # Calculate target number of missing values
  total_cells <- n * M
  n_missing_target <- round(total_cells * prop_missing)
  
  # Distribute missingness as evenly as possible across rows
  missing_per_row <- rep(0, n)
  
  for (miss in 1:n_missing_target) {
    # Find rows with least missing (under the max)
    eligible_rows <- which(missing_per_row < max_missing_per_row)
    
    if (length(eligible_rows) == 0) break
    
    # Choose row with fewest missing values so far
    min_missing <- min(missing_per_row[eligible_rows])
    candidates <- eligible_rows[missing_per_row[eligible_rows] == min_missing]
    i <- sample(candidates, 1)
    
    # Choose a random column that isn't already missing
    available_cols <- which(!is.na(Y_missing[i, ]))
    if (length(available_cols) == 0) next
    
    j <- sample(available_cols, 1)
    
    # Make it missing
    Y_missing[i, j] <- NA
    missing_per_row[i] <- missing_per_row[i] + 1
  }
  
  return(Y_missing)
}
