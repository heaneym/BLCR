#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;


NumericVector internal_helper_function(NumericVector x) {
  // Some internal processing
  return x;
}


//Function to calculate the proportions so to estimate the probabilites gamma

NumericVector calculate_proportions(IntegerVector row, int G) {
  // Create a vector of zeros with length G
  NumericVector proportions(G, 0.0);
  
  // Count occurrences of each group
  for(int group = 1; group <= G; group++) {
    proportions[group-1] = sum(row == group) / (double)row.size();
  }
  
  return proportions;
}


//Function to calculate gamma itself


NumericMatrix gamma_update_cpp(NumericVector mu, 
                               arma::mat Sigma, 
                               int n, 
                               int G, 
                               int n_samples) {
  // Generate multivariate normal samples
  arma::mat vec_epsilon_samples = arma::mvnrnd(
    arma::zeros(n*G), 
    arma::kron(arma::eye(n, n), Sigma), 
    n_samples
  );
  
  // Create 3D array to store linear predictors
  arma::cube linear_predictors(n, G, n_samples);
  
  // Convert mu to arma::vec for correct reshaping
  arma::vec mu_arma = as<arma::vec>(mu);
  
  // Reshape and add mu
  for(int sample = 0; sample < n_samples; sample++) {
    arma::mat epsilon_sample = arma::reshape(
      vec_epsilon_samples.row(sample), n, G
    );
    
    // Add mu to epsilon samples
    linear_predictors.slice(sample) = 
      epsilon_sample + arma::reshape(mu_arma, n, G);
  }
  
  // Find max group for each individual and sample
  arma::imat chosen_group(n, n_samples);
  for(int sample = 0; sample < n_samples; sample++) {
    for(int i = 0; i < n; i++) {
      arma::rowvec row = linear_predictors.slice(sample).row(i);
      chosen_group(i, sample) = row.index_max() + 1;
    }
  }
  
  // Calculate proportions
  NumericMatrix gamma_mat(n, G);
  for(int i = 0; i < n; i++) {
    IntegerVector row = wrap(chosen_group.row(i));
    gamma_mat(i, _) = calculate_proportions(row, G);
  }
  
  return gamma_mat;
}
