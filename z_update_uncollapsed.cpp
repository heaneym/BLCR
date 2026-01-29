// #include <RcppArmadillo.h>
// using namespace Rcpp;
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List z_update_uncollapsed(arma::mat log_logit_probs, arma::cube log_theta, arma::umat Y, int G, int n, int M, IntegerVector K) {
//   
//   arma::mat sum_log_theta(G, n, arma::fill::zeros);
//   
//   for (int j = 0; j < M; ++j) {
//     for (int i = 0; i < n; ++i) {
//       int y_ij = Y(i, j) - 1;  // Convert to 0-based index
//       for (int g = 0; g < G; ++g) {
//         sum_log_theta(g, i) += log_theta(g, j, y_ij);
//       }
//     }
//   }
//   
//   arma::mat log_w = log_logit_probs + sum_log_theta.t();
//   arma::mat w = exp(log_w);
//   arma::mat z(n, G, arma::fill::zeros);
//   
//   for (int i = 0; i < n; ++i) {
//     arma::vec probs = w.row(i).t();
//     probs /= arma::sum(probs);
//     
//     double u = R::runif(0, 1);
//     double cum_prob = 0.0;
//     
//     for (int g = 0; g < G; ++g) {
//       cum_prob += probs(g);
//       if (u <= cum_prob) {
//         z(i, g) = 1;
//         break;
//       }
//     }
//   }
//   
//   return List::create(Named("w") = w, Named("z") = z);
// }

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List z_update_uncollapsed(arma::mat log_logit_probs, arma::cube log_theta, 
                           arma::mat Y, arma::umat Y_missing_indicator,
                           int G, int n, int M, IntegerVector K) {
  
  arma::mat sum_log_theta(G, n, arma::fill::zeros);
  
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < n; ++i) {
      // Only use observed values
      if (Y_missing_indicator(i, j) == 0) {  // 0 = observed, 1 = missing
        int y_ij = Y(i, j) - 1;  
        for (int g = 0; g < G; ++g) {
          sum_log_theta(g, i) += log_theta(g, j, y_ij);
        }
      }
    }
  }
  
  arma::mat log_w = log_logit_probs + sum_log_theta.t();
  arma::mat w = exp(log_w);
  arma::mat z(n, G, arma::fill::zeros);
  
  for (int i = 0; i < n; ++i) {
    arma::vec probs = w.row(i).t();
    probs /= arma::sum(probs);
    
    double u = R::runif(0, 1);
    double cum_prob = 0.0;
    
    for (int g = 0; g < G; ++g) {
      cum_prob += probs(g);
      if (u <= cum_prob) {
        z(i, g) = 1;
        break;
      }
    }
  }
  
  return List::create(Named("w") = w, Named("z") = z);
}

