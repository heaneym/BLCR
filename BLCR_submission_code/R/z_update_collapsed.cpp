
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <Rmath.h>
using namespace Rcpp;


IntegerVector sample_multinom_rmath(NumericVector probs) {
  int k = probs.size();
  IntegerVector out(k);
  R::rmultinom(1, probs.begin(), k, out.begin());
  return out;
} 


List z_update_collapsed(IntegerMatrix z, 
                        IntegerVector nu,
                        NumericMatrix mu,
                        IntegerVector K,
                        double theta_hyperparam,
                        IntegerVector N_gjk,
                        IntegerVector N_g,
                        NumericVector Y_indicator,
                        IntegerMatrix Y_missing_indicator,  // NEW
                        int n, int G, int M,
                        NumericMatrix omega,
                        NumericMatrix C,
                        int p,
                        NumericMatrix beta,
                        NumericVector beta_cov_inv,
                        NumericVector beta_cov_inv_chol,
                        NumericVector beta_prior_mean,
                        NumericMatrix beta_prior_cov_inv,
                        IntegerVector gamma,
                        NumericMatrix X_current) {
  

  const int maxK = *std::max_element(K.begin(), K.end());
  IntegerVector dimN_gjk = IntegerVector::create(G, M, maxK);
  IntegerVector dimY_ind = IntegerVector::create(n, M, maxK);
  

  std::vector<int> which_item_var;
  for(int j = 0; j < M; ++j) {
    if(nu[j] == 1) which_item_var.push_back(j);
  }
  const int K_current = which_item_var.size();
  

  NumericMatrix w(n, G);
  IntegerVector N_gjk_flat = clone(N_gjk);
  IntegerVector N_g_current = clone(N_g);
  

  IntegerMatrix N_g_obs_j(G, M);
  for(int g = 0; g < G; ++g) {
    for(int j = 0; j < M; ++j) {
      int count = 0;
      for(int i = 0; i < n; ++i) {
        if(z(i, g) == 1 && Y_missing_indicator(i, j) == 0) {
          count++;
        }
      }
      N_g_obs_j(g, j) = count;
    }
  }
  

  for(int i = 0; i < n; ++i) {
    IntegerVector curr_z_i = z(i, _);
    

    IntegerMatrix mask(M, maxK);
    for(int j = 0; j < M; ++j) {
      for(int k = 0; k < maxK; ++k) {
        mask(j, k) = (Y_indicator[i + j*n + k*n*M] > 0.5) && 
          (Y_missing_indicator(i, j) == 0);
      }
    }
    
    for(int g = 0; g < G; ++g) {
      IntegerVector N_g_temp = clone(N_g_current);
      IntegerMatrix N_g_obs_j_temp = clone(N_g_obs_j);
      
      for(int gg = 0; gg < G; ++gg) {
        N_g_temp[gg] -= curr_z_i[gg];
        if(gg == g) N_g_temp[gg] += 1;
      }
      
      for(int j = 0; j < M; ++j) {
        if(Y_missing_indicator(i, j) == 0) {  
          for(int gg = 0; gg < G; ++gg) {
            N_g_obs_j_temp(gg, j) -= curr_z_i[gg];
            if(gg == g) N_g_obs_j_temp(gg, j) += 1;
          }
        }
      }
      
      double term2 = 0.0, term3 = 0.0;
      IntegerVector N_gjk_temp = clone(N_gjk_flat);
      
      for(int j : which_item_var) {
        for(int k = 0; k < K[j]; ++k) {
          if(mask(j, k)) {  
            const int offset = j*G + k*G*M;
            for(int gg = 0; gg < G; ++gg) {
              const int idx = gg + offset;
              int delta = (gg == g) ? (1 - curr_z_i[gg]) 
                : (-curr_z_i[gg]);
              N_gjk_temp[idx] += delta;
              term2 += R::lgammafn(N_gjk_temp[idx] + theta_hyperparam);
            }
          }
        }
        
        for(int gg = 0; gg < G; ++gg) {
          term3 += R::lgammafn(N_g_obs_j_temp(gg, j) + K[j] * theta_hyperparam);
        }
      }
      
      w(i, g) = mu(i, g) + term2 - term3;
    }
    
    NumericVector log_probs = w(i, _);
    double max_log = *std::max_element(log_probs.begin(), log_probs.end());
    NumericVector probs = exp(log_probs - max_log);
    probs = probs / sum(probs);
    
    int k = probs.size();
    IntegerVector new_z(k);
    R::rmultinom(1, probs.begin(), k, new_z.begin());
    
    z(i, _) = new_z;
    for(int gg = 0; gg < G; ++gg) {
      N_g_current[gg] += (new_z[gg] - curr_z_i[gg]);
    }
    
    for(int j = 0; j < M; ++j) {
      if(Y_missing_indicator(i, j) == 0) {  
        for(int gg = 0; gg < G; ++gg) {
          N_g_obs_j(gg, j) += (new_z[gg] - curr_z_i[gg]);
        }
      }
    }
    
    for(int j : which_item_var) {
      for(int k = 0; k < K[j]; ++k) {
        if(mask(j, k)) {
          const int offset = j*G + k*G*M;
          for(int gg = 0; gg < G; ++gg) {
            const int idx = gg + offset;
            N_gjk_flat[idx] += (new_z[gg] - curr_z_i[gg]);
          }
        }
      }
    }
  }
  
  NumericMatrix kappa(n, G);
  for(int i = 0; i < n; ++i) {
    for(int g = 0; g < G; ++g) {
      kappa(i, g) = z(i, g) - 0.5;
    }
  }
  
  N_gjk_flat.attr("dim") = dimN_gjk;
  return List::create(
    Named("z") = z,
    Named("w") = w,
    Named("kappa") = kappa,
    Named("N_g") = N_g_current,
    Named("N_gjk") = N_gjk_flat
  );
}





