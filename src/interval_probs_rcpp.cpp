#include <Rcpp.h>
#include "pdist_rcpp.h" 
using namespace Rcpp;

// Per-subject screening-interval probabilities: the probability that the event
// falls in each interval, weighted by the detection probabilities in pobs_vec.
// Returns the normalised probabilities used to sample the latent interval, and
// their per-subject sums, which weight the collapsed prevalence augmentation.
// pdist_rcpp() is declared in pdist_rcpp.h.

// [[Rcpp::export]]
List interval_probs_rcpp(NumericVector pobs_vec,
                        List v_obs,
                        NumericVector v_obs_l,
                        NumericVector v_obs_r,
                        NumericMatrix cur_dist_par_t,
                        std::string dist) {
  
  int n = v_obs.size();
  NumericVector m(n);
  for(int i = 0; i < n; ++i) {
    NumericVector temp = v_obs[i];
    m[i] = temp.size() - 1;
  }
  
  // Calculate 'par'
  int total_rows = sum(m);
  int n_par = cur_dist_par_t.ncol();
  NumericMatrix par(total_rows, n_par);
  int row_idx = 0;
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m[i]; ++j) {
      for(int k = 0; k < n_par; ++k) {
        par(row_idx, k) = cur_dist_par_t(i, k);
      }
      ++row_idx;
    }
  }
  
  // Calculate Fxl, Fxr, and pobs_
  NumericVector Fxl = pdist_rcpp(v_obs_l, par, dist);
  NumericVector Fxr = pdist_rcpp(v_obs_r, par, dist);
  NumericVector pobs_ = (Fxr - Fxl) * pobs_vec;
  
  // Split pobs_
  List pobs_split(n);
  int start_idx = 0;
  for(int i = 0; i < n; ++i) {
    NumericVector temp(m[i]);
    std::copy(pobs_.begin() + start_idx, pobs_.begin() + start_idx + (int)m[i], temp.begin());
    
    pobs_split[i] = temp;
    start_idx += m[i];
  }
  
  // Calculate pobs_norm
  List pobs_norm(n);
  NumericVector interval_sums(n);
  for(int i = 0; i < n; ++i) {
    NumericVector temp = pobs_split[i];
    double sum_temp = sum(temp);
    interval_sums[i] = sum_temp;
    pobs_norm[i] = temp / sum_temp;
  }
  
  return List::create(_["pobs_norm"] = pobs_norm, _["interval_sums"] = interval_sums);
}
