#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix look_up_mat_rcpp(List v_obs, IntegerVector interval_indices) {
  int n = interval_indices.size();

  if (v_obs.size() != n) {
    stop("`v_obs` and `interval_indices` must have the same length.");
  }

  NumericMatrix b(n, 2);
  
  for(int i = 0; i < n; ++i) {
    NumericVector subject_times = as<NumericVector>(v_obs[i]);
    int interval_index = interval_indices[i];

    int max_index = subject_times.size() - 1;
    if (IntegerVector::is_na(interval_index) || interval_index < 1 || interval_index > max_index) {
      stop("`interval_indices[%d]` must be between 1 and length(v_obs[[%d]]) - 1 (%d).",
           i + 1, i + 1, max_index);
    }

    b(i, 0) = subject_times[interval_index - 1];
    b(i, 1) = subject_times[interval_index];
  }
  
  return b;
}
