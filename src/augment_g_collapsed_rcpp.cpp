#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector augment_g_collapsed_rcpp(NumericVector interval_sums, List v_obs, double kappa, NumericVector prob_g, NumericVector r, LogicalVector g_fixed) {
  int n = v_obs.size();
  LogicalVector cens(n);
  NumericVector m(n);
  NumericVector p0(n, NumericVector::get_na()), p1(n, NumericVector::get_na());
  NumericVector g(n, NumericVector::get_na());
  NumericVector prob_nonprev = 1.0 - prob_g;

  for(int i = 0; i < n; ++i) {
    NumericVector currentVec = v_obs[i];
    cens[i] = std::isinf(currentVec[currentVec.size() - 1]);
    m[i] = currentVec.size();
  }

  // Calculate p0 and p1 for non-fixed g
  for(int i = 0; i < n; ++i) {
    if(!g_fixed[i]) {
      p0[i] = prob_nonprev[i] * interval_sums[i];
      p1[i] = prob_g[i] * std::pow(kappa, 1 - cens[i]) * std::pow(1 - kappa, m[i] - 2 + r[i]);
    }
  }

  // Normalize probabilities
  // for(int i = 0; i < n; ++i) {
  //   if(!g_fixed[i]) {
  //     psum[i] = p0[i] + p1[i];
  //     p1[i] /= psum[i];
  //   }
  // }
  // Normalize probabilities (robust)
  for (int i = 0; i < n; ++i) {
    if (!g_fixed[i]) {
      double weight_sum = p0[i] + p1[i];
      
      if (!R_finite(p0[i]) || p0[i] < 0.0 ||
          !R_finite(p1[i]) || p1[i] < 0.0 ||
          !R_finite(weight_sum) || weight_sum <= 0.0) {
          stop(
            "Collapsed prevalence weights for subject %d must be finite, "
            "nonnegative, and have a positive sum.",
            i + 1
          );
      }
      
      p1[i] /= weight_sum;
    }
  }

  // Sample new g values for non-fixed g
  for(int i = 0; i < n; ++i) {
    if(!g_fixed[i]) {
      g[i] = R::rbinom(1, p1[i]);
    } else {
      g[i] = 1;
    }
  }

  return g;
}


