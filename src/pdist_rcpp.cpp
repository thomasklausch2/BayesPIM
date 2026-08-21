#include <Rcpp.h>
#include <cmath>
#include <limits>
#include "pdist_rcpp.h"
using namespace Rcpp;

namespace {

// Prentice generalized-gamma CDF.  For Q != 0,
//
//   z     = (log(x) - mu) / sigma
//   shape = 1 / Q^2
//   u     = shape * exp(Q * z)
//
// and the CDF is the lower gamma tail for Q > 0 and the upper gamma
// tail for Q < 0.  The Q = 0 limit is lognormal.
double pgengamma_prentice(double x, double mu, double sigma, double Q) {
  if(R_IsNA(x) || R_IsNaN(x) ||
     R_IsNA(mu) || R_IsNaN(mu) ||
     R_IsNA(sigma) || R_IsNaN(sigma) ||
     R_IsNA(Q) || R_IsNaN(Q) ||
     !R_finite(mu) ||
     !R_finite(sigma) || sigma <= 0.0 ||
     !R_finite(Q)) {
    return NA_REAL;
  }

  if(x <= 0.0) {
    return 0.0;
  }
  if(!R_finite(x)) {
    return 1.0;
  }

  const double abs_Q = std::abs(Q);
  const double sqrt_epsilon =
    std::sqrt(std::numeric_limits<double>::epsilon());

  // At and extremely close to Q = 0, the gamma shape 1 / Q^2 is so
  // large that the gamma representation loses floating-point accuracy.
  if(abs_Q < sqrt_epsilon) {
    return R::plnorm(x, mu, sigma, 1, 0);
  }

  const double z = (std::log(x) - mu) / sigma;
  const double shape = 1.0 / (Q * Q);
  const double u = shape * std::exp(Q * z);

  return R::pgamma(u, shape, 1.0, Q > 0.0, 0);
}

} // namespace

// Internal C++ helper: declared in pdist_rcpp.h and called from interval_probs_rcpp.cpp.
// Deliberately not exported to R -- no R-level caller exists.
NumericVector pdist_rcpp(NumericVector q, NumericMatrix par, std::string dist) {
  int n = q.size();
  NumericVector result(n);
  
  if(dist == "exp") {
    for(int i = 0; i < n; ++i) {
      result[i] = R::pexp(q[i], par(i, 0), 1, 0); // q, rate, lower_tail, log_p
    }
  } else if(dist == "gamma") {
    for(int i = 0; i < n; ++i) {
      result[i] = R::pgamma(q[i], par(i, 0), 1.0 / par(i, 1), 1, 0); // q, shape, scale, lower_tail, log_p
    }
  } else if(dist == "weibull") {
    for(int i = 0; i < n; ++i) {
      result[i] = R::pweibull(q[i], par(i, 1), par(i, 0), 1, 0); // q, shape, scale, lower_tail, log_p
    }
  } else if(dist == "loglog") {
    for(int i = 0; i < n; ++i) {
      double q_i = q[i];
      double lambda = par(i, 0);
      double gamma = par(i, 1);

      if(R_IsNA(q_i) || R_IsNaN(q_i) ||
         R_IsNA(lambda) || R_IsNaN(lambda) ||
         R_IsNA(gamma) || R_IsNaN(gamma) ||
         !R_finite(lambda) || lambda <= 0.0 ||
         !R_finite(gamma) || gamma <= 0.0) {
        result[i] = NA_REAL;
      } else if(q_i <= 0.0) {
        result[i] = 0.0;
      } else if(!R_finite(q_i)) {
        result[i] = 1.0;
      } else {
        double z = gamma * (std::log(q_i) - std::log(lambda));
        result[i] = R::plogis(z, 0.0, 1.0, 1, 0);
      }
    }
  } else if(dist == "lognormal") {
    for(int i = 0; i < n; ++i) {
      result[i] = R::plnorm(q[i], par(i, 0), par(i, 1), 1, 0); // q, meanlog, sdlog, lower_tail, log_p
    }
  } else if(dist == "gengamma") {
    if(par.ncol() < 3) {
      stop("The generalized gamma distribution requires mu, sigma, and Q");
    }
    for(int i = 0; i < n; ++i) {
      result[i] = pgengamma_prentice(
        q[i], par(i, 0), par(i, 1), par(i, 2)
      );
    }
  } else {
    stop("Unknown distribution");
  }
  
  return result;
}
