#include <cmath>
#include "RcppArmadillo.h"

using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double dmvnorm (const vec & x, const arma::vec & mu, const mat & S);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double dmvnorm_rcpp (const arma::vec & x, const arma::vec & mu,
		     const arma::mat & S) {
  return dmvnorm(x,mu,S);
}

// Compute the log-probability density of the multivariate normal
// distribution with mean mu and covariance matrix S.
double dmvnorm (const vec & x, const arma::vec & mu, const mat & S) {
  double n = (double) x.n_elem;
  mat    L = chol(S,"lower");
  double d = norm(solve(L,x - mu),2);
  return -(n*log(2*M_PI) + d*d)/2 - sum(log(L.diag()));
}
