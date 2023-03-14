#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------

// Compute the softmax of x, and return the result in x. Guard against
// underflow or overflow by adjusting the entries of x so that the
// largest value is zero.
void softmax (vec& x) {
  x -= max(x);
  x  = exp(x);
  x /= sum(x);
}

// Compute the log-probability density of the multivariate normal
// distribution with zero mean and covariance matrix S, omitting terms
// that do not depend on x or S.
double ldmvnorm (const vec& x, const mat& S) {
  mat    L = chol(S,"lower");
  double d = norm(solve(L,x),2);
  return -(d*d)/2 - sum(log(L.diag()));
}

// Compute the difference of two multivariate normal log-densities,
//
//   ldmvnorm(x,S0 + S) - ldmvnorm(x,S)
//
// where S_chol is the right-hand Cholesky factor of S, and
// SplusS0_chol is the right-hand Cholesky factor of S + S0 (also an
// upper triangular matrix).
double ldmvnormdiff (const vec& x, const mat& S_chol,
                     const mat& SplusS0_chol) {
  return (chol2ldet(S_chol) -
          chol2ldet(SplusS0_chol) +
          dot(x, backsolve(S_chol, forwardsolve(S_chol, x))) - 
          dot(x, backsolve(SplusS0_chol, forwardsolve(SplusS0_chol, x))))/2;
}

// Compute the log determinant from Cholesky decomposed matrix
double chol2ldet (const mat& R) {
  return 2*sum(log(R.diag()));
}

// Function to compute some terms of the ELBO
void compute_ELBO_terms (double& var_part_tr_wERSS, double& neg_KL, const vec& x_j,
                         const mat& rbar_j, double logbf, const mat& mu1, const mat& S1, 
                         double xtx, const mat& Vinv){
  
  var_part_tr_wERSS += (as_scalar(accu(Vinv % S1))*xtx);
  
  neg_KL += (logbf + as_scalar(0.5*(-2*accu((Vinv*trans(rbar_j)) % trans(x_j*trans(mu1)))+
    accu(Vinv % (S1+mu1*trans(mu1)))*xtx)));
}

// // [[Rcpp::plugins("cpp11")]]
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List compute_ELBO_terms_rcpp (double var_part_tr_wERSS_init, double neg_KL_init, double x_j,
//                               const mat& rbar_j, double logbf, const mat& mu1, const mat& S1, 
//                               double xtx, const mat& Vinv) {
//   
//   double var_part_tr_wERSS = var_part_tr_wERSS_init;
//   double neg_KL = neg_KL_init;
//   
//   compute_ELBO_terms(var_part_tr_wERSS, neg_KL, x_j, rbar_j, logbf, mu1, S1, xtx, Vinv);
//   
//   return List::create(Named("var_part_tr_wERSS") = var_part_tr_wERSS,
//                       Named("neg_KL")            = neg_KL);
// }

// Function to compute the variance part of the ERSS
void compute_var_part_ERSS (mat& var_part_ERSS, const mat& S1, double xtx){
  
  var_part_ERSS += (S1*xtx);
}

// // [[Rcpp::plugins("cpp11")]]
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// mat compute_var_part_ERSS_rcpp (mat& var_part_ERSS_init, const mat& S1, double xtx) {
//   
//   mat var_part_ERSS = var_part_ERSS_init;
//   
//   compute_var_part_ERSS(var_part_ERSS, S1, xtx);
//   
//   return var_part_ERSS;
// }

// Function to compute some terms of the ELBO with summary data
void compute_ELBO_rss_terms (double& var_part_tr_wERSS, double& neg_KL, const mat& XtRbar_j,
                             double logbf, const mat& mu1, const mat& S1, 
                             double xtx, const mat& Vinv){
  
  var_part_tr_wERSS += (as_scalar(accu(Vinv % S1))*xtx);
  
  neg_KL += (logbf + as_scalar(0.5*accu(Vinv % (-mu1 * XtRbar_j - XtRbar_j * trans(mu1) + mu1 * trans(mu1) * xtx + S1*xtx))));
}

// // [[Rcpp::plugins("cpp11")]]
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List compute_ELBO_terms_rcpp (double var_part_tr_wERSS_init, double neg_KL_init, double x_j,
//                               const mat& rbar_j, double logbf, const mat& mu1, const mat& S1, 
//                               double xtx, const mat& Vinv) {
//   
//   double var_part_tr_wERSS = var_part_tr_wERSS_init;
//   double neg_KL = neg_KL_init;
//   
//   compute_ELBO_terms(var_part_tr_wERSS, neg_KL, x_j, rbar_j, logbf, mu1, S1, xtx, Vinv);
//   
//   return List::create(Named("var_part_tr_wERSS") = var_part_tr_wERSS,
//                       Named("neg_KL")            = neg_KL);
// }
