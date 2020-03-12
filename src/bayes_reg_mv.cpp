#include "bayes_reg_mv.h"
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Bayesian multivariate simple regression with normal prior with scaled x.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List bayes_mvr_ridge_scaled_X_rcpp (const arma::vec& b, const arma::mat& S0,
                                    const arma::mat& S, const arma::mat& S1,
                                    const arma::mat& SplusS0_chol,
                                    const arma::mat& S_chol) {
  unsigned int r = b.n_elem;
  vec    mu1(r);
  double logbf = bayes_mvr_ridge_scaled_X(b, S0, S, S1, SplusS0_chol, S_chol,
                                          mu1);
  return List::create(Named("mu1")   = mu1,
                      Named("logbf") = logbf);
}


// Bayesian multivariate simple regression with normal prior with centered x.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List bayes_mvr_ridge_centered_X_rcpp (const arma::mat& V, const arma::vec& b, const arma::mat& S,
                                      const arma::mat& S0, double xtx,
                                      const arma::mat& V_chol, const arma::mat& S_chol,
                                      const arma::mat& U0, const arma::vec& d,
                                      const arma::mat& Q) {
  unsigned int r = b.n_elem;
  vec    mu1(r);
  mat    S1(r,r);
  double logbf = bayes_mvr_ridge_centered_X(V, b, S, S0, xtx, V_chol,
                                            S_chol, U0, d, Q, mu1, S1);
  return List::create(Named("mu1")   = mu1,
                      Named("S1")   = S1,
                      Named("logbf") = logbf);
}


// // Bayesian multivariate simple regression with mixture-of-normals prior with scaled x.
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List bayes_mvr_mix_scaled_X_rcpp (const arma::vec& x, const arma::mat& Y,
//                                   const arma::vec& w0, const arma::cube& S0,
//                                   const arma::mat& S, const arma::cube& S1,
//                                   const arma::cube& SplusS0_chol,
//                                   const arma::mat& S_chol) {
//   unsigned int r = Y.n_cols;
//   unsigned int k = w0.n_elem;
//   vec mu1_mix(r);
//   mat S1_mix(r,r);
//   vec w1(k);
//   double logbf_mix = bayes_mvr_mix_scaled_X(x,Y,w0,S0,S,S1,SplusS0_chol,
//                                             S_chol, mu1_mix, S1_mix, w1);
//   return List::create(Named("mu1")   = mu1_mix,
//                       Named("S1")    = S1_mix,
//                       Named("w1")    = w1,
//                       Named("logbf") = logbf_mix);
// }
// 
// 
// // Bayesian multivariate simple regression with mixture-of-normals prior with centered x.
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List bayes_mvr_mix_centered_X_rcpp (const arma::vec& x, const arma::mat& Y,
//                                     const arma::mat& V, const arma::vec& w0,
//                                     const arma::cube& S0, double xtx, 
//                                     const arma::mat& V_chol,
//                                     const arma::cube& U0, const arma::mat& d, 
//                                     const arma::cube& Q) {
//   unsigned int r = Y.n_cols;
//   unsigned int k = w0.n_elem;
//   vec mu1_mix(r);
//   mat S1_mix(r,r);
//   vec w1(k);
//   double logbf_mix = bayes_mvr_mix_centered_X(x, Y, V, w0, S0, xtx, V_chol,
//                                               U0, d, Q, mu1_mix, S1_mix, w1);
//   return List::create(Named("mu1")   = mu1_mix,
//                       Named("S1")    = S1_mix,
//                       Named("w1")    = w1,
//                       Named("logbf") = logbf_mix);
// }


// Perform Bayesian multivariate simple regression with normal prior with scaled x.
double bayes_mvr_ridge_scaled_X (const vec& b, const mat& S0, const mat& S,
                                 const mat& S1, const mat& SplusS0_chol,
                                 const mat& S_chol, vec& mu1) {
  
  // Compute the posterior mean (mu1) aassuming a multivariate 
  // normal prior with zero mean and covariance S0.
  mu1 = S1 * backsolve(S_chol, forwardsolve(S_chol, b));
  
  // Compute the log-Bayes factor. This should give the same result as:
  //
  //   ldmvnorm(bhat,S0 + S) - ldmvnorm(bhat,S)
  //
  return ldmvnormdiff(b,S_chol,SplusS0_chol);
}


// Perform Bayesian multivariate simple regression with normal prior with centered x.
double bayes_mvr_ridge_centered_X (const mat& V, const vec& b, const mat& S, 
                                   const mat& S0, double xtx,
                                   const mat& V_chol, const mat& S_chol,
                                   const mat& U0, const vec& d, const mat& Q,
                                   vec& mu1, mat& S1) {
  
  // Compute the posterior mean (mu1) and covariance (S1) assuming a
  // multivariate normal prior with zero mean and covariance S0.
  mat D  = diagmat(1/(1 + xtx * d));
  mat U1 = U0 * Q * D * trans(Q);
  S1  = trans(V_chol) * U1 * V_chol;
  mu1 = S1 * backsolve(S_chol, forwardsolve(S_chol, b));
  
  // Compute the log-Bayes factor.
  // return ldmvnorm(bhat,S0 + S) - ldmvnorm(bhat,S)
  return ldmvnormdiff(b,S_chol,chol(S + S0, "upper"));
}


// Perform Bayesian multivariate simple regression with mixture-of-normals prior with scaled x.
double bayes_mvr_mix_scaled_X (const vec& x, const mat& Y, const vec& w0,
                               const cube& S0, const mat& S, const cube& S1,
                               const cube& SplusS0_chol, const mat& S_chol,
                               vec& mu1_mix, mat& S1_mix, vec& w1) {
  unsigned int k = w0.n_elem;
  unsigned int r = Y.n_cols;
  unsigned int n = Y.n_rows;
  
  mat mu1mix(r,k);
  vec logbfmix(k);
  vec mu1(r);
  
  // Compute the least-squares estimate.
  vec b = trans(Y)*x/(n-1);
  
  // Compute the quantities separately for each mixture component.
  for (unsigned int i = 0; i < k; i++) {
    logbfmix(i) = bayes_mvr_ridge_scaled_X(b, S0.slice(i), S, S1.slice(i),
             SplusS0_chol.slice(i), S_chol, mu1);
    mu1mix.col(i) = mu1;
  }
  
  // Compute the posterior assignment probabilities for the latent
  // indicator variable.
  logbfmix += log(w0);
  w1 = logbfmix;
  softmax(w1);
  
  // Compute the posterior mean (mu1) and covariance (S1_mix) of the
  // regression coefficients.
  S1_mix.fill(0);
  mu1_mix.fill(0);
  for (unsigned int i = 0; i < k; i++) {
    b    = mu1mix.col(i);
    mu1_mix += w1(i) * b;
    S1_mix  += w1(i) * (S1.slice(i) + b * trans(b));
  }
  S1_mix -= mu1_mix * trans(mu1_mix);
  
  // Compute the log-Bayes factor as a linear combination of the
  // individual Bayes factors for each mixture component.
  double u = max(logbfmix);
  return u + log(sum(exp(logbfmix - u)));
}


// Perform Bayesian multivariate simple regression with mixture-of-normals prior with centered x.
double bayes_mvr_mix_centered_X (const vec& x, const mat& Y, const mat& V,
                                 const vec& w0, const cube& S0, double xtx, 
                                 const mat& V_chol, const cube& U0,
                                 const mat& d, const cube& Q, vec& mu1_mix,
                                 mat& S1_mix, vec& w1) {
  unsigned int k = w0.n_elem;
  unsigned int r = Y.n_cols;
  
  mat  mu1mix(r,k);
  cube S1mix(r,r,k);
  vec  logbfmix(k);
  vec  mu1(r);
  mat  S1(r,r);
  
  // Compute the least-squares estimate.
  vec b = trans(Y)*x/xtx;
  mat S = V/xtx;
  
  // Compute quantities needed for bayes_mvr_ridge_centered_X()
  mat S_chol = V_chol/sqrt(xtx);
  
  // Compute the quantities separately for each mixture component.
  for (unsigned int i = 0; i < k; i++) {
    logbfmix(i) = bayes_mvr_ridge_centered_X(V, b, S, S0.slice(i), xtx, V_chol,
             S_chol, U0.slice(i), d.col(i),
             Q.slice(i), mu1, S1);
    mu1mix.col(i)  = mu1;
    S1mix.slice(i) = S1;
  }
  
  // Compute the posterior assignment probabilities for the latent
  // indicator variable.
  logbfmix += log(w0);
  w1 = logbfmix;
  softmax(w1);
  
  // Compute the posterior mean (mu1) and covariance (S1_mix) of the
  // regression coefficients.
  S1_mix.fill(0);
  mu1_mix.fill(0);
  for (unsigned int i = 0; i < k; i++) {
    b    = mu1mix.col(i);
    mu1_mix += w1(i) * b;
    S1_mix  += w1(i) * (S1mix.slice(i) + b * trans(b));
  }
  S1_mix -= mu1_mix * trans(mu1_mix);
  
  // Compute the log-Bayes factor as a linear combination of the
  // individual Bayes factors for each mixture component.
  double u = max(logbfmix);
  return u + log(sum(exp(logbfmix - u)));
}
