#include "bayes_reg_mv.h"
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// CLASS DEFINITIONS
// -----------------
// This class is used to implement the multithreaded computation in 
// bayes_mvr_mix_standardized_X.
struct bayes_mvr_mix_standardized_X_worker : public RcppParallel::Worker {
  const vec&  b;
  const mat&  S;
  const mat&  S_chol;
  const cube& S0;
  const cube& S1;
  const cube& SplusS0_chol;
  vec&        logbfmix;
  mat&        mu1mix;
  
  // This is used to create a bayes_mvr_mix_worker object.
  bayes_mvr_mix_standardized_X_worker (const vec& b, const mat& S,
				       const mat& S_chol, const cube& S0,
				       const cube& S1,
				       const cube& SplusS0_chol,
				       vec& logbfmix, mat& mu1mix) :
    b(b), S(S), S_chol(S_chol), S0(S0), S1(S1), SplusS0_chol(SplusS0_chol),
    logbfmix(logbfmix), mu1mix(mu1mix) { };
  
  // This function performs the Bayesian multivariate regression
  // calculations for a given range of prior mixture components.
  void operator() (std::size_t begin, std::size_t end) {
    unsigned int r = b.n_elem;
    vec mu1(r);
    for (unsigned int i = begin; i < end; i++) {
      logbfmix(i) = bayes_mvr_ridge_standardized_X(b,S0.slice(i),S,S1.slice(i),
  						   SplusS0_chol.slice(i),
						   S_chol,mu1);
      mu1mix.col(i) = mu1;
    }
  }
};

// This class is used to implement the multithreaded computation in 
// bayes_mvr_mix_centered_X.
struct bayes_mvr_mix_centered_X_worker : public RcppParallel::Worker {
  const vec&  b;
  double      xtx;
  const mat&  V;
  const mat&  Vinv;
  const mat&  V_chol;
  const mat&  S;
  const mat&  S_chol;
  const mat&  d;
  const cube& S0;
  const cube& QtimesV_chol;
  vec&        logbfmix;
  mat&        mu1mix;
  cube&       S1mix;
  
  // This is used to create a bayes_mvr_mix_worker object.
  bayes_mvr_mix_centered_X_worker (const vec& b, double xtx, const mat& V,
				   const mat& Vinv, const mat& V_chol,
				   const mat& S, const mat& S_chol,
				   const mat& d, const cube& S0,
				   const cube& QtimesV_chol, vec& logbfmix,
				   mat& mu1mix, cube& S1mix) :
    b(b), xtx(xtx), V(V), Vinv(Vinv), V_chol(V_chol), S(S), S_chol(S_chol),
    d(d), S0(S0), QtimesV_chol(QtimesV_chol), logbfmix(logbfmix),
    mu1mix(mu1mix), S1mix(S1mix) { };
  
  // This function performs the Bayesian multivariate regression
  // calculations for a given range of prior mixture components.
  void operator() (std::size_t begin, std::size_t end) {
    unsigned int r = b.n_elem;
    vec mu1(r);
    mat S1(r,r);
    for (unsigned int i = begin; i < end; i++) {
      logbfmix(i) = bayes_mvr_ridge_centered_X(V,b,S,S0.slice(i),xtx,Vinv,  
					       V_chol,S_chol,d.col(i),
					       QtimesV_chol.slice(i),mu1,S1);
      mu1mix.col(i)  = mu1;
      S1mix.slice(i) = S1;
    }
  }
};

// FUNCTION DEFINITIONS
// --------------------
// // Bayesian multivariate simple regression with normal prior with standardized x.
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List bayes_mvr_ridge_standardized_X_rcpp (const arma::vec& b, const arma::mat& S0,
//                                     const arma::mat& S, const arma::mat& S1,
//                                     const arma::mat& SplusS0_chol,
//                                     const arma::mat& S_chol) {
//   unsigned int r = b.n_elem;
//   vec    mu1(r);
//   double logbf = bayes_mvr_ridge_standardized_X(b, S0, S, S1, SplusS0_chol, S_chol,
//                                           mu1);
//   return List::create(Named("mu1")   = mu1,
//                       Named("logbf") = logbf);
// }
// 
// 
// // Bayesian multivariate simple regression with normal prior with centered x.
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List bayes_mvr_ridge_centered_X_rcpp (const arma::mat& V, const arma::vec& b, const arma::mat& S,
//                                       const arma::mat& S0, double xtx, const arma::mat& Vinv,
//                                       const arma::mat& V_chol, const arma::mat& S_chol,
//                                       const arma::vec& d, const arma::mat& QtimesV_chol) {
//   unsigned int r = b.n_elem;
//   vec    mu1(r);
//   mat    S1(r,r);
//   double logbf = bayes_mvr_ridge_centered_X(V, b, S, S0, xtx, Vinv, V_chol,
//                                             S_chol, d, QtimesV_chol, mu1, S1);
//   return List::create(Named("mu1")   = mu1,
//                       Named("S1")   = S1,
//                       Named("logbf") = logbf);
// }
// 
// 
// // Bayesian multivariate simple regression with mixture-of-normals prior with standardized x.
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List bayes_mvr_mix_standardized_X_rcpp (const arma::vec& x, const arma::mat& Y,
//                                   const arma::vec& w0, const arma::cube& S0,
//                                   const arma::mat& S, const arma::cube& S1,
//                                   const arma::cube& SplusS0_chol,
//                                   const arma::mat& S_chol) {
//   unsigned int r = Y.n_cols;
//   unsigned int k = w0.n_elem;
//   vec mu1_mix(r);
//   mat S1_mix(r,r);
//   vec w1(k);
//   double logbf_mix = bayes_mvr_mix_standardized_X(x,Y,w0,S0,S,S1,SplusS0_chol,
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
//                                     const arma::mat& Vinv, const arma::mat& V_chol,
//                                     const arma::mat& d, const arma::cube& QtimesV_chol) {
//   unsigned int r = Y.n_cols;
//   unsigned int k = w0.n_elem;
//   vec mu1_mix(r);
//   mat S1_mix(r,r);
//   vec w1(k);
//   double logbf_mix = bayes_mvr_mix_centered_X(x, Y, V, w0, S0, xtx, Vinv, V_chol,
//                                               d, QtimesV_chol, mu1_mix, S1_mix, w1);
//   return List::create(Named("mu1")   = mu1_mix,
//                       Named("S1")    = S1_mix,
//                       Named("w1")    = w1,
//                       Named("logbf") = logbf_mix);
// }

// Perform Bayesian multivariate simple regression with normal prior
// with standardized x.
double bayes_mvr_ridge_standardized_X (const vec& b, const mat& S0, const mat& S,
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
                                   const mat& S0, double xtx, const mat& Vinv,
                                   const mat& V_chol, const mat& S_chol,
                                   const vec& d, const mat& QtimesV_chol,
                                   vec& mu1, mat& S1) {
  
  // Compute the posterior mean (mu1) and covariance (S1) assuming a
  // multivariate normal prior with zero mean and covariance S0.
  vec dx = d/(1 + xtx * d);
  mat A = sqrt(dx) % QtimesV_chol.each_col();
  S1  = trans(A) * A;
  mu1 = trans(A) * (A * (Vinv * (xtx*b)));
  
  // Compute the log-Bayes factor.
  // return ldmvnorm(bhat,S0 + S) - ldmvnorm(bhat,S)
  return ldmvnormdiff(b,S_chol,chol(S+S0, "upper"));
}

// Perform Bayesian multivariate simple regression with
// mixture-of-normals prior with standardized x.
double bayes_mvr_mix_standardized_X (const vec& x, const mat& Y, const vec& w0,
				     const cube& S0, const mat& S,
				     const cube& S1, const cube& SplusS0_chol,
				     const mat& S_chol, double eps,
                                     unsigned int nthreads, vec& mu1_mix,
				     mat& S1_mix, vec& w1) {
  unsigned int k = w0.n_elem;
  unsigned int r = Y.n_cols;
  unsigned int n = Y.n_rows;
  
  mat mu1mix(r,k);
  vec logbfmix(k);
  vec mu1(r);
  
  // Compute the least-squares estimate.
  vec b = trans(Y)*x/(n-1);
  
  // Compute the quantities separately for each mixture component.
  if (nthreads > 1) {
    bayes_mvr_mix_standardized_X_worker worker(b,S,S_chol,S0,S1,SplusS0_chol,
					       logbfmix,mu1mix);
    parallelFor(0,k,worker);
  } else {
    for (unsigned int i = 0; i < k; i++) {
      logbfmix(i) = bayes_mvr_ridge_standardized_X(b,S0.slice(i),S,S1.slice(i),
  						   SplusS0_chol.slice(i),
						   S_chol,mu1);
      mu1mix.col(i) = mu1;
    }
  }
  
  // Compute the posterior assignment probabilities for the latent
  // indicator variable.
  logbfmix += log(w0 + eps);
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

// Perform Bayesian multivariate simple regression with
// mixture-of-normals prior with centered x.
double bayes_mvr_mix_centered_X (const vec& x, const mat& Y, const mat& V,
                                 const vec& w0, const cube& S0, double xtx, 
                                 const mat& Vinv, const mat& V_chol,
				 const mat& d, const cube& QtimesV_chol,
				 double eps, unsigned int nthreads,
				 vec& mu1_mix, mat& S1_mix, vec& w1) {
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
  if (nthreads > 1) {
    bayes_mvr_mix_centered_X_worker worker(b,xtx,V,Vinv,V_chol,S,S_chol,d,
					   S0,QtimesV_chol,logbfmix,mu1mix,
					   S1mix);
    parallelFor(0,k,worker);
  } else {
    for (unsigned int i = 0; i < k; i++) {
      logbfmix(i) = bayes_mvr_ridge_centered_X(V,b,S,S0.slice(i),xtx,Vinv,  
					       V_chol,S_chol,d.col(i),
					       QtimesV_chol.slice(i),mu1,S1);
      mu1mix.col(i)  = mu1;
      S1mix.slice(i) = S1;
    }
  }
  
  // Compute the posterior assignment probabilities for the latent
  // indicator variable.
  logbfmix += log(w0 + eps);
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
