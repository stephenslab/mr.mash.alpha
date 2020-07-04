// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include "bayes_reg_mv.h"
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// TYPE DEFINITIONS
// ----------------

// A list of precomputed quantities that are invariant to any updates
// to the mr-mash model parameters.
struct mr_mash_precomputed_quantities {
  const mat  S;
  const mat  S_chol;
  const cube S1;
  const cube SplusS0_chol;
  const mat  V_chol;
  const mat  d;
  const cube QtimesV_chol;
  const vec  xtx;
  
  // This is used to create a mr_mash_precomputed_quantities object.
  mr_mash_precomputed_quantities (const mat& S, const mat& S_chol,
                                  const cube& S1, const cube& SplusS0_chol,
                                  const mat& V_chol, const mat& d, const cube& QtimesV_chol,
                                  const vec& xtx) :
    S(S), S_chol(S_chol), S1(S1), SplusS0_chol(SplusS0_chol),
    V_chol(V_chol), d(d), QtimesV_chol(QtimesV_chol), xtx(xtx) { };
};


// FUNCTION DECLARATIONS
// ---------------------
// Compute logbf for each variable
vec compute_logbf (const mat& X, const mat& Y, const mat& V,
                   const mat& Vinv, const vec& w0, const cube& S0,
                   const mr_mash_precomputed_quantities& precomp_quants,
                   bool standardize, double eps, unsigned int nthreads);

// FUNCTION DEFINITIONS
// --------------------
// Compute logbf for each variable
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec compute_logbf_rcpp (const arma::mat& X, const arma::mat& Y, const arma::mat& V,
                              const arma::mat& Vinv, const arma::vec& w0, const arma::cube& S0,
                              const List& precomp_quants_list,
                              bool standardize, double eps, unsigned int nthreads) {
  
  mr_mash_precomputed_quantities precomp_quants
  (as<mat>(precomp_quants_list["S"]),
   as<mat>(precomp_quants_list["S_chol"]),
   as<cube>(precomp_quants_list["S1"]),
   as<cube>(precomp_quants_list["SplusS0_chol"]),
   as<mat>(precomp_quants_list["V_chol"]),
   as<mat>(precomp_quants_list["d"]),
   as<cube>(precomp_quants_list["QtimesV_chol"]),
   as<vec>(precomp_quants_list["xtx"]));
  
  
  return compute_logbf(X, Y, V, Vinv, w0, S0, precomp_quants, standardize, eps, nthreads);
}

vec compute_logbf (const mat& X, const mat& Y, const mat& V,
                   const mat& Vinv, const vec& w0, const cube& S0,
                   const mr_mash_precomputed_quantities& precomp_quants,
                   bool standardize, double eps, unsigned int nthreads){
  
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int r = Y.n_cols;
  unsigned int k = S0.n_slices;
  vec x(n);
  vec mu1_mix(r);
  mat S1_mix(r,r);
  vec w1_mix(k);
  double logbf_mix;
  vec logbf_mix_all(p);
  
  // Repeat for each predictor.
  for (unsigned int j = 0; j < p; j++) {
    x = X.col(j);
    // Update the posterior quantities for the jth
    // predictor.
    if (standardize)
      logbf_mix = bayes_mvr_mix_standardized_X(x, Y, w0, S0, precomp_quants.S,
                                               precomp_quants.S1, precomp_quants.SplusS0_chol,
                                               precomp_quants.S_chol, eps, nthreads,
                                               mu1_mix, S1_mix, w1_mix);
    else {
      double xtx_j = precomp_quants.xtx(j);
      logbf_mix = bayes_mvr_mix_centered_X(x, Y, V, w0, S0, xtx_j, Vinv,
                                           precomp_quants.V_chol, precomp_quants.d, 
                                           precomp_quants.QtimesV_chol, eps, nthreads,
                                           mu1_mix, S1_mix, w1_mix);
    }
    
    logbf_mix_all(j) = logbf_mix;
  }
  
  return logbf_mix_all;
}
