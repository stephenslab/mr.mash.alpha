// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include "bayes_reg_mv.h"

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
  const cube U0;
  const mat  d;
  const cube Q;
  const vec  xtx;

  // This is used to create a mr_mash_precomputed_quantities object.
  mr_mash_precomputed_quantities (const mat& S, const mat& S_chol,
				  const cube& S1, const cube& SplusS0_chol,
				  const mat& V_chol, const cube& U0,
				  const mat& d, const cube& Q,
				  const vec& xtx) :
    S(S), S_chol(S_chol), S1(S1), SplusS0_chol(SplusS0_chol),
    V_chol(V_chol), U0(U0), d(d), Q(Q), xtx(xtx) { };
};


// FUNCTION DECLARATIONS
// ---------------------
void inner_loop_general (const mat& X, mat& Rbar, mat& mu1, const mat& V,
			 const vec& w0, const cube& S0,
			 const mr_mash_precomputed_quantities& precomp_quants,
			 bool standardize, cube& S1, mat& w1);


// FUNCTION DEFINITIONS
// --------------------

//Inner loop
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List inner_loop_general_rcpp (const arma::mat& X, arma::mat& Rbar, mat& mu1,
                              const arma::mat& V, const arma::vec& w0,
                              const arma::cube& S0,
                              const List& precomp_quants_list,
                              bool standardize) {
  unsigned int r = Rbar.n_cols;
  unsigned int p = X.n_cols;
  unsigned int k = w0.n_elem;
  cube S1(r,r,p);
  mat  w1(p,k);
  mat  mu1_new  = mu1;
  mat  Rbar_new = Rbar;
  mr_mash_precomputed_quantities precomp_quants
    (as<mat>(precomp_quants_list["S"]),
     as<mat>(precomp_quants_list["S_chol"]),
     as<cube>(precomp_quants_list["S1"]),
     as<cube>(precomp_quants_list["SplusS0_chol"]),
     as<mat>(precomp_quants_list["V_chol"]),
     as<cube>(precomp_quants_list["U0"]),
     as<mat>(precomp_quants_list["d"]),
     as<cube>(precomp_quants_list["Q"]),
     as<vec>(precomp_quants_list["xtx"]));
  inner_loop_general(X, Rbar_new, mu1_new, V, w0, S0, precomp_quants,
                     standardize, S1, w1);
  return List::create(Named("rbar") = Rbar_new,
                      Named("mu1")  = mu1_new,
                      Named("S1")   = S1,
                      Named("w1")   = w1);
}


// Perform the inner loop
void inner_loop_general (const mat& X, mat& Rbar, mat& mu1, const mat& V,
                         const vec& w0, const cube& S0,
                         const mr_mash_precomputed_quantities& precomp_quants,
                         bool standardize, cube& S1, mat& w1) {
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int r = Rbar.n_cols;
  unsigned int k = w0.n_elem;
  vec x(n);
  mat Rbar_j(n,r);
  vec mu1_j(r);
  vec mu1_mix(r);
  mat S1_mix(r,r);
  vec w1_mix(k);
  
  // Repeat for each predictor.
  for (unsigned int j = 0; j < p; j++) {
    x = X.col(j);
    mu1_j = trans(mu1.row(j));
    
    // Disregard the ith predictor in the expected residuals.
    Rbar_j = Rbar + (x * trans(mu1_j));
    
    // Update the posterior quantities for the jth
    // predictor.
    if (standardize)
      bayes_mvr_mix_scaled_X(x, Rbar_j, w0, S0, precomp_quants.S,
                             precomp_quants.S1, precomp_quants.SplusS0_chol,
                             precomp_quants.S_chol, mu1_mix, S1_mix, w1_mix);
    else {
      double xtx_j = precomp_quants.xtx(j);
      bayes_mvr_mix_centered_X(x, Rbar_j, V, w0, S0, xtx_j,
                               precomp_quants.V_chol, precomp_quants.U0,
                               precomp_quants.d, precomp_quants.Q,
                               mu1_mix, S1_mix, w1_mix);
    }
    
    mu1.row(j)  = trans(mu1_mix);
    S1.slice(j) = S1_mix;
    w1.row(j)   = trans(w1_mix);
    
    // Update the expected residuals.
    Rbar = Rbar_j - (x * trans(mu1_mix));
  }
}
