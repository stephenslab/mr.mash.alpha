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

// Loop to compute mixsqp update
arma::mat compute_mixsqp_update_loop (const arma::mat& X, const arma::mat& Rbar, const arma::mat& V, const arma::cube& S0,
                                      const arma::mat& mu1, const arma::mat& Vinv, const mr_mash_precomputed_quantities& precomp_quants,
                                      bool standardize, const vec& update_order);

// FUNCTION DEFINITIONS
// --------------------

// Loop to compute mixsqp update in Rcpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat compute_mixsqp_update_loop_rcpp (const arma::mat& X, const arma::mat& Rbar, const arma::mat& V, const arma::cube& S0,
                                           const arma::mat& mu1, const arma::mat& Vinv, const List& precomp_quants_list,
                                           bool standardize, const arma::vec& update_order) {
  
  mr_mash_precomputed_quantities precomp_quants
  (as<mat>(precomp_quants_list["S"]),
   as<mat>(precomp_quants_list["S_chol"]),
   as<cube>(precomp_quants_list["S1"]),
   as<cube>(precomp_quants_list["SplusS0_chol"]),
   as<mat>(precomp_quants_list["V_chol"]),
   as<mat>(precomp_quants_list["d"]),
   as<cube>(precomp_quants_list["QtimesV_chol"]),
   as<vec>(precomp_quants_list["xtx"]));
  
  
  return compute_mixsqp_update_loop(X, Rbar, V, S0, mu1, Vinv, precomp_quants, standardize, update_order);
}

// Loop to compute mixsqp update in Rcpp
mat compute_mixsqp_update_loop (const mat& X, const mat& Rbar, const mat& V, const cube& S0,
                                const mat& mu1, const mat& Vinv, const mr_mash_precomputed_quantities& precomp_quants,
                                bool standardize, const vec& update_order) {
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int r = Rbar.n_cols;
  unsigned int k = S0.n_slices;
  vec x(n);
  mat Rbar_j(n,r);
  vec mu1_j(r);
  vec b(r);
  mat S(r,r);
  mat S_chol(r,r);
  vec mu1_ridge(r);
  mat S1_ridge(r,r);
  double xtx_j;
  double logbf;
  mat L(p,k);
  
  // Repeat for each predictor.
  for (unsigned int j : update_order) {
    x = X.col(j);
    mu1_j = trans(mu1.row(j));
    
    // Disregard the ith predictor in the expected residuals.
    Rbar_j = Rbar + (x * trans(mu1_j));
    
    if(standardize){
      // Compute the least-squares estimate.
      b = trans(Rbar_j)*x/(n-1);
    } else {
      xtx_j = precomp_quants.xtx(j);
      // Compute the least-squares estimate.
      b = trans(Rbar_j)*x/xtx_j;
      S = V/xtx_j;
      
      // Compute quantities needed for bayes_mvr_ridge_centered_X()
      S_chol = precomp_quants.V_chol/sqrt(xtx_j);
    }
    
    // Repeat for each mixture component.
    for (unsigned int i = 0; i < k; i++) {
      // Update the posterior quantities for the jth
      // predictor.
      if (standardize){
        logbf = bayes_mvr_ridge_standardized_X(b, S0.slice(i), precomp_quants.S,
                                               precomp_quants.S1.slice(i), precomp_quants.SplusS0_chol.slice(i),
                                               precomp_quants.S_chol, mu1_ridge);
      }else {
        logbf = bayes_mvr_ridge_centered_X(V, b, S, S0.slice(i), xtx_j, Vinv,
                                           precomp_quants.V_chol, S_chol, precomp_quants.d.col(i), 
                                           precomp_quants.QtimesV_chol.slice(i), mu1_ridge, S1_ridge);
      }
      L(j,i) = logbf;
    }
  }
  return L;
}

