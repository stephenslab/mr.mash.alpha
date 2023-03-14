#include "bayes_reg_mv.h"
#include "misc.h"
#include <algorithm>
#include <cmath>

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

// Inner loop
void inner_loop_general_rss (unsigned int n, const mat& XtY, mat& XtXmu1, mat& mu1, const mat& V,
                             const mat& Vinv, const vec& w0, const cube& S0,
                             const mr_mash_precomputed_quantities& precomp_quants,
                             bool standardize, bool compute_ELBO, bool update_V,
                             const vec& update_order, double eps, unsigned int nthreads,
                             cube& S1, mat& w1, double& var_part_tr_wERSS, 
                             double& neg_KL, mat& var_part_ERSS);


// FUNCTION DEFINITIONS
// --------------------

// Inner loop
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
List inner_loop_general_rss_rcpp (unsigned int n, const arma::mat& XtY, arma::mat& XtXmu1, arma::mat& mu1,
                                  const arma::mat& V, const arma::mat& Vinv, const arma::vec& w0,
                                  const arma::cube& S0, const List& precomp_quants_list,
                                  bool standardize, bool compute_ELBO, bool update_V,
                                  const arma::vec& update_order, double eps, unsigned int nthreads) {
  unsigned int r = mu1.n_cols;
  unsigned int p = mu1.n_rows;
  unsigned int k = w0.n_elem;
  cube S1(r,r,p);
  mat  w1(p,k);
  mat  mu1_new  = mu1;
  mat  XtXmu1_new = XtXmu1;
  double var_part_tr_wERSS;
  double neg_KL;
  mat var_part_ERSS(r,r);
  mr_mash_precomputed_quantities precomp_quants
    (as<mat>(precomp_quants_list["S"]),
     as<mat>(precomp_quants_list["S_chol"]),
     as<cube>(precomp_quants_list["S1"]),
     as<cube>(precomp_quants_list["SplusS0_chol"]),
     as<mat>(precomp_quants_list["V_chol"]),
     as<mat>(precomp_quants_list["d"]),
     as<cube>(precomp_quants_list["QtimesV_chol"]),
     as<vec>(precomp_quants_list["xtx"]));
  inner_loop_general_rss(n, XtY, XtXmu1_new, mu1_new, V, Vinv, w0, S0, precomp_quants,
                         standardize, compute_ELBO, update_V, update_order, eps,
		                     nthreads, S1, w1, var_part_tr_wERSS, neg_KL, var_part_ERSS);
  return List::create(Named("mu1")                = mu1_new,
                      Named("S1")                 = S1,
                      Named("w1")                 = w1,
                      Named("var_part_tr_wERSS")  = var_part_tr_wERSS,
                      Named("neg_KL")             = neg_KL,
                      Named("var_part_ERSS")      = var_part_ERSS);
}

// Perform the inner loop
void inner_loop_general_rss (unsigned int n, mat& XtY, mat& XtXmu1, mat& mu1, const mat& V,
                             const mat& Vinv, const vec& w0, const cube& S0,
                             const mr_mash_precomputed_quantities& precomp_quants,
                             bool standardize, bool compute_ELBO, bool update_V, 
                             const vec& update_order, double eps, unsigned int nthreads,
                             cube& S1, mat& w1, double& var_part_tr_wERSS, 
                             double& neg_KL, mat& var_part_ERSS) {
  unsigned int p = mu1.n_rows;
  unsigned int r = mu1.n_cols;
  unsigned int k = w0.n_elem;
  vec x(n);
  mat XtRbar_j(p,r);
  vec mu1_j(r);
  vec mu1_mix(r);
  mat S1_mix(r,r);
  vec w1_mix(k);
  double logbf_mix;
  double xtx;
  
  // Initialize ELBO parameters
  var_part_tr_wERSS = 0;
  neg_KL = 0;
  
  // Initialize V parameters
  var_part_ERSS.zeros(r,r);
  
  // Repeat for each predictor.
  for (unsigned int j : update_order) {
    
    if (standardize)
      double xtx_j = n-1;
    else
      double xtx_j = precomp_quants(j);
    
    mu1_j = trans(mu1.row(j));
    XtY_j = trans(XtY.row(j));
    XtXmu1_j = trans(XtXmu1.row(j));
    
    // Disregard the ith predictor in the expected residuals.
    XtRbar_j = XtY_j - XtXmu1_j + xtx * mu1_j;
    
    // Update the posterior quantities for the jth
    // predictor.
    if (standardize)
      logbf_mix = bayes_mvr_mix_standardized_X_rss(n, XtRbar_j, w0, S0, precomp_quants.S,
					                                         precomp_quants.S1,
					                                         precomp_quants.SplusS0_chol,
					                                         precomp_quants.S_chol, eps, nthreads,
					                                         mu1_mix, S1_mix, w1_mix);
    else
      logbf_mix = bayes_mvr_mix_centered_X_rss(XtRbar_j, V, w0, S0, xtx_j, Vinv,
					                                     precomp_quants.V_chol, precomp_quants.d, 
					                                     precomp_quants.QtimesV_chol, eps, nthreads,
					                                     mu1_mix, S1_mix, w1_mix);
    
    mu1.row(j)  = trans(mu1_mix);
    S1.slice(j) = S1_mix;
    w1.row(j)   = trans(w1_mix);
    
    // Compute ELBO parameters
    if (compute_ELBO)
      compute_ELBO_rss_terms(var_part_tr_wERSS, neg_KL, x, Rbar_j, logbf_mix, mu1_mix, S1_mix, xtx_j, Vinv);
    
    // Compute V parameters
    if (update_V)
      compute_var_part_ERSS(var_part_ERSS, S1_mix, xtx_j);
  }
}
