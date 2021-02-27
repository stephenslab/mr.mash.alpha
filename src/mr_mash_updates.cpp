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
void inner_loop_general (const mat& X, mat& Rbar, mat& mu1, const mat& V,
                         const mat& Vinv, const vec& w0, const cube& S0,
                         const mr_mash_precomputed_quantities& precomp_quants,
                         bool standardize, bool compute_ELBO, bool update_V,
                         const vec& update_order, double eps, unsigned int nthreads,
                         cube& S1, mat& w1, double& var_part_tr_wERSS, 
                         double& neg_KL, mat& var_part_ERSS);


// Missing Y imputation
void impute_missing_Y (mat& Y, const mat& mu, const mat& Vinv, 
                       const cube& miss, const cube& non_miss,
                       mat& Y_cov, double& sum_neg_ent_Y_miss);


// FUNCTION DEFINITIONS
// --------------------

// Inner loop
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
List inner_loop_general_rcpp (const arma::mat& X, arma::mat& Rbar, arma::mat& mu1,
                              const arma::mat& V, const arma::mat& Vinv, const arma::vec& w0,
                              const arma::cube& S0, const List& precomp_quants_list,
                              bool standardize, bool compute_ELBO, bool update_V,
                              const arma::vec& update_order, double eps, unsigned int nthreads) {
  unsigned int r = Rbar.n_cols;
  unsigned int p = X.n_cols;
  unsigned int k = w0.n_elem;
  cube S1(r,r,p);
  mat  w1(p,k);
  mat  mu1_new  = mu1;
  mat  Rbar_new = Rbar;
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
  inner_loop_general(X, Rbar_new, mu1_new, V, Vinv, w0, S0, precomp_quants,
                     standardize, compute_ELBO, update_V, update_order, eps,
		     nthreads, S1, w1, var_part_tr_wERSS, neg_KL, var_part_ERSS);
  return List::create(Named("Rbar")               = Rbar_new,
                      Named("mu1")                = mu1_new,
                      Named("S1")                 = S1,
                      Named("w1")                 = w1,
                      Named("var_part_tr_wERSS")  = var_part_tr_wERSS,
                      Named("neg_KL")             = neg_KL,
                      Named("var_part_ERSS")      = var_part_ERSS);
}

// Perform the inner loop
void inner_loop_general (const mat& X, mat& Rbar, mat& mu1, const mat& V,
                         const mat& Vinv, const vec& w0, const cube& S0,
                         const mr_mash_precomputed_quantities& precomp_quants,
                         bool standardize, bool compute_ELBO, bool update_V, 
                         const vec& update_order, double eps, unsigned int nthreads,
                         cube& S1, mat& w1, double& var_part_tr_wERSS, 
                         double& neg_KL, mat& var_part_ERSS) {
  unsigned int n = X.n_rows;
  // unsigned int p = X.n_cols;
  unsigned int r = Rbar.n_cols;
  unsigned int k = w0.n_elem;
  vec x(n);
  mat Rbar_j(n,r);
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
    x = X.col(j);
    mu1_j = trans(mu1.row(j));
    
    // Disregard the ith predictor in the expected residuals.
    Rbar_j = Rbar + (x * trans(mu1_j));
    
    // Update the posterior quantities for the jth
    // predictor.
    if (standardize)
      logbf_mix = bayes_mvr_mix_standardized_X(x, Rbar_j, w0, S0, precomp_quants.S,
					       precomp_quants.S1,
					       precomp_quants.SplusS0_chol,
					       precomp_quants.S_chol, eps, nthreads,
					       mu1_mix, S1_mix, w1_mix);
    else {
      double xtx_j = precomp_quants.xtx(j);
      logbf_mix = bayes_mvr_mix_centered_X(x, Rbar_j, V, w0, S0, xtx_j, Vinv,
					   precomp_quants.V_chol, precomp_quants.d, 
					   precomp_quants.QtimesV_chol, eps, nthreads,
					   mu1_mix, S1_mix, w1_mix);
    }
    
    mu1.row(j)  = trans(mu1_mix);
    S1.slice(j) = S1_mix;
    w1.row(j)   = trans(w1_mix);
    
    // Compute ELBO parameters
    if (compute_ELBO){
      if (standardize){
        xtx = n-1;
      } else {
        xtx = precomp_quants.xtx(j);
      }
      compute_ELBO_terms(var_part_tr_wERSS, neg_KL, x, Rbar_j, logbf_mix, mu1_mix, S1_mix, xtx, Vinv);
    }
    
    // Compute V parameters
    if (update_V){
      if (standardize){
        xtx = n-1;
      } else {
        xtx = precomp_quants.xtx(j);
      }
      compute_var_part_ERSS(var_part_ERSS, S1_mix, xtx);
    }
    
    // Update the expected residuals.
    Rbar = Rbar_j - (x * trans(mu1_mix));
  }
}






// Impute missing Y
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List impute_missing_Y_rcpp (arma::mat& Y, const arma::mat& mu, const arma::mat& Vinv, 
                            const arma::cube& miss, const arma::cube& non_miss) {
  unsigned int r = Y.n_cols;
  mat Y_cov(r,r);
  double sum_neg_ent_Y_miss;

  impute_missing_Y(Y, mu, Vinv, miss, non_miss, Y_cov, sum_neg_ent_Y_miss);
  return List::create(Named("Y")                  = Y,
                      Named("Y_cov")              = Y_cov,
                      Named("sum_neg_ent_Y_miss") = sum_neg_ent_Y_miss);
}


// Perform missing Y imputation
void impute_missing_Y (mat& Y, const mat& mu, const mat& Vinv, 
                       const cube& miss, const cube& non_miss,
                       mat& Y_cov, double& sum_neg_ent_Y_miss){
  
  unsigned int n = Y.n_rows;
  unsigned int r = Y.n_cols;
  
  Y_cov.zeros();
  sum_neg_ent_Y_miss = 0;
  
  for (unsigned int i = 0; i < n; i++) {
    // unsigned int z = non_miss.slice(i).n_elem
    // unsigned int x = miss.slice(i).n_elem
    
    // vec non_miss_i(z);
    // vec miss_i(x);
    // mat Vinv_mo(x,z);
    // mat Vinv_mm(x,x);
    
    vec non_miss_i = non_miss.slice(i);
    vec miss_i = miss.slice(i);
    mat Vinv_mo = Vinv(miss_i, non_miss_i);
    mat Vinv_mm = Vinv(miss_i, miss_i);
    
    if(std::any_of(miss_i.std::begin(), miss_i.std::end(), [](bool v) { return v; })){
      mat Y_cov_i.zeros(r,r);
      mat Vinv_mm_chol = chol(Vinv_mm);
      mat Y_cov_mm = inv_sympd(Vinv_mm);
      mat Y_cov_i(miss_i, miss_i) = Y_cov_mm;
      
      Y_cov += Y_cov_i;
      
      Y(i, miss_i) = mu(i, miss_i) - Y_cov_mm * Vinv_mo * (Y(i, non_miss_i) - mu(i, non_miss_i));
      
      sum_neg_ent_Y_miss += (0.5 * (Vinv_mm_chol.n_cols * log((1/(2*M_PI*exp(1)))) + chol2ldet(Vinv_mm_chol)));
    }
  }
}
