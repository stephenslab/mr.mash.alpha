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
                       const mat& miss, const mat& non_miss,
                       mat& Y_cov, double& sum_neg_ent_Y_miss);

// Inner loop rss
void inner_loop_general_rss (unsigned int n, const mat& XtX, const mat& XtY, mat& XtRbar, 
                             mat& mu1, const mat& V, const mat& Vinv, const vec& w0, 
                             const cube& S0, const mr_mash_precomputed_quantities& precomp_quants,
                             bool standardize, bool compute_ELBO, bool update_V,
                             const vec& update_order, double eps, unsigned int nthreads,
                             cube& S1, mat& w1, double& var_part_tr_wERSS, 
                             double& neg_KL, mat& var_part_ERSS);

// Inner loop rss sprse
void inner_loop_general_rss_sparse (unsigned int n, const sp_mat& XtX, const mat& XtY, mat& XtRbar, 
                                    mat& mu1, const mat& V, const mat& Vinv, const vec& w0, 
                                    const cube& S0, const mr_mash_precomputed_quantities& precomp_quants,
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
    
    // Disregard the jth predictor in the expected residuals.
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
                            const arma::mat& miss, const arma::mat& non_miss) {
  unsigned int r = Y.n_cols;
  mat Y_cov(r,r, fill::zeros);
  double sum_neg_ent_Y_miss = 0;

  impute_missing_Y(Y, mu, Vinv, miss, non_miss, Y_cov, sum_neg_ent_Y_miss);
  return List::create(Named("Y")                  = Y,
                      Named("Y_cov")              = Y_cov,
                      Named("sum_neg_ent_Y_miss") = sum_neg_ent_Y_miss);
}


// Perform missing Y imputation
void impute_missing_Y (mat& Y, const mat& mu, const mat& Vinv, 
                       const mat& miss, const mat& non_miss,
                       mat& Y_cov, double& sum_neg_ent_Y_miss){
  
  unsigned int n = Y.n_rows;
  unsigned int r = Y.n_cols;
  
  for (unsigned int i = 0; i < n; i++) {
    vec non_miss_i = non_miss.col(i);
    vec miss_i = miss.col(i);
    uvec non_miss_i_idx = find(non_miss_i);
    uvec miss_i_idx = find(miss_i);
    
    mat Vinv_mo = Vinv.submat(miss_i_idx, non_miss_i_idx);
    mat Vinv_mm = Vinv.submat(miss_i_idx, miss_i_idx);
    
    if(std::any_of(miss_i.begin(), miss_i.end(), [](bool v) { return v; })){
      // Compute covariance
      mat Y_cov_i(r,r,fill::zeros);
      mat Vinv_mm_chol = chol(Vinv_mm);
      unsigned int s = Vinv_mm_chol.n_rows;
      // mat Y_cov_mm = inv_sympd(Vinv_mm);
      mat Y_cov_mm = solve(trimatu(Vinv_mm_chol), solve(trimatl(trans(Vinv_mm_chol)),eye(s,s)));
      Y_cov_i.submat(miss_i_idx, miss_i_idx) = Y_cov_mm;
      
      Y_cov += Y_cov_i;
      
      // Compute mean
      rowvec Y_i = Y.row(i);
      rowvec mu_i_miss = mu.row(i);
      mu_i_miss = mu_i_miss.elem(miss_i_idx);
      rowvec mu_i_non_miss = mu.row(i);
      mu_i_non_miss = trans(mu_i_non_miss.elem(non_miss_i_idx));
      rowvec Y_i_non_miss = Y.row(i);
      Y_i_non_miss = trans(Y_i_non_miss.elem(non_miss_i_idx));
      Y_i.elem(miss_i_idx) = mu_i_miss - Y_cov_mm * Vinv_mo * trans(Y_i_non_miss - mu_i_non_miss);
      
      Y.row(i) = Y_i;
      
      sum_neg_ent_Y_miss += (0.5 * (s * log((1/(2*M_PI*exp(1)))) + chol2ldet(Vinv_mm_chol)));
    }
  }
}

// Inner loop rss
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
List inner_loop_general_rss_rcpp (unsigned int n, const arma::mat& XtX, const arma::mat& XtY, 
                                  arma::mat& XtRbar, arma::mat& mu1, const arma::mat& V, 
                                  const arma::mat& Vinv, const arma::vec& w0,
                                  const arma::cube& S0, const List& precomp_quants_list,
                                  bool standardize, bool compute_ELBO, bool update_V,
                                  const arma::vec& update_order, double eps, unsigned int nthreads) {
  unsigned int r = mu1.n_cols;
  unsigned int p = mu1.n_rows;
  unsigned int k = w0.n_elem;
  cube S1(r,r,p);
  mat  w1(p,k);
  mat  mu1_new  = mu1;
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
  inner_loop_general_rss(n, XtX, XtY, XtRbar, mu1_new, V, Vinv, w0, S0, precomp_quants,
                         standardize, compute_ELBO, update_V, update_order, eps,
                         nthreads, S1, w1, var_part_tr_wERSS, neg_KL, var_part_ERSS);
  return List::create(Named("mu1")                = mu1_new,
                      Named("S1")                 = S1,
                      Named("w1")                 = w1,
                      Named("var_part_tr_wERSS")  = var_part_tr_wERSS,
                      Named("neg_KL")             = neg_KL,
                      Named("var_part_ERSS")      = var_part_ERSS);
}

// Perform the inner loop rss
void inner_loop_general_rss (unsigned int n, const mat& XtX, const mat& XtY, mat& XtRbar, 
                             mat& mu1, const mat& V, const mat& Vinv, const vec& w0, 
                             const cube& S0, const mr_mash_precomputed_quantities& precomp_quants,
                             bool standardize, bool compute_ELBO, bool update_V, 
                             const vec& update_order, double eps, unsigned int nthreads,
                             cube& S1, mat& w1, double& var_part_tr_wERSS, 
                             double& neg_KL, mat& var_part_ERSS) {
  unsigned int p = mu1.n_rows;
  unsigned int r = mu1.n_cols;
  unsigned int k = w0.n_elem;
  vec Xtx(p);
  vec XtRbar_j(r);
  vec mu1_j(r);
  vec mu1_mix(r);
  mat S1_mix(r,r);
  vec w1_mix(k);
  double logbf_mix;
  double xtx_j;
  
  // Initialize ELBO parameters
  var_part_tr_wERSS = 0;
  neg_KL = 0;
  
  // Initialize V parameters
  var_part_ERSS.zeros(r,r);
  
  // Repeat for each predictor.
  for (unsigned int j : update_order) {
    
    if (standardize)
      xtx_j = n-1;
    else
      xtx_j = precomp_quants.xtx(j);
    
    Xtx = XtX.col(j);
    mu1_j = trans(mu1.row(j));
    
    // Disregard the jth predictor in the expected residuals.
    XtRbar += (Xtx * trans(mu1_j));
    XtRbar_j = trans(XtRbar.row(j));
    
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
      compute_ELBO_rss_terms(var_part_tr_wERSS, neg_KL, XtRbar_j, logbf_mix, mu1_mix, S1_mix, xtx_j, Vinv);
    
    // Compute V parameters
    if (update_V)
      compute_var_part_ERSS(var_part_ERSS, S1_mix, xtx_j);
    
    // Update the expected residuals.
    XtRbar -= (Xtx * trans(mu1_mix));
  }
}


// Inner loop rss sparse
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
List inner_loop_general_rss_sparse_rcpp (unsigned int n, const arma::sp_mat& XtX, const arma::mat& XtY, 
                                          arma::mat& XtRbar, arma::mat& mu1, const arma::mat& V, 
                                          const arma::mat& Vinv, const arma::vec& w0,
                                          const arma::cube& S0, const List& precomp_quants_list,
                                          bool standardize, bool compute_ELBO, bool update_V,
                                          const arma::vec& update_order, double eps, unsigned int nthreads) {
  unsigned int r = mu1.n_cols;
  unsigned int p = mu1.n_rows;
  unsigned int k = w0.n_elem;
  cube S1(r,r,p);
  mat  w1(p,k);
  mat  mu1_new  = mu1;
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
  inner_loop_general_rss_sparse(n, XtX, XtY, XtRbar, mu1_new, V, Vinv, w0, S0, precomp_quants,
                                standardize, compute_ELBO, update_V, update_order, eps,
                                nthreads, S1, w1, var_part_tr_wERSS, neg_KL, var_part_ERSS);
  return List::create(Named("mu1")                = mu1_new,
                      Named("S1")                 = S1,
                      Named("w1")                 = w1,
                      Named("var_part_tr_wERSS")  = var_part_tr_wERSS,
                      Named("neg_KL")             = neg_KL,
                      Named("var_part_ERSS")      = var_part_ERSS);
}

// Perform the inner loop rss
void inner_loop_general_rss_sparse (unsigned int n, const sp_mat& XtX, const mat& XtY, mat& XtRbar, 
                                    mat& mu1, const mat& V, const mat& Vinv, const vec& w0, 
                                    const cube& S0, const mr_mash_precomputed_quantities& precomp_quants,
                                    bool standardize, bool compute_ELBO, bool update_V, 
                                    const vec& update_order, double eps, unsigned int nthreads,
                                    cube& S1, mat& w1, double& var_part_tr_wERSS, 
                                    double& neg_KL, mat& var_part_ERSS) {
  unsigned int p = mu1.n_rows;
  unsigned int r = mu1.n_cols;
  unsigned int k = w0.n_elem;
  vec Xtx(p);
  vec XtRbar_j(r);
  vec mu1_j(r);
  vec mu1_mix(r);
  mat S1_mix(r,r);
  vec w1_mix(k);
  double logbf_mix;
  double xtx_j;
  
  // Initialize ELBO parameters
  var_part_tr_wERSS = 0;
  neg_KL = 0;
  
  // Initialize V parameters
  var_part_ERSS.zeros(r,r);
  
  // Repeat for each predictor.
  for (unsigned int j : update_order) {
    
    if (standardize)
      xtx_j = n-1;
    else
      xtx_j = precomp_quants.xtx(j);
    
    Xtx = XtX.col(j);
    mu1_j = trans(mu1.row(j));
    
    // Disregard the jth predictor in the expected residuals.
    XtRbar += (Xtx * trans(mu1_j));
    XtRbar_j = trans(XtRbar.row(j));
    
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
      compute_ELBO_rss_terms(var_part_tr_wERSS, neg_KL, XtRbar_j, logbf_mix, mu1_mix, S1_mix, xtx_j, Vinv);
    
    // Compute V parameters
    if (update_V)
      compute_var_part_ERSS(var_part_ERSS, S1_mix, xtx_j);
    
    // Update the expected residuals.
    XtRbar -= (Xtx * trans(mu1_mix));
  }
}
