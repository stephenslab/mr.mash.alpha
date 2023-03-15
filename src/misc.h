#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Solve for x in U'x = b by forward substitution.
inline arma::vec forwardsolve (const arma::mat& U, const arma::vec& b) {
  return arma::solve(arma::trimatl(arma::trans(U)),b);  
}

// Solve for x in Ux = b by back substitution.
inline arma::vec backsolve (const arma::mat& U, const arma::vec& b) {
  return arma::solve(arma::trimatu(U),b);
}


// FUNCTION DECLARATIONS
// ---------------------
void softmax (arma::vec& x);

double ldmvnorm (const arma::vec& x, const arma::mat& S);

double ldmvnormdiff (const arma::vec& x, const arma::mat& S_chol,
                     const arma::mat& SplusS0_chol);

double chol2ldet (const arma::mat& R);

void compute_ELBO_terms (double& var_part_tr_wERSS, double& neg_KL, const arma::vec& x_j,
                         const arma::mat& rbar_j, double logbf, const arma::mat& mu1, const arma::mat& S1, 
                         double xtx, const arma::mat& Vinv);

void compute_ELBO_rss_terms (double& var_part_tr_wERSS, double& neg_KL,
                             const arma::mat& XtRbar_j, double logbf, const arma::mat& mu1, const arma::mat& S1, 
                             double xtx, const arma::mat& Vinv);

void compute_var_part_ERSS (arma::mat& var_part_ERSS, const arma::mat& S1, double xtx);


#endif