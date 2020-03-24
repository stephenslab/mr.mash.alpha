// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// FUNCTION DECLARATIONS
// ---------------------

// 
void scale (arma::mat& M, const arma::vec& a, const arma::vec& b);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat scale_rcpp(const arma::mat& M, const arma::vec& a, const arma::vec& b){
  mat M_scaled = M;
  scale(M_scaled, a, b);
  
  return M_scaled;
}

// Scale a matrix
void scale (arma::mat& M, const arma::vec& a, const arma::vec& b) {
  unsigned int n = M.n_rows;
  unsigned int p = M.n_cols;
  
  for (unsigned int j = 0; j < p; j++)
    for (unsigned int i = 0; i < n; i++)
      M(i,j) = (M(i,j) - a(j)) / b(j);
}
