// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include <stdexcept>


using namespace Rcpp;
using namespace arma;


// FUNCTION DECLARATIONS
// ---------------------

// 
void scale (arma::mat& M, const arma::vec& a, const arma::vec& b);

// FUNCTION DEFINITIONS
// --------------------
// Scale a matrix providing column means and sds externally
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat scale_rcpp(const arma::mat& M, const arma::vec& a, const arma::vec& b){
  mat M_scaled = M;
  scale(M_scaled, a, b);
  
  return M_scaled;
}

void scale (arma::mat& M, const arma::vec& a, const arma::vec& b) {
  unsigned int n = M.n_rows;
  unsigned int p = M.n_cols;
  
  for (unsigned int j = 0; j < p; j++)
    for (unsigned int i = 0; i < n; i++)
      M(i,j) = (M(i,j) - a(j)) / b(j);
}


// Scale a matrix
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List scale2_rcpp (const arma::mat& M, bool scale, bool na_rm) {
  unsigned int n = M.n_rows;
  unsigned int p = M.n_cols;
  
  vec M_j(p);
  vec a(p);
  vec b(p);
  mat M1(n,p);
  
  for (unsigned int j = 0; j < p; j++){
    M_j = M.col(j);
    if(na_rm)
      //Keep only finite values
      M_j = M_j.elem(find_finite(M_j));
    a(j) = mean(M_j);
    if(scale){
      b(j) = stddev(M_j, 0);
      // Exit if a column has standard deviation of 0
      if(b(j)==0){
        throw std::runtime_error("Column "+std::to_string(j+1)+" has 0 standard deviation");
      }
    } else {
      b(j) = 1;
    }
    
    for (unsigned int i = 0; i < n; i++){    
      M1(i,j) = (M(i,j) - a(j)) / b(j);
    }
  }
  
  return List::create(Named("M")     = M1,
                      Named("means") = a,
                      Named("sds")   = b);
}
