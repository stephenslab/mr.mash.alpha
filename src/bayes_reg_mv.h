#ifndef INCLUDE_BAYES_REG_MV
#define INCLUDE_BAYES_REG_MV

#include <RcppArmadillo.h>


// FUNCTION DECLARATIONS
// ---------------------
double bayes_mvr_ridge_standardized_X (const arma::vec& b, const arma::mat& S0, const arma::mat& S,
                                 const arma::mat& S1, const arma::mat& SplusS0_chol,
                                 const arma::mat& S_chol, arma::vec& mu1);

double bayes_mvr_ridge_centered_X (const arma::mat& V, const arma::vec& b, const arma::mat& S, 
                                   const arma::mat& S0, double xtx, const arma::mat& Vinv,
                                   const arma::mat& V_chol, const arma::mat& S_chol,
                                   const arma::vec& d, const arma::mat& QtimesV_chol,
                                   arma::vec& mu1, arma::mat& S1);

double bayes_mvr_mix_standardized_X (const arma::vec& x, const arma::mat& Y, const arma::vec& w0,
                               const arma::cube& S0, const arma::mat& S, const arma::cube& S1,
                               const arma::cube& SplusS0_chol, const arma::mat& S_chol, double eps,
                               arma::vec& mu1_mix, arma::mat& S1_mix, arma::vec& w1);

double bayes_mvr_mix_centered_X (const arma::vec& x, const arma::mat& Y, const arma::mat& V,
                                 const arma::vec& w0, const arma::cube& S0, double xtx, 
                                 const arma::mat& Vinv, const arma::mat& V_chol,
                                 const arma::mat& d, const arma::cube& QtimesV_chol, double eps,
                                 arma::vec& mu1_mix, arma::mat& S1_mix, arma::vec& w1);

#endif