#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------

// Compute the softmax of x, and return the result in x. Guard against
// underflow or overflow by adjusting the entries of x so that the
// largest value is zero.
void softmax (vec& x) {
  x -= max(x);
  x  = exp(x);
  x /= sum(x);
}

// Compute the log-probability density of the multivariate normal
// distribution with zero mean and covariance matrix S, omitting terms
// that do not depend on x or S.
double ldmvnorm (const vec& x, const mat& S) {
  mat    L = chol(S,"lower");
  double d = norm(solve(L,x),2);
  return -(d*d)/2 - sum(log(L.diag()));
}

// Compute the difference of two multivariate normal log-densities,
//
//   ldmvnorm(x,S0 + S) - ldmvnorm(x,S)
//
// where S_chol is the right-hand Cholesky factor of S, and
// SplusS0_chol is the right-hand Cholesky factor of S + S0 (also an
// upper triangular matrix).
double ldmvnormdiff (const vec& x, const mat& S_chol,
                     const mat& SplusS0_chol) {
  return (chol2ldet(S_chol) -
          chol2ldet(SplusS0_chol) +
          dot(x, backsolve(S_chol, forwardsolve(S_chol, x))) - 
          dot(x, backsolve(SplusS0_chol, forwardsolve(SplusS0_chol, x))))/2;
}

// Compute the log determinant from Cholesky decomposed matrix
double chol2ldet (const mat& R) {
  return 2*sum(log(R.diag()));
}
