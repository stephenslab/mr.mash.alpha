###Function to simulate data to test mr.mash
## X ~ N_r(0, Gamma); B ~ N_r(0, Sigma);
## Y ~ MN_nxr(XB, V)
## where Gamma=Gamma_cor*Gamma_scale, Sigma=Sigma_cor*Sigma_scale 
## and V=V_cor*V_scale. Gamma_cor, Sigma_cor, V_cor are *positive*
## correlation matrices. N.B. The diagonal elements of V are actually
## adjusted to have so that Y has the desired per-response PVE.
##
## Inputs are:
## n=number of sample; p=number of variables; 
## p_causal=number of causal variables; r=number of responses;
## intercepts=intercept value for each response; pve=per-response pve;
## Sigma_cor_offdiag=value of the off-diagonal elements of Sigma;
## Sigma_scale=scaling factor for Sigma;
## Gamma_cor_offdiag=value of the off-diagonal elements of Sigma;
## Gamma_scale=scaling factor for Sigma;
## V_cor_offdiag=value of the off-diagonal elements of V;
## V_scale=scaling factor for Sigma.
##
#' @importFrom mvtnorm rmvnorm
#' @importFrom MBSP matrix.normal
#' @importFrom matrixStats colVars
#' 
simulate_mr_mash_data <- function(n, p, p_causal, r, intercepts=rep(1, r),
                                  pve=0.2, Sigma_cor_offdiag=1, Sigma_scale=0.8,
                                  Gamma_cor_offdiag=0, Gamma_scale=0.8,
                                  V_cor_offdiag=0.25, V_offdiag_scale=0.8){
  
  ##Simulate true effects from N_r(0, Sigma) where Sigma is a given covariance matrix across traits
  Sigma <- create_cov_canonical(r, singletons=FALSE, hetgrid=Sigma_cor_offdiag)[[1]]*Sigma_scale
  B_causal <- rmvnorm(n=p_causal, mean=rep(0, r), sigma=Sigma)
  B <- matrix(0, ncol=r, nrow=p)
  causal_variables <- sample(x=(1:p), size=p_causal)
  B[causal_variables, ] <- B_causal
  
  ##Simulate X from N_r(0, Gamma) where Gamma is a given covariance matrix across variables
  Gamma <- create_cov_canonical(p, singletons=FALSE, hetgrid=Gamma_cor_offdiag)[[1]]*Gamma_scale
  X <- rmvnorm(n=n, mean=rep(0, p), sigma=Gamma)
  X <- scale_fast2(X, scale=FALSE)$M
  
  ##Compute G and its variance
  G <- X%*%B
  Var_G <- colVars(G)
  
  ##Compute residual covariance
  Var_E <- ((1/rep(pve, r))-1)*Var_G
  V <- create_cov_canonical(r, singletons=FALSE, hetgrid=V_cor_offdiag)[[1]]*V_offdiag_scale
  diag(V) <- Var_E
  
  ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
  Y <- matrix.normal(G + matrix(intercepts, n, r, byrow=TRUE), diag(n), V)
  
  return(list(X=X, Y=Y, B=B, causal_variables=sort(causal_variables), V=V, Sigma=Sigma, Gamma=Gamma))
}