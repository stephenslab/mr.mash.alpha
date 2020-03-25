###Function to simulate data to test mr.mash
## X ~ N_r(0, Gamma_scale*Gamma); B ~ N_r(0, Sigma_scale*Sigma);
## Y ~ MN_nxr(XB, V_scale*V)
## where Gamma, Sigma and V are correlation matrices.
## Inputs are:
## n=number of sample; p=number of variables; 
## p_causal=number of causal variables; r=number of responses;
## Sigma_offdiag=value of the off-diagonal elements of Sigma;
## Sigma_scale=scaling factor for Sigma;
## Gamma_offdiag=value of the off-diagonal elements of Sigma;
## Gamma_scale=scaling factor for Sigma;
## V_offdiag=value of the off-diagonal elements of V;
## V_scale=scaling factor for Sigma.
##
#' @importFrom mvtnorm rmvnorm
#' 
simulate_mr_mash_data <- function(n, p, p_causal, r, 
                                  Sigma_offdiag=1, Sigma_scale=0.8,
                                  Gamma_offdiag=0, Gamma_scale=0.8,
                                  V_offdiag=0.25, V_scale=0.8){
  ##Compute residual covariance
  V <- create_cov_canonical(r, singletons=FALSE, hetgrid=V_offdiag)
  V <- V_scale * V[[1]]
  
  ##Simulate true effects from N_r(0, Sigma) where Sigma is a given covariance matrix across traits
  Sigma <- create_cov_canonical(r, singletons=FALSE, hetgrid=Sigma_offdiag)
  B_causal <- rmvnorm(n=p_causal, mean=rep(0, r), sigma=Sigma_scale*Sigma[[1]])
  B <- matrix(0, ncol=r, nrow=p)
  B[1:p_causal, ] <- B_causal
  
  ##Simulate X from N_r(0, Gamma) where Gamma is a given covariance matrix across variables
  Gamma <- create_cov_canonical(p, singletons=FALSE, hetgrid=Gamma_offdiag)
  X <- rmvnorm(n=n, mean=rep(0, p), sigma=Gamma_scale*Gamma[[1]])
  X <- scale_fast2(X, scale=FALSE)$M
  
  ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
  Y <- sim_mvr(X, B, V)
  
  return(list(X=X, Y=Y, B=B, V=V))
}