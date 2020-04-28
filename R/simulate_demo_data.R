###Function to simulate data to test mr.mash
## X ~ N_r(0, Gamma); B ~ N_r(0, Sigma);
## Y ~ MN_nxr(XB, V)
## where Gamma, Sigma, and V are chosen to achieve the desired correlations
## and PVE. 
##
## Inputs are:
## n=number of sample; p=number of variables; 
## p_causal=number of causal variables; r=number of responses;
## intercepts=intercept value for each response; pve=per-response pve;
## B_cor=positive correlation [0, 1] between causal effects;
## B_scale=diagonal value for Sigma;
## X_cor=positive correlation [0, 1] between predictors;
## X_scale=diagonal value for Gamma;
## V_cor=positive correlation [0, 1] between residuals;
##
#' @importFrom mvtnorm rmvnorm
#' @importFrom MBSP matrix.normal
#' @importFrom matrixStats colVars
#' 
simulate_mr_mash_data <- function(n, p, p_causal, r, intercepts=rep(1, r),
                                  pve=0.2, B_cor=0, B_scale=0.8,
                                  X_cor=0.5, X_scale=0.8, V_cor=0.25){
  
  if(length(intercepts)!=r)
    stop("intercepts must be of length equalt to r.")
  
  ##Simulate true effects from N_r(0, Sigma) where Sigma is a given covariance matrix across traits
  Sigma_offdiag <- B_scale*B_cor
  Sigma <- matrix(Sigma_offdiag, nrow=r, ncol=r)
  diag(Sigma) <- B_scale
  B_causal <- rmvnorm(n=p_causal, mean=rep(0, r), sigma=Sigma)
  B <- matrix(0, ncol=r, nrow=p)
  causal_variables <- sample(x=(1:p), size=p_causal)
  B[causal_variables, ] <- B_causal
  
  ##Simulate X from N_r(0, Gamma) where Gamma is a given covariance matrix across variables
  Gamma_offdiag <- X_scale*X_cor
  Gamma <- matrix(Gamma_offdiag, nrow=p, ncol=p)
  diag(Gamma) <- X_scale
  X <- rmvnorm(n=n, mean=rep(0, p), sigma=Gamma)
  X <- scale_fast2(X, scale=FALSE)$M
  
  ##Compute G and its variance
  G <- X%*%B
  Var_G <- colVars(G)
  
  ##Compute residual covariance
  Var_E <- ((1/pve)-1)*Var_G
  D <- diag(x=sqrt(Var_E))
  V_cor_mat <- matrix(V_cor, nrow=r, ncol=r)
  diag(V_cor_mat) <- 1
  V <- D %*% V_cor_mat %*% D
  
  ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
  Y <- matrix.normal(G + matrix(intercepts, n, r, byrow=TRUE), diag(n), V)
  
  return(list(X=X, Y=Y, B=B, causal_variables=sort(causal_variables), V=V, Sigma=Sigma, Gamma=Gamma, intercepts=intercepts))
}







### Older version only kept until some analyses are re-run
# simulate_mr_mash_data <- function(n, p, p_causal, r, intercepts=rep(1, r),
#                                   pve=0.2, Sigma_cor_offdiag=1, Sigma_scale=0.8,
#                                   Gamma_cor_offdiag=0, Gamma_scale=0.8,
#                                   V_cor_offdiag=0.25, V_offdiag_scale=0.8){
#   
#   ##Simulate true effects from N_r(0, Sigma) where Sigma is a given covariance matrix across traits
#   Sigma <- create_cov_canonical(r, singletons=FALSE, hetgrid=Sigma_cor_offdiag)[[1]]*Sigma_scale
#   B_causal <- rmvnorm(n=p_causal, mean=rep(0, r), sigma=Sigma)
#   B <- matrix(0, ncol=r, nrow=p)
#   causal_variables <- sample(x=(1:p), size=p_causal)
#   B[causal_variables, ] <- B_causal
#   
#   ##Simulate X from N_r(0, Gamma) where Gamma is a given covariance matrix across variables
#   Gamma <- create_cov_canonical(p, singletons=FALSE, hetgrid=Gamma_cor_offdiag)[[1]]*Gamma_scale
#   X <- rmvnorm(n=n, mean=rep(0, p), sigma=Gamma)
#   X <- scale_fast2(X, scale=FALSE)$M
#   
#   ##Compute G and its variance
#   G <- X%*%B
#   Var_G <- colVars(G)
#   
#   ##Compute residual covariance
#   Var_E <- ((1/rep(pve, r))-1)*Var_G
#   V <- create_cov_canonical(r, singletons=FALSE, hetgrid=V_cor_offdiag)[[1]]*V_offdiag_scale
#   diag(V) <- Var_E
#   
#   ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
#   Y <- matrix.normal(G + matrix(intercepts, n, r, byrow=TRUE), diag(n), V)
#   
#   return(list(X=X, Y=Y, B=B, causal_variables=sort(causal_variables), V=V, Sigma=Sigma, Gamma=Gamma))
# }
