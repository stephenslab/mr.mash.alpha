#' @title Simulate data to test \code{mr.mash}.
#' 
#' @description Function to simulate data from \eqn{MN_{nxr}(XB, I, V)}, where \eqn{X \sim N_p(0, Gamma)},
#'   \eqn{B \sim \sum_k w_k N_r(0, Sigma_k)}, with \eqn{Gamma}, \eqn{w_k}, \eqn{Sigma_k}, and \eqn{V} defined by the user.
#'   
#' @param n scalar indicating the number of samples.
#' 
#' @param p scalar indicating the number of variables.
#' 
#' @param p_causal scalar indicating the number of causal variables.
#' 
#' @param r scalar indicating the number of responses.
#' 
#' @param r_causal a list of numeric vectors (one element for each mixture component) 
#'   indicating in which responses the causal variables have an effect.
#' 
#' @param intercepts numeric vector of intercept for each response.
#' 
#' @param pve per-response proportion of variance explained by the causal variables.
#' 
#' @param B_cor scalar or numeric vector (one element for each mixture component)
#'   with positive correlation [0, 1] between causal effects.   
#'   
#' @param B_scale scalar or numeric vector (one element for each mixture component) with the 
#'   diagonal value for Sigma_k;
#'   
#' @param w scalar or numeric vector (one element for each mixture component) with mixture 
#'   proportions associated to each mixture component.
#'   
#' @param X_cor scalar indicating the positive correlation [0, 1] between variables.
#' 
#' @param X_scale scalar indicating the diagonal value for Gamma.
#' 
#' @param V_cor scalar indicating the positive correlation [0, 1] between residuals
#'
#' 
#' @return A list with some or all of the
#' following elements:
#' 
#' \item{X}{n x p matrix of variables.}
#' 
#' \item{Y}{n x r matrix of responses.}
#' 
#' \item{B}{p x r matrix of effects.}
#' 
#' \item{V}{r x r residual covariance matrix among responses.}
#'   
#' \item{Sigma}{list of r x r covariance matrices among the effects.} 
#' 
#' \item{Gamma}{p x p covariance matrix among the variables.}
#' 
#' \item{intercepts}{r-vector of intercept for each response.}
#' 
#' \item{causal_responses}{a list of numeric vectors of indexes indicating which responses have
#'   causal effects for each mixture component.}
#'   
#' \item{causal_variables}{p_causal-vector of indexes indicating which variables are causal.}  
#' 
#' \item{causal_vars_to_mixture_comps}{p_causal-vector of indexes indicating from which
#'   mixture components each causal effect comes.}
#'   
#' @importFrom mvtnorm rmvnorm
#' @importFrom matrixStats colVars
#' 
#' @export
#' 
#' 
#' @examples
#' set.seed(1)
#' dat <- simulate_mr_mash_data(n=50, p=40, p_causal=20, r=5,
#'                              r_causal=list(1:2, 3:4), intercepts=rep(1, 5),
#'                              pve=0.2, B_cor=c(0, 1), B_scale=c(0.5, 1),
#'                              w=c(0.5, 0.5), X_cor=0.5, X_scale=1, 
#'                              V_cor=0)
#'                              
#'                              
simulate_mr_mash_data <- function(n, p, p_causal, r, r_causal=list(1:r), intercepts=rep(1, r),
                                  pve=0.2, B_cor=1, B_scale=1, w=1,
                                  X_cor=0, X_scale=1, V_cor=0){
  ##Check that the inputs are correct
  if(length(intercepts)!=r)
    stop("intercepts must be of length equal to r.")
  if(any(sapply(r_causal, length)>r))
    stop("r_causal cannot be greater than r.")
  if(!(length(B_cor)==length(B_scale) & length(w)==length(B_cor) & length(B_scale)==length(w)))
    stop("B_cor, B_scale, and w must have the same length.")
  if(abs(sum(w) - 1) > 1e-10)
    stop("Elements of w must sum to 1.")
  if(length(pve)!=1 & length(pve)!=r)
    stop("pve must be of length equal to 1 or r.")
  
  ##Get number of mixture components
  K <- length(w)
  
  ##Simulate true effects from N_r(0, Sigma) or \sum_K w_k N_r(0, Sigma_k) where Sigma and Sigma_k are given 
  ##covariance matrices across traits and w_k is the mixture proportion associated to Sigma_k
  Sigma <- vector("list", K)
  for(i in 1:K){
    r_mix_length <- length(r_causal[[i]])
    Sigma_offdiag <- B_scale[i]*B_cor[i]
    Sigma[[i]] <- matrix(Sigma_offdiag, nrow=r_mix_length, ncol=r_mix_length)
    diag(Sigma[[i]]) <- B_scale[i]
  }

  #Sample effects from a mixture of MVN distributions or a single MVN distribution
  B_causal <- matrix(0, nrow=p_causal, ncol=r)
  if(K>1){
    mixcomps <- sample(x=1:K, size=p_causal, prob=w, replace=TRUE)
    for(j in 1:p_causal){
      comp_to_use <- mixcomps[j]
      r_causal_mix <- r_causal[[comp_to_use]]
      B_causal[j, r_causal_mix] <- rmvnorm(n=1, mean=rep(0, length(r_causal_mix)), sigma=Sigma[[comp_to_use]])
    }
  } else {
    r_causal_length <- length(r_causal[[1]])
    r_causal_index <- r_causal[[1]]
    B_causal[, r_causal_index] <- rmvnorm(n=p_causal, mean=rep(0, r_causal_length), sigma=Sigma[[1]])
  }
  B <- matrix(0, ncol=r, nrow=p)
  causal_variables <- sample(x=(1:p), size=p_causal)
  B[causal_variables, ] <- B_causal
  
  ##Simulate X from N_r(0, Gamma) where Gamma is a given covariance matrix across variables
  if(X_cor != 0){
    Gamma_offdiag <- X_scale*X_cor
    Gamma <- matrix(Gamma_offdiag, nrow=p, ncol=p)
    diag(Gamma) <- X_scale
    X <- rmvnorm(n=n, mean=rep(0, p), sigma=Gamma)
  } else {
    X <- sapply(1:p, sample_norm, n=n, m=0, s2=X_scale)
    Gamma <- paste0("I_", p)
  }
  X <- scale_fast2(X, scale=FALSE)$M
  
  ##Compute G and its variance
  G <- X%*%B
  Var_G <- colVars(G)
  
  ##Compute residual covariance
  Var_E <- ((1/pve)-1)*Var_G
  Var_E[which(Var_E<=.Machine$double.eps)] <- 1
  D <- diag(x=sqrt(Var_E))
  V_cor_mat <- matrix(V_cor, nrow=r, ncol=r)
  diag(V_cor_mat) <- 1
  V <- D %*% V_cor_mat %*% D
  
  ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
  Y <- matrix_normal_indep_rows(M=(G + matrix(intercepts, n, r, byrow=TRUE)), V=V)
  
  ##Compile output
  causal_responses <- r_causal
  names(causal_responses) <- paste0("Component", 1:K)
  names(Sigma) <- paste0("Component", 1:K)
  out <- list(X=X, Y=Y, B=B, V=V, Sigma=Sigma, Gamma=Gamma,
              intercepts=intercepts, causal_responses=causal_responses)
  if(K>1){
    if(p_causal>1){
      causal_variables_mixcomps <- cbind(causal_variables, mixcomps)
      causal_variables_mixcomps <- causal_variables_mixcomps[order(causal_variables_mixcomps[, 1]), ]
      out$causal_variables <- causal_variables_mixcomps[, 1]
      out$causal_vars_to_mixture_comps <- causal_variables_mixcomps[, 2]
    } else {
      out$causal_variables <- causal_variables
      out$causal_vars_to_mixture_comps <- mixcomps
    }
  } else {
    out$causal_variables <- sort(causal_variables)
    out$causal_vars_to_mixture_comps <- rep(1, p_causal)
  }
  
  return(out)
}


#' @importFrom stats rnorm
sample_norm <- function(i, n, m, s2){
  x <- rnorm(n=n, mean=m, sd=sqrt(s2))
  
  return(x)
}

#' @importFrom stats rnorm
matrix_normal_indep_rows = function(M, V, seed){
  a <- nrow(M)
  b <- ncol(M)
  
  # Draw Z from MN(O, I, I)
  Z <- matrix(rnorm(a*b,0,1), a, b)
  
  # Cholesky decomposition of V
  L2 <- chol(V)
  
  # Return draw from MN(M,I,V)
  return(M + Z %*% L2)
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
#   Y <- matrix_normal(G + matrix(intercepts, n, r, byrow=TRUE), diag(n), V)
#   
#   return(list(X=X, Y=Y, B=B, causal_variables=sort(causal_variables), V=V, Sigma=Sigma, Gamma=Gamma))
# }
