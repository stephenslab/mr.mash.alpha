#' @title Simulate data to test \code{mr.mash}.
#' 
#' @description Function to simulate data \eqn{Y \sim MN_{nxr}(XB, I, V)}, where \eqn{X \sim N_p(0, Gamma)},
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
#' @param r_causal a list whose elelments ar numeric vector(s) indicating in which responses
#'  the causal variables have an effect (defined as a block).
#'  
#' @param wb scalar or numeric vector (one element for each block of causal responses)
#'  of probabilities of causal variables coming from each of the blocks.
#' 
#' @param intercepts numeric vector of intercept for each response.
#' 
#' @param pve per-response proportion of variance explained by the causal variables.
#' 
#' @param B_cor list whose elements (one element for each block of causal responses) are scalars or numeric vectors
#'   (one element for each mixture components) of positive correlation [0, 1] between causal effects.   
#'   
#' @param B_scale list whose elements (one element for each block of causal responses) are scalars or numeric vectors
#'   (one element for each mixture components) with the diagonal value for Sigma_k.
#'   
#' @param w list whose elements (one element for each block of causal responses) are scalars or numeric vectors
#'   (one element for each mixture components) with mixture proportions.
#'   
#' @param X_cor scalar indicating the positive correlation [0, 1] between variables.
#' 
#' @param X_scale scalar indicating the diagonal value for Gamma.
#' 
#' @param V_cor scalar indicating the positive correlation [0, 1] between residuals
#' 
#' @return A list with the following elements:
#' 
#' \item{X}{n x p matrix of variables.}
#' 
#' \item{Y}{n x r matrix of responses.}
#' 
#' \item{B}{p x r matrix of effects.}
#' 
#' \item{V}{r x r residual covariance matrix among responses.}
#'   
#' \item{Sigma}{list of r x r covariance matrices among the effects for each causal variable.} 
#' 
#' \item{Gamma}{p x p covariance matrix among the variables.}
#' 
#' \item{intercepts}{r-vector of intercept for each response.}
#' 
#' \item{causal_responses}{a list of numeric vectors of indexes indicating which responses have
#'   causal effects for each causal variable.}
#'   
#' \item{causal_variables}{p_causal-vector of indexes indicating which variables are causal.}  
#'   
#' @importFrom mvtnorm rmvnorm
#' @importFrom MBSP matrix.normal
#' @importFrom matrixStats colVars
#' 
#' @export
#' 
#' @examples
#' 
#' dat <- simulate_mr_mash_data(n=50, p=40, p_causal=20, r=5, r_causal=list(1:2, 3:5), wb=c(0.5, 0.5), intercepts=rep(1, 5),
#'                              pve=0.2, B_cor=list(c(0,1), c(0,0.5,1)), B_scale=list(c(0.5,0.8), c(0.7,1.5,2)), 
#'                              w=list(c(0.5,0.5), c(0.33,0.33,0.34)),
#'                              X_cor=0.5, X_scale=1, V_cor=0)
#'                              
#'                              
simulate_mr_mash_data <- function(n, p, p_causal, r, r_causal=list(1:r), wb=1, intercepts=rep(1, r),
                                  pve=0.2, B_cor=list(0), B_scale=list(1), w=list(1),
                                  X_cor=0.5, X_scale=1, V_cor=0){
  ##Check that the inputs are correct
  if(length(intercepts)!=r)
    stop("intercepts must be of length equal to r.")
  if(any(sapply(r_causal, length)>r))
    stop("r_causal cannot be greater than r.")
  if(!(length(B_cor)==length(B_scale) & length(w)==length(B_cor) & length(w)==length(r_causal)))
    stop("B_cor, B_scale, w, and r_causal must have the same length.")
  if(!(any(sapply(B_cor, length)==sapply(B_scale, length) & sapply(w, length)==sapply(B_cor, length))))
    stop("Corresponding elements of B_cor, B_scale, and w must have the same length.")
  if(any(sapply(w, sum) - 1 > 1e-10))
    stop("Each element of w must sum to 1.")
  if(abs(sum(wb) - 1) > 1e-10)
    stop("Elements of wb must sum to 1.")
  if(length(r_causal)!=length(wb))
    stop("r_causal and wb must have the same length.")
  
  ##Get number of blocks
  blocks <- length(wb)
  
  ##Simulate true effects from N_r(0, Sigma) or \sum_K w_k N_r(0, Sigma_k) where Sigma and Sigma_k are given 
  ##covariance matrices across traits and w_k is the mixture proportion associated to Sigma_k
  B_causal <- matrix(0, nrow=p_causal, ncol=r)
  Sigma_all <- vector("list", p_causal)
  blocks_all <- rep(as.numeric(NA), p_causal)
  mixcomp_all <- rep(as.numeric(NA), p_causal)
  
  for(j in 1:p_causal){
    #Sample block
    block_to_use <- sample(x=1:blocks, size=1, prob=wb, replace=TRUE)
    #Sample mixture component
    K <- length(B_cor[[block_to_use]])
    mixcomp_to_use <- sample(x=1:K, size=1, prob=w[[block_to_use]], replace=TRUE)
    #Build covariance matrix corresponding to the sampled block and component
    r_mix_length <- length(r_causal[[block_to_use]])
    Sigma_offdiag <- B_scale[[block_to_use]][mixcomp_to_use]*B_cor[[block_to_use]][mixcomp_to_use]
    Sigma <- matrix(Sigma_offdiag, nrow=r_mix_length, ncol=r_mix_length)
    diag(Sigma) <- B_scale[[block_to_use]][mixcomp_to_use]
    #Sample effects from a MVN with the given covariance matrix
    r_causal_mix <- r_causal[[block_to_use]]
    B_causal[j, r_causal_mix] <- rmvnorm(n=1, mean=rep(0, length(r_causal_mix)), sigma=Sigma)
    #Store infor for the current causal variable
    Sigma_all[[j]] <- Sigma
    blocks_all[j] <- block_to_use
    mixcomp_all[j] <- mixcomp_to_use
  }
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
  Var_E[which(Var_E<=.Machine$double.eps)] <- 1
  D <- diag(x=sqrt(Var_E))
  V_cor_mat <- matrix(V_cor, nrow=r, ncol=r)
  diag(V_cor_mat) <- 1
  V <- D %*% V_cor_mat %*% D
  
  ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
  Y <- matrix.normal(G + matrix(intercepts, n, r, byrow=TRUE), diag(n), V)
  
  ##Compile output
  causal_variables_blocks_mixcomps <- cbind(causal_variables, blocks_all, mixcomp_all, 1:p_causal)
  causal_variables_blocks_mixcomps <- causal_variables_blocks_mixcomps[order(causal_variables_blocks_mixcomps[, 1]), ]
  out <- list(X=X, Y=Y, B=B, V=V, intercepts=intercepts, causal_variables=causal_variables_blocks_mixcomps[,1], 
              causal_responses=r_causal[causal_variables_blocks_mixcomps[, 2]], 
              Sigma=Sigma_all[causal_variables_blocks_mixcomps[, 4]], Gamma=Gamma)
  
  return(out)
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
