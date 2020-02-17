# Return the sigmoid of the elements of x. The sigmoid function is
# also known as the logistic link function. It is the inverse of
# logit(x).
sigmoid <- function (x) {
  if (x > -500)
    y <- 1/(1 + exp(-x))
  else
    y <- 0
  return(y)
}

# Compute the softmax of vector x in a more numerically prudent manner
# that avoids overflow or underflow.
softmax <- function (x) {
  y <- exp(x - max(x))
  return(y/sum(y))
}

# Return the dot product of vectors x and y.
dot <- function (x,y)
  sum(x*y)

# Return the quadratic norm (2-norm) of vector x.
norm2 <- function (x)
  sqrt(dot(x,x))

###Compute the trace of a square matrix
tr <- function(x)
  sum(diag(x))

# Should be the same as mvtnorm::dmvnorm(x,mu,S,log = TRUE)
#
#' @importFrom Rcpp evalCpp
#' @useDynLib mr.mash.alpha
#' 
dmvnorm <- function (x, mu, S)
  dmvnorm_rcpp(x,mu,S)

###Function to simulate from MN distribution
#
#' @importFrom MBSP matrix.normal
#' 
sim_mvr <- function (X, B, V) {
  
  # Get the number of samples (n) and conditions (m).
  n <- nrow(X)
  R <- ncol(B)
  
  # Simulate the responses, Y.
  M <- X%*%B
  U <- diag(n)
  Y <- matrix.normal(M, U, V)
  
  # Output the simulated responses.
  return(Y)
}

###Function to compute canonical covariance matrices scaled by a grid 
compute_cov_canonical <- function(ntraits, singletons, hetgrid, grid, zeromat=T){
  S <- mmbr:::create_cov_canonical(ntraits, singletons, hetgrid)
  U <- list()
  t <- 0
  for(i in 1:length(S)){
    for(j in 1:length(grid)){
      t <- t+1
      U[[t]] <- S[[i]]*grid[j]
    }
  }
  
  names(U) <- paste0("S0_", seq(1, length(U)))
  
  if(zeromat){
    zero_mat <- matrix(0, ntraits, ntraits)
    zero_mat[upper.tri(zero_mat)] <- 1e-10
    zero_mat[lower.tri(zero_mat)] <- 1e-10
    U[[paste0("S0_", length(U)+1)]] <- zero_mat
  }
  
  return(U)
}

###Update mixture weights
update_weights <- function(x){
  w <- colSums(x)
  w <- w/sum(w)
  return(w)
}

##Compute ELBO from intermediate components
compute_ELBO_fun <- function(rbar, V, Vinv, ldetV, var_part_ERSS, neg_KL){
  n <- nrow(rbar)
  R <- ncol(rbar)
  ERSS <- tr(Vinv%*%(crossprod(rbar))) + var_part_ERSS
  ELBO <- -log(n)/2 - (n*R)/2*log(2*pi) - n/2 * ldetV - 0.5*ERSS + neg_KL
  
  return(ELBO)
}

###Compute log-determinant from Cholesky decomposition
chol2ldet <- function(R){
  logdet <- log(prod(diag(R)))*2
  
  return(logdet)
}

###Compute quantities needed when using scaled X
precompute_quants_scaled_X <- function(n, V, S0){
  ###Quantities that don't depend on S0
  R <- chol(V)
  S <- V/(n-1)
  S_chol <- R/sqrt(n-1)
  ldetS_chol <- chol2ldet(S_chol)
  
  ###Quantities that depend on S0
  SplusS0_chol <- list()
  S1 <- list()
  ldetSplusS0_chol <- c()
  for(i in 1:length(S0)){
    SplusS0_chol[[i]] <- chol(S+S0[[i]])
    ldetSplusS0_chol[i] <- chol2ldet(SplusS0_chol[[i]])
    S1[[i]] <- S0[[i]]%*%backsolve(SplusS0_chol[[i]], forwardsolve(t(SplusS0_chol[[i]]), S))
  }
  
  return(list(V_chol=R, S=S, S1=S1, S_chol=S_chol, SplusS0_chol=SplusS0_chol, 
              ldetS_chol=ldetS_chol, ldetSplusS0_chol=ldetSplusS0_chol))
}

###Compute quantities needed when using transformed X
precompute_quants_transformed_X <- function(X, V, S0){
  ###Quantities that don't depend on S0
  xtx <- diag(crossprod(X))
  R <- chol(V)
  Rtinv <- solve(t(R))
  Rinv <- solve(R)
  
  ###Quantities that depend on S0
  U0 <- list()
  d <- list()
  Q <- list()
  for(i in 1:length(S0)){
    U0[[i]]  <- Rtinv %*% S0[[i]] %*% Rinv
    out <- eigen(U0[[i]])
    d[[i]]   <- out$values
    Q[[i]]   <- out$vectors   
  }
  
  return(list(xtx=xtx, V_chol=R, U0=U0, d=d, Q=Q))
}