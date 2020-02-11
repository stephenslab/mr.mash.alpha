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
compute_ELBO_fun <- function(rbar, V, Vinv, var_part_ERSS, neg_KL){
  n <- nrow(rbar)
  R <- ncol(rbar)
  ERSS <- tr(Vinv%*%(crossprod(rbar))) + var_part_ERSS
  ELBO <- -log(n)/2 - (n*R)/2*log(2*pi) - n/2 * as.numeric(determinant(V, logarithm = TRUE)$modulus) - 0.5*ERSS + neg_KL
  
  return(ELBO)
}

###Update variational parameters, expected residuals, and ELBO components
inner_loop <- function(X, rbar, mu, V, Vinv, w0, S0){
  ###Create variables to store quantities
  R <- ncol(rbar)
  p <- ncol(X)
  K <- length(S0)
  mu1   <- mu
  S1    <- array(0, c(R, R, p))
  w1    <- matrix(0, nrow=p, ncol=K)
  
  if(!is.null(Vinv)){
    ##Initialize ELBO parameters
    var_part_ERSS <- 0
    neg_KL <- 0
  }

  ##Loop through the variables
  for(j in 1:p){
    
    #Remove j-th effect from expected residuals 
    rbar_j <- rbar + outer(X[, j], mu1[j, ])
    
    #Run Bayesian SLR
    bfit <- bayes_mvr_mix(X[, j], rbar_j, V, w0, S0)
    
    #Update variational parameters
    mu1[j, ]         <- bfit$mu1
    S1[, , j]        <- bfit$S1
    w1[j, ]          <- bfit$w1
    
    #Compute ELBO params
    if(!is.null(Vinv)){
      xtx <- sum(X[, j]^2)
      var_part_ERSS <- var_part_ERSS + (tr(Vinv%*%bfit$S1)*xtx)
      mu1_mat <- matrix(bfit$mu1, ncol=1)
      neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, rbar_j)%*%tcrossprod(matrix(X[, j], ncol=1), mu1_mat))+
                                             tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*(xtx)))
    }
    
    #Update expected residuals
    rbar <- rbar_j - outer(X[, j], mu1[j, ])
  }
  
  ###Return output
  if(!is.null(Vinv)){
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1, var_part_ERSS=var_part_ERSS, neg_KL=neg_KL))
  } else {
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1))
  }
}

###Perform one iteration of the outer loop
mr_mash_update <- function(Y, X, mu1_t, w1_t, V, Vinv, w0, S0, update_w0, compute_ELBO){
  ##Compute expected residuals
  rbar <- Y - X%*%mu1_t
  
  #Update w0 if requested
  if(update_w0 && !is.null(w1_t)){
    w0 <- update_weights(w1_t)
  }
  
  ##Update variational parameters, expected residuals, and ELBO components
  if(compute_ELBO){
    updates <- inner_loop(X=X, rbar=rbar, mu=mu1_t, V=V, Vinv=Vinv, w0=w0, S0=S0) 
  } else {
    updates <- inner_loop(X=X, rbar=rbar, mu=mu1_t, V=V, Vinv=NULL, w0=w0, S0=S0)
  }
  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  rbar    <- updates$rbar
  
  if(compute_ELBO){
    ##Compute ELBO
    var_part_ERSS <- updates$var_part_ERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, var_part_ERSS=var_part_ERSS, neg_KL=neg_KL)
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, ELBO=ELBO))
  } else {
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t))
  }
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
  
  ###Quantities that depend on S0
  SplusS0_chol <- list()
  S1 <- list()
  for(i in 1:length(S0)){
    SplusS0_chol[[i]] <- chol(S+S0[[i]])
    S1[[i]] <- S0[[i]]%*%backsolve(SplusS0_chol[[i]], forwardsolve(t(SplusS0_chol[[i]]), S))
  }
  
  return(list(R=R, S=S, S1=S1, S_chol=S_chol, SplusS0_chol=SplusS0_chol))
}

###Update variational parameters, expected residuals, and ELBO components with scaled X
inner_loop_scaled_X <- function(X, rbar, mu, Vinv, w0, S0, S, S1, SplusS0_chol, S_chol){
  ###Create variables to store quantities
  R <- ncol(rbar)
  p <- ncol(X)
  K <- length(S0)
  mu1c   <- mu
  S1c    <- array(0, c(R, R, p))
  w1c    <- matrix(0, nrow=p, ncol=K)
  
  if(!is.null(Vinv)){
    ##Initialize ELBO parameters
    var_part_ERSS <- 0
    neg_KL <- 0
  }

  ##Loop through the variables
  for(j in 1:p){
    
    #Remove j-th effect from expected residuals 
    rbar_j <- rbar + outer(X[, j], mu1[j, ])
    
    #Run Bayesian SLR
    bfit <- bayes_mvr_mix_scaled_X(X[, j], rbar_j, w0, S0, S, S1, SplusS0_chol, S_chol)
    
    #Update variational parameters
    mu1c[j, ]         <- bfit$mu1
    S1c[, , j]        <- bfit$S1
    w1c[j, ]          <- bfit$w1
    
    #Compute ELBO params
    if(!is.null(Vinv)){
      xtx <- sum(X[, j]^2)
      var_part_ERSS <- var_part_ERSS + (tr(Vinv%*%bfit$S1)*xtx)
      mu1_mat <- matrix(bfit$mu1, ncol=1)
      neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, rbar_j)%*%tcrossprod(matrix(X[, j], ncol=1), mu1_mat))+
                                             tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*(xtx)))
    }
    
    #Update expected residuals
    rbar <- rbar_j - outer(X[, j], mu1c[j, ])
  }
  
  ###Return output
  if(!is.null(Vinv)){
    return(list(rbar=rbar, mu1=mu1c, S1=S1c, w1=w1c, var_part_ERSS=var_part_ERSS, neg_KL=neg_KL))
  } else {
    return(list(rbar=rbar, mu1=mu1c, S1=S1c, w1=w1c))
  }
}

###Perform one iteration of the outer loop
mr_mash_update_scaled_X <- function(Y, X, mu1_t, w1_t, V, Vinv, w0, S0, S, S1, 
                                    SplusS0_chol, S_chol, update_w0, compute_ELBO){
  ##Compute expected residuals
  rbar <- Y - X%*%mu1_t
  
  #Update w0 if requested
  if(update_w0 && !is.null(w1_t)){
    w0 <- update_weights(w1_t)
  }
  
  ##Update variational parameters, expected residuals, and ELBO components
  if(compute_ELBO){
    updates <- inner_loop_scaled_X(X=X, rbar=rbar, mu=mu1_t, Vinv=Vinv, w0=w0, S0=S0, 
                          S=S, S1=S1, SplusS0_chol=SplusS0_chol, S_chol=S_chol) 
  } else {
    updates <- inner_loop_scaled_X(X=X, rbar=rbar, mu=mu1_t, Vinv=NULL, w0=w0, S0=S0,
                          S=S, S1=S1, SplusS0_chol=SplusS0_chol, S_chol=S_chol)
  }
  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  rbar    <- updates$rbar
  
  if(compute_ELBO){
    ##Compute ELBO
    var_part_ERSS <- updates$var_part_ERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, var_part_ERSS=var_part_ERSS, neg_KL=neg_KL)
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, ELBO=ELBO))
  } else {
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t))
  }
}

