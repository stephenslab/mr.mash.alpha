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
compute_cov_canonical <- function(ntraits, singletons, hetgrid, grid){
  S <- mmbr:::create_cov_canonical(ntraits, singletons, hetgrid)
  U <- list()
  t <- 0
  for(i in 1:length(S)){
    for(j in 1:length(grid)){
      t <- t+1
      U[[t]] <- S[[i]]*grid[j]
    }
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
  ELBO <- -log(n)/2 -log(n*R)/2 -n/2 * as.numeric(determinant(V, logarithm = TRUE)$modulus) - 0.5*ERSS + neg_KL
  
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
