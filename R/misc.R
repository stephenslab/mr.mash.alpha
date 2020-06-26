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

# Add b[i] to each column A[,i].
addtocols <- function (A, b)
  t(t(A) + b)

###Function to simulate from MN distribution
#
#' @importFrom MBSP matrix.normal
#' 
sim_mvr <- function (X, B, V) {
  
  # Get the number of samples (n) and conditions (m).
  n <- nrow(X)
  r <- ncol(B)
  
  # Simulate the responses, Y.
  M <- X%*%B
  U <- diag(n)
  Y <- matrix.normal(M, U, V)
  
  # Output the simulated responses.
  return(Y)
}

###Compute log-determinant from Cholesky decomposition
chol2ldet <- function(R){
  logdet <- 2*sum(log(diag(R)))
  return(logdet)
}

###Similar to base::simplify2array but returns appropriate output when r=1
simplify2array_custom <- function (x, higher = TRUE) {
  common.len <- unique(lengths(x))
  if (common.len >= 1L){
    n <- length(x)
    r <- unlist(x, recursive = FALSE, use.names = FALSE)
    if (higher && length(c.dim <- unique(lapply(x, dim))) == 
        1 && is.numeric(c.dim <- c.dim[[1L]]) && prod(d <- c(c.dim, n)) == length(r)) {
      iN1 <- is.null(n1 <- dimnames(x[[1L]]))
      n2 <- names(x)
      dnam <- if (!(iN1 && is.null(n2))) 
        c(if (iN1) rep.int(list(n1), length(c.dim)) else n1, 
          list(n2))
      array(r, dim = d, dimnames = dnam)
    }
    else if (prod(d <- c(common.len, n)) == length(r)) 
      array(r, dim = d, dimnames = if (!(is.null(n1 <- names(x[[1L]])) & 
                                         is.null(n2 <- names(x)))) 
        list(n1, n2))
    else x
  } else {
    x
  }
}

###Add small number e to diagonal elements to a matrix 
makePD <- function(S0, e){
  S0_PD <- S0+(diag(nrow(S0))*e)
  
  return(S0_PD)
}

###Precompute quantities in any case
precompute_quants <- function(X, V, S0, standardize, version){
  if(standardize){
    n <- nrow(X)
    xtx <- n-1
    
    ###Quantities that don't depend on S0
    R <- chol(V)
    S <- V/xtx
    S_chol <- R/sqrt(xtx)
    
    ###Quantities that depend on S0
    SplusS0_chol <- list()
    S1 <- list()
    for(i in 1:length(S0)){
      SplusS0_chol[[i]] <- chol(S+S0[[i]])
      S1[[i]] <- S0[[i]]%*%backsolve(SplusS0_chol[[i]], forwardsolve(t(SplusS0_chol[[i]]), S))
    }
    
    if(version=="R"){
      return(list(V_chol=R, S=S, S1=S1, S_chol=S_chol, SplusS0_chol=SplusS0_chol))      
    } else if(version=="Rcpp"){
      xtx <- c(0, 0) ##Vector
      d <- matrix(0, nrow=1, ncol=1)
      QtimesR <- array(0, c(1, 1, 1))
      
      return(list(V_chol=R, S=S, S1=simplify2array_custom(S1), S_chol=S_chol, SplusS0_chol=simplify2array_custom(SplusS0_chol), 
                  xtx=xtx, d=d, QtimesV_chol=QtimesR))
    }
    
  } else {
    ###Quantities that don't depend on S0
    R <- chol(V)
    #Rtinv <- solve(t(R))
    #Rinv <- solve(R)
    Rtinv <- forwardsolve(t(R), diag(nrow(R)))
    Rinv <- backsolve(R, diag(nrow(R)))
    
    ###Quantities that depend on S0
    d <- list()
    QtimesR <- list()
    for(i in 1:length(S0)){
      U0  <- Rtinv %*% S0[[i]] %*% Rinv
      out <- eigen(U0)
      d[[i]]   <- out$values
      QtimesR[[i]]   <- crossprod(out$vectors, R)   
    }
    
    if(version=="R"){
      return(list(V_chol=R, d=d, QtimesV_chol=QtimesR))
    } else if(version=="Rcpp"){
      S <- matrix(0, nrow=1, ncol=1)
      S1 <- array(0, c(1, 1, 1))
      S_chol <- matrix(0, nrow=1, ncol=1)
      SplusS0_chol <- array(0, c(1, 1, 1))
      
      return(list(V_chol=R, d=simplify2array_custom(d), QtimesV_chol=simplify2array_custom(QtimesR), 
                  S=S, S1=S1, S_chol=S_chol, SplusS0_chol=SplusS0_chol))
    }
  }
}

###Filter out quantities corresponding to components that are dropped by w0_threshold
filter_precomputed_quants <- function(precomp_quants, to_keep, standardize, version){
  if(standardize){
    if(version=="R"){
      precomp_quants$SplusS0_chol <- precomp_quants$SplusS0_chol[to_keep]
      precomp_quants$S1 <- precomp_quants$S1[to_keep]
    } else if(version=="Rcpp"){
      precomp_quants$SplusS0_chol <- precomp_quants$SplusS0_chol[, , to_keep]
      precomp_quants$S1 <- precomp_quants$S1[, , to_keep]
    }
  } else {
    if(version=="R"){
      precomp_quants$d <- precomp_quants$d[to_keep]
      precomp_quants$QtimesV_chol <- precomp_quants$QtimesV_chol[to_keep]
    } else if(version=="Rcpp"){
      precomp_quants$d <- precomp_quants$d[, to_keep]
      precomp_quants$QtimesV_chol <- precomp_quants$QtimesV_chol[, , to_keep]
    }
  }
  
  return(precomp_quants)
}

###Compute variance part of the ERSS
compute_var_part_ERSS <- function(var_part_ERSS, bfit, xtx){
  var_part_ERSS <- var_part_ERSS + (bfit$S1*xtx)
  
  return(var_part_ERSS)
}

###Rescale posterior mean and covariance of the regression coefficients when standardizing X
rescale_post_mean_covar <- function(mu1, S1, sx){
  p <- nrow(mu1)
  r <- ncol(mu1)
  
  mu1_orig <- mu1/sx
  
  S1_orig <- array(0, c(r, r, p))
  for(j in 1:p){
    S1_orig[, , j] <- S1[, , j]/sx[j]^2
  }
  
  return(list(mu1_orig=mu1_orig, S1_orig=S1_orig))
}

###Faster version of rescale_post_mean_covar()
rescale_post_mean_covar_fast <- function(mu1, S1, sx){
  rescale_post_mean_covar_rcpp(mu1, S1, sx)
}
  

###Scale a matrix (similar to but faster than base::scale())
#' @importFrom matrixStats colSds colMeans2
#' 
scale_fast <- function(M, scale=TRUE, na.rm=TRUE){
  ##Check whether M is a matrix. If not, coerce into it. 
  if(!is.matrix(M))
    M <- as.matrix(M)
  
  ##Store dimnames
  col_names <- colnames(M)
  row_names <- rownames(M)
  
  ###Compute column means and sds
  a <- colMeans2(M, na.rm=na.rm)
  names(a) <- col_names
  if(scale){  
    b <- colSds(M, na.rm=na.rm)
    if(any(b==0))
      stop("Some column(s) have 0 standard deviation")
  } else{
    b <- rep(1, ncol(M))
  }
  names(b) <- col_names
  
  ###Scale
  M <- scale_rcpp(M, a, b)
  
  ###Attach dimension names
  colnames(M) <- col_names
  rownames(M) <- row_names
  
  return(list(M=M, means=a, sds=b))
}

###Scale a matrix (similar to the above but does not use R to compute means and sds)
scale_fast2 <- function(M, scale=TRUE, na.rm=TRUE){
  ##Check whether M is a matrix. If not, coerce into it. 
  if(!is.matrix(M))
    M <- as.matrix(M)
  
  ##Store dimnames
  col_names <- colnames(M)
  row_names <- rownames(M)
  
  ###Scale
  out <- scale2_rcpp(M, scale=scale, na_rm=na.rm)
  means <- drop(out$means)
  sds <- drop(out$sds)
  M <- out$M
  rm(out)
  
  ###Attach dimension names
  colnames(M) <- col_names
  rownames(M) <- row_names
  names(means) <- col_names
  names(sds) <- col_names
  
  return(list(M=M, means=means, sds=sds))
}

###Compute initial estimate of V
compute_V_init <- function(X, Y, B){
  R <- Y - X%*%B
  V <- cov(R)
  
  return(V)
}

