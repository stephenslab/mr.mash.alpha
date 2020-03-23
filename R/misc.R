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
  R <- ncol(B)
  
  # Simulate the responses, Y.
  M <- X%*%B
  U <- diag(n)
  Y <- matrix.normal(M, U, V)
  
  # Output the simulated responses.
  return(Y)
}

###Function to compute canonical covariance matrices scaled by a grid 
compute_cov_canonical <- function(ntraits, singletons, hetgrid, grid, zeromat=TRUE){
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

###Compute log-determinant from Cholesky decomposition
chol2ldet <- function(R){
  logdet <- 2*sum(log(diag(R)))
  return(logdet)
}

###Similar to base::simplify2array but returns appropriate output when R=1
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
