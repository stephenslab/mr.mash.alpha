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

###Update mixture weights
update_weights_em <- function(x){
  w <- colSums(x)
  w <- w/sum(w)
  return(w)
}

###Compute intermediate components of the ELBO
compute_ELBO_terms <- function(var_part_tr_wERSS, neg_KL, x_j, rbar_j, bfit, xtx, Vinv){
  mu1_mat <- matrix(bfit$mu1, ncol=1)
  # var_part_tr_wERSS <- var_part_tr_wERSS + (tr(Vinv%*%bfit$S1)*xtx)
  # neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, rbar_j)%*%tcrossprod(matrix(x_j, ncol=1), mu1_mat))+
  #                                        tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*xtx))
  ##Equivalent to the above but more efficient
  var_part_tr_wERSS <- var_part_tr_wERSS + (sum(Vinv*bfit$S1)*xtx)
  neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*sum(tcrossprod(Vinv, rbar_j)*t(tcrossprod(matrix(x_j, ncol=1), mu1_mat)))+
                                         sum(Vinv*(bfit$S1+tcrossprod(mu1_mat)))*xtx))
  
  
  return(list(var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
}

###Compute ELBO from intermediate components
compute_ELBO_fun <- function(rbar, V, Vinv, ldetV, var_part_tr_wERSS, neg_KL){
  n <- nrow(rbar)
  R <- ncol(rbar)
  # tr_wERSS <- tr(Vinv%*%(crossprod(rbar))) + var_part_tr_wERSS
  tr_wERSS <- sum(Vinv*(crossprod(rbar))) + var_part_tr_wERSS
  ELBO <- -log(n)/2 - (n*R)/2*log(2*pi) - n/2 * ldetV - 0.5*tr_wERSS + neg_KL
  
  return(ELBO)
}

###Compute log-determinant from Cholesky decomposition
chol2ldet <- function(R){
  logdet <- 2*sum(log(diag(R)))
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

###Compute quantities needed when using centered X
precompute_quants_centered_X <- function(X, V, S0){
  ###Quantities that don't depend on S0
  #xtx <- diag(crossprod(X))
  xtx <- colSums(X^2)
  R <- chol(V)
  #Rtinv <- solve(t(R))
  #Rinv <- solve(R)
  Rtinv <- forwardsolve(t(R), diag(nrow(R)))
  Rinv <- backsolve(R, diag(nrow(R)))
  
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

###Precompute quantities in any case
precompute_quants <- function(n, X, V, S0, standardize, version){
  if(standardize){
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
    
    if(version=="R"){
      return(list(V_chol=R, S=S, S1=S1, S_chol=S_chol, SplusS0_chol=SplusS0_chol))      
    } else if(version=="Rcpp"){
      xtx <- c(0, 0) ##Vector
      U0 <- array(0, c(1, 1, 1))
      d <- matrix(0, nrow=1, ncol=1)
      Q <- array(0, c(1, 1, 1))
      
      return(list(V_chol=R, S=S, S1=simplify2array(S1), S_chol=S_chol, SplusS0_chol=simplify2array(SplusS0_chol), 
                  xtx=xtx, U0=U0, d=d, Q=Q))
    }
    
  } else {
    ###Quantities that don't depend on S0
    #xtx <- diag(crossprod(X))
    xtx <- colSums(X^2)
    R <- chol(V)
    #Rtinv <- solve(t(R))
    #Rinv <- solve(R)
    Rtinv <- forwardsolve(t(R), diag(nrow(R)))
    Rinv <- backsolve(R, diag(nrow(R)))
    
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
    
    if(version=="R"){
      return(list(xtx=xtx, V_chol=R, U0=U0, d=d, Q=Q))
    } else if(version=="Rcpp"){
      S <- matrix(0, nrow=1, ncol=1)
      S1 <- array(0, c(1, 1, 1))
      S_chol <- matrix(0, nrow=1, ncol=1)
      SplusS0_chol <- array(0, c(1, 1, 1))
      
      return(list(xtx=xtx, V_chol=R, U0=simplify2array(U0), d=simplify2array(d), Q=simplify2array(Q), 
                  S=S, S1=S1, S_chol=S_chol, SplusS0_chol=SplusS0_chol))
    }
  }
}

###Update mixture weights with mixsqp
#' @importFrom mixsqp mixsqp
#' 
compute_mixsqp_update <- function (X, Y, V, S0, mu1_t, precomp_quants, standardize) {
  
  # Get the number of predictors (p), the number of mixture
  # components in the prior (k), and the number of samples (n).
  p <- ncol(X)
  K <- length(S0)
  n <- nrow(Y)
  
  # Compute the p x k matrix of log-likelihoods conditional on each
  # prior mixture component.
  L <- matrix(0, p, K)
  rbar <- Y - X%*%mu1_t
  
  if(standardize){
    for(j in 1:p){
      # Compute the least-squares estimate.
      b <- drop(X[, j] %*% Y)/(n-1)
      
      #Remove j-th effect from expected residuals 
      rbar_j <- rbar + outer(X[, j], mu1_t[j, ])
      
      for(k in 1:K){
        L[j, k] <- bayes_mvr_ridge_scaled_X(b, precomp_quants$S0[[k]], precomp_quants$S, 
                                            precomp_quants$S1[[k]], precomp_quants$SplusS0_chol[[k]], precomp_quants$S_chol)$logbf
      }
    }
  } else {
    for(j in 1:p){
      # Compute the least-squares estimate and covariance.
      b <- drop(X[, j] %*% Y)/precomp_quants$xtx[j]
      S <- V/precomp_quants$xtx[j]
      
      # Compute quantities needed for bayes_mvr_ridge_centered_X()
      S_chol <- precomp_quants$V_chol/sqrt(precomp_quants$xtx[j])
      
      #Remove j-th effect from expected residuals 
      rbar_j <- rbar + outer(X[, j], mu1_t[j, ])
      
      for(k in 1:K){
        L[j, k] <- bayes_mvr_ridge_centered_X(V, b, S, S0[[k]], precomp_quants$xtx[j], 
                                              precomp_quants$V_chol, S_chol, precomp_quants$U0[[k]], precomp_quants$d[[k]], 
                                              precomp_quants$Q[[k]])$logbf
      }
    }
  }
  
  out <- mixsqp(L,log = TRUE,control = list(verbose = FALSE))
  if (out$status != "converged to optimal solution")
    warning("mixsqp did not converge to optimal solution")
  
  # Return the updated mixture weights ("w0") and the number of
  # mix-SQP iterations performed.
  return(list(w0=out$x, numiter=nrow(out$progress)))
}

# Perform backtracking line search to identify a step size for the
# mixture weights update that increases the ELBO.
backtracking_line_search <- function (X, Y, V, Vinv, ldetV, S0, mu1_t, w0em, w0mixsqp,
                                      precomp_quants, standardize, compute_ELBO, update_V, 
                                      version, stepsize.reduce, stepsize.min) {
  
  # Compute the objective (ELBO) at the current iterate.
  ##Update variational parameters, expected residuals, and ELBO components
  rbar <- Y - X%*%mu1_t
  if(version=="R"){
    updates <- inner_loop_general(X=X, rbar=rbar, mu=mu1_t, V=V, Vinv=Vinv, w0=w0em, S0=S0, 
                                  precomp_quants=precomp_quants, standardize=standardize, 
                                  compute_ELBO=compute_ELBO, update_V=update_V)
  } else if(version=="Rcpp"){
    updates <- inner_loop_general_rcpp_wrapper(X=X, Rbar=rbar, mu1=mu1_t, V=V, Vinv=Vinv, w0=w0em, 
                                               S0=simplify2array(S0), precomp_quants=precomp_quants, 
                                               standardize=standardize, compute_ELBO=compute_ELBO, 
                                               update_V=update_V)
  }
  f <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=updates$var_part_tr_wERSS, neg_KL=updates$neg_KL)
    
  
  # Perform backtracking line search to identify a step size that
  # increases the ELBO.
  a <- 1
  while (TRUE) {
    w0new <- a*w0mixsqp + (1 - a)*w0em
    if(version=="R"){
      updates <- inner_loop_general(X=X, rbar=rbar, mu=mu1_t, V=V, Vinv=Vinv, w0=w0new, S0=S0, 
                                    precomp_quants=precomp_quants, standardize=standardize,
                                    compute_ELBO=compute_ELBO, update_V=update_V) 
    } else if(version=="Rcpp"){
      updates <- inner_loop_general_rcpp_wrapper(X=X, Rbar=rbar, mu1=mu1_t, V=V, Vinv=Vinv, w0=w0new, 
                                                 S0=simplify2array(S0), precomp_quants=precomp_quants, 
                                                 standardize=standardize, compute_ELBO=compute_ELBO, 
                                                 update_V=update_V)
    }
    fnew <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=updates$var_part_tr_wERSS, neg_KL=updates$neg_KL)
    
    # Check whether the new candidate increases the ELBO.
    if (fnew > f)
      break
    
    # If we cannot decrease the step size further, terminate the
    # backtracking line search, and set the step size to be the
    # minimum step size.
    else if (a * stepsize.reduce < stepsize.min) {
      
      # We need to terminate backtracking line search because we have
      # arrived at the smallest allowable step size.
      a     <- 0
      w0new <- w0em
      break
    }
    
    # The new candidate does not increase the ELBO, so we need to try
    # again with a smaller step size.
    a = a * stepsize.reduce
  }
  
  # Return the updated mixture weights ("w0") and the step size
  # determined by the backtracking line search ("a").
  return(list(w0 = w0new,a = a))
}



# Update the mixture weights with mix-SQP, following by backtracking
# line search to ensure that the ELBO does not decrease.
update_weights_mixsqp <- function (X, Y, mu1_t, V, Vinv, ldetV, w0em, S0,
                                   precomp_quants, update_w0, standardize,
                                   compute_ELBO=TRUE, update_V=FALSE, version,
                                   stepsize.reduce = 0.5, stepsize.min = 1e-8) {
  
  # Compute the mix-SQP update for the mixture weights. Note that this
  # update is not guaranteed to increase the ELBO.
  out1 <- compute_mixsqp_update(X, Y, V, S0, mu1_t, precomp_quants, standardize)
  
  # Perform backtracking line search to identify a step size that
  # increases the ELBO.
  out2 <- backtracking_line_search(X, Y, V, Vinv, ldetV, S0, mu1_t, w0em, out1$w0,
                                   precomp_quants, standardize, compute_ELBO, update_V, 
                                   version, stepsize.reduce, stepsize.min)
  
  # Return the updated mixture weights ("w0"), the number of mix-SQP
  # iterations performed ("numiter"), and the step size determined by the
  # backtracking line search ("a").
  return(list(w0      = out2$w0,
              numiter = out1$numiter,
              a       = out2$a))
}

###Compute variance part of the ERSS
compute_var_part_ERSS <- function(var_part_ERSS, bfit, xtx){
  var_part_ERSS <- var_part_ERSS + (bfit$S1*xtx)
  
  return(var_part_ERSS)
}

###Update V
update_V_fun <- function(Y, X, mu1_t, var_part_ERSS){
  rbar <- Y - X%*%mu1_t
  ERSS <- crossprod(rbar) + var_part_ERSS
  n <- nrow(X)
  V <- ERSS/n
  
  return(V)
}
