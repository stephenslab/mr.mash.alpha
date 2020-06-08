# Should be the same as mvtnorm::dmvnorm(x,mu,S,log = TRUE)
#
# #' @importFrom Rcpp evalCpp
# #' @useDynLib mr.mash.alpha
#' 
dmvnorm <- function (x, mu, S)
  dmvnorm_rcpp(x,mu,S)

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
  d <- list()
  QtimesR <- list()
  for(i in 1:length(S0)){
    U0 <- Rtinv %*% S0[[i]] %*% Rinv
    out <- eigen(U0)
    d[[i]]   <- out$values
    QtimesR[[i]]   <- crossprod(out$vectors, R)   
  }
  
  return(list(xtx=xtx, V_chol=R, d=d, QtimesV_chol=QtimesR))
}




# Bayesian multivariate regression with spike-and-slab prior 
#
# The outputs are: b, the least-squares estimate of the regression
# coefficients; S, the covariance of b; mu1, the posterior mean of the
# regression coefficients given that the coefficients are not all
# zero; S1, the posterior covariance of the regression coefficients
# given that the coefficients are not all zero; and p1, the posterior
# probability that the coefficients are not all zero; logbf,
# the log-Bayes factor.
#
# Input argument p0 specifies the prior probability that the
# coefficients are not all zero.
bayes_mvr_spike_slab <- function (x, Y, V, S0, p0) {
  
  # Compute the least-squares estimate and its covariance.
  f <- bayes_mvr_ridge(x, Y, V, S0)
  
  # Compute the posterior probability that the coefficient is nonzero.
  p1 <- sigmoid(log(p0/(1 - p0)) + f$logbf)
  
  # Return the least-squares estimate (b, S), the posterior mean and
  # standard deviation (mu1, S1), the log-Bayes factor (logbf), and the
  # posterior inclusion probability (p1).
  return(list(b = f$b,S = f$S,mu1 = f$mu1,S1 = f$S1,logbf = f$logbf,p1 = p1))
}

# Bayesian multivariate regression with mixture-of-normals prior
# (mixture weights w0 and covariance matrices S0) using MASH
# TO DO: Move MashInitializer$new outside the function when using in mr.mash
# because it only needs to be done once!                
#
# The outputs are: the log-Bayes factor (logbf), the posterior
# assignment probabilities (w1), the posterior mean of the
# coefficients given that all the coefficients are not nonzero (mu1),
# and the posterior covariance of the coefficients given that all the
# coefficients are not zero (S1).
bayes_mvr_mash <- function(x, Y, V, w0, S0){
  if(!is.matrix(x)){x <- matrix(x, ncol=1)}
  data <- mmbr:::DenseData$new(x, Y)
  data$standardize(FALSE, FALSE)
  mash_init <- mmbr:::MashInitializer$new(S0, grid=1, prior_weights=w0, null_weight=0, top_mixtures=-1)
  B <- mmbr:::MashRegression$new(1, V, mash_init)
  B$fit(data, save_var=TRUE)
  
  return(list(mu1=drop(B$posterior_b1), S1=drop(B$posterior_variance), w1=drop(B$mixture_posterior_weights[, -1]), logbf=B$lbf))
}


###Update variational parameters, expected residuals, and ELBO components
inner_loop <- function(X, rbar, mu, V, Vinv, w0, S0, xtx, V_chol, d, QtimesV_chol){
  ###Create variables to store quantities
  R <- ncol(rbar)
  p <- ncol(X)
  K <- length(S0)
  mu1   <- mu
  S1    <- array(0, c(R, R, p))
  w1    <- matrix(0, nrow=p, ncol=K)
  
  if(!is.null(Vinv)){
    ##Initialize ELBO parameters
    var_part_tr_wERSS <- 0
    neg_KL <- 0
  }
  
  ##Loop through the variables
  for(j in 1:p){
    
    #Remove j-th effect from expected residuals 
    rbar_j <- rbar + outer(X[, j], mu1[j, ])
    
    #Run Bayesian SLR
    bfit <- bayes_mvr_mix_centered_X(X[, j], rbar_j, V, w0, S0, xtx[j], V_chol, d, QtimesV_chol)
    
    #Update variational parameters
    mu1[j, ]         <- bfit$mu1
    S1[, , j]        <- bfit$S1
    w1[j, ]          <- bfit$w1
    
    #Compute ELBO params
    if(!is.null(Vinv)){
      var_part_tr_wERSS <- var_part_tr_wERSS + (tr(Vinv%*%bfit$S1)*xtx[j])
      mu1_mat <- matrix(bfit$mu1, ncol=1)
      neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, rbar_j)%*%tcrossprod(matrix(X[, j], ncol=1), mu1_mat))+
                                             tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*(xtx[j])))
    }
    
    #Update expected residuals
    rbar <- rbar_j - outer(X[, j], mu1[j, ])
  }
  
  ###Return output
  if(!is.null(Vinv)){
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
  } else {
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1))
  }
}

###Perform one iteration of the outer loop
mr_mash_update <- function(Y, X, mu1_t, w1_t, V, Vinv, ldetV, w0, S0, 
                           xtx, V_chol, d, QtimesV_chol, update_w0, compute_ELBO){
  ##Compute expected residuals
  rbar <- Y - X%*%mu1_t
  
  #Update w0 if requested
  if(update_w0 && !is.null(w1_t)){
    w0 <- update_weights(w1_t)
  }
  
  ##Update variational parameters, expected residuals, and ELBO components
  updates <- inner_loop(X=X, rbar=rbar, mu=mu1_t, V=V, Vinv=Vinv, w0=w0, S0=S0, xtx=xtx, V_chol=V_chol, d=d, QtimesV_chol=QtimesV_chol) 
  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  rbar    <- updates$rbar
  
  if(compute_ELBO){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, ELBO=ELBO))
  } else {
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t))
  }
}

###Update variational parameters, expected residuals, and ELBO components with scaled X
inner_loop_scaled_X <- function(X, rbar, mu, Vinv, w0, S0, S, S1, SplusS0_chol, S_chol){
  ###Create variables to store quantities
  n <- nrow(rbar)
  R <- ncol(rbar)
  p <- ncol(X)
  K <- length(S0)
  mu1c   <- mu
  S1c    <- array(0, c(R, R, p))
  w1c    <- matrix(0, nrow=p, ncol=K)
  
  if(!is.null(Vinv)){
    ##Initialize ELBO parameters
    var_part_tr_wERSS <- 0
    neg_KL <- 0
  }
  
  ##Loop through the variables
  for(j in 1:p){
    
    #Remove j-th effect from expected residuals 
    rbar_j <- rbar + outer(X[, j], mu1c[j, ])
    
    #Run Bayesian SLR
    bfit <- bayes_mvr_mix_scaled_X(X[, j], rbar_j, w0, S0, S, S1, SplusS0_chol, S_chol)
    
    #Update variational parameters
    mu1c[j, ]         <- bfit$mu1
    S1c[, , j]        <- bfit$S1
    w1c[j, ]          <- bfit$w1
    
    #Compute ELBO params
    if(!is.null(Vinv)){
      xtx <- n-1
      var_part_tr_wERSS <- var_part_tr_wERSS + (tr(Vinv%*%bfit$S1)*xtx)
      mu1_mat <- matrix(bfit$mu1, ncol=1)
      neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, rbar_j)%*%tcrossprod(matrix(X[, j], ncol=1), mu1_mat))+
                                             tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*(xtx)))
    }
    
    #Update expected residuals
    rbar <- rbar_j - outer(X[, j], mu1c[j, ])
  }
  
  ###Return output
  if(!is.null(Vinv)){
    return(list(rbar=rbar, mu1=mu1c, S1=S1c, w1=w1c, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
  } else {
    return(list(rbar=rbar, mu1=mu1c, S1=S1c, w1=w1c))
  }
}

###Perform one iteration of the outer loop with scaled X
mr_mash_update_scaled_X <- function(Y, X, mu1_t, w1_t, V, Vinv, ldetV, w0, S0, S, S1, 
                                    SplusS0_chol, S_chol, update_w0, compute_ELBO){
  ##Compute expected residuals
  rbar <- Y - X%*%mu1_t
  
  #Update w0 if requested
  if(update_w0 && !is.null(w1_t)){
    w0 <- update_weights(w1_t)
  }
  
  ##Update variational parameters, expected residuals, and ELBO components
  updates <- inner_loop_scaled_X(X=X, rbar=rbar, mu=mu1_t, Vinv=Vinv, w0=w0, S0=S0, 
                                 S=S, S1=S1, SplusS0_chol=SplusS0_chol, S_chol=S_chol) 
  
  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  rbar    <- updates$rbar
  
  if(compute_ELBO){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, ELBO=ELBO))
  } else {
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t))
  }
}

