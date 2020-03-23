# Bayesian multivariate regression with Normal prior 
#
# The outputs are: b, the least-squares estimate of the regression
# coefficients; S, the covariance of b; mu1, the posterior mean of the
# regression coefficients; S1, the posterior covariance of the
# regression coefficients; logbf, the log-Bayes factor.
bayes_mvr_ridge <- function (x, Y, V, S0) {
  
  # Compute the least-squares estimate and its covariance.
  b <- drop(x %*% Y)/sum(x^2)
  S <- V/sum(x^2)
  
  # Compute the log-Bayes factor.
  # logbf <- mvtnorm::dmvnorm(x=b, sigma=(S+S0), log=T) - mvtnorm::dmvnorm(x=b, sigma=S, log=T)  ##Slow
  # logbf <- (log(prod(abs(Re(diag(qr(S)$qr))))) +
  #             - log(prod(abs(Re(diag(qr(S0+S)$qr)))))
  #           + dot(b,solve(S,b)) - dot(b,solve(S0 + S,b)))/2  ##Not as fast as with determinant() but more stable
  logbf <- (as.numeric(determinant(S)$modulus) +
              - as.numeric(determinant(S0 + S)$modulus)
            + dot(b,solve(S,b)) - dot(b,solve(S0 + S,b)))/2
  
  # Compute the posterior mean and covariance assuming a multivariate
  # normal prior with zero mean and covariance S0.
  # r   <- ncol(Y)
  # I   <- diag(r)
  # S1  <- solve(solve(S0) + solve(S))
  # S1 <- S0%*%solve(S+S0)%*%S
  #Avoid inverting matrices
  SplusS0_chol <- chol(S+S0)
  S1 <- S0%*%backsolve(SplusS0_chol, forwardsolve(t(SplusS0_chol), S))
  # mu1 <- solve(S %*% solve(S0) + I,b)
  # mu1 <- drop(S1%*%solve(S)%*%b)
  #Avoid inverting matrices
  S_chol <- chol(S)
  mu1 <- drop(S1%*%backsolve(S_chol, forwardsolve(t(S_chol), b)))
  
  # Return the least-squares estimate and covariance (b, S), the posterior mean and
  # covariance (mu1, S1), and the log-Bayes factor (logbf)
  return(list(b = b,S = S,mu1 = mu1,S1 = S1,logbf = logbf))
}


# Bayesian multivariate regression with Normal prior with standardized X
#
# The outputs are: mu1, the posterior mean of the
# regression coefficients; logbf, the log-Bayes factor.
bayes_mvr_ridge_standardized_X <- function (b, S0, S, S1, SplusS0_chol, S_chol) {
  
  # Compute the log-Bayes factor.
  logbf <- (chol2ldet(S_chol) - chol2ldet(SplusS0_chol) +
              dot(b,backsolve(S_chol, forwardsolve(t(S_chol), b))) - 
              dot(b,backsolve(SplusS0_chol, forwardsolve(t(SplusS0_chol), b))))/2
  
  # Compute the posterior mean assuming a multivariate
  # normal prior with zero mean and covariance S0.
  mu1 <- drop(S1%*%backsolve(S_chol, forwardsolve(t(S_chol), b)))
  
  # Return the posterior mean
  # (mu1), and the log-Bayes factor (logbf)
  return(list(mu1 = mu1, logbf = logbf))
}


# Bayesian multivariate regression with Normal prior with transformed 
# xtilde = x%*%solve(chol(V)) [this calculation is not needed but useful
# to understand the derivation] that allows to precompute some
# quantities
#
# The outputs are: mu1, the posterior mean of the
# regression coefficients; S1, the posterior covariance of the
# regression coefficients; logbf, the log-Bayes factor.
bayes_mvr_ridge_centered_X <- function (V, b, S, S0, xtx, Vinv, V_chol, S_chol, d, QtimesV_chol) {
  
  # Compute the log-Bayes factor.
  SplusS0_chol <- chol(S+S0)
  logbf <- (chol2ldet(S_chol) - chol2ldet(SplusS0_chol) +
              dot(b,backsolve(S_chol, forwardsolve(t(S_chol), b))) -
              dot(b,backsolve(SplusS0_chol, forwardsolve(t(SplusS0_chol), b))))/2
  #logbf <- dmvnorm(b, rep(0, times=length(b)), (S+S0)) - dmvnorm(b, rep(0, times=length(b)), S)
  
  # Compute the posterior mean assuming a multivariate
  # normal prior with zero mean and covariance S0.
  dx <- d/(1 + xtx*d)
  A <- sqrt(dx)*QtimesV_chol
  S1 <- crossprod(A)
  mu1 <- drop(crossprod(A, (A %*% (Vinv %*% (xtx*b)))))
  # Return the posterior mean and covariance
  # (mu1, S1), and the log-Bayes factor (logbf)
  return(list(mu1 = mu1, S1=S1, logbf=logbf))
}


# Bayesian multivariate regression with mixture-of-normals prior
# (mixture weights w0 and covariance matrices S0) 
#
# The outputs are: the log-Bayes factor (logbf), the posterior assignment probabilities
# (w1), the posterior mean of the coefficients given that all the
# coefficients are not nonzero (mu1), and the posterior covariance of
# the coefficients given that all the coefficients are not zero (S1).
bayes_mvr_mix <- function (x, Y, V, w0, S0) {
  
  # Get the number of variables (n) and the number of mixture
  # components (k).
  r <- ncol(Y)
  K <- length(S0)
  
  # Compute the Bayes factors and posterior statistics separately for
  # each mixture component.
  # out <- vector("list",K)
  # for (k in 1:K){
  #   out[[k]] <- bayes_mvr_ridge(x,Y,V,S0[[k]])
  # }
  bayes_mvr_ridge_lapply <- function(i){
    bayes_mvr_ridge(x, Y, V, S0[[i]])
  }
  out <- lapply(1:K, bayes_mvr_ridge_lapply)
  
  # Compute the posterior assignment probabilities for the latent
  # indicator variable.
  logbf <- sapply(out,function (x) x$logbf)
  w1    <- softmax(logbf + log(w0))
  
  # Compute the posterior mean (mu1_mix) and covariance (S1_mix) of the
  # regression coefficients.
  A   <- matrix(0,r,r)
  mu1_mix <- rep(0,r)
  for (k in 1:K) {
    wk  <- w1[k]
    muk <- out[[k]]$mu1
    Sk  <- out[[k]]$S1
    mu1_mix <- mu1_mix + wk*muk
    A   <- A   + wk*(Sk + tcrossprod(muk))
  }
  S1_mix <- A - tcrossprod(mu1_mix)
  ##The following code does not work in the univariate case
  # muk <- t(sapply(out, function (x) x$mu1))
  # Sk <- lapply(out, function (x) x$S1)
  # mu1_mix <- colSums(muk*w1)
  # S1_mix <- Reduce("+", lapply(1:K, function(i){w1[i]*(Sk[[i]] + tcrossprod(muk[i, ]))})) - tcrossprod(mu1_mix)
  
  # Compute the log-Bayes factor for the mixture as a linear combination of the
  # individual BFs foreach mixture component.
  u     <- max(logbf)
  logbf_mix <- u + log(sum(w0 * exp(logbf - u)))
  
  # Return the log-Bayes factor for the mixture (logbf), the posterior assignment probabilities (w1), the
  # posterior mean of the coefficients (mu1), and the posterior
  # covariance of the coefficients (S1).
  return(list(logbf = logbf_mix,w1 = w1,mu1 = mu1_mix,S1 = S1_mix))
}

                  
# Bayesian multivariate regression with mixture-of-normals prior
# (mixture weights w0 and covariance matrices S0) with standardized X.
#
# The outputs are: the log-Bayes factor (logbf), the posterior assignment probabilities
# (w1), the posterior mean of the coefficients given that all the
# coefficients are not nonzero (mu1), and the posterior covariance of
# the coefficients given that all the coefficients are not zero (S1).
bayes_mvr_mix_standardized_X <- function (x, Y, w0, S0, S, S1, SplusS0_chol, S_chol) {
  
  
  # Get the number of conditions (r), the number of mixture
  # components (K), and the number of samples (n).
  r <- ncol(Y)
  K <- length(S0)
  n <- nrow(Y)
  
  # Compute the least-squares estimate.
  b <- drop(x %*% Y)/(n-1)
  
  # Compute the Bayes factors and posterior statistics separately for
  # each mixture component.
  # out <- vector("list",K)
  # for (k in 1:K){
  #   out[[k]] <- bayes_mvr_ridge_standardized_X(b, S0[[k]], S, S1[[k]], SplusS0_chol[[k]], S_chol)
  # }
  bayes_mvr_ridge_lapply <- function(i){
    bayes_mvr_ridge_standardized_X(b, S0[[i]], S, S1[[i]], SplusS0_chol[[i]], S_chol)
  }
  out <- lapply(1:K, bayes_mvr_ridge_lapply)
  
  # Compute the posterior assignment probabilities for the latent
  # indicator variable.
  logbf <- sapply(out,function (x) x$logbf)
  w1    <- softmax(logbf + log(w0))
  
  # Compute the posterior mean (mu1_mix) and covariance (S1_mix) of the
  # regression coefficients.
  A   <- matrix(0,r,r)
  mu1_mix <- rep(0,r)
  for (k in 1:K) {
    wk  <- w1[k]
    muk <- out[[k]]$mu1
    Sk  <- S1[[k]]
    mu1_mix <- mu1_mix + wk*muk
    A   <- A   + wk*(Sk + tcrossprod(muk))
  }
  S1_mix <- A - tcrossprod(mu1_mix)
  ##The following code does not work in the univariate case
  # muk <- t(sapply(out, function (x) x$mu1))
  # Sk <- lapply(out, function (x) x$S1)
  # mu1_mix <- colSums(muk*w1)
  # S1_mix <- Reduce("+", lapply(1:K, function(i){w1[i]*(Sk[[i]] + tcrossprod(muk[i, ]))})) - tcrossprod(mu1_mix)
  
  # Compute the log-Bayes factor for the mixture as a linear combination of the
  # individual BFs foreach mixture component.
  u     <- max(logbf)
  logbf_mix <- u + log(sum(w0 * exp(logbf - u)))
  
  # Return the log-Bayes factor for the mixture (logbf), the posterior assignment probabilities (w1), the
  # posterior mean of the coefficients (mu1), and the posterior
  # covariance of the coefficients (S1).
  return(list(logbf = logbf_mix,w1 = w1,mu1 = mu1_mix,S1 = S1_mix))
}


# Bayesian multivariate regression with mixture-of-normals prior
# (mixture weights w0 and covariance matrices S0) and centered X
#
# The outputs are: the log-Bayes factor (logbf), the posterior assignment probabilities
# (w1), the posterior mean of the coefficients given that all the
# coefficients are not nonzero (mu1), and the posterior covariance of
# the coefficients given that all the coefficients are not zero (S1).
bayes_mvr_mix_centered_X <- function (x, Y, V, w0, S0, xtx, Vinv, V_chol, d, QtimesV_chol) {
  
  # Get the number of variables (n) and the number of mixture
  # components (k).
  r <- ncol(Y)
  K <- length(S0)
  
  # Compute the least-squares estimate and covariance.
  b <- drop(x %*% Y)/xtx
  S <- V/xtx
  
  # Compute quantities needed for bayes_mvr_ridge_centered_X()
  S_chol <- V_chol/sqrt(xtx)
  
  # Compute the Bayes factors and posterior statistics separately for
  # each mixture component.
  # out <- vector("list",K)
  # for (k in 1:K){
  #   out[[k]] <- bayes_mvr_ridge_centered_X(V, b, S, S0[[k]], xtx, V_chol, d[[k]], QtimesV_chol[[k]])
  # }
  bayes_mvr_ridge_lapply <- function(i){
    bayes_mvr_ridge_centered_X(V, b, S, S0[[i]], xtx, Vinv, V_chol, S_chol, d[[i]], QtimesV_chol[[i]])
  }
  out <- lapply(1:K, bayes_mvr_ridge_lapply)
  
  # Compute the posterior assignment probabilities for the latent
  # indicator variable.
  logbf <- sapply(out,function (x) x$logbf)
  w1    <- softmax(logbf + log(w0))
  
  # Compute the posterior mean (mu1_mix) and covariance (S1_mix) of the
  # regression coefficients.
  A   <- matrix(0,r,r)
  mu1_mix <- rep(0,r)
  for (k in 1:K) {
    wk  <- w1[k]
    muk <- out[[k]]$mu1
    Sk  <- out[[k]]$S1
    mu1_mix <- mu1_mix + wk*muk
    A   <- A   + wk*(Sk + tcrossprod(muk))
  }
  S1_mix <- A - tcrossprod(mu1_mix)
  ##The following code does not work in the univariate case
  # muk <- t(sapply(out, function (x) x$mu1))
  # Sk <- lapply(out, function (x) x$S1)
  # mu1_mix <- colSums(muk*w1)
  # S1_mix <- Reduce("+", lapply(1:K, function(i){w1[i]*(Sk[[i]] + tcrossprod(muk[i, ]))})) - tcrossprod(mu1_mix)
  
  # Compute the log-Bayes factor for the mixture as a linear combination of the
  # individual BFs foreach mixture component.
  u     <- max(logbf)
  logbf_mix <- u + log(sum(w0 * exp(logbf - u)))
  
  # Return the log-Bayes factor for the mixture (logbf), the posterior assignment probabilities (w1), the
  # posterior mean of the coefficients (mu1), and the posterior
  # covariance of the coefficients (S1).
  return(list(logbf = logbf_mix,w1 = w1,mu1 = mu1_mix,S1 = S1_mix))
}