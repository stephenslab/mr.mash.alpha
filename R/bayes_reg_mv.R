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
  # R   <- ncol(Y)
  # I   <- diag(R)
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


# Bayesian multivariate regression with Normal prior with scaled X
#
# The outputs are: mu1, the posterior mean of the
# regression coefficients; logbf, the log-Bayes factor.
bayes_mvr_ridge_scaled_X <- function (x, Y, b, S0, S, S1, SplusS0_chol, S_chol, ldetSplusS0_chol, ldetS_chol) {
  
  # Compute the log-Bayes factor.
  logbf <- (ldetS_chol - ldetSplusS0_chol +
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
bayes_mvr_ridge_centered_X <- function (x, Y, V, b, S, S0, xtx, V_chol, S_chol, U0, d, Q) {
  
  # Compute the log-Bayes factor.
  SplusS0_chol <- chol(S+S0)
  logbf <- (chol2ldet(S_chol) - chol2ldet(SplusS0_chol) +
              dot(b,backsolve(S_chol, forwardsolve(t(S_chol), b))) -
              dot(b,backsolve(SplusS0_chol, forwardsolve(t(SplusS0_chol), b))))/2
  #logbf <- dmvnorm(b, rep(0, times=length(b)), (S+S0)) - dmvnorm(b, rep(0, times=length(b)), S)
  
  # Compute the posterior mean assuming a multivariate
  # normal prior with zero mean and covariance S0 (handling the univariate case).
  if(length(d)>1){
    D <- diag(1/(1 + xtx*d))
  } else {
    D <- 1/(1 + xtx*d)
  }
  U1 <- U0 %*% Q %*% tcrossprod(D, Q)
  S1 <- crossprod(V_chol, U1) %*% V_chol
  mu1 <- drop(S1%*%backsolve(S_chol, forwardsolve(t(S_chol), b)))
  
  # Return the posterior mean and covariance
  # (mu1, S1), and the log-Bayes factor (logbf)
  return(list(mu1 = mu1, S1=S1, logbf=logbf))
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
# (mixture weights w0 and covariance matrices S0) 
#
# The outputs are: the log-Bayes factor (logbf), the posterior assignment probabilities
# (w1), the posterior mean of the coefficients given that all the
# coefficients are not nonzero (mu1), and the posterior covariance of
# the coefficients given that all the coefficients are not zero (S1).
bayes_mvr_mix <- function (x, Y, V, w0, S0) {
  
  # Get the number of variables (n) and the number of mixture
  # components (k).
  R <- ncol(Y)
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
  A   <- matrix(0,R,R)
  mu1_mix <- rep(0,R)
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
# (mixture weights w0 and covariance matrices S0) with scaled X.
#
# The outputs are: the log-Bayes factor (logbf), the posterior assignment probabilities
# (w1), the posterior mean of the coefficients given that all the
# coefficients are not nonzero (mu1), and the posterior covariance of
# the coefficients given that all the coefficients are not zero (S1).
bayes_mvr_mix_scaled_X <- function (x, Y, w0, S0, S, S1, SplusS0_chol, S_chol, ldetSplusS0_chol, ldetS_chol) {
  
  
  # Get the number of conditions (R), the number of mixture
  # components (K), and the number of samples (n).
  R <- ncol(Y)
  K <- length(S0)
  n <- nrow(Y)
  
  # Compute the least-squares estimate.
  b <- drop(x %*% Y)/(n-1)
  
  # Compute the Bayes factors and posterior statistics separately for
  # each mixture component.
  # out <- vector("list",K)
  # for (k in 1:K){
  #   out[[k]] <- bayes_mvr_ridge_scaled_X(x, Y, S0[[k]], S, S1[[k]], SplusS0_chol[[k]], S_chol, ldetSplusS0_chol[k], ldetS_chol)
  # }
  bayes_mvr_ridge_lapply <- function(i){
    bayes_mvr_ridge_scaled_X(x, Y, b, S0[[i]], S, S1[[i]], SplusS0_chol[[i]], S_chol, ldetSplusS0_chol[i], ldetS_chol)
  }
  out <- lapply(1:K, bayes_mvr_ridge_lapply)
  
  # Compute the posterior assignment probabilities for the latent
  # indicator variable.
  logbf <- sapply(out,function (x) x$logbf)
  w1    <- softmax(logbf + log(w0))
  
  # Compute the posterior mean (mu1_mix) and covariance (S1_mix) of the
  # regression coefficients.
  A   <- matrix(0,R,R)
  mu1_mix <- rep(0,R)
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
# (mixture weights w0 and covariance matrices S0) and transformed X
#
# The outputs are: the log-Bayes factor (logbf), the posterior assignment probabilities
# (w1), the posterior mean of the coefficients given that all the
# coefficients are not nonzero (mu1), and the posterior covariance of
# the coefficients given that all the coefficients are not zero (S1).
bayes_mvr_mix_centered_X <- function (x, Y, V, w0, S0, xtx, V_chol, U0, d, Q) {
  
  # Get the number of variables (n) and the number of mixture
  # components (k).
  R <- ncol(Y)
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
  #   out[[k]] <- bayes_mvr_ridge_centered_X(x, Y, V, S0[[k]], xtx, V_chol, U0[[k]], d[[k]], Q[[k]])
  # }
  bayes_mvr_ridge_lapply <- function(i){
    bayes_mvr_ridge_centered_X(x, Y, V, b, S, S0[[i]], xtx, V_chol, S_chol, U0[[i]], d[[i]], Q[[i]])
  }
  out <- lapply(1:K, bayes_mvr_ridge_lapply)
  
  # Compute the posterior assignment probabilities for the latent
  # indicator variable.
  logbf <- sapply(out,function (x) x$logbf)
  w1    <- softmax(logbf + log(w0))
  
  # Compute the posterior mean (mu1_mix) and covariance (S1_mix) of the
  # regression coefficients.
  A   <- matrix(0,R,R)
  mu1_mix <- rep(0,R)
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
