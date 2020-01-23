# Bayesian multivariate regression with Normal prior 
#
# The outputs are: b, the least-squares estimate of the regression
# coefficients; S, the covariance of b; mu1, the posterior mean of the
# regression coefficients; S1, the posterior covariance of the regression coefficients;
# logbf, the log-Bayes factor.
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
  f <- bayes_ridge_mvr(x, Y, V, S0)
  
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
  # A   <- matrix(0,R,R)
  # mu1 <- rep(0,R)
  # for (k in 1:K) {
  #   wk  <- w1[k]
  #   muk <- out[[k]]$mu1
  #   Sk  <- out[[k]]$S1
  #   mu1 <- mu1 + wk*muk
  #   A   <- A   + wk*(Sk + tcrossprod(muk))
  # }
  # S1 <- A - tcrossprod(mu1)
  muk <- t(sapply(out, function (x) x$mu1))
  Sk <- lapply(out, function (x) x$S1)
  mu1_mix <- colSums(muk*w1)
  S1_mix <- Reduce("+", lapply(1:K, function(i){w1[i]*(Sk[[i]] + tcrossprod(muk[i, ]))})) - tcrossprod(mu1_mix)
  
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
#
# The outputs are: the log-Bayes factor (logbf), the posterior assignment probabilities
# (w1), the posterior mean of the coefficients given that all the
# coefficients are not nonzero (mu1), and the posterior covariance of
# the coefficients given that all the coefficients are not zero (S1).
bayes_mvr_mash <- function(x, Y, V, w0, s0){
  data <- mmbr::DenseData$new(x, Y)
  data$standardize(F, F)
  mash_init <- mmbr::MashInitializer$new(S0, grid=1, prior_weights=w0, null_weight=0, top_mixtures=-1)
  B <- mmbr::MashRegression$new(1, V, mash_init)
  B$fit(data, save_var=T)

  return(mu1=B$posterior_b1, S1=B$posterior_variance, w1=B$mixture_posterior_weights, logbf=B$lbf)
}
