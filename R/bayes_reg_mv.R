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


# Bayesian multivariate multiple regression with mixture-of-normals prior
# (mixture weights w0 and covariance matrices S0) --> mr.mash
#
# The outputs are: the posterior assignment probabilities
# (w1), the posterior mean of the coefficients given that all the
# coefficients are not nonzero (mu1), and the posterior covariance of
# the coefficients given that all the coefficients are not zero (S1)
# the intercept (intercept), and the Evidence Lower Bound (ELBO).
mr.mash <- function(Y, X, V, S0, w0, mu_init = matrix(0, nrow=ncol(X), ncol=ncol(Y)), 
                    tol=1e-8, max_iter=1e5, update_w0=T, compute_ELBO=T) {
  ###Center Y and X
  Y <- scale(Y, center=T, scale=F)
  X <- scale(X, center=T, scale=F)
 
  if(compute_ELBO){ 
    ###Compute inverse of V (needed for the ELBO)
    Vinv <- solve(V)
  }
  
  ###Initilize mu1, S1, w1, error, ELBO and iterator
  p       <- ncol(X)
  n       <- nrow(X)
  R       <- ncol(Y)
  K       <- length(S0)
  mu1_t   <- mu_init 
  S1_t    <- array(0, c(R, R, p))
  w1_t    <- matrix(0, nrow=p, ncol=K)
  err     <- matrix(Inf, nrow=p, ncol=R)
  t       <- 0
  if(compute_ELBO){
    ELBO    <- -Inf
  }
  
  ###Repeat the following until convergence
  cat("iter beta_max.diff ELBO_diff ELBO\n")
  while(any(err>tol)){
    ##Initialize ELBO parameters
    var_part_ERSS <- 0
    neg_KL <- 0
    
    if(compute_ELBO){
      ##Set last value of ELBO as ELBO0
      ELBO0 <- ELBO
    }
    
    ##Update iterator
    t <- t+1
    
    ##Exit loop if maximum number of iterations is reached
    if(t>max_iter){
      warning("Max number of iterations reached. Try increasing max_iter.")
      break
    }
    
    ##Compute expected residuals
    rbar <- Y - X%*%mu1_t
    
    ##Save current estimates.
    mu1_tminus1 <- mu1_t
    
    ##Loop through the variables
    for(j in 1:p){
      
      #Remove j-th effect from expected residuals 
      rbar_j <- rbar + outer(X[, j], mu1_t[j, ])
      
      #Run Bayesian multivariate simple linear regression with mixture-of-normals prior
      bfit <- bayes_mvr_mix(X[, j], rbar_j, V, w0, S0)
      
      #Update variational parameters
      mu1_t[j, ]           <- bfit$mu1
      S1_t[, , j]          <- bfit$S1
      w1_t[j, ]            <- bfit$w1
      
      if(compute_ELBO){
        #Compute ELBO components
        # var_part_ERSS <- var_part_ERSS + (tr(Vinv%*%bfit$S1)*(sum(X[, j]^2)))
        # neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(Vinv%*%t(rbar_j)%*%matrix(X[, j], ncol=1)%*%matrix(bfit$mu1, nrow=1, byrow=T))+
        #                                        tr(Vinv%*%(bfit$S1+tcrossprod(matrix(bfit$mu1, ncol=1))))*(sum(X[, j]^2))))
        xtx <- sum(X[, j]^2)
        mu1_mat <- matrix(bfit$mu1, ncol=1)
      
        var_part_ERSS <- var_part_ERSS + (tr(Vinv%*%bfit$S1)*xtx)
        neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, rbar_j)%*%tcrossprod(matrix(X[, j], ncol=1), mu1_mat))+
                                               tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*(xtx)))
      }
      
      #Update expected residuals
      rbar <- rbar_j - outer(X[, j], mu1_t[j, ])
    }
    
    ##Update w0 if requested
    if(update_w0){
      w0 <- colSums(w1_t)/p
    }
    
    ##Compute distance in mu1 between two successive iterations
    err <- abs(mu1_t - mu1_tminus1)
    
    if(compute_ELBO){
      ##Compute ELBO
      ERSS <- tr(Vinv%*%(crossprod(rbar))) + var_part_ERSS
      ELBO <- -log(n)/2 -log(n*R)/2 -n/2 * as.numeric(determinant(V, logarithm = TRUE)$modulus) - 0.5*ERSS + neg_KL
    
      ##Print out useful info
      cat(sprintf("%4d %0.2e %0.2e %0.20e\n", t, max(err), ELBO - ELBO0, ELBO))
    } else {
      ##Print out useful info
      cat(sprintf("%4d %0.2e\n", t, max(err)))
    }
  }
  
  ###Compute intercept
  Ybar <- attr(Y, 'scaled:center')
  Xbar <- matrix(rep(attr(X, 'scaled:center'), each=ncol(mu1_t)), ncol=ncol(mu1_t), byrow=T)
  intercept <- Ybar - colSums(Xbar * mu1_t)
  
  if(compute_ELBO){
    ###Return the posterior assignment probabilities (w1), the posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), the intercept (intercept), and the Evidence Lower Bound (ELBO).    
    return(list(mu1=mu1_t, S1=S1_t, w1=w1_t, intercept=intercept, ELBO=ELBO))
  } else {
    ###Return the posterior assignment probabilities (w1), the posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), and the intercept (intercept).    
    return(list(mu1=mu1_t, S1=S1_t, w1=w1_t, intercept=intercept))
  }
}
