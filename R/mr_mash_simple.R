# Run several iterations of the co-ordinate ascent updates for the
# mr-mash model allowing for missing values in Y.
#
# The special case of univariate regression, when Y is a vector or a
# matrix with 1 column, is also handled.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_mash_simple_missing_Y <- function (X, Y, V, S0, w0, B, numiter = 100,
                                      tol=1e-4, update_w0=TRUE, update_V=FALSE,
                                      verbose=FALSE) {
  
  r <- ncol(Y)
  n <- nrow(Y)
  
  # Get missing Y patterns (code from Yuxin)
  Y_missing <- is.na(Y)
  Y_non_missing <- !Y_missing
  missing_pattern <- unique(Y_non_missing)
  Y_missing_pattern_assign <- numeric(n)
  for(z in 1:nrow(missing_pattern)){
    idx = which(apply(Y_non_missing, 1, function(x) identical(x, missing_pattern[z,])))
    Y_missing_pattern_assign[idx] <- z
  }
  Y[Y_missing] <- 0
  
  # Compute inverse of V for each missing pattern (code from Yuxin)
  V_inv <- list()
  for(z in 1:nrow(missing_pattern)){
    if(r == 1){
      V_inv[[z]] <- missing_pattern[z,] / V
    }else{
      Vz = V[which(missing_pattern[z,]), which(missing_pattern[z,])]
      eigenVz <- eigen(Vz, symmetric = TRUE)
      dinv <- 1/(eigenVz$values)
      V_inv[[z]] <- eigenVz$vectors %*% (dinv * t(eigenVz$vectors))
    }
  }
  
  # Center X
  X <- scale(X, scale=FALSE)
  mux <- attr(X,"scaled:center")
  
  # This variable is used to keep track of the algorithm's progress.
  maxd <- rep(0,numiter)
  
  # Iterate the updates.
  for (t in 1:numiter) {
    
    # Save the current estimates of the posterior means.
    B0 <- B
    
    # M-step, if not the first iteration
    if(t!=1){
      # Update mixture weights, if requested
      if(update_w0){
        w0 <- colSums(out$W1)
        w0 <- w0/sum(w0)
      }
      
      # # Update V, if requested
      # if(update_V){
      #   R <- Y - X%*%B
      #   ERSS <- crossprod(R) + out$var_part_ERSS
      #   V <- ERSS/nrow(Y)
      # }
    }
    
    # Impute missing Y (code from Yuxin)
    for (i in 1:n){
      missing_pattern_i = !missing_pattern[Y_missing_pattern_assign[i],]
      V_mo = V[missing_pattern_i, !missing_pattern_i, drop=FALSE]
      if(any(missing_pattern_i)){
        imp_mean = B[i, missing_pattern_i] + V_mo %*% V_inv[[Y_missing_pattern_assign[i]]] %*%
          (Y[i, !missing_pattern_i] - B[i, !missing_pattern_i])
        Y[i,missing_pattern_i] = imp_mean
      }
    }
      
    # Center X
    Y <- scale(Y, scale=FALSE)
    muy <- attr(Y,"scaled:center")
    
    # E-step: Update the posterior means of the regression coefficients.
    out <- mr_mash_update_simple(X,Y,B,V,w0,S0)
    B <- out$B 
    
    # Store the largest change in the posterior means.
    delta_B <- abs(max(B - B0))
    maxd[t] <- delta_B
      
    # Print out some info, if requested
    if(verbose)
      cat("Iter:", t, "    max_deltaB:", delta_B, "\n")
      
    # Break if convergence is reached
    if(delta_B<tol)
      break
  }
    
  # Compute the intercept
  intercept <- drop(muy - mux %*% B)
    
  # Return the updated posterior means of the regression coefficicents
  # (B), the maximum change at each iteration (maxd), the prior weights,
  # and V.
  return(list(intercept=intercept,B = B,maxd = maxd,w0=w0,V=V))
}
  

# Run several iterations of the co-ordinate ascent updates for the
# mr-mash model.
#
# The special case of univariate regression, when Y is a vector or a
# matrix with 1 column, is also handled.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_mash_simple <- function (X, Y, V, S0, w0, B, numiter = 100,
                            tol=1e-4, update_w0=TRUE, update_V=TRUE,
                            verbose=FALSE) {
  
  #Center X and Y
  X <- scale(X, scale=FALSE)
  Y <- scale(Y, scale=FALSE)
  mux <- attr(X,"scaled:center")
  muy <- attr(Y,"scaled:center")
  
  # This variable is used to keep track of the algorithm's progress.
  maxd <- rep(0,numiter)

  # Iterate the updates.
  for (t in 1:numiter) {

    # Save the current estimates of the posterior means.
    B0 <- B
    
    # M-step, if not the first iteration
    if(t!=1){
      # Update mixture weights, if requested
      if(update_w0){
        w0 <- colSums(out$W1)
        w0 <- w0/sum(w0)
      }
      
      # Update V, if requested
      if(update_V){
        R <- Y - X%*%B
        ERSS <- crossprod(R) + out$var_part_ERSS
        V <- ERSS/nrow(Y)
      }
    }
      
    # E-step: Update the posterior means of the regression coefficients.
    out <- mr_mash_update_simple(X,Y,B,V,w0,S0)
    B <- out$B 
    
    # Store the largest change in the posterior means.
    delta_B <- abs(max(B - B0))
    maxd[t] <- delta_B
    
    # Print out some info, if requested
    if(verbose)
      cat("Iter:", t, "    max_deltaB:", delta_B, "\n")
    
    # Break if convergence is reached
    if(delta_B<tol)
      break
  }
  
  # Compute the intercept
  intercept <- drop(muy - mux %*% B)

  # Return the updated posterior means of the regression coefficicents
  # (B), the maximum change at each iteration (maxd), the prior weights,
  # and V.
  return(list(intercept=intercept,B = B,maxd = maxd,w0=w0,V=V))
}

# Perform a single pass of the co-ordinate ascent updates for the
# mr-mash model.
#
# The special case of univariate regression, when Y is a vector or a
# matrix with 1 column, is also handled.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_mash_update_simple <- function (X, Y, B, V, w0, S0) {

  # Make sure B is a matrix.
  B <- as.matrix(B)
  
  # Get the number of predictors, responses, and mixture components.
  p <- ncol(X)
  k <- length(w0)
  r <- ncol(Y)

  # Create matrix to store posterior assignment probabilities
  W1 <- matrix(as.numeric(NA), p, k)
  
  # Initialize var_part_ERSS
  var_part_ERSS <- matrix(0, r, r)
  
  # Compute the expected residuals.
  R <- Y - X %*% B

  # Repeat for each predictor.
  for (i in 1:p) {
    x <- X[,i]
    b <- B[i,]
    
    # Disregard the ith predictor in the expected residuals.
    R <- R + outer(x,b)

    # Update the posterior of the regression coefficients for the ith
    # predictor.
    out   <- bayes_mvr_mix_simple(x,R,V,w0,S0)
    b     <- out$mu1
    B[i,] <- b
    
    # Update the posterior assignment probabilities for the ith predictor
    W1[i,] <- out$w1
    
    # Update quantity needed to update V
    var_part_ERSS <- var_part_ERSS + out$S1*sum(x^2)
    
    # Update the expected residuals.
    R <- R - outer(x,b)
  }

  # Output the updated predictors.
  return(list(B=drop(B), W1=W1, var_part_ERSS=var_part_ERSS))
}

# Compute quantities for a basic Bayesian multivariate regression with
# a multivariate normal prior on the regression coefficients: Y = xb'
# + E, E ~ MN(0,I,V), b ~ N(0,S0). The outputs are: bhat, the
# least-squares estimate of the regression coefficients; S, the
# covariance of bhat; mu1, the posterior mean of the regression
# coefficients; S1, the posterior covariance of the regression
# coefficients; and logbf, the logarithm of the Bayes factor.
#
# The special case of univariate regression, when Y is a vector or a
# matrix with 1 column, is also handled.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
bayes_mvr_ridge_simple <- function (x, Y, V, S0) {

  # Make sure Y, V and S0 are matrices.
  Y  <- as.matrix(Y)
  V  <- as.matrix(V)
  S0 <- as.matrix(S0)
    
  # Compute the least-squares estimate of the coefficients (bhat) and
  # the covariance of the standard error (S).
  xx   <- norm2(x)^2
  bhat <- drop(x %*% Y)/xx
  S    <- V/xx
  
  # Compute the posterior mean (mu1) and covariance (S1) assuming a
  # multivariate normal prior with zero mean and covariance S0.
  r   <- ncol(Y)
  I   <- diag(r)
  S1  <- S0 %*% solve(I + solve(S) %*% S0)
  mu1 <- drop(S1 %*% solve(S,bhat))

  # Compute the log-Bayes factor.
  logbf <- ldmvnorm(bhat,S0 + S) - ldmvnorm(bhat,S)
  
  # Return the least-squares estimate (bhat) and its covariance (S), the
  # posterior mean (mu1) and covariance (S1), and the log-Bayes factor
  # (logbf).
  return(list(bhat  = bhat,
              S     = drop(S),
              mu1   = mu1,
              S1    = drop(S1),
              logbf = logbf))
}

# Compute quantities for Bayesian multivariate regression with a
# mixture-of-multivariate-normals prior on the regression
# coefficients. mu1, the posterior mean of the regression
# coefficients; S1, the posterior covariance of the regression
# coefficients; w1, the posterior "weights" for the individual
# components; and logbf, the logarithm of the Bayes factor.
#
# The special case of univariate regression, when Y is a vector or a
# matrix with 1 column, is also handled.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
bayes_mvr_mix_simple <- function (x, Y, V, w0, S0) {
    
  # Make sure Y is a matrix.
  Y <- as.matrix(Y)
  
  # Get the dimension of the response (r) and the number of mixture
  # components (k).
  r <- ncol(Y)
  k <- length(w0)

  # Compute the quantities separately for each mixture component.
  out <- vector("list",k)
  for (i in 1:k)
    out[[i]] <- bayes_mvr_ridge_simple(x,Y,V,S0[[i]])

  # Compute the posterior assignment probabilities for the latent
  # indicator variable.
  logbf <- sapply(out,"[[","logbf")
  w1    <- softmax(logbf + log(w0))

  # Compute the log-Bayes factor as a linear combination of the
  # individual Bayes factors for each mixture component.
  z     <- log(w0) + logbf
  u     <- max(z)
  logbf <- u + log(sum(exp(z - u)))
  
  # Compute the posterior mean (mu1) and covariance (S1) of the
  # regression coefficients.
  S1  <- matrix(0,r,r)
  mu1 <- rep(0,r)
  for (i in 1:k) {
    w   <- w1[i]
    mu  <- out[[i]]$mu1
    S   <- out[[i]]$S1
    mu1 <- mu1 + w*mu
    S1  <- S1 + w*(S + tcrossprod(mu))
  }
  S1 <- S1 - tcrossprod(mu1)
  
  # Return the the posterior mean (mu1) and covariance (S1), the
  # posterior assignment probabilities (w1), and the log-Bayes factor
  # (logbf).
  return(list(mu1   = mu1,
              S1    = drop(S1),
              w1    = w1,
              logbf = logbf))
}

# Return the log-density of the multivariate normal with zero mean
# and covarirance S at x.
#
#' @importFrom mvtnorm dmvnorm
ldmvnorm <- function (x, S){
  dmvnorm(x,sigma = S,log = TRUE)
}