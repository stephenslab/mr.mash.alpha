# Run several iterations of the co-ordinate ascent updates for the
# mr-mash model, allowing for missing values in Y.
#
# The special case of univariate regression, when Y is a vector or a
# matrix with 1 column, is also handled.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
mr_mash_simple <- function (X, Y, V, S0, w0, B, numiter = 100,
                            tol=1e-4, update_w0=TRUE, update_V=FALSE,
                            verbose=FALSE) {
  
  r <- ncol(Y)
  n <- nrow(Y)
  
  Y_has_missing <- any(is.na(Y))
  
  if(Y_has_missing){
    # Store missingness patterns for each individual
    miss <- vector("list", n)
    non_miss <- vector("list", n)
    for(i in 1:n){
      miss_i <- is.na(Y[i, ])
      non_miss_i <- !miss_i
      miss[[i]] <- miss_i
      non_miss[[i]] <- non_miss_i
    }
    
    # Initialize missing Ys and intercept
    intercept <- colMeans(Y, na.rm=TRUE)
    for(l in 1:r){
      Y[is.na(Y[, l]), l] <- intercept[l]
    }
  } else {
    # Center Y
    Y <- scale(Y, scale=FALSE)
    muy <- attr(Y,"scaled:center")
  }
  
  # ty <- 0
  
  # Center X
  X <- scale(X, scale=FALSE)
  mux <- attr(X,"scaled:center")
  
  # These variables is used to keep track of the algorithm's progress.
  maxd <- rep(0,numiter)
  ELBO <- rep(0,numiter)
  
  # Iterate the updates.
  for (t in 1:numiter) {
    
    # Save the current estimates of the posterior means.
    B0 <- B
    
    if(Y_has_missing){
      # Compute expected Y
      mu <- t(t(X%*%B) + intercept)
    } else {
      # Compute expected Y
      mu <- X%*%B
      
      # Set quantities needed for computing the ELBO with missing Ys to 0
      y_var <- matrix(0, r, r)
      sum_entropy_Y <- 0
    }
    
    # M-step, if not the first iteration
    if(t!=1){
      # Update mixture weights, if requested
      if(update_w0){
        w0 <- colSums(out$W1)
        w0 <- w0/sum(w0)
      }
      
      # Update V, if requested
      if(update_V){
        R <- Y - mu
        ERSS <- crossprod(R) + out$var_part_ERSS + y_var
        V <- ERSS/n
      }
    }
    
    if(Y_has_missing){
      # Impute missing Y (code adapted from Yuxin)
      y_var <- matrix(0, r, r)
      sum_entropy_Y <- 0
      Vinv <- chol2inv(chol(V)) 
      
      for (i in 1:n){
        non_miss_i <- non_miss[[i]]
        miss_i <- miss[[i]]
        Vinv_mo = Vinv[miss_i, non_miss_i, drop=FALSE]
        Vinv_mm = Vinv[miss_i, miss_i, drop=FALSE]
        if(any(miss_i)){
          # Compute variance
          y_var_i = matrix(0, r, r)
          y_var_mm <- chol2inv(chol(Vinv_mm))
          y_var_i[miss_i, miss_i] = y_var_mm
          
          y_var <- y_var + y_var_i
          
          # Compute mean
          imp_mean = mu[i, miss_i] - y_var_mm %*% Vinv_mo %*% (Y[i, non_miss_i] - mu[i, non_miss_i])
          Y[i, miss_i] = imp_mean
          
          # Compute sum of the negative entropy of Y missing
          sum_entropy_Y = sum_entropy_Y + (0.5 * as.numeric(determinant(1/(2*pi*exp(1))*Vinv_mm, logarithm = TRUE)$modulus))
        }
      }
      
      # t2 <- proc.time()
      # ty <- ty + t2["elapsed"] - t1["elapsed"]
      
      # Update the intercept
      intercept <- colMeans(Y)
      
      # E-step: Update the posterior means of the regression coefficients.
      out <- mr_mash_update_simple(X,scale(Y, scale=FALSE),B,V,w0,S0, y_var, sum_entropy_Y)
    } else {
      out <- mr_mash_update_simple(X,Y,B,V,w0,S0, y_var, sum_entropy_Y)
    }
    
    
    B <- out$B 
    
    # Store the largest change in the posterior means.
    delta_B <- abs(max(B - B0))
    maxd[t] <- delta_B
    ELBO[t] <- out$ELBO
    
    # Print out some info, if requested
    if(verbose)
      cat("Iter:", t, "    max_deltaB:", delta_B, "\n")
    
    # Break if convergence is reached
    if(delta_B<tol)
      break
  }
  
  # Compute the intercept
  if(Y_has_missing){
    intercept <- drop(intercept - mux %*% B)
  } else {
    intercept <- drop(muy - mux %*% B)
  }
  
  # Return the updated posterior means of the regression coefficicents
  # (B), the maximum change at each iteration (maxd), the prior weights,
  # and V.
  return(list(intercept = intercept,B = B,maxd = maxd,w0 = w0,
         V = V,ELBO = ELBO[1:t],Y = Y))
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
mr_mash_update_simple <- function (X, Y, B, V, w0, S0, yvar, sum_ent_Y) {
  
  # Make sure B is a matrix.
  B <- as.matrix(B)
  
  # Get the number of predictors, responses, and mixture components.
  p <- ncol(X)
  k <- length(w0)
  r <- ncol(Y)
  n <- nrow(Y)
  
  # Create matrix to store posterior assignment probabilities
  W1 <- matrix(as.numeric(NA), p, k)
  
  # Initialize quantities need to update V and ELBO
  var_part_tr_wERSS <- 0
  neg_KL <- 0
  var_part_ERSS <- matrix(0, nrow=r, ncol=r)
  
  # Compute inverse of V
  Vinv <- chol2inv(chol(V))
  
  # Compute the expected residuals.
  R <- Y - X %*% B
  
  # Repeat for each predictor.
  for (i in 1:p) {
    x <- X[,i]
    xtx <- sum(x^2)
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
    var_part_ERSS <- var_part_ERSS + out$S1*xtx
    
    # Update quantities needed for the ELBO
    b_mat <- matrix(b, ncol=1)
    var_part_tr_wERSS <- var_part_tr_wERSS + (tr(Vinv%*%out$S1)*xtx)
    neg_KL <- neg_KL + (out$logbf +0.5*(-2*tr(tcrossprod(Vinv, R)%*%tcrossprod(matrix(x, ncol=1), b_mat))+
                                          tr(Vinv%*%(out$S1+tcrossprod(b_mat)))*xtx))
    # Update the expected residuals.
    R <- R - outer(x,b)
  }
  
  # Compute the ELBO
  tr_wERSS <- tr(Vinv%*%(crossprod(R))) + var_part_tr_wERSS
  if(all(yvar==0)){
    e2 <- 0
  } else {
    e2 <- tr(Vinv%*%yvar)
  }
  ELBO <- -log(n)/2 - (n*r)/2*log(2*pi) - n/2 * as.numeric(determinant(V, logarithm = TRUE)$modulus) - 
    0.5*(tr_wERSS+e2) + neg_KL - sum_ent_Y
  
  # Output the updated predictors.
  return(list(B=drop(B), W1=W1, var_part_ERSS=var_part_ERSS, ELBO=drop(ELBO)))
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
  mvtnorm::dmvnorm(x,sigma = S,log = TRUE)
}
