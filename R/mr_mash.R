#' @title  Multiple Regression with Multivariate Adaptive Shrinkage.
#' @details Performs multivariate multiple regression with mixture-of-normals prior.
#' 
#' @param Y an NxR matrix of responses.
#' @param X an NxP matrix of covariates.
#' @param V an RxR residual covariance matrix.
#' @param S0 a list of length K containing the desired RxR prior covariance matrices 
#' on the regression coefficients.
#' @param w0 a K-vector with prior mixture weights, each associated with the 
#' respective covariance matrix in \code{S0}.   
#' @param mu_init a PxR matrix of initial estimates for the regression coefficients.
#' @param tol convergence tolerance.
#' @param max_iter maximum number of iterations for the optimization algorithm.
#' @param update_w0 if TRUE, prior weights are updated.
#' @param compute_ELBO if TRUE, ELBO is computed.
#' 
#' @return a mr.mash fit, which is a list with some or all of the following elements\cr
#' \item{mu1}{a PxR matrix of posterior means for the regression coeffcients}
#' \item{S1}{a RxRxP array of posterior covariances for the regression coeffcients}
#' \item{mu1}{a PxK matrix of posterior assignment probabilities to the mixture components}
#' \item{intercept}{an R-vector with the estimated intercepts}
#' \item{ELBO}{the Evidence Lower Bound at convergence} 
#' 
#' @export
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
